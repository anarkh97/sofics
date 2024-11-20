#include "DistrExplicitUDEIMPodProjectionNonLinDynamic.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>

#include <Driver.d/GeoSource.h>

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "VecNodeDof6Map.h"
#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"
#include "PtrPtrIterAdapter.h"

#include <utility>
#include <cstddef>
#include <algorithm>

#include <sys/time.h>

extern GeoSource *geoSource;

namespace Rom {

DistrExplicitUDEIMPodProjectionNonLinDynamic::DistrExplicitUDEIMPodProjectionNonLinDynamic(Domain *domain) :
  DistrExplicitLumpedPodProjectionNonLinDynamic(domain)
{}

void
DistrExplicitUDEIMPodProjectionNonLinDynamic::preProcess() {

  DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess();

  kelArrayCopy = new FullSquareMatrix*[decDomain->getNumSub()];

  numOfIndices = domain->solInfo().maxDeimBasisSize;

  {
   unassembledInfo.domLen = new int[MultiDomainDynam::solVecInfo().numDom];
   unassembledInfo.numDom = MultiDomainDynam::solVecInfo().numDom;
   int totLen = 0;
   for(int iSub = 0; iSub < MultiDomainDynam::solVecInfo().numDom; ++iSub) {
     unassembledInfo.domLen[iSub] = (iSub==0) ? numOfIndices : 0;
     totLen += unassembledInfo.domLen[iSub];
   }

   unassembledInfo.len = totLen;
   unassembledInfo.setMasterFlag();
  }

  buildInterpolationBasis();
  buildReducedLinearOperator();

}

void
DistrExplicitUDEIMPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {
  
  DistrVector uFint(unassembledInfo,0);

  execParal(decDomain->getNumSub(),this,&DistrExplicitUDEIMPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly, uFint,t,tIndex);

  if(!domain->solInfo().reduceFollower) 
    execParal(decDomain->getNumSub(),this,&DistrExplicitUDEIMPodProjectionNonLinDynamic::subGetFollowerForceOnly,*fExt,t,tIndex);

  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    GenMDDynamMat<double> ops;
    ops.K = K;
    decDomain->rebuildOps(ops, 0.0, 0.0, 0.0, kelArray);
  }
  
  if (domain->solInfo().filterFlags) {
    trProject(*fInt);
  }

  DistrVector fExt_reduced(solVecInfo());
  DistrVector fLin_reduced(reducedVecInfo());

  normalizedBasis_.sparseVecReduce(*fExt,fExt_reduced);
  ReducedStiffness.fullExpand(d,fLin_reduced); //(V^T*K*V)*d_reduced
  udeimBasis_.compressedVecReduce(uFint,f);
  f -= fExt_reduced;
  f += fLin_reduced;
}

void
DistrExplicitUDEIMPodProjectionNonLinDynamic::subGetKtimesU(int isub, DistrVector &d, DistrVector &f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subf(f.subData(isub), f.subLen(isub));
  StackVector subd(d.subData(isub), d.subLen(isub));
  sd->getKtimesU(subd, (double *) 0, subf, 1.0, (kelArrayCopy) ? kelArrayCopy[isub] : (FullSquareMatrix *) 0);
}

void
DistrExplicitUDEIMPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)
  
  sd->getUDEIMInternalForceOnly(unassembledElemDOFMask[iSub], *(*geomState)[iSub], eIF,
                                allCorot[iSub], kelArray[iSub], residual, fExt->subLen(iSub),
                                1.0, t, (*geomState)[iSub], melArray[iSub],kelArrayCopy[iSub]); // residual -= internal force);     

  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, -1.0); // f = -residual
} 

void
DistrExplicitUDEIMPodProjectionNonLinDynamic::subGetFollowerForceOnly(int iSub, DistrVector &f, double &t, int &tIndex) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)

  sd->getFollowerForce(*(*geomState)[iSub], eIF, allCorot[iSub], kelArray[iSub], residual, 1.0, t, NULL, false);

  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, 1.0); // f = -residual
}
 
void
DistrExplicitUDEIMPodProjectionNonLinDynamic::buildInterpolationBasis() {

 filePrint(stderr, " ... Reading Interpolation basis ...\n");

 filePrint(stderr, " ... Interplation subspace of dimension %d x %d ...\n", numOfIndices, domain->solInfo().forcePodSize);
 udeimBasis_.dimensionIs(normalizedBasis_.numVec(),unassembledInfo);//we will read in the V^T*U*(P^T*U)*P^T in R(n x m)

 int vector = 0; int ind = 0;
  for(std::vector<double>::const_iterator it = geoSource->UDEIMVecBegin(), it_end = geoSource->UDEIMVecEnd(); it != it_end; ++it){
    udeimBasis_[vector][ind] = *it; ++ind;
    if(ind == domain->solInfo().maxDeimBasisSize){//each process get a copy of the UDEIM basis
      ind = 0;
      ++vector;
    }
  }

 //if we are doing UDEIM, then we must only compute the DOFs we need in each element to prevent conflicts with neighboring elements
 columnKey.resize(decDomain->getNumSub());
 unassembledElemDOFMask.resize(decDomain->getNumSub());
 execParal(decDomain->getNumSub(),this,&DistrExplicitUDEIMPodProjectionNonLinDynamic::subBuildUnassembledMask);

 filePrint(stderr," ... Compressing Interpolation Basis ...\n");
 udeimBasis_.makeSparseBasis(columnKey); 

}

void
DistrExplicitUDEIMPodProjectionNonLinDynamic::buildReducedLinearOperator() {
 //build reduced stiffness matrix
 ReducedStiffness.dimensionIs(normalizedBasis_.numVectors(),reducedVecInfo()); //each mpi process gets a reduced linear operator

  if(domain->solInfo().ReducedStiffness){
    filePrint(stderr," ... Reading Pre-computed Reduced Linear Stiffness Matrix ...\n");

    int vector = 0; int ind = 0;
    for(std::vector<double>::const_iterator it = geoSource->RedKVecBegin(), it_end = geoSource->RedKVecEnd(); it != it_end; ++it){
      ReducedStiffness[vector][ind] = *it; ++ind;
      if(ind == normalizedBasis_.numVectors()){
        ind = 0; 
        ++vector;
      }
    }

  } else {
    filePrint(stderr," ... Constructing Reduced Linear Stiffness Matrix ...\n");
 
    for( int column = 0; column != normalizedBasis_.numVectors(); ++column){
      DistrVector columnOfKtimesV(MultiDomainDynam::solVecInfo());
      columnOfKtimesV = 0;
      //K*V
      execParal(decDomain->getNumSub(),this,&DistrExplicitUDEIMPodProjectionNonLinDynamic::subGetKtimesU, normalizedBasis_[column],columnOfKtimesV); 
      //V^T*(K*V)
      normalizedBasis_.reduce(columnOfKtimesV,ReducedStiffness[column]);
    } 
  }

}

void
DistrExplicitUDEIMPodProjectionNonLinDynamic::subBuildUnassembledMask(int iSub) {
  //create mask for unassembled DEIM 
  SubDomain *sd = decDomain->getSubDomain(iSub);
  std::map<int, std::vector<int> > &subUnassembledMaskBuf = unassembledElemDOFMask[iSub];
  std::map<int,std::vector<int> > &subColumnKey = columnKey[iSub];

  sd->createKelArray(kelArrayCopy[iSub]);

  int row = 0;
  for (GeoSource::ElemDofPairVec::const_iterator it = geoSource->elemDofBegin(),
                                                   it_end = geoSource->elemDofEnd();
                                                   it != it_end; ++it) {

     const int elemId   = it->first;
     const int packedId = sd->glToPackElem(elemId);
     const int elemDOF  = it->second;

     row += 1;

     if(packedId < 0) {continue;}
    
     //unassembled element contibution are added to residual in the order the are iterated through
     if(subColumnKey.find(packedId) != subColumnKey.end()){
       subColumnKey[packedId].push_back(row-1);
     }else{
       std::vector<int> dummyVec;
       dummyVec.push_back(row-1);
       subColumnKey.insert(std::make_pair(packedId,dummyVec));
     }
  
     if(subUnassembledMaskBuf[packedId].size() > 0){
       subUnassembledMaskBuf[packedId].push_back(it->second);       
     }else{
       std::vector<int> DOFs;
       subUnassembledMaskBuf.insert(std::make_pair(packedId,DOFs));
       subUnassembledMaskBuf[packedId].push_back(elemDOF);
     }
  }

}

} // end namespace Rom
