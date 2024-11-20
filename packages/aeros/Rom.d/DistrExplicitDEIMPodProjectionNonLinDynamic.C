#include "DistrExplicitDEIMPodProjectionNonLinDynamic.h"

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

DistrExplicitDEIMPodProjectionNonLinDynamic::DistrExplicitDEIMPodProjectionNonLinDynamic(Domain *domain) :
  DistrExplicitLumpedPodProjectionNonLinDynamic(domain)
{}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::preProcess() {

  DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess();

  kelArrayCopy = new FullSquareMatrix*[decDomain->getNumSub()];

  buildInterpolationBasis();
  buildReducedLinearOperator();
}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {

  execParal(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly,*fInt,t,tIndex);

  if(!domain->solInfo().reduceFollower) 
    execParal(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subGetFollowerForceOnly,*fExt,t,tIndex);

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
  deimBasis_.compressedVecReduce(*fInt,f);
  f -= fExt_reduced;
  f += fLin_reduced;
}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::subGetKtimesU(int isub, DistrVector &d, DistrVector &f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subf(f.subData(isub), f.subLen(isub));
  StackVector subd(d.subData(isub), d.subLen(isub));
  sd->getKtimesU(subd, (double *) 0, subf, 1.0, (kelArrayCopy) ? kelArrayCopy[isub] : (FullSquareMatrix *) 0);
}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)
  
  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    sd->getWeightedStiffAndForceOnly(packedElementWeights_[localReducedMeshId_][iSub], *(*geomState)[iSub], eIF,
                                     allCorot[iSub], kelArray[iSub], residual,
                                     1.0, t, (*geomState)[iSub], melArray[iSub]); // residual -= internal force);
  }
  else {
      sd->getWeightedInternalForceOnly(packedElementWeights_[localReducedMeshId_][iSub], *(*geomState)[iSub], eIF,
                                    allCorot[iSub], kelArray[iSub], residual,
                                    1.0, t, (*geomState)[iSub], melArray[iSub],kelArrayCopy[iSub]); // residual -= internal force);
  }

  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, -1.0); // f = -residual
} 

void
DistrExplicitDEIMPodProjectionNonLinDynamic::subGetFollowerForceOnly(int iSub, DistrVector &f, double &t, int &tIndex) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)

  sd->getFollowerForce(*(*geomState)[iSub], eIF, allCorot[iSub], kelArray[iSub], residual, 1.0, t, NULL, false);

  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, 1.0); // f = -residual
}
 
void
DistrExplicitDEIMPodProjectionNonLinDynamic::buildInterpolationBasis() {

 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, BasisId::FORCE, BasisId::POD);
 fileName = fileName + ".deim";

 DistrBasisInputFile BasisFile(fileName); //read in mass-normalized basis
 filePrint(stderr, " ... Reading Interpolation basis from file %s ...\n", fileName.c_str());
 const int interpBasisSize = domain->solInfo().maxSizePodRom ?
                             std::min(domain->solInfo().maxSizePodRom, BasisFile.stateCount()) :
                             BasisFile.stateCount();

 filePrint(stderr, " ... Interplation subspace of dimension = %d ...\n", interpBasisSize);
 deimBasis_.dimensionIs(interpBasisSize,decDomain->masterSolVecInfo());

 DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

 typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
 DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                  SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
 DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

 for (DistrVecBasis::iterator it = deimBasis_.begin(),
                              it_end = deimBasis_.end();
                              it != it_end; ++it) {
   assert(BasisFile.validCurrentState());
   
  
   BasisFile.currentStateBuffer(buffer);

   converter.vector(buffer, *it);
   BasisFile.currentStateIndexInc();
 }
 
 std::vector< std::vector<std::pair<int,int> > > maskedIndicesBuf;
 maskedIndicesBuf.resize(decDomain->getNumSub());
 execParal(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subBuildInterpolationBasis, maskedIndicesBuf); 
 
 filePrint(stderr," ... Compressing Interpolation Basis ...\n");
 DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
 for(int i=0; i<decDomain->getNumSub(); ++i) {all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();}
 deimBasis_.makeSparseBasis(maskedIndicesBuf, all_cdsa); 

}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::buildReducedLinearOperator() {
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
      execParal(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subGetKtimesU, normalizedBasis_[column],columnOfKtimesV); 
      //V^T*(K*V)
      normalizedBasis_.reduce(columnOfKtimesV,ReducedStiffness[column]);
    } 
  }

/*
 DistrVecBasis DEIMReducedStiffness(normalizedBasis_.numVectors(),reducedVecInfo());

 for( int column = 0; column != normalizedBasis_.numVectors(); ++column){
   DistrVector columnOfKtimesV(MultiDomainDynam::solVecInfo());
   columnOfKtimesV = 0;
   execParal(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subGetKtimesU, deimBasis_[column],columnOfKtimesV);
   deimBasis_.reduce(columnOfKtimesV,DEIMReducedStiffness[column]);
 }

 Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > symmetricK(ReducedStiffness.data(),ReducedStiffness.size(),ReducedStiffness.numVectors());
 Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > antiSymmetricK(DEIMReducedStiffness.data(),DEIMReducedStiffness.size(),DEIMReducedStiffness.numVectors());

 Eigen::EigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > EigSolve(symmetricK);
 std::cout << "symmetric K: \n" << symmetricK-symmetricK.transpose() << std::endl;
 std::cout << "symmetric K eigenvalues : \n" << EigSolve.eigenvalues() << std::endl;

 Eigen::EigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > EigSolve2(antiSymmetricK);
 std::cout << "anti-symmetric K: \n" << antiSymmetricK << std::endl;
 std::cout << "anti-symmetric K difference: \n" << antiSymmetricK-antiSymmetricK.transpose() << std::endl;
 std::cout << "anti-symmetric K eigenvalues : \n" << EigSolve2.eigenvalues() << std::endl;
*/
}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::subBuildInterpolationBasis(int iSub, std::vector< std::vector<std::pair<int,int> > > &maskedIndicesBuf) {

  SubDomain *sd = decDomain->getSubDomain(iSub);
  std::vector<std::pair<int,int> > &subMaskedIndicesBuf = maskedIndicesBuf[iSub];

  for (GeoSource::NodeDofPairVec::const_iterator it = geoSource->nodeDofSlotBegin(),
                                                   it_end = geoSource->nodeDofSlotEnd();
                                                   it != it_end; ++it) {
    
     const int nodeId = it->first;
     const int packedId = sd->globalToLocal(nodeId);
   
     if(packedId < 0) {continue;}

     subMaskedIndicesBuf.push_back(std::make_pair(packedId,it->second));

  }

  sd->createKelArray(kelArrayCopy[iSub]);
 
}

} // end namespace Rom
