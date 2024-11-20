#include "UDEIMPodProjectionNonLinDynamic.h"

#include <Problems.d/NonLinDynam.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>

#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "NodeDof6Buffer.h"
#include "VecNodeDof6Conversion.h"
 
#include "VecBasisFile.h"
#include "PodProjectionSolver.h"

#include <utility>

extern GeoSource *geoSource;

namespace Rom {

UDEIMPodProjectionNonLinDynamic::UDEIMPodProjectionNonLinDynamic(Domain *domain) :
  LumpedPodProjectionNonLinDynamic(domain)
{}

void
UDEIMPodProjectionNonLinDynamic::preProcess() {

  LumpedPodProjectionNonLinDynamic::preProcess();

  kelArrayCopy = new FullSquareMatrix;
  domain->createKelArray(kelArrayCopy);

  numOfIndices = domain->solInfo().maxDeimBasisSize;

  readInterpolationBasis(); 
  buildSparseInterpolationBasis();
  buildReducedLinearOperator();

}

void
UDEIMPodProjectionNonLinDynamic::reBuild(ModalGeomState &geomState, int iteration, double localDelta, double t)
{

  double Kcoef;

  if(t == domain->solInfo().initialTime) {
     Kcoef = 0;
   }
   else {
     double beta, gamma, alphaf, alpham, dt = 2*localDelta;
     getNewmarkParameters(beta, gamma, alphaf, alpham);
     Kcoef = (domain->solInfo().order == 1) ? dt*gamma : dt*dt*beta;
   }

  spm->zeroAll();
  computeRedKDyn(Kcoef);
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > LinearStiffnessMatrixMap(LinearStiffness_.data(), LinearStiffness_.vectorSize(), LinearStiffness_.numVec());
  getSolver()->addToReducedMatrix(LinearStiffnessMatrixMap,Kcoef);
#endif
  NonLinDynamic::reBuild(*geomState_Big, iteration, localDelta, t);
}

void
UDEIMPodProjectionNonLinDynamic::computeRedKDyn(double Kcoef){

 GenVecBasis<double> &projectionBasis = const_cast<GenVecBasis<double> &>(dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis());
 GenVecBasis<double> RedKDyn(projectionBasis.numVec(),projectionBasis.numVec());
 for(int column = 0; column != projectionBasis.numVec(); column++) {
   GenVector<double> columnOfKtimesV(numOfIndices);
   columnOfKtimesV = 0.0;
   //(K_dyn - K_stiff)*V
   domain->getUnassembledKtimesU(unassembledElemDOFMask, projectionBasis[column], (double *) 0, columnOfKtimesV, Kcoef, kelArray);
   udeimBasis_.compressedVecReduce(columnOfKtimesV,RedKDyn[column]);
 }
#ifdef USE_EIGEN3
 Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > RKDMap(RedKDyn.data(),RedKDyn.size(),RedKDyn.numVec());
 getSolver()->addToReducedMatrix(RKDMap);
#endif
}

double
UDEIMPodProjectionNonLinDynamic::getStiffAndForce(ModalGeomState &geomState, Vector &residual,
                                             Vector &elementInternalForce, double t, ModalGeomState *refState, bool forceOnly)
{
  Vector q_Big(NonLinDynamic::solVecInfo()),
  residual_Big(numOfIndices, 0.0);
  const GenVecBasis<double> &projectionBasis = dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis();
  projectionBasis.expand(*geomState.getModalq(), q_Big);
  geomState_Big->explicitUpdate(domain->getNodes(), q_Big);

  NonLinDynamic::getStiffAndForce(*geomState_Big, residual_Big, elementInternalForce, t, refState_Big, forceOnly);

  Vector r(solVecInfo(),0.);
  Vector r_lin(solVecInfo(),0.);

  udeimBasis_.compressedVecReduce(residual_Big,r);
  LinearStiffness_.fullExpand(*geomState.getModalq(),r_lin);
  residual += r;
  residual -= r_lin;

  return residual.norm();
}

void
UDEIMPodProjectionNonLinDynamic::getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                                             Corotator **allCorot, FullSquareMatrix *kelArray,
                                                             Vector &residual, double lambda, double time, GeomState *refState,
                                                             FullSquareMatrix *melArray, bool forceOnly) {
  
  if(forceOnly) {
    domain->getUDEIMInternalForceOnly(unassembledElemDOFMask,
                                      geomState, elementInternalForce,
                                      allCorot, kelArray,
                                      residual, NonLinDynamic::solVecInfo(), lambda, time, refState, melArray, kelArrayCopy);
  }
  else {
    domain->getUnassembledStiffAndForceOnly(unassembledElemDOFMask,
                                            geomState, elementInternalForce,
                                            allCorot, kelArray,
                                            residual, NonLinDynamic::solVecInfo(), lambda, time, refState, melArray, kelArrayCopy);
  }
}

void
UDEIMPodProjectionNonLinDynamic::readInterpolationBasis(){

 GenVecBasis<double> &projectionBasis = const_cast<GenVecBasis<double> &>(dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis());

 filePrint(stderr, " ... Reading Empirical Interpolation basis ...\n");
 filePrint(stderr, " ... Interplation subspace of dimension %d x %d ...\n", numOfIndices, domain->solInfo().forcePodSize);
 udeimBasis_.dimensionIs(projectionBasis.numVec(),numOfIndices);//we will read in the V^T*U*(P^T*U)*P^T in R(n x m)

 int vector = 0; int ind = 0;
  for(std::vector<double>::const_iterator it = geoSource->UDEIMVecBegin(), it_end = geoSource->UDEIMVecEnd(); it != it_end; ++it){
    udeimBasis_[vector][ind] = *it; ++ind;
    if(ind == domain->solInfo().maxDeimBasisSize){//each process get a copy of the UDEIM basis
      ind = 0;
      ++vector;
    }
  }

 getSolver()->EmpiricalSolver();
}

void
UDEIMPodProjectionNonLinDynamic::buildSparseInterpolationBasis(){

  //create mask for unassembled DEIM 
  std::map<int,std::vector<int> > columnKey;

  int row = 0;
  for (GeoSource::ElemDofPairVec::const_iterator it = geoSource->elemDofBegin(),
                                                   it_end = geoSource->elemDofEnd();
                                                   it != it_end; ++it) {

     const int elemId   = it->first;
     const int packedId = geoSource->glToPackElem(elemId);
     const int elemDOF  = it->second;

     row += 1;

     if(packedId < 0) {continue;}

     //unassembled element contibution are added to residual in the order the are iterated through
     if(columnKey.find(packedId) != columnKey.end()){
       columnKey[packedId].push_back(row-1);
     }else{
       std::vector<int> dummyVec;
       dummyVec.push_back(row-1);
       columnKey.insert(std::make_pair(packedId,dummyVec));
     }

     if(unassembledElemDOFMask[packedId].size() > 0){
       unassembledElemDOFMask[packedId].push_back(it->second);
     }else{
       std::vector<int> DOFs;
       unassembledElemDOFMask.insert(std::make_pair(packedId,DOFs));
       unassembledElemDOFMask[packedId].push_back(elemDOF);
     }
  }

  filePrint(stderr," ... Compressing Interpolation Basis ...\n");
  udeimBasis_.makeSparseBasis(columnKey);
}

void
UDEIMPodProjectionNonLinDynamic::buildReducedLinearOperator(){

  GenVecBasis<double> &projectionBasis = const_cast<GenVecBasis<double> &>(dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis());

  const int interpBasisSize = domain->solInfo().maxSizePodRom ?
                              std::min(domain->solInfo().maxSizePodRom, projectionBasis.numVec()) :
                              projectionBasis.size();

  //build reduced stiffness matrix
  LinearStiffness_.dimensionIs(projectionBasis.numVec(),projectionBasis.numVec());

  if(domain->solInfo().ReducedStiffness){
    filePrint(stderr," ... Reading Pre-computed Reduced Linear Stiffness Matrix ...\n");

    int vector = 0; int ind = 0;
    for(std::vector<double>::const_iterator it = geoSource->RedKVecBegin(), it_end = geoSource->RedKVecEnd(); it != it_end; ++it){
      LinearStiffness_[vector][ind] = *it; ++ind;
      if(ind == projectionBasis.numVec()){
        ind = 0;
        ++vector;
      }
    }

  } else {
    filePrint(stderr," ... Constructing Reduced Linear Stiffness Matrix ...\n");

    for( int column = 0; column != projectionBasis.numVectors(); ++column){
      GenVector<double> columnOfKtimesV(projectionBasis.size());
      columnOfKtimesV = 0.0;
      //K*V
      domain->getKtimesU(projectionBasis[column], (double *) 0, columnOfKtimesV, 1.0, kelArrayCopy);
      //V^T*(K*V)
      projectionBasis.reduce(columnOfKtimesV,LinearStiffness_[column]);
    }
  }
}

} /* end namespace Rom */
