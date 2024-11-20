#include "DEIMPodProjectionNonLinDynamic.h"

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

DEIMPodProjectionNonLinDynamic::DEIMPodProjectionNonLinDynamic(Domain *domain) :
  LumpedPodProjectionNonLinDynamic(domain)
{}

void
DEIMPodProjectionNonLinDynamic::preProcess() {

  LumpedPodProjectionNonLinDynamic::preProcess();

  kelArrayCopy = new FullSquareMatrix;
  domain->createKelArray(kelArrayCopy);

  readInterpolationBasis(); 
  buildSparseInterpolationBasis();
  buildReducedLinearOperator();

}

void
DEIMPodProjectionNonLinDynamic::reBuild(ModalGeomState &geomState, int iteration, double localDelta, double t)
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
DEIMPodProjectionNonLinDynamic::computeRedKDyn(double Kcoef){

 GenVecBasis<double> &projectionBasis = const_cast<GenVecBasis<double> &>(dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis());
 GenVecBasis<double> RedKDyn(projectionBasis.numVec(),projectionBasis.numVec());
 for(int column = 0; column != projectionBasis.numVec(); column++) {
   GenVector<double> columnOfKtimesV(projectionBasis.size());
   columnOfKtimesV = 0.0;
   //(K_dyn - K_stiff)*V
   domain->getWeightedKtimesU(packedElementWeights_[0], projectionBasis[column], (double *) 0, columnOfKtimesV, Kcoef, kelArray);
   deimBasis_.compressedVecReduce(columnOfKtimesV,RedKDyn[column]);
 }
#ifdef USE_EIGEN3
 Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > RKDMap(RedKDyn.data(),RedKDyn.size(),RedKDyn.numVec());
 getSolver()->addToReducedMatrix(RKDMap);
#endif
}

double
DEIMPodProjectionNonLinDynamic::getStiffAndForce(ModalGeomState &geomState, Vector &residual,
                                             Vector &elementInternalForce, double t, ModalGeomState *refState, bool forceOnly)
{
  Vector q_Big(NonLinDynamic::solVecInfo()),
  residual_Big(NonLinDynamic::solVecInfo(), 0.0);
  const GenVecBasis<double> &projectionBasis = dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis();
  projectionBasis.expand(*geomState.getModalq(), q_Big);
  geomState_Big->explicitUpdate(domain->getNodes(), q_Big);

  NonLinDynamic::getStiffAndForce(*geomState_Big, residual_Big, elementInternalForce, t, refState_Big, forceOnly);

  Vector r(solVecInfo(),0.);
  Vector r_lin(solVecInfo(),0.);

  deimBasis_.compressedVecReduce(residual_Big,r);
  LinearStiffness_.fullExpand(*geomState.getModalq(),r_lin);
  residual += r;
  residual -= r_lin;

  return residual.norm();
}

void
DEIMPodProjectionNonLinDynamic::getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                                             Corotator **allCorot, FullSquareMatrix *kelArray,
                                                             Vector &residual, double lambda, double time, GeomState *refState,
                                                             FullSquareMatrix *melArray, bool forceOnly) {
  
  if(forceOnly) {
    domain->getWeightedInternalForceOnly(packedElementWeights_[0],
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray, kelArrayCopy);
  }
  else {
    domain->getWeightedStiffAndForceOnly(packedElementWeights_[0],
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray, kelArrayCopy);
  }
}

void
DEIMPodProjectionNonLinDynamic::readInterpolationBasis(){
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, BasisId::FORCE, BasisId::POD);
 fileName = fileName + ".deim";

 VecNodeDof6Conversion vecNodeDof6Conversion_(*domain->getCDSA());

 BasisInputStream<6> deimBasisInput(fileName, vecNodeDof6Conversion_);

 filePrint(stderr, " ... Reading Empirical Interpolation basis from file %s ...\n", fileName.c_str());
 const int interpBasisSize = domain->solInfo().maxSizePodRom ?
                             std::min(domain->solInfo().maxSizePodRom, deimBasisInput.size()) :
                             deimBasisInput.size();

 filePrint(stderr, " ... Interplation subspace of dimension = %d ...\n", interpBasisSize);

 readVectors(deimBasisInput, deimBasis_, interpBasisSize);

 getSolver()->EmpiricalSolver();
}

void
DEIMPodProjectionNonLinDynamic::buildSparseInterpolationBasis(){

 std::vector<std::pair<int,int> > maskedIndicesBuf;
 
 {
   for (GeoSource::NodeDofPairVec::const_iterator it = geoSource->nodeDofSlotBegin(),
                                                  it_end = geoSource->nodeDofSlotEnd();
                                                  it != it_end; ++it) {

     const int nodeId = it->first;
//     const int packedId = domain->globalToLocal(nodeId);
     const int packedId = nodeId;

     if(packedId < 0) {continue;}

     maskedIndicesBuf.push_back(std::make_pair(packedId,it->second));

   }
 }

  filePrint(stderr," ... Compressing Interpolation Basis ...\n");
  deimBasis_.makeSparseBasis(maskedIndicesBuf, domain->getCDSA());
}

void
DEIMPodProjectionNonLinDynamic::buildReducedLinearOperator(){

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
