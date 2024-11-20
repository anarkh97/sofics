#include "DEIMSamplingDriver.h"

#include "ElementSamplingDriver.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "NodalRestrictionMapping.h"
#include "ConnectivityUtils.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"
#include "BasisOps.h"

#include "RenumberingUtils.h"
#include "MeshDesc.h"

#include "DistrBasisFile.h"
#include <Timers.d/StaticTimers.h>

#include <Driver.d/Domain.h>
#include <Driver.d/SysState.h>
#include <Driver.d/GeoSource.h>
#include <Math.d/DBSparseMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
#include <cmath>
#include <utility>
#include <algorithm>
#include <numeric>

namespace Rom {

DEIMSamplingDriver::DEIMSamplingDriver(Domain *domain) :
  SingleDomainDynamic(domain),
  converter(NULL)
{}

void
DEIMSamplingDriver::solve() {

  SingleDomainDynamic::preProcess();
  if(domain->solInfo().newmarkBeta == 0) {
    domain->assembleNodalInertiaTensors(melArray);
  }

  converter = new VecNodeDof6Conversion(*domain->getCDSA());
  nodeDofMap = new VecNodeDof6Map(*domain->getCDSA());  

  const int podSizeMax = domain->solInfo().maxSizePodRom; 
  bool normalized = true;

  domain->createKelArray(kelArrayCopy);  

  VecBasis forceBasis;

  if(domain->solInfo().computeForceSnap){
    std::vector<double> timeStamps;
    readInBasis(podBasis_, BasisId::STATE, BasisId::POD, podSizeMax);
    writeProjForceSnap(forceBasis,timeStamps);
    writeBasisToFile(forceBasis, timeStamps, BasisId::FORCE, BasisId::SNAPSHOTS);
  } else if(domain->solInfo().orthogForceSnap) {
    int forcePodSizeMax = domain->solInfo().forcePodSize;
    readInBasis(forceBasis, BasisId::FORCE, BasisId::SNAPSHOTS,forcePodSizeMax);
    std::vector<double> SVs;
    OrthoForceSnap(forceBasis,SVs);
    writeBasisToFile(forceBasis, SVs, BasisId::FORCE, BasisId::ROB);
  } else {
    readInBasis(podBasis_, BasisId::STATE, BasisId::POD, podSizeMax, normalized);
    std::vector<int> maskIndices;
    std::vector<double> timeStamps;
    computeInterpIndices(forceBasis, maskIndices);
    computeAndWriteDEIMBasis(forceBasis, maskIndices);
    writeSampledMesh(maskIndices);
    if(domain->solInfo().statePodRomFile.size() > 0)
      computeErrorBound();
    }

}

void
DEIMSamplingDriver::computeErrorBound(){

  std::cout << "Computing Error Bound" << std::endl;
  VecBasis forceSnapshots;
  std::vector<double> timeStamps;
  readInBasis(forceSnapshots, BasisId::STATE, BasisId::SNAPSHOTS,0);

  double squareNorm     = 0.0;
  double squareNormDiff = 0.0;

  for(int column = 0; column != forceSnapshots.numVectors(); column++){
    Vector dummy1(deimBasis.numVec());
    Vector dummy2(podBasis_.numVec());
    deimBasis.reduce(forceSnapshots[column],dummy1);
    podBasis_.reduce(forceSnapshots[column],dummy2);

    squareNorm = dummy2.squareNorm();
    dummy2 = dummy2 - dummy1;
    squareNormDiff += dummy2.squareNorm();  
  }
  
  squareNormDiff = sqrt(squareNormDiff);
  squareNorm     = sqrt(squareNorm);

  std::cout << "||V^T(F-U(P^TU)^-1(P^T)F)||_F/||V^TF||_F = " << squareNormDiff/squareNorm << std::endl;

}

void
DEIMSamplingDriver::readInBasis(VecBasis &podBasis, BasisId::Type type, BasisId::Level level, int podSizeMax, bool normalized)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 if(normalized) fileName.append(".massorthonormalized");
 BasisInputStream<6> in(fileName, *converter);
 if (podSizeMax != 0) {
   std::cout << "reading in " << podSizeMax << " vectors from " << fileName.c_str() << std::endl;
   readVectors(in, podBasis, podSizeMax);
   if(type == BasisId::FORCE){
     Vector dummyV(solVecInfo());
     double dummyD = 0;
     std::pair<double,Vector> foo;
     foo = std::make_pair(dummyD,dummyV);
     in >> foo;
     MPOSingularValue = foo.first;
   }
 } else {
   std::cout << "reading in all vectors from " << fileName.c_str() << std::endl;
   readVectors(in, podBasis);
 }
}

template <typename Scalar>
void
DEIMSamplingDriver::writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type, BasisId::Level level)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 BasisOutputStream<6> outputNormalized(fileName, *converter, false);
 filePrint(stderr, " ... Writing basis to file %s ...\n", fileName.c_str());
 for (int iVec = 0; iVec < OutputBasis.vectorCount(); ++iVec) {
   if(singularValue.size() > 0)
     outputNormalized << std::make_pair(double(singularValue[iVec]), OutputBasis[iVec]);
   else
     outputNormalized << OutputBasis[iVec];
 }
}

void
DEIMSamplingDriver::writeProjForceSnap(VecBasis &forceBasis,std::vector<double> &timeStamps) 
{
  //First read in state snapshots 
  VecBasis displac;
  std::vector<int> snapshotCounts;
  readAndProjectSnapshots(BasisId::STATE, converter->vectorSize(), podBasis_, converter,
                          snapshotCounts, timeStamps, displac);

  VecBasis *veloc_;
  if(!domain->solInfo().velocPodRomFile.empty()) {
   std::vector<double> velTimeStamps;
   std::vector<int> velSnapshotCounts;
   veloc_ = new VecBasis;
   readAndProjectSnapshots(BasisId::VELOCITY, converter->vectorSize(), podBasis_, converter,
                           velSnapshotCounts, velTimeStamps, *veloc_);
   if(velSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  } else { veloc_ = NULL;}

  VecBasis *accel_;
  if(!domain->solInfo().accelPodRomFile.empty()) {
   std::vector<double> accTimeStamps;
   std::vector<int> accSnapshotCounts;
   accel_ = new VecBasis;
   readAndProjectSnapshots(BasisId::ACCELERATION, converter->vectorSize(), podBasis_, converter,
                           accSnapshotCounts, accTimeStamps, *accel_);
   if(accSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  } else { accel_ = NULL; } 

  //Now build forcevectors
  buildForceArray(forceBasis, displac, veloc_, accel_, timeStamps, snapshotCounts);
}

void
DEIMSamplingDriver::computeInterpIndices(VecBasis &forceBasis, std::vector<int> &maskIndices) {
  //member function for determining the sampling indicies for DEIM
  //result is (U*(P^T*U)^-1) where U is the left singular vectors of force snapshots and P
  //is the column selection matrix

  int forcePodSizeMax = domain->solInfo().forcePodSize;

  readInBasis(forceBasis, BasisId::FORCE, BasisId::SNAPSHOTS,forcePodSizeMax);  
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > forceMatrix(forceBasis.data(),forceBasis.vectorSize(),forceBasis.vectorCount());
  
  std::set<int> auxilaryIndices;//container for all dofs of selected nodes. Pseudo-GNAT type hyper-reduction
 
  int maxCoeffSlot;

  {
   Eigen::Matrix<double,Eigen::Dynamic,1> firstCol(forceMatrix.rows());
   firstCol = forceMatrix.col(0);
   //take absolute value of first column
   for(int row = 0; row != firstCol.rows(); row++)
       firstCol(row) = firstCol(row)*firstCol(row);
   //find maximum component
   firstCol.maxCoeff(&maxCoeffSlot);
   if(domain->solInfo().selectFullNode)
      getFullNodeIndices(firstCol,maxCoeffSlot,maskIndices,auxilaryIndices);
   else
      maskIndices.push_back(maxCoeffSlot);
  }

  //start loop to compute mask indicies 
  for(int i = 1; i < domain->solInfo().forcePodSize; ++i){ //loop starts at 1 i.e. the 2nd column
    filePrint(stderr,"\r %4.2f%% complete", double(i)/double(domain->solInfo().forcePodSize)*100.);

    //allocate space for P^T*U and P^T*u_i
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Umasked(maskIndices.size(),maskIndices.size());
    Eigen::Matrix<double,Eigen::Dynamic,1>              u_i_masked(maskIndices.size());

    for(int j = 0; j < maskIndices.size(); ++j) {//select proper rows and columns of force basis
      Umasked.row(j) = forceMatrix.block(maskIndices[j],0,1,maskIndices.size()); //(P^T*U) is square
      u_i_masked(j) = forceMatrix(maskIndices[j],maskIndices.size());//mask next column over
    }

    Eigen::FullPivLU< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > luOfUmasked(Umasked); //invert row reduced basis

    Eigen::Matrix<double,Eigen::Dynamic,1> residual(forceBasis.vectorSize());

    if(luOfUmasked.isInvertible())
      residual = forceMatrix.col(maskIndices.size()) - forceMatrix.leftCols(maskIndices.size())*luOfUmasked.inverse()*u_i_masked;
    else
      throw std::runtime_error("... Matrix Not Invertible ...");

    //take absolute value of residual vector
    for(int row = 0; row != residual.rows(); row++)
      residual(row) = residual(row)*residual(row);

    residual.maxCoeff(&maxCoeffSlot); //find indice of maximum component

    if(domain->solInfo().selectFullNode)
      getFullNodeIndices(residual,maxCoeffSlot,maskIndices,auxilaryIndices);
    else
      maskIndices.push_back(maxCoeffSlot);
     
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  if(domain->solInfo().selectFullNode){
    for(int i = 0; i != maskIndices.size(); i++)
      auxilaryIndices.erase(auxilaryIndices.erase(maskIndices[i]));//remove redundant auxilary indices from set

   std::copy(auxilaryIndices.begin(),auxilaryIndices.end(),std::back_inserter(maskIndices));//append auxilarry indices to end to selected indices list
  }

  Eigen::Map< Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> > indSol(maskIndices.data(),maskIndices.size(),1);
  std::cout << "selected indices:" << indSol.transpose() << std::endl;
#endif
}

#ifdef USE_EIGEN3
void
DEIMSamplingDriver::getFullNodeIndices(Eigen::Matrix<double,Eigen::Dynamic,1> res,int MaxCoeff,std::vector<int> &container,std::set<int> &auxilaryIndices){
  //this member function adds all indices for a selected node in order of largest residual contribution
  std::vector<int> nodalIndices;

  int selectedNode = nodeDofMap->nodeDof(MaxCoeff).nodeRank;
  
  nodeDofMap->locations(selectedNode,std::back_inserter(nodalIndices));

  std::map<double,int> myMap;
  for(int i = 0; i != nodalIndices.size(); i++){
    myMap.insert(std::make_pair(res(nodalIndices[i]),nodalIndices[i]));
  }

  container.push_back(MaxCoeff); 

  myMap.erase(myMap.find(res(MaxCoeff)));
  for(std::map<double,int>::reverse_iterator rit = myMap.rbegin(); rit != myMap.rend(); rit++)
    auxilaryIndices.insert(rit->second);
}
#endif

void
DEIMSamplingDriver::computeAndWriteDEIMBasis(VecBasis &forceBasis, std::vector<int> &maskIndices)
{
  //member function for computing and writing the DEIM basis P*(P^T*U)^-T*U^T*V
  //where V is the mass-orthogonal POD basis, U is the left Singular vectors of the force snapshots
  //and P column selection matrix derived form the sampled indicies computed above
  int maxDeimBasisSize = domain->solInfo().maxDeimBasisSize;//if this is less than the number of indices, then we solve a least squares problem
  if(maxDeimBasisSize == 0)
    maxDeimBasisSize = maskIndices.size();

  deimBasis.dimensionIs(podBasis_.vectorCount(),podBasis_.vectorInfo());
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > podMap(podBasis_.data(),podBasis_.size(), podBasis_.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > forceMap(forceBasis.data(),forceBasis.size(), forceBasis.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > deimMap(deimBasis.data(),deimBasis.size(), deimBasis.vectorCount());

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> compressedDBTranspose(podBasis_.numVectors(),maskIndices.size());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> rowReducedFM(maskIndices.size(),maxDeimBasisSize);

  //initialize (P^T*U)
  for(int row = 0; row != maskIndices.size(); row++)
    rowReducedFM.row(row) = forceMap.block(maskIndices[row],0,1,maxDeimBasisSize);

  //compute W^T = V^T*U*(P^T*U)^-1
  Eigen::JacobiSVD< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > SVDOfUmasked(rowReducedFM,Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix<double,Eigen::Dynamic,1> invSVs(SVDOfUmasked.nonzeroSingularValues());

  for(int i = 0; i != invSVs.rows(); i++) invSVs(i) = 1.0/SVDOfUmasked.singularValues()(i);

  std::cout << "condition Number of (P^T*U) = " << SVDOfUmasked.singularValues()(0)/SVDOfUmasked.singularValues()(SVDOfUmasked.nonzeroSingularValues()-1) << std::endl;
  std::cout << "||(P^T*U)^-1)||_2 = " << invSVs(SVDOfUmasked.nonzeroSingularValues()-1) << std::endl;
  std::cout << "Sigma_m+1 = " << MPOSingularValue << std::endl;
  std::cout << "E(f) = " << invSVs(SVDOfUmasked.nonzeroSingularValues()-1)*MPOSingularValue << std::endl;

  compressedDBTranspose = podMap.transpose()*forceMap.leftCols(maxDeimBasisSize)*SVDOfUmasked.matrixV()*invSVs.asDiagonal()*SVDOfUmasked.matrixU().transpose();
  //we are computing the transpose of the basis

  std::cout << "compressed Basis" << std::endl;
  std::cout << compressedDBTranspose.transpose() << std::endl;

  //initialize deim basis to all zeros
  for(int i = 0; i != deimMap.rows(); i++)
    for(int j = 0; j != deimMap.cols(); j++)
      deimMap(i,j) = 0.0;

  //expand rows W^T*P^T
  std::cout << " num of rows = " << compressedDBTranspose.rows() << " num of cols = " << compressedDBTranspose.cols() << std::endl;
  for(int row = 0; row != maskIndices.size(); row++) {
    deimMap.row(maskIndices[row]) = compressedDBTranspose.col(row);//use .col member to get rows of transposed basis
  }

  std::vector<double> dummySVs; 
  writeBasisToFile(deimBasis, dummySVs, BasisId::FORCE, BasisId::ROB);
#endif
}

void
DEIMSamplingDriver::writeSampledMesh(std::vector<int> &maskIndices) {

  //get nodes & dofs belonging to sampled indices (global numbering)
  std::set<int> selectedNodeSet; //To ensure unicity, since several locations (i.e. vector indices) can correspond to one node
  std::vector<std::pair<int,int> > compressedNodeKey; //need separate vector container for repeated nodes with different selected dofs
  for (std::vector<int>::const_iterator it = maskIndices.begin(); it != maskIndices.end(); ++it) {
    selectedNodeSet.insert(nodeDofMap->nodeDof(*it).nodeRank);
    int selectedNode = nodeDofMap->nodeDof(*it).nodeRank;
    int selectedNodeDof = std::log10(nodeDofMap->nodeDof(*it).dofId)/std::log10(2.0); //this function call is returning 2^dof for some reason, call log to normalize
    compressedNodeKey.push_back(std::make_pair(selectedNode,selectedNodeDof));
  }

//  std::sort(compressedNodeKey.begin(),compressedNodeKey.end());

  //print nodes to screen
  filePrint(stderr,"selected Nodes:");
  for(std::set<int>::iterator it = selectedNodeSet.begin(); it != selectedNodeSet.end(); it++)
    filePrint(stderr," %d", *it);

  filePrint(stderr,"\n");

  filePrint(stderr,"number of selectedNodeSet = %d, number of Indices = %d, number of dofs = %d\n", selectedNodeSet.size(),maskIndices.size(), compressedNodeKey.size());

  {
   DofSetArray *cdsa = domain->getCDSA();
   for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++)
     std::cout << "node: " << it->first+1 << " dof1 = " << cdsa->firstdof(it->first)+1 << " slot = " << it->second << std::endl;
  }

  // compute the reduced forces (constant only)
  Vector constForceFull(solVecInfo());
  getConstForce(constForceFull);
  Vector constForceRed(podBasis_.vectorCount());
  reduce(podBasis_, constForceFull,  constForceRed);

  // compute the reduced initial conditions
  Vector d0Full(SingleDomainDynamic::solVecInfo()),
         v0Full(SingleDomainDynamic::solVecInfo());
  Vector tmp(SingleDomainDynamic::solVecInfo());
  SysState<Vector> inState(d0Full, v0Full, tmp, tmp);
  SingleDomainDynamic::getInitState(inState);
  Vector d0Red(podBasis_.vectorCount()),
         v0Red(podBasis_.vectorCount());
  bool reduce_idis = (d0Full.norm() != 0),
       reduce_ivel = (v0Full.norm() != 0);
  if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) {
    AllOps<double> allOps;
    if(reduce_idis || reduce_ivel) {
      //std::cerr << "building mass matrix\n";
      allOps.M = domain->constructDBSparseMatrix<double>();
      domain->makeSparseOps(allOps, 0, 0, 0);
    }
    if(reduce_idis) {
      allOps.M->mult(d0Full, tmp);
      reduce(podBasis_, tmp, d0Red);
    }
    if(reduce_ivel) {
      allOps.M->mult(v0Full, tmp);
      reduce(podBasis_, tmp, v0Red);
    }
  }
  else {
    if(reduce_idis) reduce(podBasis_, d0Full, d0Red);
    if(reduce_ivel) reduce(podBasis_, v0Full, v0Red);
  }

  //fill element container with all elements 
  std::vector<int> packedToInput(elementCount());
  Elemset &inputElemSet = *(geoSource->getElemSet());
  {
   for (int iElem = 0, iElemEnd = inputElemSet.size(); iElem != iElemEnd; ++iElem) {
     Element *elem = inputElemSet[iElem];
     if (elem) {
       const int iPackElem = domain->glToPackElem(iElem);
       if(iPackElem >= 0) {
         assert(iPackElem < packedToInput.size());
         packedToInput[iPackElem] = iElem;
       }
     }
   }
  }
 
  // Determine mapping between elements and nodes
  std::unique_ptr<Connectivity> elemToNode(new Connectivity(inputElemSet.asSet()));
  std::unique_ptr<Connectivity> nodeToElem(elemToNode->alloc_reverse());
 

   //get elements belonging to sampledNodes
   std::vector<int> sampleElemRanks;
   std::vector<int> sampleNodeIds(selectedNodeSet.begin(),selectedNodeSet.end()); //fill vector container with sampled nodes
   std::sort(sampleNodeIds.begin(),sampleNodeIds.end());//sort nodes in ascending order
   connections(*nodeToElem, sampleNodeIds.begin(), sampleNodeIds.end(), std::back_inserter(sampleElemRanks));//get all adjecent elements

   //fill sampled elements container
   std::vector<int> sampleElemIds;
   sampleElemIds.reserve(sampleElemRanks.size());
   std::map<int, double> weights;
   for (std::vector<int>::const_iterator it = sampleElemRanks.begin(), it_end = sampleElemRanks.end(); it != it_end; ++it) {
     const int elemRank = *it;
     weights.insert(std::make_pair(elemRank, 1.0));
     sampleElemIds.push_back(elemRank);
   }
  std::cout << std::endl;
  //construct element weight solution vector
  std::vector<double> solution(elementCount());
  for(int i = 0; i != elementCount(); i++)
    if(weights[packedToInput[i]])
     solution[i] = weights[packedToInput[i]];
    else 
     solution[i] = 0.0;
 
  //output Full Weights for compatibility with full mesh hyperreduction
  {
   // 1. write the ATTRIBUTE data
   outputFullWeights(solution, packedToInput);

   // 2. append the SNSLOT data
   const std::string fileName = domain->solInfo().reducedMeshFile;
   const std::ios_base::openmode mode = std::ios_base::out | std::ios_base::app;
   std::ofstream weightOut(fileName.c_str(), mode);
   weightOut.precision(std::numeric_limits<double>::digits10+1);

   weightOut << "*\n";
   weightOut << "SNSLOT\n";
   for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++)
     weightOut << it->first + 1 << " " << it->second << std::endl; 
  }
 
  //construct and print renumbered mesh
  const FileNameInfo fileInfo;
  const MeshRenumbering meshRenumbering(sampleElemIds.begin(), sampleElemIds.end(), *elemToNode, true);
  const MeshDesc reducedMesh(domain, geoSource, meshRenumbering, weights);
  outputMeshFile(fileInfo, reducedMesh, podBasis_.vectorCount());

  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
  meshOut.precision(std::numeric_limits<double>::digits10+1);
  if(domain->solInfo().reduceFollower) meshOut << "EXTFOL\n";
 
  // reduced stiffness
  meshOut << "*\nREDSTIFF\n"  ;
  
  for(int column=0; column<podBasis_.vectorCount(); ++column){
    Vector columnOfKtimesV(solVecInfo());
    Vector columnOfRedK(podBasis_.vectorCount());
    columnOfKtimesV.zero();
    domain->getKtimesU(podBasis_[column], bcx, columnOfKtimesV, 1.0, kelArrayCopy);
    reduce(podBasis_, columnOfKtimesV, columnOfRedK);
    for(int i=0; i<podBasis_.vectorCount(); ++i){
      meshOut << columnOfRedK[i] << std::endl;}
  }

  // output the reduced forces
  meshOut << "*\nFORCES\nMODAL\n";
  for(int i=0; i<podBasis_.vectorCount(); ++i)
    meshOut << i+1 << " "  << constForceRed[i] << std::endl; 

  // output the reduced initial conditions
  if(reduce_idis) {
    meshOut << "*\nIDISPLACEMENTS\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << d0Red[i] << std::endl;
  }
  if(reduce_ivel) {
    meshOut << "*\nIVELOCITIES\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << v0Red[i] << std::endl;
  }


  #ifdef USE_EIGEN3
  // build and output compressed basis and mask operator P
  podBasis_.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain->getCDSA());
  deimBasis.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain->getCDSA());
  {
    std::string PODfilename = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    std::string DEIMfilename = BasisFileId(fileInfo, BasisId::FORCE, BasisId::ROB);
    PODfilename.append(".reduced");
    DEIMfilename.erase(DEIMfilename.end()-5,DEIMfilename.end());
    DEIMfilename.append(".reduced.deim");
    if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) {
      PODfilename.append(".massorthonormalized");
    }
    filePrint(stderr," ... Writing compressed POD/DEIM basis to file %s and %s ...\n", PODfilename.c_str(), DEIMfilename.c_str());
    DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
    ConstrainedDSA reduced_cdsa(reduced_dsa, reducedMesh.dirichletBConds().size(), const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]));
    VecNodeDof6Conversion converter(reduced_cdsa);
    BasisOutputStream<6> outputPOD(PODfilename, converter, false);
    BasisOutputStream<6> outputDEIM(DEIMfilename, converter, false);

    for (int iVec = 0; iVec < podBasis_.vectorCount(); ++iVec) {
      outputPOD  << podBasis_.compressedBasis().col(iVec);
      outputDEIM << deimBasis.compressedBasis().col(iVec);      
    }

    std::map<int,int> nodeRenum(meshRenumbering.nodeRenumbering());
    std::map<int,int> elemRenum(meshRenumbering.elemRenumbering());

    meshOut << "*\nSNSLOT\n";
    for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++){
     //output assembled indices in form of a node plus a node dof and an element with and element dof (to keep neighbor elements for adding to the same dof)
     meshOut << nodeRenum[it->first] + 1 << " " << it->second << std::endl;
    }

  }
#endif 

}

void
DEIMSamplingDriver::buildForceArray(VecBasis &forceBasis, const VecBasis &displac, const VecBasis *veloc,
                                     const VecBasis *accel, std::vector<double> timeStamps_, std::vector<int> snapshotCounts_)
{//this memeber function is for converting state snapshots to force snapshots in the absence of precollected force snapshots from model I
  //most of the code is copied from assembleTrainingData in ElementSamplingDriver.C
  std::vector<double>::iterator timeStampIt = timeStamps_.begin();

  forceBasis.dimensionIs(displac.vectorCount(), displac.vectorInfo());

  int iSnap = 0;
  double gamma  = domain->solInfo().newmarkGamma;
  double alphaf = domain->solInfo().newmarkAlphaF;

  double TotalSum = 0;
  double NonlinearSum = 0;
  double LinearSum = 0;
  int counter = 0;

  for(int i = 0; i < snapshotCounts_.size(); i++) {
    GenVector<double> FNLint(solVecInfo());
    GenVector<double> FLint(solVecInfo());
    for(int jSnap = 0; jSnap != snapshotCounts_[i]; ++iSnap, ++jSnap){
      filePrint(stderr,"\r %4.2f%% complete, time = %f", double(iSnap)/double(std::accumulate(snapshotCounts_.begin(),snapshotCounts_.end(),0))*100.,*timeStampIt);

      geomState->explicitUpdate(domain->getNodes(), displac[iSnap]);
      if(veloc){ geomState->setVelocity((*veloc)[iSnap]);} //just set the velocity at the nodes
      if(accel){ geomState->setAcceleration((*accel)[iSnap]);} //just set the acceleration at the nodes

      //set up for internal force vector
      Vector dummy(solVecInfo()); 
      dummy = displac[iSnap];
      FLint.zero();
      SingleDomainDynamic::getInternalForce( dummy, FNLint, *timeStampIt, jSnap);
      domain->getKtimesU(dummy, bcx, FLint, 1.0, kelArrayCopy);
   
      //set vector in force snapshot container 
      //remove linear component from force snapshots
      forceBasis[iSnap] = (FNLint-FLint); 
      TotalSum     += FNLint.norm();
      NonlinearSum += forceBasis[iSnap].norm();
      LinearSum    += FLint.norm();
     
      counter += 1;
      timeStampIt++;
    }
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  std::cout << "Avg.     Total Norm = " << TotalSum/double(counter) << std::endl;
  std::cout << "Avg. Nonlinear Norm = " << NonlinearSum/double(counter) << std::endl;
  std::cout << "Avg.    Linear Norm = " << LinearSum/double(counter) << std::endl;

}

void DEIMSamplingDriver::OrthoForceSnap(VecBasis &forceBasis,std::vector<double> &SVs)
{
#ifdef USE_EIGEN3
  std::cout << "... Orthogonalizting Snapshots ..." << std::endl;
  SVs.resize(forceBasis.numVectors());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > SingularValueMap(SVs.data(),forceBasis.numVectors());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ForceMap(forceBasis.data(), forceBasis.size(), forceBasis.numVectors());
  Eigen::JacobiSVD<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ForceSVD(ForceMap, Eigen::ComputeThinU);
  ForceMap = ForceSVD.matrixU();
  SingularValueMap = ForceSVD.singularValues();
#endif
}

int DEIMSamplingDriver::elementCount() const {
  return domain->numElements();
}

//member functions for fascilitating the reading of snapshot data
int 
DEIMSamplingDriver::snapSize(BasisId::Type type, std::vector<int> &snapshotCounts)
{
  // compute the number of snapshots that will be used for a given skipping strategy
  FileNameInfo fileInfo;
  snapshotCounts.clear();
  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    DistrBasisInputFile in(fileName);
    const double N = in.stateCount() - skipOffSet;
    const int singleSnapshotCount = (N > 0) ? 1+(N-1)/skipFactor : 0;
    snapshotCounts.push_back(singleSnapshotCount);
  }
  // return the total
  return std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0);
}

void 
DEIMSamplingDriver::readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                                              const VecNodeDof6Conversion *vecDofConversion,
                                              std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config)
{
  const int snapshotCount = snapSize(type, snapshotCounts);
  filePrint(stderr, " ... Reading in %d %s Snapshots ...\n", snapshotCount, toString(type).c_str());

  config.dimensionIs(snapshotCount, vectorSize);
  timeStamps.clear();
  timeStamps.reserve(snapshotCount);

  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  const int podVectorCount = podBasis.vectorCount();
  Vector snapshot(vectorSize);
  Vector podComponents(podVectorCount);
  const FileNameInfo fileInfo;

  int offset = 0;
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
    BasisInputStream<6> in(fileName, *vecDofConversion);

    int count = 0;
    int skipCounter = skipFactor - skipOffSet;
    while(count < snapshotCounts[i]) {
      std::pair<double, double *> data;
      data.second = snapshot.data();
      in >> data;
      assert(in);
      if(skipCounter == skipFactor) {
    //    expand(podBasis, reduce(podBasis, snapshot, podComponents), config[offset+count]);
        config[offset+count] = snapshot;
        timeStamps.push_back(data.first);
        skipCounter = 1;
        ++count;
        filePrint(stderr, "\r ... timeStamp = %7.2e, %4.2f%% complete ...", data.first, double(count)/snapshotCounts[i]*100);
      }
      else {
        ++skipCounter;
      }
    }

    filePrint(stderr,"\n");
    offset += snapshotCounts[i];
  }

  assert(timeStamps.size() == snapshotCount);
}

} /* end namespace Rom */

Rom::DriverInterface *deimSamplingDriverNew(Domain *domain) {
  return new Rom::DEIMSamplingDriver(domain);
}
