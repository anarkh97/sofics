#include "ConstraintSamplingDriver.h"

#ifdef USE_EIGEN3
#include "ElementSamplingDriver.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "NodalRestrictionMapping.h"
#include "ConnectivityUtils.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"
#include "BasisOps.h"

#include "NonnegativeMatrixFactorization.h"
#include <Element.d/MpcElement.d/MpcElement.h>

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

#include <Eigen/Dense>
#include <cmath>
#include <utility>
#include <algorithm>
#include <numeric>

namespace Rom {

ConstraintSamplingDriver::ConstraintSamplingDriver(Domain *domain) :
  SingleDomainDynamic(domain),
  converter6(NULL),
  converter1(NULL)
{}

void
ConstraintSamplingDriver::solve() {

  // preprocess data structures
  SingleDomainDynamic::preProcess();
  if(domain->solInfo().newmarkBeta == 0) {
    domain->assembleNodalInertiaTensors(melArray);
  }

  // read in dual and primal snapshots
  
  //writeSampledMesh(expandedSolution); // only output part of reduced mesh that contains contact specification

}

void
ConstraintSamplingDriver::readInBasis(VecBasis &basis, BasisId::Type type, BasisId::Level level, bool vectorQuant, int i, int podSizeMax, bool normalized)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level, i);
 if(normalized) fileName.append(".massorthonormalized");

 if(vectorQuant){ // read in vector valued data
   BasisInputStream<6> in(fileName, *converter6);
   if (podSizeMax != 0) {
     std::cout << " ... reading in " << podSizeMax << " vectors from " << fileName.c_str() << " ... " << std::endl;
     readVectors(in, basis, podSizeMax);
   } else {
     std::cout << " ... reading in all vectors from " << fileName.c_str() << " ... " << std::endl;
     readVectors(in, basis);
   }
 } else {       // read in scalar valued data
     BasisInputStream<1> in(fileName, *converter1);
     if (podSizeMax != 0) {
       std::cout << " ... reading in " << podSizeMax << " vectors from " << fileName.c_str() << " ... " << std::endl;
       readVectors(in, basis, podSizeMax);
     } else {
       std::cout << " ... reading in all vectors from " << fileName.c_str() << " ... " << std::endl;
       readVectors(in, basis);
     }
 }

}

void
ConstraintSamplingDriver::setSolverFlags() {

  solver_.relativeToleranceIs(domain->solInfo().tolPodRom);
  solver_.verboseFlagIs(true);
  solver_.scalingFlagIs(domain->solInfo().useScalingSpnnls);
  solver_.centerFlagIs(domain->solInfo().useCenterSpnnls);
  solver_.reverseFlagIs(domain->solInfo().useReverseOrder);
  solver_.projectFlagIs(domain->solInfo().projectSolution);
  solver_.positivityIs(domain->solInfo().positiveElements);
  solver_.solverTypeIs(domain->solInfo().solverTypeSpnnls);
  solver_.maxSizeRatioIs(domain->solInfo().maxSizeSpnnls);
  solver_.maxIterRatioIs(domain->solInfo().maxIterSpnnls);
  solver_.maxNumElemsIs(domain->solInfo().maxElemSpnnls);

}

template <typename Scalar>
void
ConstraintSamplingDriver::writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type, BasisId::Level level, bool vectorQuant)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 if(vectorQuant){  // write vector valued data
   BasisOutputStream<6> outputNormalized(fileName, *converter6, false);
   filePrint(stderr, " ... Writing basis to file %s ...\n", fileName.c_str());
   if(singularValue.size() > 0){
     for (int iVec = 0; iVec < singularValue.size(); ++iVec) {
       outputNormalized << std::make_pair(double(singularValue[iVec]), OutputBasis[iVec]);
     }
   } else {
     for (int iVec = 0; iVec < OutputBasis.numVec(); ++iVec) {
       outputNormalized << OutputBasis[iVec];
     }
   }
 } else {         // write scalar valued data
   BasisOutputStream<1> outputNormalized(fileName, *converter1, false);
   filePrint(stderr, " ... Writing basis to file %s ...\n", fileName.c_str());
   if(singularValue.size() > 0){
     for (int iVec = 0; iVec < singularValue.size(); ++iVec) {
       outputNormalized << std::make_pair(double(singularValue[iVec]), OutputBasis[iVec]);
     }
   } else {
     for (int iVec = 0; iVec < OutputBasis.numVec(); ++iVec) {
       outputNormalized << OutputBasis[iVec];
     }
   }
 }
}

void
ConstraintSamplingDriver::expandSolution(std::vector<double> &solution, std::vector<double> &expandedSolution){

  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dualMap(dualBasis.data(),dualBasis.vectorSize(),dualBasis.numVectors());

  int counter = 0;
  for(int row = 0; row < dualMap.rows(); row++){
    if(dualMap.row(row).norm() > 1e-16){
      expandedSolution[row] = solution[counter]; // could be zero
      counter++;
    } else {
      expandedSolution[row] = 0.0;
    }
  }

}

void
ConstraintSamplingDriver::writeWeightedDualBasis(std::vector<double> &solution, std::vector<double> &expandedSolution){

  // solution vector should have some sparsity, remove zero elements for constructing compressed dual basis
  std::vector<double> compressedSolution; 
  int nnz = 0;
  for (int i = 0; i < solution.size(); i++){
    if(solution[i] >1e-16){ // keep non-zero elements
      compressedSolution.push_back(solution[i]);
      nnz++;
    }
  }

  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > solutionMap(compressedSolution.data(),compressedSolution.size());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > expSolutionMap(expandedSolution.data(),expandedSolution.size());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dualMap(dualBasis.data(),dualBasis.vectorSize(),dualBasis.numVectors());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> compressedDB(nnz,dualBasis.numVectors());
  
  int counter = 0;
  for (int row = 0; row < dualMap.rows(); row++){
    if(dualMap.row(row).norm() > 1e-16 && expandedSolution[row] > 1e-16){ // only keep non-zeros rows of dual basis that have a non-zero weight
      compressedDB.row(counter) = dualMap.row(row);
      counter++;
    }
  }
 
  // output compressed dual basis weighted by compressed solution vector
  {
    FileNameInfo fileInfo;
    std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD); // filename is set under OUTPUT command robdataf
    fileName += ".ecsw.compressed";
    VecNodeDof1Conversion converterReduced(nnz);
    BasisOutputStream<1> output(fileName, converterReduced, false);
    for(int col = 0; col != compressedDB.cols(); ++col){
      Eigen::Matrix<double,Eigen::Dynamic,1> buffer(compressedDB.rows());
      buffer = compressedDB.col(col); // eigen is stored column major, use buffer to stream data correctly
      buffer.array() *= solutionMap.array(); // weight the current column and then output it
      output << std::make_pair(double(col)+1.0,buffer.data());
    }
  }

  // output full dual basis weighted by expanded solution vector
  {
    FileNameInfo fileInfo;
    std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD); // filename is set under OUTPUT command robdataf
    fileName += ".ecsw";
    VecNodeDof1Conversion converterReduced(dualMap.rows());
    BasisOutputStream<1> output(fileName, converterReduced, false);
    for(int col = 0; col != dualMap.cols(); ++col){
      Eigen::Matrix<double,Eigen::Dynamic,1> buffer(dualMap.rows());
      buffer = dualMap.col(col); // eigen is stored column major, use buffer to stream data correctly
      buffer.array() *= expSolutionMap.array(); // weight the current column and then output it 
      output << std::make_pair(double(col)+1.0,buffer.data());
    }
  }

}

void
ConstraintSamplingDriver::writeSampledMesh(std::vector<double> &solution) {
// keep only constraint functions that correspond to selected indices

  FileNameInfo fileInfo;
  const std::ios_base::openmode mode = std::ios_base::out;
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), mode); //file name specified under OUTPUT samplmsh
  meshOut << "TOPOLOGY" << std::endl;
  meshOut << "*constraints needed for DEIM" << std::endl; 

  Elemset &constraintEset = geoSource->getPackedEsetConstraintElementIeq();

  // need indices to be ordered so that multiplication is done correctly
  std::set<int> indiceSet;
  int ind = 0;
  for(std::vector<double>::iterator it = solution.begin(); it !=  solution.end(); ++it, ind++){
    if(*it > 1e-16) {
      indiceSet.insert(ind);
    }
  }

  int counter = 1; 
  for(std::set<int>::iterator it = indiceSet.begin(); it !=  indiceSet.end(); ++it ){ // loop over the selected indices
    int numLMPC = 0;
    for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) { //loop over constraint functions to find which element the current indice belongs to
      Element *ele = constraintEset[i];
      int n = ele->getNumMPCs(); // number of constraints specified by this element
      int nn = ele->numNodes();  // number of nodes in this constraint element
      int *nvec = ele->allNodes();  // vector of node number -1 
      int tnum = ele->getElementType(); // element type
      int glnum = ele->getGlNum();
      for(int j = 0; j < n; ++j) {
        if(*it == numLMPC){
          meshOut << domain->numElements() + counter << " "; // output new element number
          meshOut << tnum << " ";                            // then output element type
          for( int w = 0;  w < nn - 1; ++w){
             meshOut << nvec[w] + 1 << " ";                  // lastly the element nodes
          }
          meshOut << std::endl;
          counter++;
        }     
        numLMPC++; 
      }
    }
  }
  meshOut << "ATTRIBUTES" << std::endl;
  meshOut << domain->numElements() + 1 << " " << domain->numElements() + counter - 1 << " 2" << std::endl;

}

void
ConstraintSamplingDriver::buildConstraintArray(const VecBasis &dualBasis, const VecBasis &displac)
{//this memeber function is for converting state snapshots to constraint function snapshots in the absence of precollected snapshots from model I

  // set up constraint functions
  if(geoSource->getNumConstraintElementsIeq() > 0) {

// if lagrange multiplier snapshots are included, use them to filter out stateshots
    int nsnap = 0;
    VecBasis *lagrangeMP = NULL;
    if(domain->solInfo().statePodRomFile.size() > 1){
      lagrangeMP = new VecBasis;
      std::cout << " ... Selective Snapshot Procedure ...\n";
      readInBasis(*lagrangeMP, BasisId::STATE, BasisId::SNAPSHOTS, false, 1);
      for(int col = 0; col < lagrangeMP->numVectors(); col++){ // count the number of non-zero columns to set up solver problem size
        if((*lagrangeMP)[col].norm() > 1e-16) {
          nsnap++;
        }
      }
      std::cout << " ... there are " << nsnap << " non-zero vectors ... " << std::endl;
    } else {
      nsnap = displac.numVectors();
    }

    if(lagrangeMP){
      if(lagrangeMP->numVectors() != displac.numVectors()){
        std::cout << " ... Lagrange snapshots do not correspond to state snapshots ... " << std::endl;
        std::cout << "lagrangeMP->numVectors() = " << lagrangeMP->numVectors() << std::endl;
        std::cout << "displac.numVectors() = " << displac.numVectors() << std::endl;
      }
    }

    // set up useful maps
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dualMap(dualBasis.data(),dualBasis.vectorSize(),dualBasis.numVectors());

    // allocate space for constraint snapshot
    Eigen::Matrix<double,Eigen::Dynamic,1> targetSnap(dualBasis.vectorSize());
    Eigen::Matrix<double,Eigen::Dynamic,1> reduceTarg(dualBasis.numVectors());    

    // set up target and element contribution containers
    int nnzrows = NNZRows(dualBasis);

    solver_.problemSizeIs(dualBasis.numVectors()*nsnap,nnzrows);
    double *trainingTarget = solver_.rhsBuffer();

    std::cout << " ... Initializing Solver ..." << std::endl;
    int iSnap = 0;
    for(int jSnap = 0; jSnap != displac.numVectors(); ++jSnap){ // loop over snapshots
      bool notZeroColumn = true;
      if(lagrangeMP) { // if we have read in the lagrange multipliers, see if the current column has any active constraints
        if((*lagrangeMP)[jSnap].norm() < 1e-16){
          notZeroColumn = false;
        }
      }

      if(notZeroColumn){
        filePrint(stderr,"\r %4.2f%% complete", double(jSnap)/double(displac.numVectors())*100.);
        geomState->explicitUpdate(domain->getNodes(), displac[jSnap]);

        Elemset &eset = geoSource->getPackedEsetConstraintElementIeq();
        ResizeArray<LMPCons *> lmpc(0);
        int numLMPC = 0;
        for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) { //loop over constraint functions
          Element *ele = eset[i];
          static_cast<MpcElement*>(ele)->update(geomState, *geomState, domain->getNodes(), double(jSnap));
          int n = ele->getNumMPCs();
          LMPCons **l = ele->getMPCs();
          for(int j = 0; j < n; ++j) {
            lmpc[numLMPC++] = l[j];
          }
          delete [] l;
        }
        targetSnap.setZero();
        int eleCounter = 0;
        for(int i=0; i<numLMPC; ++i) { // get right hand side which is Constraint evaluated at current state
          if(dualMap.row(i).norm() > 1e-16) {
            targetSnap[i] = -1.0*lmpc[i]->rhs.r_value;
            reduceTarg.setZero();
            reduceTarg = dualMap.row(i)*targetSnap[i];
            for(int q = 0; q < dualMap.cols(); q++){
              solver_.matrixEntry(q+iSnap*dualMap.cols(),eleCounter) = reduceTarg[q]; // set elem contribution for current snapshot
            }
            eleCounter++;
          }
        }
        reduceTarg.setZero();
        reduceTarg = dualMap.transpose()*targetSnap;
        for( int vec = 0; vec < dualMap.cols(); vec++) { // dump full snapshot
           *trainingTarget = reduceTarg[vec];
            trainingTarget++;
        }
        iSnap++;
      }
    }
    filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  } else {

    std::cout << "No constraints specified\n" << std::endl; 

  }

}

int
ConstraintSamplingDriver::NNZRows(const VecBasis &basis) {

  // zero out rows of constraint snapshots that correspond to zeros of the dual basis
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > basisMap(basis.data(), basis.size(), basis.numVectors());

  int nnz = 0;
  std::cout << " ... Filtering Constraint Snapshots ... " << std::endl;
  for (int row = 0; row < basisMap.rows(); ++row)
    if(basisMap.row(row).norm() > 1e-16) // row is all zeros
      nnz++;

  return nnz; 
}


//member functions for fascilitating the reading of snapshot data
int 
ConstraintSamplingDriver::snapSize(BasisId::Type type, std::vector<int> &snapshotCounts)
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
ConstraintSamplingDriver::readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
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
        expand(podBasis, reduce(podBasis, snapshot, podComponents), config[offset+count]);
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
#endif

Rom::DriverInterface *constraintSamplingDriverNew(Domain *domain) {
#ifdef USE_EIGEN3
  return new Rom::ConstraintSamplingDriver(domain);
#else
  std::cerr << " *** ERROR: ConstraintSamplingDriver requires AERO-S configured with the Eigen library. Exiting...\n";
  exit(-1);
#endif
}
