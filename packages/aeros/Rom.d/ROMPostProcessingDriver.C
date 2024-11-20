#include "ROMPostProcessingDriver.h"

#include <Driver.d/Domain.h>
#include <Problems.d/DynamDescr.h>

#include "VecBasis.h"
#include "VecBasisOps.h"
#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "VecBasisFile.h"

#include "MasterMapping.h"
#include "VecNodeDof6Conversion.h"
#include "NodeDof6Buffer.h"

#include <Math.d/Vector.h>

#include "PtrPtrIterAdapter.h"
#include "PodProjectionSolver.h"

#include <algorithm>
#include <memory>
#include <fstream>
#include <iostream>
#include <cassert>
#include <numeric>

class Rbm;
 
namespace Rom {

ROMPostProcessingDriver::ROMPostProcessingDriver(Domain *domain_) :
NonLinDynamic(domain_),
normalizedBasis_(),
curState(NULL), fullDispBuffer(NULL), fullVelBuffer(NULL), fullAccBuffer(NULL),
fullVel2Buffer(NULL), fullDummyBuffer(NULL)
{}

ROMPostProcessingDriver::~ROMPostProcessingDriver()
{
 if(curState) delete curState;
 if(fullDispBuffer) delete fullDispBuffer;
 if(fullVelBuffer) delete fullVelBuffer;
 if(fullAccBuffer) delete fullAccBuffer;
 if(fullVel2Buffer) delete fullVel2Buffer;
 if(fullDummyBuffer) delete fullDummyBuffer;
}

void
ROMPostProcessingDriver::preProcess()
{
  NonLinDynamic::preProcess();
  bufferReducedFiles();

  VecNodeDof6Conversion converter(*domain->getCDSA());
  FileNameInfo fileInfo;

  normalizedBasis_.dimensionIs(projectionSubspaceSize, solVecInfo());

  for(int j=0; j<domain->solInfo().readInROBorModes.size(); ++j)
  {
    // Load projection basis
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD, j);
    if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) {
      if(j==0) filePrint(stderr, " ... Using Mass-normalized Basis    ...\n");
      fileName.append(".normalized");
    }
    BasisInputStream<6> projectionBasisInput(fileName, converter);
    if(verboseFlag) filePrint(stderr, " ... Reading basis from file %s ...\n", fileName.c_str());

    const int projectionSubspaceSize = domain->solInfo().localBasisSize[j] ?
                                       std::min(domain->solInfo().localBasisSize[j], projectionBasisInput.size()) :
                                       projectionBasisInput.size();

    readVectors(projectionBasisInput, normalizedBasis_, 
                std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0),
                projectionSubspaceSize,
                std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.begin()+j, 0));
  }

  filePrint(stderr, " ... Proj. Subspace Dimension = %-3d ...\n", projectionSubspaceSize);
  if(domain->solInfo().readInROBorModes.size() > 1)
    filePrint(stderr, " ... Number of Local Bases = %-3d    ...\n", domain->solInfo().readInROBorModes.size());

  if(domain->solInfo().sensitivity) {
 
    if(!domain->solInfo().readInAdjointROB.empty()) {
      int adjointSubspaceMaxSize = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(),
                                                   domain->solInfo().maxSizeAdjointBasis.end(), 0); 
      adjointBasis_.dimensionIs(adjointSubspaceMaxSize, solVecInfo());

      for(int j=0; j<domain->solInfo().readInAdjointROB.size(); ++j) {
        std::string &fileName = domain->solInfo().readInAdjointROB[j];
        if(verboseFlag) filePrint(stderr, " ... Reading adjoint basis from file %s ...\n", fileName.c_str());
        BasisInputStream<6> adjointBasisInput(fileName, converter);

        const int &adjointSubspaceSize = domain->solInfo().maxSizeAdjointBasis[j];
        if(adjointSubspaceSize > adjointBasisInput.size()) {
          std::cerr << 
            "ERROR: selected size of adjoint basis is larger than the available number of basis vectors in file " 
            << fileName << std::endl;
          exit(-1);
        }

        readVectors(adjointBasisInput, adjointBasis_,
                    std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(), domain->solInfo().maxSizeAdjointBasis.end(), 0),
                    adjointSubspaceSize,
                    std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(), domain->solInfo().maxSizeAdjointBasis.begin()+j, 0));
        filePrint(stderr, " ... Adjoint Subspace #%d Dim. = %-3d ...\n", j+1, adjointSubspaceSize);
      }

      PodProjectionSolver *ppsolver = dynamic_cast<PodProjectionSolver*>(solver);
      if(ppsolver) {
        ppsolver->projectionBasisIs(adjointBasis_);
        ppsolver->fullSolutionIs(true);
        ppsolver->setLocalBasis(0,0);
      }
    }

    preProcessSA();
  }
  
  VECsize = normalizedBasis_.size();

  //initialize containers for full coordinates 
  fullDispBuffer  = new GenVector<double>(NonLinDynamic::solVecInfo()); fullDispBuffer->zero();
  fullVelBuffer   = new GenVector<double>(NonLinDynamic::solVecInfo()); fullVelBuffer->zero();
  fullAccBuffer   = new GenVector<double>(NonLinDynamic::solVecInfo()); fullAccBuffer->zero();
  fullVel2Buffer  = new GenVector<double>(NonLinDynamic::solVecInfo()); fullVel2Buffer->zero();
  fullDummyBuffer = new GenVector<double>(NonLinDynamic::solVecInfo()); fullDummyBuffer->zero();

  //initialize system state vector container for use in multi domain dynamic post processor
  curState = new SysState<GenVector<double> >( *fullDispBuffer, *fullVelBuffer, *fullAccBuffer, *fullVel2Buffer);
}  //end preProcessing

void 
ROMPostProcessingDriver::bufferReducedFiles()
{
  // get output information needed for parsing reduced data
  numConversionFiles = domain->solInfo().numRODFile;
  int skipTime = domain->solInfo().skipPodRom;

  // loop over reduced coordinate files
  for(int i = 0; i < numConversionFiles; i++) {
    //have all threads parse the reduced coordinate input file
    //there should be plenty of memory per node since projectionSubspaceSize is small
    std::ifstream reducedCoordFile(domain->solInfo().RODConversionFiles[i].c_str());
    if(reducedCoordFile.is_open()) {
      if(skipTime > 1) filePrint(stderr, " ... Skipping every %3d snapshots   ...\n", skipTime);

#ifdef ANDROID
      float time, dummyVar; // XXX DEBUG ANDROID
#else
      double time, dummyVar;
#endif
      int datatype, podsize, skipCounter;
      skipCounter = skipTime; // need to include t0 
      reducedCoordFile>>datatype; reducedCoordFile>>podsize;
      DataType.push_back(std::make_pair(datatype, podsize));
      switch(DataType[i].first) {
        case 0 : { // read reduced acceleration data
          filePrint(stderr, " ... Buffering Reduced Accelerations...\n");
          std::vector<double> timestamps;
          while(reducedCoordFile>>time) {
            if(skipCounter == skipTime) {
              skipCounter = 1;
              timestamps.push_back(time);
              for(int j = 0; j < podsize; j++) {
                reducedCoordFile>>dummyVar;
                reducedAccBuffer.push_back(dummyVar);
              }
              filePrint(stderr,"\r Timestamp = %f", time);
            }
            else {
              for(int j = 0; j < podsize; j++) 
                reducedCoordFile>>dummyVar;
              skipCounter += 1;
            }
          }
          TimeStamps.push_back(timestamps);
          filePrint(stderr,"\n");
        } break;
        case 1 : { // read reduced displacement data
          filePrint(stderr, " ... Buffering Reduced Displacements...\n");
          std::vector<double> timestamps;
          while(reducedCoordFile>>time) {
            if(skipCounter == skipTime) {
              skipCounter = 1;
              timestamps.push_back(time);
              for(int j = 0; j < podsize; j++) {
                reducedCoordFile>>dummyVar;
                reducedDispBuffer.push_back(dummyVar);
              }
              filePrint(stderr,"\r Timestamp = %f", time);
            }
            else {
              for(int j = 0; j < podsize; j++) 
                reducedCoordFile>>dummyVar;
              skipCounter += 1;
            }
          }
          TimeStamps.push_back(timestamps);
          filePrint(stderr,"\n");
        } break;
        case 2 : { // read reduced velocity data
          filePrint(stderr, " ... Buffering Reduced Velocities   ...\n");
          std::vector<double> timestamps;
          while(reducedCoordFile>>time) {
            if(skipCounter == skipTime){
              skipCounter = 1;
              timestamps.push_back(time);
              for(int j = 0; j < podsize; j++) {
                reducedCoordFile>>dummyVar;
                reducedVelBuffer.push_back(dummyVar);
              }
              filePrint(stderr,"\r Timestamp = %f", time);
            }
            else {
              for(int j = 0; j < podsize; j++)
                reducedCoordFile>>dummyVar;
              skipCounter += 1;
            }
          }
          TimeStamps.push_back(timestamps);
          filePrint(stderr,"\n");
        } break;
        default :
          filePrint(stderr, "\n... ROD conversion only supports Acceleration, Displacement, and Velocity ...\n");
      }
    }
    else {
      filePrint(stderr,"\nFailure to open file \n");
    }

    if(i != 0) { 
      if(DataType[i].second != DataType[i-1].second) {
        filePrint(stderr,"\n *** WARNING: Incompatible Input files %f %f\n", DataType[i-1].second, DataType[i].second);
      }
    }

  } //end loop over input files, finished reading reduced data

  if(DataType.size() > 0) projectionSubspaceSize = DataType[0].second;
}

void
ROMPostProcessingDriver::solve()
{
  preProcess();

  GeomState *geomState = new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes(), 
                                       &domain->getElementSet(),domain->getNodalTemperatures());
  Vector fext(solVecInfo()), constantForce(solVecInfo());
  getConstForce(constantForce);

  int counter = 0; //TODO: make this portion more general so it doesn't depend on the assumption
                   //that all files have matching timestamps
  if(TimeStamps.size() > 0) {
    for(std::vector<double>::iterator it = TimeStamps[0].begin(); it != TimeStamps[0].end(); it++) {

      // load current state for output 
      for(int i = 0; i < numConversionFiles; i++) {
        GenVector<double> buffer(projectionSubspaceSize);
        switch(DataType[i].first) {
          case 0 :
            for(int j = 0; j < projectionSubspaceSize; j++) 
              buffer[j] = reducedAccBuffer[counter*projectionSubspaceSize+j];
            normalizedBasis_.expand(buffer, *fullAccBuffer);
            break;
          case 1 :
            for(int j = 0; j < projectionSubspaceSize; j++)
              buffer[j] = reducedDispBuffer[counter*projectionSubspaceSize+j];
            normalizedBasis_.expand(buffer, *fullDispBuffer);
            break;
          case 2 :
            if(counter != 0) *fullVel2Buffer = *fullVelBuffer;
            for(int j = 0; j < projectionSubspaceSize; j++) 
              buffer[j] = reducedVelBuffer[counter*projectionSubspaceSize+j];
            normalizedBasis_.expand(buffer, *fullVelBuffer);
            break;
          default :
            break; 
        }
      }

      geomState->explicitUpdate(domain->getNodes(), *fullDispBuffer);
      geomState->transform(*fullVelBuffer, 0, true); // transform angular velocity from the 1st time derivative
                                                     // of the (unscaled) total rotation vector to convected
      geomState->transform(*fullAccBuffer, 4, true); // transform angular acceleration from the 2nd time derivative
                                                     // of the (unscaled) total rotation vector to convected
      domain->updateStates(geomState,*geomState,allCorot,*it);
      getExternalForce(fext, constantForce, counter, *it, geomState, *fullDummyBuffer, *fullDummyBuffer,
                       domain->solInfo().getTimeStep()/2);
      dynamOutput(geomState, *fullVelBuffer, *fullVel2Buffer, *it, counter-1, fext, *fullDummyBuffer, 
                  *fullAccBuffer, geomState);

      filePrint(stderr,"\r ... ROM Conversion Loop: t = %9.3e, %3d%% complete ...",
                *it, int(*it/(TimeStamps[0].back())*100));

      counter += 1;
    } //end of loop over time stamps
    filePrint(stderr,"\n");
  }
#ifdef USE_EIGEN3
  if(domain->solInfo().sensitivity) { // XXX should be inside loop ?
    Vector elementInternalForce(domain->maxNumDOF(), 0.0);
    Vector residual(domain->numUncon(), 0.0);
    residual = fext;
    domain->getStiffAndForce(*geomState, elementInternalForce, allCorot, kelArray, residual,
                             1.0, TimeStamps[0].back(), geomState, (Vector*) NULL, 
                             ((melArray) ? melArray : NULL));
    reBuild(*geomState, counter, domain->solInfo().getTimeStep(), TimeStamps[0].back());
    allSens->residual = new Eigen::Matrix<double,Eigen::Dynamic,1>(residual.size());
    *(allSens->residual) = Eigen::Map<Eigen::VectorXd>(residual.data(),residual.size());    
    sensitivityAnalysis(geomState, geomState);
  }
#endif

  delete geomState;
}

} //end namespace Rom

Rom::DriverInterface *ROMPostProcessingDriverNew(Domain *D) {
  return new Rom::ROMPostProcessingDriver(D);
}
