#include "DistrSnapshotClusteringDriver.h"

#include <Utils.d/DistHelper.h>

#if defined (USE_SCALAPACK) && defined (USE_EIGEN3)
#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrSnapshotClusteringSolver.h"
#include "DistrSvdOrthogonalization.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <algorithm>
#include <memory>
#include <cassert>
#include <sstream>

namespace Rom {

DistrSnapshotClusteringDriver::DistrSnapshotClusteringDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

//Non-member functions
//===============
template<int DOFS_PER_NODE>
void readIntoSolver(DistrSnapshotClusteringSolver &solver, DistrNodeDofBuffer<DOFS_PER_NODE> &inputBuffer, int& nodeCount,
                    DistrVecNodeDofConversion<DOFS_PER_NODE> &converter, const DistrInfo &vectorSize, BasisId::Type workload,
                    BasisId::Level type, int numEntries, int& solverCol, std::vector<int> &basisStateOffset,
                    int skipFactor=1)
{
  FileNameInfo fileInfo;
  basisStateOffset.clear(); basisStateOffset.push_back(0);
  for(int i = 0; i < numEntries; i++) {
    std::string fileName = BasisFileId(fileInfo,workload,type,i);
    DistrBasisInputFileTemplate<DOFS_PER_NODE> inputFile(fileName);
    filePrint(stderr, " ... Reading in Snapshot file: %s ...\n", fileName.c_str());
    nodeCount = inputFile.nodeCount();
    int basisStateCount = 1+(inputFile.stateCount()-1)/skipFactor;
    basisStateOffset.push_back(basisStateOffset[i]+basisStateCount);
    {
      int count = 0;
      int skipCounter = skipFactor;
      while (count < basisStateCount) {
        assert(inputFile.validCurrentState());
        double curTimeStamp = inputFile.currentStateHeaderValue();
        inputFile.currentStateBuffer(inputBuffer);
        if (skipCounter >= skipFactor) {
          double *vecBuffer = solver.matrixColBuffer(solverCol);
          GenStackDistVector<double> vec(vectorSize, vecBuffer);
          solver.setColTimeStamp(solverCol,curTimeStamp);

          if(DOFS_PER_NODE == 1) {
            converter.vector(inputBuffer, vec);
          } else {
            converter.unpaddedMasterVector(inputBuffer, vec);
          }
          std::fill(vecBuffer + vectorSize.totLen(), vecBuffer + solver.localRows(), 0.0);
          skipCounter = 1;
          ++solverCol;
          ++count;
        } else {
          ++skipCounter;
        }
        inputFile.currentStateIndexInc();
      }
    }
  }
}

template<int DOFS_PER_NODE>
void writeOutofSolver(DistrSnapshotClusteringSolver &solver, DistrNodeDofBuffer<DOFS_PER_NODE> &outputBuffer, int nodeCount,
                      DistrVecNodeDofConversion<DOFS_PER_NODE> &converter, const DistrInfo &vectorSize, BasisId::Type workload,
                      BasisId::Level type, int numClusters, std::vector<int> &basisStateOffset, Communicator *comm)
{
  FileNameInfo fileInfo;
  // loop over each cluster 
  for(int i=0; i<numClusters; ++i) {
    // get number of snapshots in that cluster
    int clusterDim = solver.clusterColCount(i);
    std::string fileName = BasisFileId(fileInfo, workload, type);
    std::ostringstream ss;
    ss << ".cluster" << i+1 << ".of." << numClusters;
    fileName.append(ss.str());
    // open output file for clustered snapshots named same as input file with new file extension
    DistrBasisOutputFile outputFile(fileName,
                                    nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                    comm, false, DOFS_PER_NODE);
    filePrint(stderr, " ... Writing %d clustered snapshots to file %s ...\n", clusterDim, fileName.c_str());
    // loop over each snapshot in current cluster
    for (int iVec = 0; iVec < clusterDim; ++iVec) {
      double * const vecBuffer = const_cast<double *>(solver.clusterColBuffer(i,iVec));
      double timeStampBuffer = solver.getClusterColTimeStamp(i,iVec);
      const GenStackDistVector<double> vec(vectorSize, vecBuffer);
      if(DOFS_PER_NODE == 1)
        converter.paddedNodeDof6(vec, outputBuffer);
      else
        converter.unpaddedNodeDof6(vec, outputBuffer);
      outputFile.stateAdd(outputBuffer, timeStampBuffer);
    }

    // this is writes a file used for hyper reduction that tells AEROS which parameter each snapshot is associated with
    if(workload == BasisId::STATE || (domain->solInfo().solverTypeCluster == 2 && workload == BasisId::VELOCITY)) {
      std::string fileName2 = fileName+".sources";
      std::ofstream sources(fileName2.c_str());
      for (int iVec = 0; iVec < clusterDim; ++iVec) {
        int j = solver.clusterCol(i,iVec);
        int k;
        for(k=1; k<basisStateOffset.size(); ++k) if(j < basisStateOffset[k]) break;
        sources << k << " " << j+1-basisStateOffset[k-1] << std::endl;
      }

      // open output file for centroids of clusters
      fileName.append(".centroid");
      DistrBasisOutputFile outputFile2(fileName,
                                       nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                       comm, false);
      filePrint(stderr, " ... Writing snapshot cluster centroid to file %s ...\n", fileName.c_str());
      double * const vecBuffer = const_cast<double *>(solver.clusterCentroidBuffer(i));
      const GenStackDistVector<double> vec(vectorSize, vecBuffer);
      converter.unpaddedNodeDof6(vec, outputBuffer);
      outputFile2.stateAdd(outputBuffer, 1.0);
    }
  }
}

void
DistrSnapshotClusteringDriver::solve() {
  
  MultiDomainDynam::preProcess();

  const int blockSize = domain->solInfo().svdBlockSize;

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  for(int i=0; i<decDomain->getNumSub(); ++i) decDomain->getSubDomain(i)->computeMasterFlag(decDomain->getMpcToSub());
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(),
                                       decDomain->getAllSubDomains() + decDomain->getNumSub());

  FileNameInfo fileInfo;
  const BasisId::Type workload = BasisId::STATE;

  const int skipFactor = domain->solInfo().skipPodRom;
  int stateCount = 0;
  int nodeCount = 0;
  int snapBasisStateCount = 0;
  if(!domain->solInfo().snapfiPodRom.empty()) {
    for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
      std::string fileName = BasisFileId(fileInfo, workload, BasisId::SNAPSHOTS, i);
      DistrBasisInputFileTemplate<6> inputFile(fileName);
      snapBasisStateCount += 1+(inputFile.stateCount()-1)/skipFactor;
    }
  }
  else {
    filePrint(stderr, " *** ERROR: no files provided\n");
    exit(-1);
  }

  DistrInfo distrInfo; decDomain->makeNonOverlappingDistrInfo(distrInfo);

  int localLength = distrInfo.totLen();
  int numClusters = domain->solInfo().clustering;
  const int globalProbSize = domain->getCDSA()->size();

  DistrSnapshotClusteringSolver solver(comm_, globalProbSize, snapBasisStateCount, localLength, numClusters, blockSize);
  solver.setNNLSTolerance(domain->solInfo().tolPodRom);
 
  int solverCol = 0;
  std::vector<int> basisStateOffset;
  DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  readIntoSolver<6>(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::STATE, BasisId::SNAPSHOTS,
                    domain->solInfo().snapfiPodRom.size(), solverCol, basisStateOffset, domain->solInfo().skipPodRom); // read in snapshots

  filePrint(stderr, " ... Partitioning snapshots into %d clusters ...\n", numClusters);
  solver.solverTypeIs(domain->solInfo().solverTypeCluster);
  solver.solve(); // this calls the clustering routine
  numClusters = solver.getNumClusters();  

  DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
  writeOutofSolver<6>(solver, outputBuffer, nodeCount, converter, distrInfo, BasisId::STATE, BasisId::SNAPSHOTS,
                      numClusters, basisStateOffset, comm_);

  if(!domain_->solInfo().velocPodRomFile.empty()) {
    solverCol = 0;
    readIntoSolver<6>(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS,
                   domain->solInfo().velocPodRomFile.size(), solverCol, basisStateOffset, domain->solInfo().skipPodRom); // read in velocity snapshots
    if(domain->solInfo().solverTypeCluster == 2) solver.recomputeCentroids(); // for clustering on velocity snapshots, but still need displacement centroids
    writeOutofSolver<6>(solver, outputBuffer, nodeCount, converter, distrInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS,
                     numClusters, basisStateOffset, comm_);
  }

  if(!domain_->solInfo().accelPodRomFile.empty()) {
    solverCol = 0;
    readIntoSolver<6>(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS,
                   domain->solInfo().accelPodRomFile.size(), solverCol, basisStateOffset, domain->solInfo().skipPodRom); // read in acceleration snapshots
    writeOutofSolver<6>(solver, outputBuffer, nodeCount, converter, distrInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS,
                     numClusters, basisStateOffset, comm_);
  }


  // read in dual snapshots
  if(!domain_->solInfo().dsvPodRomFile.empty()){
    // huge pain in the ass
    solverCol = 0;

    DistrTrivialMasterMapping masterMapping1(SubDomIt(decDomain->getAllSubDomains()),
                                          SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()),
                                          domain->getNumCTC(), blockSize, comm_, decDomain->getNumSub());

    DistrVecNodeDof1Conversion converter1(decDomain->getAllSubDomains(),
                                          decDomain->getAllSubDomains() + decDomain->getNumSub(),
                                          domain->getNumCTC(), blockSize, comm_, decDomain->getNumSub());

    DistrNodeDof1Buffer inputBuffer1(masterMapping1.localNodeBegin(), masterMapping1.localNodeEnd());

    DistrNodeDof1Buffer outputBuffer1(masterMapping1.masterNodeBegin(), masterMapping1.masterNodeEnd()); 

    DistrInfo distrInfo1;
    decDomain->makeBlockCyclicDistrInfo(distrInfo1, domain->getNumCTC(), blockSize);

    readIntoSolver<1>(solver, inputBuffer1, nodeCount, converter1, distrInfo1, BasisId::DUALSTATE, BasisId::SNAPSHOTS,
                   domain->solInfo().dsvPodRomFile.size(), solverCol, basisStateOffset, domain->solInfo().skipPodRom); // read in dual variable snapshots
    writeOutofSolver(solver, outputBuffer1, nodeCount, converter1, distrInfo1, BasisId::DUALSTATE, BasisId::SNAPSHOTS,
                     numClusters, basisStateOffset, comm_);
  }

}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrSnapshotClusteringDriverNew(Domain *domain) {
  return new Rom::DistrSnapshotClusteringDriver(domain, structCom);
}

#else

Rom::DriverInterface *distrSnapshotClusteringDriverNew(Domain *domain) {
  filePrint(stderr, " *** ERROR: requested driver requires ScaLAPACK and Eigen libraries\n");
  exit(-1);
  return 0;
}

#endif
