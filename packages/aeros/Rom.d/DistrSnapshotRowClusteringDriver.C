#include "DistrSnapshotRowClusteringDriver.h"

#include <Utils.d/DistHelper.h>

#if defined (USE_SCALAPACK) && defined (USE_EIGEN3)
#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrSnapshotRowClusteringSolver.h"
#include "DistrSvdOrthogonalization.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <algorithm>
#include <memory>
#include <cassert>
#include <sstream>

namespace Rom {

DistrSnapshotRowClusteringDriver::DistrSnapshotRowClusteringDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

//Non-member functions
//===============
template<int DOFS_PER_NODE>
void readIntoSolver(DistrSnapshotRowClusteringSolver &solver, DistrNodeDofBuffer<DOFS_PER_NODE> &inputBuffer, int& nodeCount,
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
void writeOutofSolver(DistrSnapshotRowClusteringSolver &solver, DistrNodeDofBuffer<DOFS_PER_NODE> &outputBuffer, int nodeCount,
                      DistrVecNodeDofConversion<DOFS_PER_NODE> &converter, const DistrInfo &vectorSize, BasisId::Type workload,
                      BasisId::Level type, int numClusters, std::vector<int> &basisStateOffset, Communicator *comm)
{
  FileNameInfo fileInfo;
  // loop over each cluster 
  for(int i=0; i<numClusters; ++i) {
    // get number of snapshots
    int clusterDim = solver.colCount();
    std::string fileName = BasisFileId(fileInfo, workload, type);
    std::ostringstream ss;
    ss << ".cluster" << i+1 << ".of." << numClusters;
    fileName.append(ss.str());
    // open output file for clustered snapshots named same as input file with new file extension
    DistrBasisOutputFile outputFile(fileName,
                                    nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                    comm, false, DOFS_PER_NODE);
    filePrint(stderr, " ... Writing %d row-clustered snapshots to file %s ...\n", clusterDim, fileName.c_str());
    // loop over each snapshot
    for (int iVec = 0; iVec < clusterDim; ++iVec) {
      double * const vecBuffer = const_cast<double *>(solver.matrixColBuffer(iVec));
      double timeStampBuffer = solver.getColTimeStamp(iVec);
      const GenStackDistVector<double> vec(vectorSize, vecBuffer);
      GenDistrVector<double> vec2(vectorSize);
      vec2.zero();
      for (int j = 0; j < solver.clusterLocalRowCount(i); ++j) {
         int k = solver.clusterLocalRow(i,j);
         vec2[k] = vec[k];
      }
      if(DOFS_PER_NODE == 1)
        converter.paddedNodeDof6(vec2, outputBuffer);
      else
        converter.unpaddedNodeDof6(vec2, outputBuffer);
      outputFile.stateAdd(outputBuffer, timeStampBuffer);
    }
  }
}

void
DistrSnapshotRowClusteringDriver::solve() {
  
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
  int numClusters = domain->solInfo().rowClustering;
  const int globalProbSize = domain->getCDSA()->size();

  DistrSnapshotRowClusteringSolver solver(comm_, globalProbSize, snapBasisStateCount, localLength, numClusters, blockSize);
  solver.setNNLSTolerance(domain->solInfo().tolPodRom);
 
  int solverCol = 0;
  std::vector<int> basisStateOffset;
  DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  readIntoSolver<6>(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::STATE, BasisId::SNAPSHOTS,
                    domain->solInfo().snapfiPodRom.size(), solverCol, basisStateOffset, domain->solInfo().skipPodRom); // read in snapshots

  filePrint(stderr, " ... Partitioning snapshots into %d row clusters ...\n", numClusters);
  solver.solverTypeIs(domain->solInfo().solverTypeCluster);
  solver.solve(); // this calls the clustering routine
  numClusters = solver.getNumClusters();  

  DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
  writeOutofSolver<6>(solver, outputBuffer, nodeCount, converter, distrInfo, BasisId::STATE, BasisId::SNAPSHOTS,
                      numClusters, basisStateOffset, comm_);
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrSnapshotRowClusteringDriverNew(Domain *domain) {
  return new Rom::DistrSnapshotRowClusteringDriver(domain, structCom);
}

#else

Rom::DriverInterface *distrSnapshotRowClusteringDriverNew(Domain *domain) {
  filePrint(stderr, " *** ERROR: requested driver requires ScaLAPACK and Eigen libraries\n");
  exit(-1);
  return 0;
}

#endif
