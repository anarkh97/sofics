#include "DistrPositiveDualBasisDriver.h"

#include <Utils.d/DistHelper.h>

#if defined (USE_SCALAPACK) && defined (USE_EIGEN3)
#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrNonnegativeMatrixFactorization.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <algorithm>
#include <memory>
#include <cassert>
#include <sstream>

namespace Rom {

DistrPositiveDualBasisDriver::DistrPositiveDualBasisDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

//Non-member functions
//===============
void readIntoSolver(DistrNonnegativeMatrixFactorization &solver, DistrNodeDof1Buffer &inputBuffer, int& nodeCount,
                    DistrVecNodeDof1Conversion &converter, const DistrInfo &vectorSize,
                    BasisId::Level type, int numEntries, int& solverCol, int skipFactor=1)
{
  const BasisId::Type workload = BasisId::DUALSTATE;
  FileNameInfo fileInfo;
  for(int i = 0; i < numEntries; i++) {
    std::string fileName = BasisFileId(fileInfo,workload,type,i);
    DistrBasisInputFileTemplate<1> inputFile(fileName);
    filePrint(stderr, " ... Reading in Snapshot file: %s ...\n", fileName.c_str());
    nodeCount = inputFile.nodeCount();
    int basisStateCount = 1+(inputFile.stateCount()-1)/skipFactor;
    {
      int count = 0;
      int skipCounter = skipFactor;
      while (count < basisStateCount) {
        assert(inputFile.validCurrentState());
        inputFile.currentStateBuffer(inputBuffer);

        if (skipCounter >= skipFactor) {
          double *vecBuffer = solver.matrixColBuffer(solverCol);
          GenStackDistVector<double> vec(vectorSize, vecBuffer);

          converter.vector(inputBuffer, vec);
          //converter.paddedMasterVector(inputBuffer, vec);
          //std::fill(vecBuffer + vectorSize.totLen(), vecBuffer + solver.localRows(), 0.0);
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

void
DistrPositiveDualBasisDriver::solve() {
  
  MultiDomainDynam::preProcess();

  const int blockSize = domain->solInfo().svdBlockSize;

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  for(int i=0; i<decDomain->getNumSub(); ++i) decDomain->getSubDomain(i)->computeMasterFlag(decDomain->getMpcToSub());
  DistrTrivialMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                          SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()),
                                          domain->getNumCTC(), blockSize, comm_, decDomain->getNumSub());

  DistrVecNodeDof1Conversion converter(decDomain->getAllSubDomains(),
                                       decDomain->getAllSubDomains() + decDomain->getNumSub(),
                                       domain->getNumCTC(), blockSize, comm_, decDomain->getNumSub());

  FileNameInfo fileInfo;
  const BasisId::Type workload = BasisId::DUALSTATE;

  const int skipFactor = domain->solInfo().skipPodRom;
  int stateCount = 0;
  int nodeCount = 0;
  int snapBasisStateCount = 0;
  if(!domain->solInfo().dsvPodRomFile.empty()) {
    for(int i = 0; i < domain->solInfo().dsvPodRomFile.size(); i++) {
      std::string fileName = BasisFileId(fileInfo, workload, BasisId::SNAPSHOTS, i);
      DistrBasisInputFileTemplate<1> inputFile(fileName);
      snapBasisStateCount += 1+(inputFile.stateCount()-1)/skipFactor;
    }
  }
  else {
    filePrint(stderr, " *** ERROR: no files provided\n");
    exit(-1);
  }

  DistrInfo distrInfo;
  decDomain->makeBlockCyclicDistrInfo(distrInfo, domain->getNumCTC(), blockSize);

  const int localLength = distrInfo.totLen();
  int orthoBasisDim = domain->solInfo().maxSizePodRom;
  const int globalProbSize = domain->getNumCTC();
 
  DistrNonnegativeMatrixFactorization solver(comm_, globalProbSize, snapBasisStateCount, localLength, orthoBasisDim, blockSize,
                                             domain->solInfo().nmfMaxIter, domain->solInfo().nmfTol, domain->solInfo().use_nmf,
                                             domain->solInfo().nmfNumSub, domain->solInfo().nmfPqnNumInnerIter, domain->solInfo().nmfPqnAlpha);
 
  if(domain->solInfo().romEnergy > 1e-16)
    solver.setEnergy(domain->solInfo().romEnergy);

  int solverCol = 0;
  DistrNodeDof1Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  readIntoSolver(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::SNAPSHOTS,
                 domain->solInfo().dsvPodRomFile.size(), solverCol, domain->solInfo().skipPodRom); // read in snapshots

  filePrint(stderr, " ... Computation of a positive basis of size %d ...\n", orthoBasisDim);
  solver.solve();

  if(domain->solInfo().romEnergy > 1e-16)
    orthoBasisDim = solver.basisSize();

  std::string fileName = BasisFileId(fileInfo, workload, BasisId::POD);
  std::ostringstream ss;
  ss << orthoBasisDim;
  //fileName.append(ss.str());
  DistrNodeDof1Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
  DistrBasisOutputFile outputFile(fileName,
                                  nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                  comm_, false, 1);
  filePrint(stderr, " ... Writing positive basis to file %s ...\n", fileName.c_str());
  for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
    double * const vecBuffer = const_cast<double *>(solver.basisColBuffer(iVec));
    const GenStackDistVector<double> vec(distrInfo, vecBuffer);
    converter.paddedNodeDof6(vec, outputBuffer);
    outputFile.stateAdd(outputBuffer, 1.0);
  }
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrPositiveDualBasisDriverNew(Domain *domain) {
  return new Rom::DistrPositiveDualBasisDriver(domain, structCom);
}

#else

Rom::DriverInterface *distrPositiveDualBasisDriverNew(Domain *domain) {
  filePrint(stderr, " *** ERROR: requested driver requires ScaLAPACK and Eigen libraries\n");
  exit(-1);
  return 0;
}

#endif
