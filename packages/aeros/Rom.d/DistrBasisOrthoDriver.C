#include "DistrBasisOrthoDriver.h"

#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrSvdOrthogonalization.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"
#include "XPostInputFile.h"
#include "RobCodec.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>

namespace Rom {

DistrBasisOrthoDriver::DistrBasisOrthoDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

//Non-member functions
//===============
void readIntoSolver(DistrSvdOrthogonalization &solver, DistrNodeDof6Buffer &inputBuffer, int& nodeCount,
                    DistrVecNodeDof6Conversion &converter, MDDynamMat *dynOps, const DistrInfo &vectorSize,
                    BasisId::Level type, int numEntries, int& solverCol, int skipFactor=1)
{
  const BasisId::Type workload = BasisId::STATE;
  FileNameInfo fileInfo;
  for(int i = 0; i < numEntries; i++){
    DistrBasisInputFile inputFile(BasisFileId(fileInfo,workload,type,i));
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

          converter.paddedMasterVector(inputBuffer, vec);
          // TODO Multiply by weighting factor if given in input file
          if(geoSource->getMRatio() == 0 && domain->solInfo().normalize == 1) {
            dynOps->dynMat->squareRootMult(vec);
          }
          if(type== BasisId::ROB) vec *= inputFile.currentStateHeaderValue(); // multiply by singular values for robfiles
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

//Member functions
//================
const DistrInfo&
DistrBasisOrthoDriver::vectorSize() const
{
  return decDomain->masterSolVecInfo();
}

void
DistrBasisOrthoDriver::solve() {
  
  MultiDomainDynam::preProcess();

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
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
  int robBasisStateCount = 0;
  if(!domain->solInfo().snapfiPodRom.empty()) {
    for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
      BasisFileId basisFileId(fileInfo, workload, BasisId::SNAPSHOTS, i);
      std::string fileName = basisFileId.name();
      if(!basisFileId.isBinary()) {
        filePrint(stderr," ... Convert ASCII file to binary   ...\n");
        convert_rob<Rom::XPostInputFile, Rom::BasisBinaryOutputFile>(fileName, fileName+".bin");
        fileName = domain->solInfo().snapfiPodRom[i] = fileName+".bin";
      }
      DistrBasisInputFile inputFile(fileName);
      snapBasisStateCount += 1+(inputFile.stateCount()-1)/skipFactor;
    }
  }
  if(!domain->solInfo().robfi.empty()) {
    for(int i = 0; i < domain->solInfo().robfi.size(); i++) {
      BasisFileId basisFileId(fileInfo, workload, BasisId::ROB, i);
      std::string fileName = basisFileId.name();
      if(!basisFileId.isBinary()) {
        filePrint(stderr," ... Convert ASCII file to binary   ...\n");
        convert_rob<Rom::XPostInputFile, Rom::BasisBinaryOutputFile>(fileName, fileName+".bin");
        fileName = domain->solInfo().robfi[i] = fileName+".bin";
      }
      DistrBasisInputFile inputFile(fileName);
      robBasisStateCount += inputFile.stateCount();
    }
  }

  DistrSvdOrthogonalization solver(comm_, comm_->numCPUs(), 1);
  
  const int blockSize = domain->solInfo().svdBlockSize; // default: 64
  {
    solver.blockSizeIs(blockSize);
  }

  const int localLength = decDomain->solVecInfo().totLen();
  const int masterLength = decDomain->masterSolVecInfo().masterLen();
  {
    const int maxLocalLength = comm_->globalMax(localLength);
    // compute upper bound for global problem size such that solver.localRows() is 
    // always greater than or equal to localLength. see numroc.f (ScaLAPACK)
    const int globalProbSize = ((maxLocalLength/blockSize+1)*solver.rowCpus()+1)*blockSize+1;
    solver.problemSizeIs(globalProbSize, snapBasisStateCount+robBasisStateCount);
    assert(solver.localRows() >= localLength);
  }
  const int singularValueCount = std::min(masterLength, snapBasisStateCount+robBasisStateCount);

  double mratio = geoSource->getMRatio();
  if(mratio != 0 && domain->solInfo().normalize == 1) {
    filePrint(stderr, " *** ERROR: \"mnorma 1\" is only available for LUMPED mass matrices when a decomposition is used.\n");
    exit(-1);
  }
  if(domain->solInfo().normalize <= 0) domain->solInfo().solvercntl->type = static_cast<SolverSelection>(-1); // no solver required
  // Assembling mass matrix
  MDDynamMat *dynOps = MultiDomainDynam::buildOps(1.0, 0.0, 0.0);
  assert(dynOps->M);
  if(domain->solInfo().normalize == 1) assert(dynOps->dynMat);
 
  if(domain->solInfo().snapfiPodRom.empty() && domain->solInfo().robfi.empty()) {
    filePrint(stderr, " *** ERROR: no files provided\n");
    exit(-1);
  }
   
  int solverCol = 0;
  DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  readIntoSolver(solver, inputBuffer, nodeCount, converter, dynOps, decDomain->solVecInfo(), BasisId::SNAPSHOTS,
                 domain->solInfo().snapfiPodRom.size(), solverCol, domain->solInfo().skipPodRom); // read in snapshots

  readIntoSolver(solver, inputBuffer, nodeCount, converter, dynOps, decDomain->solVecInfo(), BasisId::ROB,
                 domain->solInfo().robfi.size(), solverCol); // read in robs

  if(domain->solInfo().robcSolve) solver.solve();

  int podVectorCount = domain_->solInfo().maxSizePodRom ?
                       std::min(domain_->solInfo().maxSizePodRom, singularValueCount) :
                       singularValueCount;

  // Compute and output the truncation error
  {
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    fileName += ".truncation_error.txt";
    std::vector<double> toto(podVectorCount+1);
    toto[podVectorCount] = 0;
    for (int iVec = podVectorCount-1; iVec >= 0; --iVec) {
      toto[iVec] = toto[iVec+1]+solver.singularValue(iVec); // running sum
    }
    std::ofstream out;
    if(comm_->myID() == 0) out.open(fileName.c_str());
    bool reset = true;
    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      double energy = toto[iVec]/toto[0];
      if(energy < domain->solInfo().romEnergy && reset) {
        podVectorCount = iVec+1;
        reset = false;
      }
      if(comm_->myID() == 0) out << iVec+1 << " " << solver.singularValue(iVec) << " " << energy << std::endl;
    }
  }

  // Output solution
  std::string fileName = BasisFileId(fileInfo, workload, BasisId::POD);
  fileName.append(".orthonormalized");
  {
    DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
    DistrBasisOutputFile outputFile(fileName,
                                    nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                    comm_, false);

    if(domain->solInfo().normalize <= 0)
      filePrint(stderr, " ... Writing orthonormal basis of size %d to file %s ...\n", podVectorCount, fileName.c_str());
    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      double * const vecBuffer = (domain->solInfo().robcSolve) ? const_cast<double *>(solver.basisColBuffer(iVec))
                                                               : const_cast<double *>(solver.matrixColBuffer(iVec));
      const GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);
      converter.paddedNodeDof6(vec, outputBuffer);
      outputFile.stateAdd(outputBuffer, solver.singularValue(iVec));
    }
  }
  comm_->sync(); // this is important, do not delete.

  // Read back in output file to perform renormalization
  DistrVecBasis basis;
  {
    DistrBasisInputFile inputFile(fileName);
    DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
    basis.dimensionIs(podVectorCount, decDomain->masterSolVecInfo()); 
    int i = 0;
    for(DistrVecBasis::iterator it = basis.begin(),
        it_end = basis.end();
        it != it_end; ++it) {
      assert(inputFile.validCurrentState());

      inputFile.currentStateBuffer(inputBuffer);
      converter.vector(inputBuffer, *it);

      inputFile.currentStateIndexInc();
    }
  }

  DistrVecBasis normalizedBasis;
  if(domain->solInfo().normalize == 0) {
    // Old method: renormalize the orthonormal basis
    renormalized_basis(*dynOps->M, basis, normalizedBasis);
  }
  if(domain->solInfo().normalize == 1) {
    // New method: multiply by inverse square root of the mass matrix
    if(mratio == 0) { // lumped
      for(int col = 0; col < podVectorCount; col++) {
        dynOps->dynMat->inverseSquareRootMult(basis[col]);
      }
    }
    normalizedBasis = basis;
  }

  // Output the renormalized basis as separate file
  if(domain->solInfo().normalize >= 0) {
    std::string fileName = BasisFileId(fileInfo, workload, BasisId::POD);
    fileName.append(".massorthonormalized");
    DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
    DistrBasisOutputFile outputNormalizedFile(fileName, nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(), comm_, false);
    filePrint(stderr, " ... Writing mass-normalized basis of size %d to file %s ...\n", podVectorCount, fileName.c_str());
    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      converter.paddedNodeDof6(normalizedBasis[iVec], outputBuffer);
      outputNormalizedFile.stateAdd(outputBuffer, solver.singularValue(iVec));
    }
  }

  // Compute and output identity normalized basis if using new method
  if(domain->solInfo().normalize == 1) {
    MGSVectors(normalizedBasis);
    DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
    DistrBasisOutputFile outputOrthoNormalFile(fileName, nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(), comm_, false);
    filePrint(stderr, " ... Writing orthonormal basis of size %d to file %s ...\n", podVectorCount, fileName.c_str());
    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      converter.paddedNodeDof6(normalizedBasis[iVec], outputBuffer);
      outputOrthoNormalFile.stateAdd(outputBuffer, solver.singularValue(iVec));
    }
  }
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *domain) {
  return new Rom::DistrBasisOrthoDriver(domain, structCom);
}
