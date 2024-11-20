#include "PositiveDualBasisDriver.h"

#include "NonnegativeMatrixFactorization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#include <utility>
#include <algorithm>
#include <sstream>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

PositiveDualBasisDriver::PositiveDualBasisDriver(Domain *domain) :
  SingleDomainDynamic(domain)
{}

//Non-member functions
//====================
void readIntoSolver(NonnegativeMatrixFactorization &solver, VecNodeDof1Conversion &converter, BasisId::Level fileType,
                    int numEntries, int vectorSize, BasisId::Type type, int &colCounter, int skipTime=1)
{
  FileNameInfo fileInfo; 
  for(int i = 0 ; i < numEntries; i++) {
    std::string fileName = BasisFileId(fileInfo, type, fileType, i);
    BasisInputStream<1> input(fileName, converter);
    filePrint(stderr, " ... Reading in Snapshot file: %s ...\n", fileName.c_str());
    int skip = 1;
    for (int iCol = 0; iCol < input.size(); ++iCol) {
      if(skip == skipTime) {
        double *buffer = solver.matrixCol(colCounter);
        input >> buffer;
        assert(input);
        colCounter++;
        skip = 1;
      } else {
        SimpleBuffer<double> dummyVec;
        dummyVec.sizeIs(input.vectorSize());  
        double *dummyBuffer = dummyVec.array();
        input >> dummyBuffer;
        assert(input);
        ++skip;
      }
    }
  }
}

//Member functions
//====================

void
PositiveDualBasisDriver::solve() {
  // process model geometry
  SingleDomainDynamic::preProcess();

  // initialize solver
  NonnegativeMatrixFactorization solver(domain->solInfo().maxSizePodRom, domain->solInfo().use_nmf);
  solver.maxIterIs(domain->solInfo().nmfMaxIter);
  solver.toleranceIs(domain->solInfo().nmfTol);
  solver.numRandInitIs(domain->solInfo().nmfRandInit);
  solver.nmfcAlphaIs(domain->solInfo().nmfcAlpha);
  solver.nmfcBetaIs(domain->solInfo().nmfcBeta);
  solver.nmfcGammaIs(domain->solInfo().nmfcGamma);

  std::vector<BasisId::Type> workload;
  workload.push_back(BasisId::DUALSTATE);

  // check that some data has been provided
  if(domain->solInfo().dsvPodRomFile.empty() && domain->solInfo().muvPodRomFile.empty()  && domain->solInfo().robfi.empty()) {
    std::cerr << "*** ERROR: no files provided\n";
    exit(-1);
  }


  int vectorSize = 0; // size of vectors
  int sizeSnap = 0; // number of state snapshots
  int skipTime = domain->solInfo().skipPodRom;

  // see how many vectors will be read in and initialize converter structure
  FileNameInfo fileInfo;
  VecNodeDof1Conversion *converter; 
  BasisId::Type type; 
  int numFiles = 0;
  if(domain->solInfo().dsvPodRomFile.size() > 0){ // for snapshots collected from LMPCs or Constraint function elements
    type = BasisId::DUALSTATE;
    numFiles = domain->solInfo().dsvPodRomFile.size(); 
    converter = new VecNodeDof1Conversion(domain->getNumCTC());
  } else if(domain->solInfo().muvPodRomFile.size() > 0) { // for snapshots collected from FETI-Mortar method
    type = BasisId::MUSTATE;
    numFiles = domain->solInfo().muvPodRomFile.size();
    converter = new VecNodeDof1Conversion(*(domain->getDSA()));
  }
  for(int i = 0; i < numFiles; i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    BasisInputStream<1> input(fileName, *converter);
    vectorSize = input.vectorSize();
    sizeSnap += input.size()/skipTime;
  }

  // set basis dimension and allocate space in the solver
  int maxBasisDimension = domain->solInfo().maxSizePodRom + (domain->solInfo().nmfDelROBDim)*(domain->solInfo().nmfNumROBDim-1);
  solver.matrixSizeIs(vectorSize, sizeSnap);
  solver.robSizeIs(vectorSize, maxBasisDimension);

  // read the snapshots in to the solver buffer
  int colCounter = 0;
  readIntoSolver(solver, *converter, BasisId::SNAPSHOTS, numFiles, vectorSize, type, colCounter, skipTime); 

  // solve for incrementally larger daul bases 
  for (int iBasis=0; iBasis < domain->solInfo().nmfNumROBDim; ++iBasis) {
    int orthoBasisDim = domain->solInfo().maxSizePodRom + iBasis*domain->solInfo().nmfDelROBDim;
    filePrint(stderr, " ... Computation of a positive basis of size %d ...\n", orthoBasisDim);
    solver.basisDimensionIs(orthoBasisDim);
    if (iBasis==0)
      solver.solve(0);
    else
      solver.solve(orthoBasisDim-domain->solInfo().nmfDelROBDim);

    std::string fileName = BasisFileId(fileInfo, type, BasisId::POD);
    std::ostringstream ss;
    ss << orthoBasisDim;
    fileName.append(ss.str()); 
    BasisOutputStream<1> output(fileName, *converter, false); 
    filePrint(stderr, " ... Writing positive basis to file %s ...\n", fileName.c_str());
    for (int iVec = 0; iVec < orthoBasisDim; ++iVec) { 
      output << std::make_pair(1.0, solver.robCol(iVec));
    }
  }  
}

void
PositiveDualBasisDriver::preProcess() {
  domain->preProcessing();
 
  // Build the constrained DofSetArray incorporating the boundary conditions 
  const int numdof = domain->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain->make_bc(bc.array(), bcx.array());
  domain->make_constrainedDSA(bc.array());
}

} /* end namespace Rom */

Rom::DriverInterface *positiveDualBasisDriverNew(Domain *domain) {
  return new Rom::PositiveDualBasisDriver(domain);
}
