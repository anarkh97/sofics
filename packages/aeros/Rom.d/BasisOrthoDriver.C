#include "BasisOrthoDriver.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"
#include "XPostInputFile.h"
#include "RobCodec.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Math.d/DiagMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#include <utility>
#include <algorithm>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

namespace { // anonymous

template <typename VecType>
class VectorTransform {
public:
  virtual void operator()(VecType &) const = 0;

  virtual ~VectorTransform() {}
};

template <typename VecType>
class NoOp : public VectorTransform<VecType> {
public:
  virtual void operator()(VecType &) const  { /* Nothing */ }
};

template <typename VecType>
class RefSubtraction : public VectorTransform<VecType> {
public:
  virtual void operator()(VecType &v) const;

//  explicit RefSubtraction(const Domain *);
  explicit RefSubtraction(const Vector &);
private:
  Vector ref_;
};

template <typename VecType>
void
RefSubtraction<VecType>::operator()(VecType &v) const {
  for (int i = 0; i < ref_.size(); ++i) {
    v[i] -= ref_[i];
  }
}

template <typename VecType>
RefSubtraction<VecType>::RefSubtraction(const Vector &ref) :
  ref_(ref)
{}

/*
template <typename VecType>
RefSubtraction<VecType>::RefSubtraction(const Domain *domain) :
  ref_(const_cast<Domain *>(domain)->numUncon())
{
  Vector dummy(const_cast<Domain *>(domain)->numUncon());
  const_cast<Domain *>(domain)->initDispVeloc(ref_, dummy, dummy, dummy);
}
*/

} // end anonymous namespace

BasisOrthoDriver::BasisOrthoDriver(Domain *domain) :
  SingleDomainDynamic(domain)
{}

//Non-member functions
//====================
void readIntoSolver(SvdOrthogonalization &solver, VecNodeDof6Conversion &converter, BasisId::Level fileType,
                    int numEntries, int vectorSize, std::unique_ptr<VectorTransform<double*> > &transform, BasisId::Type type,
                    int &colCounter, GenSparseMatrix<double> *fullMass, GenSolver<double> *fullMassSolver, int skipTime=1)
{
  FileNameInfo fileInfo; 
  for(int i = 0 ; i < numEntries; i++) {
    std::string fileName = BasisFileId(fileInfo, type, fileType, i);
    BasisInputStream<6> input(fileName, converter);
    if(fileType == BasisId::SNAPSHOTS) filePrint(stderr, " ... Reading in Snapshot file: %s ...\n", fileName.c_str());
    if(fileType == BasisId::ROB) filePrint(stderr, " ... Reading in ROB file: %s ...\n", fileName.c_str());
    int skip = 1;
    std::vector<double> s; int firstCol = colCounter;
    for (int iCol = 0; iCol < input.size(); ++iCol) {
      if(skip == skipTime) {
        double *buffer = solver.matrixCol(colCounter);
        std::pair<double, double *> data;
        data.second = buffer;
        input >> data;
        assert(input);
        s.push_back(data.first);
        colCounter++;
        if(geoSource->getMRatio() == 0 && domain->solInfo().normalize == 1) {
          fullMass->squareRootMult(buffer);
        }
        if(geoSource->getMRatio() != 0 && domain->solInfo().normalize == 1) {
          fullMassSolver->upperMult(buffer);
        }
        (*transform)(buffer);
        skip = 1;
      } 
      else {
        SimpleBuffer<double> dummyVec;
        dummyVec.sizeIs(input.vectorSize());  
        double *dummyBuffer = dummyVec.array();
        input >> dummyBuffer;
        assert(input);
        ++skip;
      }
    }
    // Multiply by weighting factor if specified in input file using sscali
    if((fileType == BasisId::ROB && domain->solInfo().flagrs) || (fileType == BasisId::SNAPSHOTS && domain->solInfo().flagss)) {
      for(int i = 0; i < s.size(); ++i) {
        double weight;
        if(fileType == BasisId::ROB) weight = s[i];
        else if(s.size() > 1) {
          if(i == 0) weight = std::sqrt(s[1]-s[0]);
          else if(i == s.size()-1) weight = std::sqrt(s[s.size()-1]-s[s.size()-2]);
          else weight = std::sqrt(((s[i]-s[i-1])+(s[i+1]-s[i]))/2);
        }
        double *buffer = solver.matrixCol(firstCol+i);
        for(int row = 0; row < vectorSize; row++) {
          buffer[row] *= weight;
        }
      }
    }
  }
}

//Member functions
//====================

void
BasisOrthoDriver::solve() {
  SingleDomainDynamic::preProcess();
  VecNodeDof6Conversion converter(*domain->getCDSA());
  FileNameInfo fileInfo;
  SvdOrthogonalization solver;

  std::vector<BasisId::Type> workload;
       
  if(domain->solInfo().statevectPodRom) {
	workload.push_back(BasisId::STATE);
        if(verboseFlag) fprintf(stderr," ... For State SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().residvectPodRom) {
	workload.push_back(BasisId::RESIDUAL);
	if(verboseFlag) fprintf(stderr," ... For Residual SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().jacobvectPodRom) {
        workload.push_back(BasisId::JACOBIAN);
	if(verboseFlag) fprintf(stderr," ... For Jacobian SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().forcevectPodRom) {
        workload.push_back(BasisId::FORCE);
	if(verboseFlag) fprintf(stderr," ... For Force SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().accelvectPodRom) {
        workload.push_back(BasisId::ACCELERATION);
	if(verboseFlag) fprintf(stderr," ... For Acceleration SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().velocvectPodRom) {
        workload.push_back(BasisId::VELOCITY);
        if(verboseFlag) fprintf(stderr," ... For Velocity SVD, workload size = %zd ...\n", workload.size());}
  else { workload.push_back(BasisId::STATE);
	if(verboseFlag) fprintf(stderr," ... For default SVD, workload size = %zd ...\n", workload.size());}



  // if we want to use affine subspaces (instead of linear subspaces) we need to subtract off an offset
  Vector *centroid; 
  if(domain->solInfo().subtractRefPodRom) { // read in centroid for snapshots
    std::cout << " ... Reading Basis Offset           ... " << std::endl;
    BasisInputStream<6> centIn(domain->solInfo().readInLocalBasesCent[0], converter);
    centroid = new Vector(centIn.vectorSize());
    centIn >> centroid->data();
  }
  typedef VectorTransform<double *> VecTrans;
  std::unique_ptr<VecTrans> transform(domain->solInfo().subtractRefPodRom ?
                                    static_cast<VecTrans *>(new RefSubtraction<double *>(*centroid)) :
                                    static_cast<VecTrans *>(new NoOp<double *>));

  double mratio = geoSource->getMRatio();
  // Assemble mass matrix, and factor if necessary
  AllOps<double> allOps;
  if(mratio != 0) {
#ifdef USE_EIGEN3
    allOps.M = domain->constructEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >();
#else
    allOps.M = domain->constructDBSparseMatrix<double>();
#endif
  }
  else {
    allOps.M = new DiagMatrix(domain->getCDSA());
  }
  domain->makeSparseOps<double>(allOps, 0.0, 1.0, 0.0, (GenSparseMatrix<double>*) NULL, kelArray, melArray);

  GenSparseMatrix<double> *fullMass = allOps.M;
  GenSolver<double> *fullMassSolver;
  if(mratio != 0 && domain->solInfo().normalize == 1) { 
    fullMassSolver = dynamic_cast<GenSolver<double>*>(fullMass);
    if(fullMassSolver) {
      filePrint(stderr, " ... Factoring mass matrix          ...\n");
      fullMassSolver->factor();
    }
    else {
      std::cerr << "*** ERROR: Cannot factorize mass matrix.\n";
      exit(-1);
    }
  }
  int vectorSize = 0; // size of vectors
  int sizeSnap   = 0; // number of state snapshots
  int sizeROB    = 0;
  int skipTime = domain->solInfo().skipPodRom;
  if(domain->solInfo().snapfiPodRom.empty() && domain->solInfo().robfi.empty()) {
    std::cerr << "*** ERROR: no files provided\n";
    exit(-1);
  }
 
  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;
    // Loop over snapshots
    for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
      BasisFileId basisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
      std::string fileName = basisFileId.name();
      if(!basisFileId.isBinary()) {
        filePrint(stderr," ... Convert ASCII file to binary   ...\n");
        convert_rob<Rom::XPostInputFile, Rom::BasisBinaryOutputFile>(fileName, fileName+".bin");
        fileName = domain->solInfo().snapfiPodRom[i] = fileName+".bin";
      }
      BasisInputStream<6> input(fileName, converter);
      vectorSize = input.vectorSize();
      sizeSnap += input.size()/skipTime;
    }

    // Loop over rob files 
    for(int i = 0; i < domain->solInfo().robfi.size(); i++) {
      BasisFileId basisFileId(fileInfo,type,BasisId::ROB, i);
      std::string fileName = basisFileId.name();
      if(!basisFileId.isBinary()) {
        filePrint(stderr," ... Convert ASCII file to binary   ...\n");
        convert_rob<Rom::XPostInputFile, Rom::BasisBinaryOutputFile>(fileName, fileName+".bin");
        fileName = domain->solInfo().robfi[i] = fileName+".bin";
      }
      BasisInputStream<6> input(fileName, converter);
      vectorSize = input.vectorSize();
      sizeROB += input.size();
    }
  }
  solver.matrixSizeIs(vectorSize, sizeSnap+sizeROB);

  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;
    int colCounter = 0;
    readIntoSolver(solver, converter, BasisId::SNAPSHOTS, domain->solInfo().snapfiPodRom.size(), vectorSize, transform, type, colCounter, fullMass, fullMassSolver, skipTime); // read in snapshots
    readIntoSolver(solver, converter, BasisId::ROB, domain->solInfo().robfi.size(), vectorSize, transform, type, colCounter, fullMass, fullMassSolver); // read in ROB
    
    if(domain->solInfo().robcSolve) solver.solve();

    int orthoBasisDim = domain->solInfo().maxSizePodRom ?
                              std::min(domain->solInfo().maxSizePodRom, solver.singularValueCount()) :
                              solver.singularValueCount();

    // Compute and output the truncation error
    {
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
      fileName += ".truncation_error.txt";
      std::vector<double> toto(orthoBasisDim+1);
      toto[orthoBasisDim] = 0;
      for (int iVec = orthoBasisDim-1; iVec >= 0; --iVec) {
        toto[iVec] = toto[iVec+1]+solver.singularValue(iVec); // running sum
      }
      std::ofstream out(fileName.c_str());
      bool reset = true;
      int truncatedDim;
      for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
        double energy = toto[iVec]/toto[0];
        if(energy < domain->solInfo().romEnergy && reset) {
          truncatedDim = iVec+1;
          reset = false;
        }
        out << iVec+1 << " " << solver.singularValue(iVec) << " " << energy << std::endl;
      }
      if(!reset) orthoBasisDim = truncatedDim;
    }

    // Output solution
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    fileName.append(".orthonormalized");
    BasisOutputStream<6> output(fileName, converter, false);

    if(domain->solInfo().normalize <= 0) // old method for lumped: outputs identity normalized basis
      filePrint(stderr, " ... Writing orthonormal basis of size %d to file %s ...\n", orthoBasisDim, fileName.c_str());
    for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
      output << std::make_pair(solver.singularValue(iVec), solver.matrixCol(iVec));
    }
    if(domain->solInfo().subtractRefPodRom){ // add offset to data
      output << std::make_pair(1.0, centroid->data());
      orthoBasisDim++; 
    }

    // Read back in output file to renormalize basis
    VecBasis basis;
    BasisInputStream<6> in(fileName, converter);
    readVectors(in, basis);
    if(domain->solInfo().subtractRefPodRom) MGSVectors(basis.data(), basis.numVec(), basis.size()); // orthonormalize offset with respect to basis

    VecBasis normalizedBasis;
    if(domain->solInfo().normalize == 0) {
      // Old method: renormalize the orthonormal basis
      renormalized_basis(*fullMass, basis, normalizedBasis);
    }
    else if(domain->solInfo().normalize == 1) {
      // New method: multiply by inverse square root or cholesky factor of the mass matrix
      if(mratio == 0) {
        for(int col = 0; col < orthoBasisDim; col ++ ) {
          fullMass->inverseSquareRootMult(basis[col].data());
        }
      }
      if(mratio != 0) {
        for(int col = 0; col < orthoBasisDim; col++){
          fullMassSolver->backward(basis[col].data());
        }
      }
      normalizedBasis = basis;
    }
  
    // Output the renormalized basis as separate file
    if(domain->solInfo().normalize >= 0) {
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
      fileName.append(".massorthonormalized");
      BasisOutputStream<6> outputNormalized(fileName, converter, false); 
      filePrint(stderr, " ... Writing mass-normalized basis of size %d to file %s ...\n", orthoBasisDim, fileName.c_str());
      for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
        outputNormalized << std::make_pair(solver.singularValue(iVec), normalizedBasis[iVec]);
      }
    }
  
    // Compute and output orthonormal basis if using new method
    if(domain->solInfo().normalize == 1) {
      MGSVectors(normalizedBasis.data(), normalizedBasis.numVec(), normalizedBasis.size());
      BasisOutputStream<6> outputIdentityNormalized(fileName, converter, false); 
      filePrint(stderr, " ... Writing orthonormal basis of size %d to file %s ...\n", orthoBasisDim, fileName.c_str());
      for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
        outputIdentityNormalized << std::make_pair(solver.singularValue(iVec), normalizedBasis[iVec]);
      }
    }

  }
  if(domain->solInfo().subtractRefPodRom) delete centroid;
}

void
BasisOrthoDriver::preProcess() {
  domain->preProcessing();
 
  // Build the constrained DofSetArray incorporating the boundary conditions 
  const int numdof = domain->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain->make_bc(bc.array(), bcx.array());
  domain->make_constrainedDSA(bc.array());
}

} /* end namespace Rom */

Rom::DriverInterface *basisOrthoDriverNew(Domain *domain) {
  return new Rom::BasisOrthoDriver(domain);
}
