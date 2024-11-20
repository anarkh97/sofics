#include <iostream>
#include <algorithm>
#include <complex>
#include <Math.d/Vector.h>
#include <Driver.d/Communicator.h>
#include <unsupported/Eigen/SparseExtra>
#include <Timers.d/GetTime.h>

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::GenEiSparseMatrix(const Connectivity *cn, const DofSetArray *dsa, const int *rCN, bool _selfadjoint)
: SparseData(dsa,rCN,cn,int(!_selfadjoint)),
  selfadjoint(_selfadjoint),
  nnz(xunonz[numUncon]-1),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz.data(), rowu.data(), unonz),
  M_copy(NULL)
{
  for(int k=0; k < numUncon; k++)
	std::sort(rowu.begin() + xunonz[k]-1, rowu.begin() + xunonz[k+1]-1);
 
  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::GenEiSparseMatrix(const Connectivity *cn, const DofSetArray *dsa,
                                                         const DofSetArray *c_dsa, bool _selfadjoint)
: SparseData(dsa,c_dsa,cn,int(!_selfadjoint)),
  selfadjoint(_selfadjoint),
  nnz(xunonz[numUncon]-1),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz.data(), rowu.data(), unonz),
  M_copy(NULL)
{
  for(int k=0; k < numUncon; k++)
	std::sort(rowu.begin() + xunonz[k]-1, rowu.begin() + xunonz[k+1]-1);

  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::GenEiSparseMatrix(const Connectivity *cn, const EqNumberer *eqNums, bool _selfadjoint)
: SparseData(eqNums,cn,(int*)NULL,0,1),
  selfadjoint(_selfadjoint),
  nnz(xunonz[numUncon]-1),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz.data(), rowu.data(), unonz),
  M_copy(NULL)
{
  for(int k=0; k < numUncon; k++)
	std::sort(rowu.data() + xunonz[k]-1, rowu.data() + xunonz[k+1]-1);

  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::~GenEiSparseMatrix()
{
	delete [] unonz;
	delete M_copy;
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::zeroAll()
{
  for(int i=0; i<nnz; ++i) unonz[i] = 0.;
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::print()
{
  std::cerr << std::setprecision(15) << M << std::endl;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::printSparse(const std::string& filename)
{
  // export to ascii matrix market format
  int sym = (selfadjoint) ? Eigen::Symmetric : 0x30;
  saveMarket(M, filename, sym);
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::add(const FullSquareMatrix &kel, const int *dofs)
{
  int k,l;
  for(int i = 0; i < kel.dim(); ++i) {
	if(dofs[i] < 0 || (k = unconstrNum[dofs[i]]) < 0) continue; // Skip undefined/constrained dofs
	for(int j = 0; j < kel.dim(); ++j) {
	  if(selfadjoint && dofs[i] > dofs[j]) continue; // Work with upper symmetric half
	  if(dofs[j] < 0 || (l = unconstrNum[dofs[j]]) < 0) continue;  // Skip undefined/constrained dofs
	  for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
		if(rowu[m] == k) {
		  unonz[m] += kel[i][j];
		  break;
		}
	  }
	}
  }
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::addCoef(int k, int l, Scalar val)
{
  for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
	if(rowu[m] == k) {
	  unonz[m] += val;
	  break;
	}
  }
}

template<class Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::add(const GenAssembledFullM<Scalar> &kel, const int *dofs)
{
  // this function is used to assemble Kcc and requires dofs to be in constrained numbering
  int i, j, m, mstart, mstop, ri, rj;
  for(i = 0; i < kel.numRow(); ++i) {       // Loop over rows.
	if((ri = dofs[i]) == -1) continue;      // Skip constrained dofs
	for(j = 0; j < kel.numCol(); ++j) {     // Loop over columns.
	  if((rj = dofs[j]) == -1) continue;    // Skip constrained dofs
	  if(rj < ri) continue;
	  mstart = xunonz[rj];
	  mstop  = xunonz[rj+1];
	  for(m = mstart; m < mstop; ++m) {
		if(rowu[m] == ri) {
		  unonz[m] += kel[i][j];
		  break;
		}
	  }
	}
  }
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  int k,l;
  for(int i = 0; i < kel.dim(); ++i) {
	if(dofs[i] < 0 || (k = unconstrNum[dofs[i]]) < 0) continue; // Skip undefined/constrained dofs
	for(int j = 0; j < kel.dim(); ++j) {
	  if(selfadjoint && dofs[i] > dofs[j]) continue; // Work with upper symmetric half
	  if(dofs[j] < 0 || (l = unconstrNum[dofs[j]]) < 0) continue;  // Skip undefined/constrained dofs
	  for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
		if(rowu[m] == k) {
		  ScalarTypes::addComplex(unonz[m], complex<double>(0,kel[i][j]));
		  break;
		}
	  }
	}
  }
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::add(const FullSquareMatrixC &kel, const int *dofs)
{
  int k,l;
  for(int i = 0; i < kel.dim(); ++i) {
	if(dofs[i] < 0 || (k = unconstrNum[dofs[i]]) < 0) continue; // Skip undefined/constrained dofs
	for(int j = 0; j < kel.dim(); ++j) {
	  if(selfadjoint && dofs[i] > dofs[j]) continue; // Work with upper symmetric half
	  if(dofs[j] < 0 || (l = unconstrNum[dofs[j]]) < 0) continue;  // Skip undefined/constrained dofs
	  for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
		if(rowu[m] == k) {
		  ScalarTypes::addComplex(unonz[m], kel[i][j]);
		  break;
		}
	  }
	}
  }
}

template<typename Scalar, typename SolverClass> 
double
GenEiSparseMatrix<Scalar,SolverClass>::getMemoryUsed() const
{
  return sizeof(Scalar)*nnz/(1024.0*1024.0);
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::mult(const Scalar *_rhs, Scalar *_result) const
{
  Eigen::Map< const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > result(_result,numUncon,1);
  if(selfadjoint)
	result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
	result = M*rhs;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::mult(const GenVector<Scalar> &_rhs, Scalar *_result) const
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result,numUncon,1);
  if(selfadjoint)
	result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
	result = M*rhs;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::transposeMult(const GenVector<Scalar> &_rhs, GenVector<Scalar> &_result) const
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result.data(),numUncon,1);
  if(selfadjoint)
	result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
	result = M.adjoint()*rhs;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::transposeMult(const Scalar *_rhs, Scalar *_result) const
{
  Eigen::Map< const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > result(_result,numUncon,1);
  if(selfadjoint)
	result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
	result = M.adjoint()*rhs;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::multAdd(const Scalar *_rhs, Scalar *_result) const
{
  Eigen::Map< const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > result(_result,numUncon,1);
  if(selfadjoint)
	result += M.template selfadjointView<Eigen::Upper>()*rhs;
  else
	result += M*rhs;
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::mult(const GenVector<Scalar> &_rhs, GenVector<Scalar> &_result) const
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result.data(),numUncon,1);
  if(selfadjoint)
	result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
	result = M*rhs;
}

template<typename Scalar, typename SolverClass> 
Scalar
GenEiSparseMatrix<Scalar,SolverClass>::diag(int dof) const
{
  std::cerr << "GenEiSparseMatrix<Scalar,SolverClass>::diag is not implemented\n";
  return Scalar();
}

template<typename Scalar, typename SolverClass> 
Scalar &
GenEiSparseMatrix<Scalar,SolverClass>::diag( int dof )
{
  static Scalar defaultValue;
  std::cerr << "GenEiSparseMatrix<Scalar,SolverClass>::diag is not implemented\n";
  return defaultValue;
}

template<typename Scalar, typename SolverClass> 
long
GenEiSparseMatrix<Scalar,SolverClass>::size() const
{
  return (numUncon) ? nnz : 0;
}

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
  communicator->globalSum(nnz, unonz);
#endif
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::factor()
{
  solver.compute(M);
  if(solver.info() != Eigen::Success) std::cerr << "sparse factor failed\n";
}

#ifdef EIGEN_UMFPACK_SUPPORT
template<>
void
GenEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > >::factor();

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<complex<double> > > >::factor();
#endif

template<typename Scalar, typename SolverClass>
void 
GenEiSparseMatrix<Scalar,SolverClass>::solve(const Scalar *_rhs, Scalar *_solution)
{
	this->solveTime -= getTime();
	Eigen::Map< const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
	Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > solution(_solution,numUncon,1);
	solution = solver.solve(rhs);
	if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
	this->solveTime += getTime();
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::solve(const GenVector<Scalar> &_rhs, GenVector<Scalar> &_solution)
{
	this->solveTime -= getTime();
	Eigen::Map< const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1);
	Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > solution(_solution.data(),numUncon,1);
	solution = solver.solve(rhs);
	if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
	this->solveTime += getTime();
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::reSolve(Scalar *_rhs)
{
  this->solveTime -= getTime();
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  rhs = solver.solve(rhs);
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
  this->solveTime += getTime();
}

#ifdef EIGEN_UMFPACK_SUPPORT
template<>
void
GenEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > >::reSolve(double* _rhs);

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<complex<double> > > >::reSolve(complex<double>* _rhs);
#endif

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::reSolve(GenVector<Scalar> &_rhs)
{
  this->solveTime -= getTime();
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1);
  rhs = solver.solve(rhs);
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
  this->solveTime += getTime();
}

#ifdef EIGEN_UMFPACK_SUPPORT
template<>
void
GenEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > >::reSolve(GenVector<double> &_rhs);

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<complex<double> > > >::reSolve(GenVector<complex<double> > &_rhs);
#endif

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::upperMult(Scalar* _rhs)
{
  std::cerr << " *** ERROR: GenEiSparseMatrix::upperMult is not implemented\n";
  exit(-1);
}

template<>
void
GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >::upperMult(double* _rhs);

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::SimplicialLLT<Eigen::SparseMatrix<complex<double> >,Eigen::Upper> >::upperMult(complex<double>* _rhs);

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::backward(Scalar* _rhs)
{
  std::cerr << " *** ERROR: GenEiSparseMatrix::backward is not implemented\n";
  exit(-1);
}

template<>
void
GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >::backward(double* _rhs);

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::SimplicialLLT<Eigen::SparseMatrix<complex<double> >,Eigen::Upper> >::backward(complex<double>* _rhs);
