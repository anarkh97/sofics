#ifdef USE_EIGEN3
#include <Math.d/EiSparseMatrix.h>

#ifdef EIGEN_UMFPACK_SUPPORT
template<>
void
GenEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > >::factor()
{
  // workaround issue with umfpack and mapped sparse matrix
  if(M_copy) delete M_copy;
  M_copy = new Eigen::SparseMatrix<double>(M);
  solver.compute(*M_copy);
  if(solver.info() != Eigen::Success) std::cerr << "sparse factor failed\n";
}

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<complex<double> > > >::factor()
{
  // workaround issue with umfpack and mapped sparse matrix
  if(M_copy) delete M_copy;
  M_copy = new Eigen::SparseMatrix<complex<double> >(M);
  solver.compute(*M_copy);
  if(solver.info() != Eigen::Success) std::cerr << "sparse factor failed\n";
}

template<>
void
GenEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > >::reSolve(double* _rhs)
{
  // umfpack does not support in-place solve
  solveTime -= getTime();
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> sol(numUncon);
  sol = solver.solve(rhs);
  rhs = sol;
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
  solveTime += getTime();
}

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<complex<double> > > >::reSolve(complex<double>* _rhs)
{
  // umfpack does not support in-place solve
  solveTime -= getTime();
  Eigen::Map< Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> sol(numUncon);
  sol = solver.solve(rhs);
  rhs = sol;
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
  solveTime += getTime();
}

template<>
void
GenEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > >::reSolve(GenVector<double> &_rhs)
{
  // umfpack does not support in-place solve
  solveTime -= getTime();
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> sol(numUncon);
  sol = solver.solve(rhs);
  rhs = sol;
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
  solveTime += getTime();
}

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<complex<double> > > >::reSolve(GenVector<complex<double> > &_rhs)
{
  // umfpack does not support in-place solve
  solveTime -= getTime();
  Eigen::Map< Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1);
  Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> sol(numUncon);
  sol = solver.solve(rhs);
  rhs = sol;
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
  solveTime += getTime();
}
#endif

template<>
void
GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >::upperMult(double* _rhs)
{
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  if(solver.permutationP().size() > 0) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> prhs = solver.permutationP()*rhs;
    rhs = solver.matrixU()*prhs ;
  } 
  else 
    rhs = (solver.matrixU()*rhs).eval();
}

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::SimplicialLLT<Eigen::SparseMatrix<complex<double> >,Eigen::Upper> >::upperMult(complex<double>* _rhs)
{
  Eigen::Map< Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  if(solver.permutationP().size() > 0) {
    Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> prhs = solver.permutationP()*rhs;
    rhs = solver.matrixU()*prhs ;
  }
  else
    rhs = (solver.matrixU()*rhs).eval();
}

template<>
void
GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >::backward(double* _rhs)
{
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  solver.matrixU().solveInPlace(rhs);
  if(solver.permutationP().size() > 0) rhs = (solver.permutationPinv()*rhs).eval();
}

template<>
void
GenEiSparseMatrix<complex<double>,Eigen::SimplicialLLT<Eigen::SparseMatrix<complex<double> >,Eigen::Upper> >::backward(complex<double>* _rhs)
{
  Eigen::Map< Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  solver.matrixU().solveInPlace(rhs);
  if(solver.permutationP().size() > 0) rhs = (solver.permutationPinv()*rhs).eval();
}
#endif
