#include "CholeskyUtils.h"

#include <Math.d/FullSquareMatrix.h>
//#include <Math.d/PseudoInverse.h>
//#include <Eigen/Dense>

#include <Utils.d/linkfc.h>

#include <complex>
#include <stdexcept>

extern "C" {
  // Cholesky factorization
  void _FORTRAN(dpotrf)(const char *uplo, const int *n,
                        double *a, const int *lda, int *info);
 
  // Solve with existing Cholesky factorization
  void _FORTRAN(dpotrs)(const char *uplo, const int *n, const int *nhrs,
                        const double *a, const int *lda,
                        double *b, const int *ldb, int *info);
 
  // Inverse triangular matrix
  void _FORTRAN(dtrtri)(const char *uplo, const char *diag, const int *n,
                        double *a, const int *lda, int *info);

  // LDLT factorization for symmetric indefinite matrix
  void _FORTRAN(dsytrf)(const char *uplo, const int *n,
                        double *a, const int *lda, int *ipiv, double *work, 
                        int *lwork, int *info);

  // Solve with existing LDLT factorization
  void _FORTRAN(dsytrs)(const char *uplo, const int *n, const int *nhrs,
                        const double *a, const int *lda, const int *ipiv,
                        double *b, const int *ldb, int *info);
}

namespace Rom {

template <>
const GenFullSquareMatrix<double> &
cholesky_factor_upper(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> m_copy = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(m.data(),basisDim,basisDim);

  int info;
  _FORTRAN(dpotrf)("L", &basisDim, m.data(), &basisDim, &info);

  //Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dec(m_copy);
  //std::cerr << "here is the matrix m:\n" << m_copy << std::endl;
  //std::cerr << "here are the eigenvalues of m:\n" << dec.eigenvalues() << std::endl;

  if (info) {
    throw std::runtime_error("Error in dpotrf (Cholesky factorization)");
  };

  return m;
}

template <>
const GenFullSquareMatrix<double> &
cholesky_factor_lower(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  int info;
  _FORTRAN(dpotrf)("U", &basisDim, m.data(), &basisDim, &info);

  if (info) {
    throw std::runtime_error("Error in dpotrf (Cholesky factorization)");
  };

  return m;
}

template <>
const double *
cholesky_solve_upper(const GenFullSquareMatrix<double> &m, double *v) {
  const int basisDim = m.dim();
  const int INT_ONE = 1;
  
  int info;
  _FORTRAN(dpotrs)("L", &basisDim, &INT_ONE, m.data(), &basisDim, v, &basisDim, &info);
  
  if (info) {
    throw std::runtime_error("Error in dpotrs (Cholesky solver)");
  };

  return v;
}

template <>
const double *
cholesky_solve_lower(const GenFullSquareMatrix<double> &m, double *v) {
  const int basisDim = m.dim();
  const int INT_ONE = 1;
  
  int info;
  _FORTRAN(dpotrs)("U", &basisDim, &INT_ONE, m.data(), &basisDim, v, &basisDim, &info);
  
  if (info) {
    throw std::runtime_error("Error in dpotrs (Cholesky solver)");
  };
  
  return v;
}

template <>
const GenFullSquareMatrix<double> &
inverse_triangular_upper(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  int info;
  _FORTRAN(dtrtri)("L", "N", &basisDim, m.data(), &basisDim, &info);

  if (info) {
    throw std::runtime_error("Error in dtrtri (Triangular matrix inversion)");
  };

  return m;
}

template <>
const GenFullSquareMatrix<double> &
inverse_triangular_lower(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  int info;
  _FORTRAN(dtrtri)("U", "N", &basisDim, m.data(), &basisDim, &info);

  if (info) {
    throw std::runtime_error("Error in dtrtri (Triangular matrix inversion)");
  };

  return m;
}

template <>
const GenFullSquareMatrix<double> &
ldlt_factor_upper(GenFullSquareMatrix<double> &m, int* ipiv) {
  const int basisDim = m.dim();

  int lwork = 100*basisDim;
  double *work = new double[lwork];
  int info;
  _FORTRAN(dsytrf)("L", &basisDim, m.data(), &basisDim, ipiv, work, &lwork, &info);
  delete [] work;

  if (info) {
    throw std::runtime_error("Error in dsytrf (LDLT factorization)");
  };

  return m;
}

template <>
const GenFullSquareMatrix<double> &
ldlt_factor_lower(GenFullSquareMatrix<double> &m, int *ipiv) {
  const int basisDim = m.dim();

  int lwork = 100*basisDim;
  double *work = new double[lwork];
  int info;
  _FORTRAN(dsytrf)("U", &basisDim, m.data(), &basisDim, ipiv, work, &lwork, &info);
  delete [] work;

  if (info) {
    throw std::runtime_error("Error in dsytrf (LDLT factorization)");
  };

  return m;
}

template <>
const double *
ldlt_solve_upper(const GenFullSquareMatrix<double> &m, double *v, const int *ipiv) {
  const int basisDim = m.dim();
  const int INT_ONE = 1;
  
  int info;
  _FORTRAN(dsytrs)("L", &basisDim, &INT_ONE, m.data(), &basisDim, ipiv, v, &basisDim, &info);
  
  if (info) {
    throw std::runtime_error("Error in dsytrs (LDLT solver)");
  };
  
  return v;
}

template <>
const double *
ldlt_solve_lower(const GenFullSquareMatrix<double> &m, double *v, const int *ipiv) {
  const int basisDim = m.dim();
  const int INT_ONE = 1;
  
  int info;
  _FORTRAN(dsytrs)("U", &basisDim, &INT_ONE, m.data(), &basisDim, ipiv, v, &basisDim, &info);
  
  if (info) {
    throw std::runtime_error("Error in dsytrs (LDLT solver)");
  };
  
  return v;
}

//-------------------------------------------------------------------

template <>
const GenFullSquareMatrix<complex<double> > &
cholesky_factor_upper(GenFullSquareMatrix<complex<double> > &) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<complex<double> > &
cholesky_factor_lower(GenFullSquareMatrix<complex<double> > &) {
  throw std::logic_error("Not implemented");
}

template <>
const complex<double> *
cholesky_solve_upper(const GenFullSquareMatrix<complex<double> > &, complex<double> *) {
  throw std::logic_error("Not implemented");
}

template <>
const complex<double> *
cholesky_solve_lower(const GenFullSquareMatrix<complex<double> > &, complex<double> *) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<complex<double> > &
inverse_triangular_upper(GenFullSquareMatrix<complex<double> > &) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<complex<double> > &
inverse_triangular_lower(GenFullSquareMatrix<complex<double> > &) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<complex<double> > &
ldlt_factor_upper(GenFullSquareMatrix<complex<double> > &, int*) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<complex<double> > &
ldlt_factor_lower(GenFullSquareMatrix<complex<double> > &, int*) {
  throw std::logic_error("Not implemented");
}

template <>
const complex<double> *
ldlt_solve_upper(const GenFullSquareMatrix<complex<double> > &, complex<double> *, const int*) {
  throw std::logic_error("Not implemented");
}

template <>
const complex<double> *
ldlt_solve_lower(const GenFullSquareMatrix<complex<double> > &, complex<double> *, const int*) {
  throw std::logic_error("Not implemented");
}

} /* end namespace Rom */
