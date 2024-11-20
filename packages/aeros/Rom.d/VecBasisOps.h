#ifndef ROM_VECBASISOPS_H
#define ROM_VECBASISOPS_H

#include "VecBasis.h"
#include "CholeskyUtils.h"

#include <Math.d/Vector.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>


namespace Rom {

// Returns the matrix-basis product
// (basis and result must refer to different objects)
template <typename Scalar>
const GenVecBasis<Scalar, GenVector> &
mult(const GenSparseMatrix<Scalar> &matrix, const GenVecBasis<Scalar, GenVector> &basis, GenVecBasis<Scalar, GenVector> &result) {
  assert(&basis != &result);

  typedef GenVecBasis<Scalar, GenVector> BasisType;
  typedef typename BasisType::iterator VecIt;

  result.dimensionIs(basis.vectorCount(), basis.vectorInfo());
  
  for (VecIt it = const_cast<BasisType &>(basis).begin(),
             it_end = const_cast<BasisType &>(basis).end(),
             jt = result.begin();
       it != it_end;
       ++it, ++jt) {
    const_cast<GenSparseMatrix<Scalar> &>(matrix).mult(*it, *jt);
  }

  return result;
}

// Returns the matrix^T-basis product
// (basis and result must refer to different objects)
template <typename Scalar>
const GenVecBasis<Scalar, GenVector> &
transposeMult(const GenSparseMatrix<Scalar> &matrix, const GenVecBasis<Scalar, GenVector> &basis, GenVecBasis<Scalar, GenVector> &result) {
  assert(&basis != &result);

  typedef GenVecBasis<Scalar, GenVector> BasisType;
  typedef typename BasisType::iterator VecIt;

  result.dimensionIs(basis.vectorCount(), basis.vectorInfo());

  for (VecIt it = const_cast<BasisType &>(basis).begin(),
             it_end = const_cast<BasisType &>(basis).end(),
             jt = result.begin();
       it != it_end;
       ++it, ++jt) {
    const_cast<GenSparseMatrix<Scalar> &>(matrix).transposeMult(*it, *jt);
  }

  return result;
}

// Returns the renormalized basis Phi with respect to the metric M (assumed symmetric positive semidefinite)
// (Phi, M) -> Phi * R^{-T} where (Phi^T * M * Phi) = R * R^T
// Distributed consistency: 1) requirements: Phi fully consistent, M does not assemble 2) guarantees: result fully consistent
template <typename Scalar>
const GenVecBasis<Scalar, GenVector> &
renormalized_basis(const GenSparseMatrix<Scalar> &metric, const GenVecBasis<Scalar, GenVector> &basis, GenVecBasis<Scalar, GenVector> &result) {
  // result <- M * Phi
  mult(metric, basis, result);

  // Build the lower-triangular part of the normal matrix
  // normalMatrix <- lower( (M * Phi)^T * Phi )
  const int vecCount = result.vectorCount();
  GenFullSquareMatrix<Scalar> normalMatrix(vecCount);

  for (int row = 0; row < vecCount; ++row) {
    const GenVector<Scalar> &dual = result[row];
    for (int col = 0; col <= row; ++col) {
      // sum-consistent * fully consistent
      normalMatrix[row][col] = dual * basis[col];
    }
  }

  cholesky_factor_lower(normalMatrix); // normalMatrix <- R where A = R * R^T
  inverse_triangular_lower(normalMatrix); // normalMatrix <- R^{-1}

  // result <- Phi * R^{-T}
  for (int row = 0; row < vecCount; ++row) {
    GenVector<Scalar> &target = result[row];
    int col = 0;
    target.linC(basis[col], normalMatrix[row][col]);
    for (col = 1; col <= row; ++col) {
      target.linAdd(normalMatrix[row][col], basis[col]);
    }
  }

  return result;
}

// Calculates the reduced stiffness matrix  K_red = Phi^T * K * Phi with Phi as the mass-normalized basis
template <typename Scalar>
void calculateReducedStiffness(const GenSparseMatrix<Scalar> &K, const GenVecBasis<Scalar, GenVector> &basis, GenFullSquareMatrix<Scalar> &K_reduced) {
  // K^T * Phi
  VecBasis product; // used as a buffer for intermediate steps
  transposeMult(K, basis, product);
  // calculate transpose of product multiplied with basis  (K^T * Phi)^T * Phi = Phi^T * K * Phi
  const int vecCount = product.vectorCount();

  GenFullSquareMatrix<Scalar> normalMatrix(vecCount);
  for(int i = 0; i < vecCount; i++){
    for(int j = 0 ; j < vecCount; j++){
      normalMatrix[i][j] = product[i] * basis[j];
    }
  }
  K_reduced.copy(normalMatrix);
}

// Modified Gram-Schmidt Algorithm
template <typename Scalar>
void
MGSVectors(Scalar *d, int numVec, int lengthVec, bool RowMajor = false) {
#ifdef USE_EIGEN3
  filePrint(stderr," ... Gram-Schmidt Algorithm: orthogonalizing vectors ...\n");
  // initialize eigen matrix class with pointer to vectors
  Eigen::Map< Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> > matrix(NULL,0,0);
  if (RowMajor)
    new (&matrix) Eigen::Map< Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> >(d,lengthVec,numVec);
  else
    new (&matrix) Eigen::Map< Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::ColMajor> >(d,lengthVec,numVec);

  // loop over all vectors
  for(int i = 0; i != numVec; ++i) {  
    filePrint(stderr,"\r %5.2f%% complete", double(i)/double(numVec)*100.); 

    // normalize the current vector
    matrix.col(i).normalize();

    // loop over all other vectors and subtract off projection 
    for(int j = i+1; j != numVec; ++j) { 
      Scalar vecProj = matrix.col(i).dot(matrix.col(j));
      matrix.col(j) -= vecProj*matrix.col(i);
    }
  }
  filePrint(stderr,"\r %5.2f%% complete\n", 100.);
#else
  filePrint(stderr, " *** ERROR: MGSVectors requires Eigen 3 library\n");
  exit(-1);
#endif
}

template <typename Scalar>
void
PrintData(Scalar *d, int numVec, int lengthVec, bool RowMajor = false) {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> > matrix(NULL,0,0);
  if (RowMajor)
    new (&matrix) Eigen::Map< Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> >(d,lengthVec,numVec);
  else
    new (&matrix) Eigen::Map< Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::ColMajor> >(d,lengthVec,numVec);

  std::cout << matrix.transpose()*matrix << std::endl;
#endif
}


} // end namespace Rom

#endif /* ROM_VECBASISOPS_H */
