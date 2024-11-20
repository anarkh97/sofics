#ifndef ROM_DISTRVECBASISOPS_H
#define ROM_DISTRVECBASISOPS_H

#include "DistrVecBasis.h"
#include "CholeskyUtils.h"

#include <Feti.d/DistrVector.h>
#include <Paral.d/SubDOp.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/DistHelper.h>

#include <cassert>

namespace Rom {

// Computes the transposed matrix-vector product using the masterflag of the basis
template <typename Scalar, template <typename> class GenVecType>
const GenVecType<Scalar> &
vector_components(const GenVecBasis<Scalar, GenDistrVector> &basis,
                  const GenDistrVector<Scalar> &vec,
                  GenVecType<Scalar> &result) {
  assert(basis.vectorSize() == vec.size());
  assert(basis.vectorCount() == result.size());

  for (int iVec = 0; iVec < basis.vectorCount(); ++iVec) {
    typedef GenDistrVector<Scalar> DVec;
    result[iVec] = const_cast<DVec &>(basis[iVec]) * const_cast<DVec &>(vec);
  }

  return result;
}

// Computes the transposed matrix-vector product using the masterflag of the vector
template <typename Scalar, template <typename> class GenVecType>
const GenVecType<Scalar> &
vector_components_vector_masterflag(const GenVecBasis<Scalar, GenDistrVector> &basis,
                                    const GenDistrVector<Scalar> &vec,
                                    GenVecType<Scalar> &result) {
  assert(basis.vectorSize() == vec.size());
  assert(basis.vectorCount() == result.size());

  for (int iVec = 0; iVec < basis.vectorCount(); ++iVec) {
    typedef GenDistrVector<Scalar> DVec;
    result[iVec] = const_cast<DVec &>(vec) * const_cast<DVec &>(basis[iVec]);
  }

  return result;
}

// Computes the matrix-vector product
template <typename Scalar, template <typename> class GenVecType>
const GenDistrVector<Scalar> &
assembled_vector(const GenVecBasis<Scalar, GenDistrVector> &basis,
                 const GenVecType<Scalar> &components,
                 GenDistrVector<Scalar> &result) {
  assert(basis.vectorSize() == result.size());
  assert(basis.vectorCount() == components.size());
  result.zero();

  for (int iVec = 0; iVec < basis.vectorCount(); ++iVec) {
    GenDistrVector<Scalar> &vec = const_cast<GenDistrVector<Scalar> &>(basis[iVec]);
    result.linAdd(components[iVec], vec);
  }

  return result;
}


// Returns the matrix-basis product
// (basis and result must refer to different objects)
template <typename Scalar>
const GenVecBasis<Scalar, GenDistrVector> &
mult(const GenSubDOp<Scalar> &matrix, const GenVecBasis<Scalar, GenDistrVector> &basis, GenVecBasis<Scalar, GenDistrVector> &result) {
  assert(&basis != &result);

  typedef GenVecBasis<Scalar, GenDistrVector> BasisType;
  typedef typename BasisType::iterator VecIt;

  result.dimensionIs(basis.vectorCount(), basis.vectorInfo());
  
  for (VecIt it = const_cast<BasisType &>(basis).begin(),
             it_end = const_cast<BasisType &>(basis).end(),
             jt = result.begin();
       it != it_end;
       ++it, ++jt) {
    const_cast<GenSubDOp<Scalar> &>(matrix).mult(*it, *jt);
  }

  return result;
}

// Returns the matrix^T-basis product
// (basis and result must refer to different objects)
template <typename Scalar>
const GenVecBasis<Scalar, GenDistrVector> &
transposeMult(const GenSubDOp<Scalar> &matrix, const GenVecBasis<Scalar, GenDistrVector> &basis, GenVecBasis<Scalar, GenDistrVector> &result) {
  assert(&basis != &result);

  typedef GenVecBasis<Scalar, GenDistrVector> BasisType;
  typedef typename BasisType::iterator VecIt;

  result.dimensionIs(basis.vectorCount(), basis.vectorInfo());
  
  for (VecIt     it = const_cast<BasisType &>(basis).begin(),
             it_end = const_cast<BasisType &>(basis).end(),
                 jt = result.begin();
       it != it_end;
       ++it, ++jt) {
    const_cast<GenSubDOp<Scalar> &>(matrix).transposeMult(*it, *jt);
  }

  return result;
}

// Returns the renormalized basis Phi with respect to the metric M (assumed symmetric positive semidefinite)
// (Phi, M) -> Phi * R^{-T} where (Phi^T * M * Phi) = R * R^T
// Distributed consistency: 1) requirements: Phi fully consistent, M does not assemble 2) guarantees: result fully consistent
template <typename Scalar>
const GenVecBasis<Scalar, GenDistrVector> &
renormalized_basis(const GenSubDOp<Scalar> &metric, const GenVecBasis<Scalar, GenDistrVector> &basis, GenVecBasis<Scalar, GenDistrVector> &result) {
  // result <- M * Phi
  mult(metric, basis, result);

  // Build the lower-triangular part of the normal matrix
  // normalMatrix <- lower( (M * Phi)^T * Phi )
  const int vecCount = result.vectorCount();
  GenFullSquareMatrix<Scalar> normalMatrix(vecCount);
  for (int row = 0; row < vecCount; ++row) {
    const GenDistrVector<Scalar> &dual = result[row];
    for (int col = 0; col <= row; ++col) {
      // sum-consistent * fully consistent
      normalMatrix[row][col] = dot_ignore_master_flag(dual, basis[col]);
    }
  }
  cholesky_factor_lower(normalMatrix); // normalMatrix <- R where A = R * R^T
  inverse_triangular_lower(normalMatrix); // normalMatrix <- R^{-1}

  // result <- Phi * R^{-T}
  GenVecBasis<Scalar, GenDistrVector> &basis_fix = const_cast<GenVecBasis<Scalar, GenDistrVector> &>(basis);
  for (int row = 0; row < vecCount; ++row) {
    GenDistrVector<Scalar> &target = result[row];
    int col = 0;
    target.linC(basis_fix[col], normalMatrix[row][col]);
    for (col = 1; col <= row; ++col) {
      target.linAdd(normalMatrix[row][col], basis_fix[col]);
    }
  }

  return result;
}

// Calculates the reduced stiffness matrix  K_red = Phi^T * K * Phi with Phi as the mass-normalized basis
template <typename Scalar>
void calculateReducedStiffness(const GenSubDOp<Scalar> &K, const GenVecBasis<Scalar, GenDistrVector> &basis, GenFullSquareMatrix<Scalar> &K_reduced, bool sym = false) {
  // K^T * Phi
  DistrVecBasis product; // used as a buffer for intermediate steps
  transposeMult(K, basis, product);
  // calculate transpose of product multiplied with basis  (K^T * Phi)^T * Phi = Phi^T * K * Phi
  const int vecCount = product.vectorCount();

  GenFullSquareMatrix<Scalar> normalMatrix(vecCount); normalMatrix.zero();
  for(int i = 0; i < vecCount; i++){
    if(sym) { // if symmetric, compute only upper half
      for(int j = i ; j < vecCount; j++){
        normalMatrix[i][j] = dot_ignore_master_flag(product[i],basis[j]);
      }
    } else { // if not, compute who matrix
      for(int j = 0 ; j < vecCount; j++){
        normalMatrix[i][j] = dot_ignore_master_flag(product[i],basis[j]);
      }
    }
  }
  K_reduced.copy(normalMatrix);
}

template <typename Scalar>
const GenVecBasis<Scalar, GenDistrVector> &
MGSVectors(const GenVecBasis<Scalar, GenDistrVector> &basis) {
  filePrint(stderr," ... Distributed Modified Gram-Schmidt: orthogonalizing vectors ...\n");

  int numVectors = basis.numVec();

  for (int i = 0; i < numVectors; i++) {
    filePrint(stderr,"\r %5.2f%% complete", double(i)/double(numVectors)*100.);

    GenDistrVector<Scalar> v(basis.vectorInfo(),basis[i].data(),false);
    Scalar normCol = v.norm();
    
    v /= normCol;
    for(int j = i+1; j < numVectors; ++j) {
      GenDistrVector<Scalar> q(basis.vectorInfo(),basis[j].data(),false);
      Scalar vecProj = v*q;
      q -= vecProj*v;
    }
  }
  filePrint(stderr,"\r %5.2f%% complete\n", 100.);
  return basis;
}

} // end namespace Rom

#endif /* ROM_DISTRVECBASISOPS_H */
