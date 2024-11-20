//
// Created by Michel Lesoinne on 2018-12-14.
//

#ifndef FEM_CHOLMODIMP_H
#define FEM_CHOLMODIMP_H

#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>

#include "cholmod.h"

class Connectivity;
class EqNumberer;

class CholmodImp {
public:
	CholmodImp(const SparseData &upperTriangularStructure, bool isComplex);
	~CholmodImp();

	template <typename Scalar>
	void setData(const GenDBSparseMatrix<Scalar> &K);
	void factorize();

	template <typename Scalar>
	void solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution);

	template <typename Scalar>
	void reSolve(Scalar *rhs);

	long memorySize() const;
private:
	/// \brief The sparse matrix structure for Cholmod. Always allocated.
	cholmod_sparse *A;
	cholmod_factor *L = nullptr;
	cholmod_dense* b = nullptr;
	size_t n;
	size_t nnz;
	bool isComplex;

};


#endif //FEM_CHOLMODIMP_H
