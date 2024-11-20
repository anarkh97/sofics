//
// Created by Michel Lesoinne on 2018-12-14.
//

#include <complex>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include "CholmodImp.h"

namespace {
cholmod_common *getCholmodCommon() {
	static cholmod_common c;
	static bool is_initialized = [] {
		cholmod_start(&c);
		return true;
	} ();
	return &c;
}
}

CholmodImp::CholmodImp(const SparseData &upperTriangularStructure, bool isComplex)
	: isComplex(isComplex)
{
	n = upperTriangularStructure.numCol();
	nnz = upperTriangularStructure.nnz();
	auto c = getCholmodCommon();
	A = cholmod_allocate_sparse(
		n,
		n,
		nnz,
		false,
		true,
		1, // upper triangular
		isComplex ? CHOLMOD_COMPLEX : CHOLMOD_REAL,
		c
	);
	// Create the structure of A
	auto colptr = static_cast<int*>(A->p);
	auto rowind = static_cast<int*>(A->i);
	auto val = static_cast<double*>(A->x);

	auto &colPointers = upperTriangularStructure.colPointers();
	auto &rowIndices = upperTriangularStructure.rowIndices();
	const auto basis = colPointers[0]; // 0 or 1 depending on the SparseData implementation.
	for (int i = 0; i < colPointers.size(); ++i)
		colptr[i] = colPointers[i] - basis;
	for (int i = 0; i < rowIndices.size(); ++i)
	rowind[i] = rowIndices[i] - basis;

	L = cholmod_analyze(A, c);
//	cholmod_factorize(A, L, c);
}

CholmodImp::~CholmodImp()
{
	auto c = getCholmodCommon();
	if (A != nullptr)
		cholmod_free_sparse (&A, c);
	if (L != nullptr)
		cholmod_free_factor (&L, c) ; /* free matrices */

}

long CholmodImp::memorySize() const
{
	return 0;
}

void CholmodImp::factorize()
{
	auto c = getCholmodCommon();
	if (!L)
		L = cholmod_analyze(A, c);
	cholmod_factorize(A, L, c);
}

template<typename Scalar>
void CholmodImp::setData(const GenDBSparseMatrix <Scalar> &K)
{
	auto colptr = static_cast<int*>(A->p);
	auto rowind = static_cast<int*>(A->i);
	auto val = static_cast<Scalar*>(A->x);
	auto &colPointers = K.colPointers();
	auto &rowIndices = K.rowIndices();
	const auto basis = colPointers[0]; // 0 or 1 depending on the SparseData implementation.

	const auto &Kvalues = K.values();
	for (size_t i = 0 ; i < Kvalues.size(); ++i)
		val[i] = Kvalues[i];
}

template<typename Scalar>
void CholmodImp::solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution)
{
	auto c = getCholmodCommon();
	if (b == nullptr)
		b = cholmod_allocate_dense(n, 1, n, isComplex ? CHOLMOD_COMPLEX : CHOLMOD_REAL, c) ;
	auto *chol_data = static_cast<Scalar *>(b->x);
	for (int i = 0; i < rhs.size(); ++i)
		chol_data[i] = rhs[i];
	cholmod_dense* res = cholmod_solve(CHOLMOD_A, L, b, c) ;
	chol_data = static_cast<Scalar *>(res->x);
	for (int i = 0; i < rhs.size(); ++i)
		solution[i] = chol_data[i];
	cholmod_free_dense(&res, c) ;
}

template<typename Scalar>
void CholmodImp::reSolve(Scalar *rhs)
{
	auto c = getCholmodCommon();
	if (b == nullptr)
		b = cholmod_allocate_dense(n, 1, n, isComplex ? CHOLMOD_COMPLEX : CHOLMOD_REAL, c) ;
	auto *chol_data = static_cast<Scalar *>(b->x);
	for (int i = 0; i < n; ++i)
		chol_data[i] = rhs[i];
	cholmod_dense* res = cholmod_solve(CHOLMOD_A, L, b, c) ;
	chol_data = static_cast<Scalar *>(res->x);
	for (int i = 0; i < n; ++i)
		rhs[i] = chol_data[i];
	cholmod_free_dense(&res, c) ;
}

template void CholmodImp::setData<double>(const GenDBSparseMatrix <double> &K);
template void CholmodImp::setData(const GenDBSparseMatrix <std::complex<double>> &K);

template void  CholmodImp::solve(const GenVector<double> &rhs, GenVector<double> &solution);
template void  CholmodImp::solve(const GenVector<std::complex<double>> &rhs, GenVector<std::complex<double>> &solution);

template void  CholmodImp::reSolve(double *x);
template void  CholmodImp::reSolve(std::complex<double> *x);
