//
// Created by Michel Lesoinne on 2018-12-14.
//

#include <complex>
#include <memory>

#include "CholmodSolver.h"
#include <Math.d/DBSparseMatrix.h>


#ifndef WITH_CHOLMOD
template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getCholmod(const Connectivity *cn, const EqNumberer *_dsa)
{
	std::cerr << "Cholmod is not available in this version of the code. Crash will ensue." << std::endl;
	return {nullptr, nullptr};
}

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getCholmod(const Connectivity *cn, const DofSetArray &dsa, const ConstrainedDSA &cdsa)
{
	std::cerr << "Cholmod is not available in this version of the code. Crash will ensue." << std::endl;
	return {nullptr, nullptr};
}

#else // WITH_CHOLMOD
#include "CholmodImp.h"

class CholmodImp;

template <typename Scalar>
class CholmodSolver :
	public GenDBSparseMatrix<Scalar>, public GenSolver<Scalar> {
public:
	CholmodSolver(const Connectivity *cn, const EqNumberer *_dsa);
	CholmodSolver(const Connectivity *cn, const DofSetArray &dsa, const ConstrainedDSA &cdsa);

	~CholmodSolver() = default;

	int neqs() const override;

	void unify(FSCommunicator *communicator) override {
		GenDBSparseMatrix<Scalar>::unify(communicator);
	}

	void solve(const Scalar *rhs, Scalar *solution) override;

	void reSolve(Scalar *rhs) override;

	long size() const override;

	void factor() override;

	void parallelFactor() override;
private:
	CholmodImp &getImplementation();

	std::unique_ptr<CholmodImp> impl;
};
template <typename Scalar>
CholmodSolver<Scalar>::CholmodSolver(const Connectivity *cn, const EqNumberer *_dsa)  :
	GenDBSparseMatrix<Scalar>(cn, _dsa)
{
}

template <typename Scalar>
CholmodSolver<Scalar>::CholmodSolver(const Connectivity *cn, const DofSetArray &dsa, const ConstrainedDSA &cdsa)  :
	GenDBSparseMatrix<Scalar>(cn, &dsa, &cdsa)
{
}


template<typename Scalar>
int CholmodSolver<Scalar>::neqs() const
{
	return GenDBSparseMatrix<Scalar>::neqs();
}

template<typename Scalar>
long CholmodSolver<Scalar>::size() const
{
	return GenDBSparseMatrix<Scalar>::size()
	       + (impl ? impl->memorySize() : 0);
}

template<typename Scalar>
void CholmodSolver<Scalar>::solve(const Scalar *rhs, Scalar *solution)
{
	GenStackVector<Scalar> b(const_cast<Scalar *>(rhs), neqs()), x(solution, neqs());
	impl->solve(b, x);
}

template<typename Scalar>
void CholmodSolver<Scalar>::factor()
{
	auto &cholmodImp = getImplementation();
	cholmodImp.setData(*this);
	cholmodImp.factorize();
}


template<typename Scalar>
void CholmodSolver<Scalar>::parallelFactor()
{
	factor(); // TODO make parallel.
}

template<typename Scalar>
CholmodImp &CholmodSolver<Scalar>::getImplementation()
{
	if ( !impl )
		impl = std::make_unique<CholmodImp>(*this, std::is_same<std::complex<double>, Scalar>::value);
	return *impl;
}

template<typename Scalar>
void CholmodSolver<Scalar>::reSolve(Scalar *rhs)
{
	impl->reSolve(rhs);
}

template class CholmodSolver<double>;
template class CholmodSolver<std::complex<double>>;

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getCholmod(const Connectivity *cn, const EqNumberer *dsa)
{
	auto solver = new CholmodSolver<Scalar>(cn, dsa);
	return { solver, solver};
}

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getCholmod(const Connectivity *cn, const DofSetArray &dsa, const ConstrainedDSA &cdsa)
{
	auto solver = new CholmodSolver<Scalar>(cn, dsa, cdsa);
	return { solver, solver};
}

#endif // WITH_CHOLMOD

template std::pair<GenSolver<double> *, GenSparseMatrix<double> *>
getCholmod<double>(const Connectivity *cn, const EqNumberer *_dsa);
template std::pair<GenSolver<std::complex<double>> *, GenSparseMatrix<std::complex<double>> *>
getCholmod<std::complex<double>>(const Connectivity *cn, const EqNumberer *_dsa);


template std::pair<GenSolver<double> *, GenSparseMatrix<double> *>
getCholmod<double>(const Connectivity *cn, const DofSetArray &dsa, const ConstrainedDSA &cdsa);

template std::pair<GenSolver<std::complex<double>> *, GenSparseMatrix<std::complex<double>> *>
getCholmod<std::complex<double>>(const Connectivity *cn, const DofSetArray &dsa, const ConstrainedDSA &cdsa);
