//
// Created by Michel Lesoinne on 10/31/17.
//

#ifndef FEM_FDPSOLVER_H
#define FEM_FDPSOLVER_H

#include <complex>
#include <memory>
#if __cpp_lib_experimental_propagate_const >= 201505
#include <experimental/propagate_const>
#endif
#ifdef USE_MPI
#include <mpi.h>
namespace FetiLib {
using Com = MPI_Comm;
inline auto default_com() { return MPI_COMM_WORLD; }
}
#else
namespace FetiLib {
using Com = int;
inline auto default_com() { return -1; }
}
#endif

#include <Eigen/Dense>

#include <FetiLib/DOFInfo.h>
#include "Subdomain.h"
#include <Comm.d/OpaqueHandle.h>

namespace FetiLib {

namespace tpl {


class DPSImpl;

#if __cpp_lib_experimental_propagate_const >= 201505
template <typename O>
using const_enforcing_unique_ptr = std::experimental::propagate_const<std::unique_ptr<O>>;
#else
template <typename O>
using const_enforcing_unique_ptr = std::unique_ptr<O>;
#endif

/** \brief FETI-DP solver for double or complex. */
template <typename T>
class DPSolver {
	/** \brief Constructor with an array of subdomains and a communicator.
	 *
	 * @param subdomains
	 * @param communicator
	 */
	DPSolver(std::vector<Subdomain<T>> subdomains,
			 CommunicatorHandle communicator = getWorldComm());

	/// \brief Set the relative tolerance of the solver.
	void setTolerance(double epsilon) { setOption("tolerance", epsilon); }
	/// \brief Set the wave number for a Helmholtz problem.
	void setWaveNumber(double k) { setOption("wave_number", k); }

	bool solve(VectorReference<const T> rhs, VectorReference<T> solution);
private:
	/** \brief Pointer to the opaque implementation of the solver.
	 * \details The pointer is wrapped to enforce const correctness.
	 */
	const_enforcing_unique_ptr<DPSImpl> pImpl;

	/// \brief Generic way of setting an option.
	void setOption(const char *optionName, double v);
};

}

using FDPSolver = tpl::DPSolver<double>;
using FDPHSolver = tpl::DPSolver<std::complex<double>>;

} // Namespace FetiLib

#endif //FEM_FDPSOLVER_H
