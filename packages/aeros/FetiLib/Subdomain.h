//
// Created by Michel Lesoinne on 2/1/18.
//

#ifndef FEM_SUBDOMAIN_H
#define FEM_SUBDOMAIN_H

#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <FetiLib/DOFInfo.h>
#include <FetiLib/SharedNodes.h>
#include <FetiLib/Types.h>
#include "SparseMatrix.h"

namespace FetiLib {

/** \brief Base class for actual implementations */
class SubImpl;

template <typename S> class Subdomain;

template <typename S> SubImpl &getter(const Subdomain<S> &s);


template <typename S>
class Subdomain {
public:

	/** \brief Constructor.
	 *
	 * @param subdomainIndex Global subdomain index.
	 * @param dofInfo Information relative to the DOFs.
	 * @param globalNodeIndices Local to global index mapping of the nodes involved in this subdomain.
	 * @param X Coordinates of the nodes.
	 * @param K System matrix for this subdomain.
	 * @param sharedNodes Data relative to the shared nodes.
	 */
	Subdomain(global_subdomain_index subdomainIndex,
	          DOFInfo dofInfo,
	          VectorReference<const global_node_index> globalNodeIndices,
	          VectorReference<const std::array<double,3>> X,
	          const SparseMatrix<S> &K,
	          const SharedNodes &sharedNodes
	);

	~Subdomain();

private:
	std::unique_ptr<SubImpl> subImpl;

	friend SubImpl &getter<>(const Subdomain<S> &s);
};

}


#endif //FEM_SUBDOMAIN_H
