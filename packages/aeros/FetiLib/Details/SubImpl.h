//
// Created by Michel Lesoinne on 4/3/18.
//

#ifndef FEM_SUBIMPL_H
#define FEM_SUBIMPL_H

#include <vector>
#include <Driver.d/SComm.h>
#include <FetiLib/VectorReference.h>
#include <FetiLib/SharedNodes.h>
#include <FetiLib/DOFInfo.h>
#include <FetiLib/SparseMatrix.h>
#include "FetiLib/Types.h"

namespace FetiLib {

class SubImpl {
public:
	std::vector<global_node_index> glNodes;

	SubImpl(VectorReference<const global_node_index> globalNodeIndices, const SharedNodes &sharedNodes);
};

template <typename S>
class TSubImpl : public SubImpl {
public:
	TSubImpl(global_subdomain_index subdomainIndex,
	         DOFInfo dofInfo,
	         VectorReference<const global_node_index> globalNodeIndices, VectorReference<const std::array<double,3>> X,
	         const SparseMatrix<S> &K,
	         const SharedNodes &sharedNodes);
};

template<typename S>
TSubImpl<S>::TSubImpl(global_subdomain_index subdomainIndex, DOFInfo dofInfo,
                      VectorReference<const global_node_index> globalNodeIndices,
                      VectorReference<const std::array<double, 3>> X, const SparseMatrix<S> &K,
                      const SharedNodes &sharedNodes)
		: SubImpl{globalNodeIndices, sharedNodes} {

}

}


#endif //FEM_SUBIMPL_H
