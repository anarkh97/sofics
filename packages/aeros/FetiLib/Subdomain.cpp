//
// Created by Michel Lesoinne on 2/1/18.
//

#include <complex>

#include <Utils.d/Connectivity.h>
#include "Subdomain.h"
#include "DOFInfo.h"
#include "Details/SubImpl.h"

Connectivity dofToNode(const FetiLib::DOFInfo &dofInfo) {
	using dof_t = int;
	using nt = decltype(dofInfo[0].first);
	std::vector<std::pair<size_t, nt>> dofNodeVec;
	dofNodeVec.reserve(dofInfo.size());
	for(dof_t i = 0; i < dofInfo.size(); ++i)
		dofNodeVec.emplace_back(i, dofInfo[i].first);
	return Connectivity::fromLinkRange(dofNodeVec);
}

template <typename S>
Connectivity dofToDof(const FetiLib::SparseMatrix<S> &matrix) {
	std::vector<size_t> pointers{matrix.outerIndexPtr(), matrix.outerIndexPtr()+matrix.outerSize()};
	std::vector<int> targets{matrix.innerIndexPtr(), matrix.innerIndexPtr()+matrix.innerSize()};
	return {static_cast<int>(pointers.size()-1), std::move(pointers), std::move(targets)};
}

namespace FetiLib {

template<typename S>
Subdomain<S>::Subdomain(global_subdomain_index subdomainIndex, DOFInfo dofInfo,
                        VectorReference<const global_node_index> globalNodeIndices,
                        VectorReference<const std::array<double, 3>> X, const SparseMatrix<S> &K,
                        const SharedNodes &sharedNodes)
		: subImpl(std::make_unique<TSubImpl<S>>(subdomainIndex, dofInfo, globalNodeIndices, X, K, sharedNodes))
{
}

template<typename S>
Subdomain<S>::~Subdomain() {
}


template class Subdomain<double>;
template class Subdomain<std::complex<double>>;

} // namespace FetiLib