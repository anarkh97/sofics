//
// Created by Michel Lesoinne on 3/30/18.
//
#include "DOFManip.h"

namespace FetiLib {
namespace Details {

DofSetArray buildDSA(const DOFInfo &dofInfo) {
	local_node_index maxNodeIndex = 0;
	DofSet allDofs(0);
	for(const auto &info : dofInfo) {
		maxNodeIndex = std::max(maxNodeIndex, info.first);
		allDofs |= DofSet(static_cast<int>(info.second));
	}
	DofSetArray dofSetArray(maxNodeIndex+1);
	for(const auto &info : dofInfo)
		dofSetArray[info.first] |= allDofs;
	return dofSetArray;
}

DofSetArray buildCDSA(const DOFInfo &dofInfo) {
	local_node_index maxNodeIndex = 0;
	DofSet allDofs(0);
	for(const auto &info : dofInfo) {
		maxNodeIndex = std::max(maxNodeIndex, info.first);
		allDofs |= DofSet(static_cast<int>(info.second));
	}
	DofSetArray dofSetArray(maxNodeIndex+1);
	for(const auto &info : dofInfo)
		dofSetArray[info.first] |= DofSet(static_cast<int>(info.second));
	return dofSetArray;
}

std::pair<DofSetArray, ConstrainedDSA> buildNodalDOFInfo(const DOFInfo &dofInfo) {
	local_node_index maxNodeIndex = 0;
	DofSet allDofs(0);
	for(const auto &info : dofInfo) {
		maxNodeIndex = std::max(maxNodeIndex, info.first);
		allDofs |= DofSet(static_cast<int>(info.second));
	}
	DofSetArray dofSetArray(maxNodeIndex+1);
	DofSetArray cdsaBase(maxNodeIndex+1);
	for(const auto &info : dofInfo) {
		cdsaBase[info.first] |= DofSet(static_cast<int>(info.second));
		dofSetArray[info.first] |= allDofs;
	}
	ConstrainedDSA constrainedDSA(cdsaBase, 0);
	return { std::move(dofSetArray), std::move(constrainedDSA) };
}

} // namespace Details

}// namespace FetiLib