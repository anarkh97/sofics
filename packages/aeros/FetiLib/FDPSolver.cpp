//
// Created by Michel Lesoinne on 10/31/17.
//

#include <Driver.d/SComm.h>
#include <set>
#include <Comm.d/BaseCommunicator.h>
#include <Comm.d/MPICompatTraits.h>
#include "FDPSolver.h"
#include <Feti.d/Feti.h>
#include "Math.d/CuCSparse.h"
#include "Math.d/SparseSet.h"
#include "ConcreteSub.h"

namespace FetiLib {

struct SharedDataInfo {
	std::vector<int> cpus;
	Connectivity sharedData;
};

/** \brief A reprentation of uniform segmentation the integer interval [0,n).
 * \details The size of each segment can differ only by one, being larger by 1 when the segment index
 * is less than the remainder of the division.
 */
class Segmenter {
	Segmenter(global_node_index nNodes, int nProc) : quotient(nNodes/nProc), remainder(nNodes%nProc) {}
	int getProc(global_node_index node) {
		auto pr = node/quotient;
		auto r = node % quotient;
		if(r < pr && r < remainder)
			return pr-1;
		return pr;
	}
private:
	global_node_index quotient;
	global_node_index remainder;
};

/** \brief Function creating useful connectivities from a matrix structure
 *
 */
//template <typename T>
//Connectivity getConnectivities(const Subdomain<T> &subdomain) {
//
//}

/**
 * \details makeSubs needs to equip ConcreteSub with:
 * nodeToNode connectivity.
 *
 * @tparam T
 * @param subdomains
 * @return
 */
template <typename T>
std::vector<std::unique_ptr<ConcreteSub<T>>> makeSubs(const std::vector<Subdomain<T>> &subdomains) {

}

namespace tpl {

/** \brief Obtain the set of nodes involved in a set of DOFs.
 *
 * @param dofInfo An array of pair of global node number, DofType for all the DOFs in this CPU.
 * @return The set of all global node numbers that appear in dofInfo.
 */
std::set<global_node_index> globalNodeSet(const DOFInfo &dofInfo) {
	// Could be done in one line, using transform iterators.
	std::set<global_node_index> nodes;
	for(const auto &info : dofInfo)
		nodes.insert(info.first);
	return nodes;
}

SharedDataInfo
buildSharedDataInfo(const std::set<global_node_index> &glNodes, const BaseCommunicator &communicator) {
	auto it = std::max_element(glNodes.begin(), glNodes.end());
	global_node_index largestNode = (it == glNodes.end()) ? 0 : *it;
	global_node_index numGlNodes = communicator.globalMax(largestNode+1);

    throw "Unfinished.";
}

//SComm buildSComm() {
//	return SComm{0};
//}

class OptionSet {
public:
	void addOption(const char *name, double v) {
		doubleOptions.insert({std::string{name}, v});
	}

	double getDouble(const std::string name, double defaultValue) {
		auto it = doubleOptions.find(name);
		return it == doubleOptions.end() ? defaultValue : it->second;
	}

private:
	std::map<std::string, double> doubleOptions;
};

/** \brief The base class of FetiDP(H) implementation */
class DPSImpl {
public:
	DPSImpl(Com communicator) : communicator(CommunicatorHandle(communicator)) {}

	void setOption(const char *optionName, double v) { options.addOption(optionName, v); }

protected:
	const BaseCommunicator communicator;
	OptionSet options;
};

/** \brief Implementation class for the FetiDPSolver. */
template <typename T>
class FetiDP : public DPSImpl {
	FetiDP(std::vector<Subdomain<T>> subdomains, Com communicator);
private:
	std::vector<Subdomain<T>> subdomains;
	std::unique_ptr<GenFetiDPSolver<T>> fetiSolver;
};

template<typename T>
FetiDP<T>::FetiDP(std::vector<Subdomain<T>> subdomains, Com communicator) :
		DPSImpl(communicator),
		subdomains(std::move(subdomains)){
}
namespace {

template <typename T>
std::unique_ptr<FetiSub<T>>
makeFetiSub(global_subdomain_index globSubIndex, int localSubIndex, const Subdomain <T> &subData) {
	// The subdomain as seen from the user's point of view.
	SubImpl &userSub = getter(subData);

//	global_subdomain_index globSubIndex, local_subdomain_index locSubIndex,
//	const SubImpl &subdomain, SComm *sComm

	SComm *sc = nullptr;
//	SComm *sc = new SComm(nConnect[NESubMap[subI]], connectedDomain[NESubMap[subI]],
//	                      remoteID[NESubMap[subI]], interfNode[NESubMap[subI]]);
	sc->locSubNum = localSubIndex;

	auto sub = std::make_unique<ConcreteSub<T>>(globSubIndex, localSubIndex, subData);

	sub->setSComm(sc);

	return std::move(sub);
}

using SubRange = std::pair<global_subdomain_index,global_subdomain_index >;

struct CompareRange {
	using is_transparent = void;

	bool operator()(const SubRange &a, const SubRange &b) const { return a < b; }

	bool operator()(global_subdomain_index a, const SubRange &b) const { return a < b.first; }

	bool operator()(const SubRange &a, global_subdomain_index b) const { return a.first < b; }
};

bool operator<(global_subdomain_index i, const SubRange &sr) { return i < sr.first; }

/** \brief Obtain a vector with the processors having the subdomains with given indices.
 *
 * @param indices The indices for which the processors are desired.
 * @param firstIndex Index of the first subdomain local to this process.
 * @param numSubs Number of subdomains in this processor.
 * @param comm Communicator to which all processes belong.
 * @return
 */
std::vector<int> processWithGlobalSubs(const std::vector<global_subdomain_index> &indices,
                                       global_subdomain_index firstIndex, int numSubs,
                                       const BaseCommunicator &comm)
{
	int maxNumSubs = comm.globalMax(numSubs);
	int minNumSubs = comm.globalMin(numSubs);
	if (maxNumSubs == 1 && minNumSubs == 1)
		return {indices.begin(), indices.end()};

	global_subdomain_index totalNumSubs = comm.globalSum(static_cast<global_subdomain_index>(numSubs));
	std::vector<int> result;
	result.reserve(indices.size());


	std::pair<global_subdomain_index,global_subdomain_index > subRange{ firstIndex, firstIndex+numSubs };
	// Create a window so that all other processes can find what subdomains we own.
	auto rangeWindow = comm.window(&subRange, 1);

	std::map<SubRange, int, CompareRange> procRanges{ {{0,0}, -1}, {{totalNumSubs,totalNumSubs}, comm.commSize()}};
	for(auto glSub : indices) {
		// Get an iterator to the range strictly above the desired sub
		auto beyond = procRanges.upper_bound(glSub);
		// Now form the iterator (guaranteed to exist) that has the highest range, not beyond glSub
		auto low = beyond;
		--low;
		// Loop until glSub belongs to the range pointed to by low
		while(low->first.second <= glSub) {
			// Span of subdomains between low and and beyond
			auto span = beyond->first.first - low->first.second;
			double avgSubPerProc = double(span)/(beyond->second - low->second -1);
			// Estimate where the containing subdomain falls.
			int procGuess = low->second+1 + std::floor((glSub - low->first.second)/avgSubPerProc);
			// Get the data from the guessed process
			rangeWindow.sharedLock(procGuess);
			SubRange remRange;
			rangeWindow.get(&remRange, 1, procGuess, 0);
			rangeWindow.unlock(procGuess);
			auto insIt = procRanges.insert({remRange, procGuess}).first;
			if(remRange.first > glSub)
				beyond = insIt;
			else
				low = insIt;
		}
		result.push_back(low->second);
	}

	return result;
}


class ProcFinder {
public:
	ProcFinder(SubRange myRange, const BaseCommunicator &comm)
			: myRange{myRange},
			  rangeWindow{ comm.window(&myRange, 1) } {
		auto numSubs = myRange.second-myRange.first;
		int maxNumSubs = comm.globalMax(numSubs);
		int minNumSubs = comm.globalMin(numSubs);
		isTrivial = maxNumSubs == 1 && minNumSubs == 1;
		if(!isTrivial) {
			global_subdomain_index totalNumSubs = comm.globalSum(static_cast<global_subdomain_index>(numSubs));
			procRanges.insert({{0,0}, -1});
			procRanges.insert({{totalNumSubs,totalNumSubs}, comm.commSize()});
		}
	}

	std::vector<int> processWithGlobalSubs(const std::vector<global_subdomain_index> &indices);
private:
	bool isTrivial; //!< Trivial if each process and a single subdomain.
	SubRange myRange;
	Window rangeWindow;
	std::map<SubRange, int, CompareRange> procRanges;
};

std::vector<int> FetiLib::tpl::ProcFinder::processWithGlobalSubs(const std::vector<global_subdomain_index> &indices) {
	if(isTrivial)
		return { indices.begin(), indices.end() };
	std::vector<int> result;
	result.reserve(indices.size());


	for(auto glSub : indices) {
		// Get an iterator to the range strictly above the desired sub
		auto beyond = procRanges.upper_bound(glSub);
		// Now form the iterator (guaranteed to exist) that has the highest range, not beyond glSub
		auto low = beyond;
		--low;
		// Loop until glSub belongs to the range pointed to by low
		while(low->first.second <= glSub) {
			// Span of subdomains between low and and beyond
			auto span = beyond->first.first - low->first.second;
			double avgSubPerProc = double(span)/(beyond->second - low->second -1);
			// Estimate where the containing subdomain falls.
			int procGuess = low->second+1 + std::floor((glSub - low->first.second)/avgSubPerProc);
			// Get the data from the guessed process
			rangeWindow.sharedLock(procGuess);
			SubRange remRange;
			rangeWindow.get(&remRange, 1, procGuess, 0);
			rangeWindow.unlock(procGuess);
			auto insIt = procRanges.insert({remRange, procGuess}).first;
			if(remRange.first > glSub)
				beyond = insIt;
			else
				low = insIt;
		}
		result.push_back(low->second);
	}

	return result;
}

}



template<typename T>
DPSolver<T>::DPSolver(std::vector<Subdomain<T>> subdomains, CommunicatorHandle communicator) {
	if(subdomains.size() > 1)
		throw "Multiple subs/CPU not yet supported.";

	BaseCommunicator comm{communicator};
	// Get a global numbering by adding the local number to the number of subdomains in the processes
	// with strictly lower indices in the communicator.
	global_subdomain_index firstIndex = comm.exScan(static_cast<global_node_index>(subdomains.size()));

	// Build the FetiSubs
	for(int iSub = 0; iSub < subdomains.size(); ++iSub)
		std::unique_ptr<FetiSub<T>> fetiSub{makeFetiSub(firstIndex+iSub, iSub, subdomains[iSub])};
}

template<typename T>
void DPSolver<T>::setOption(const char *optionName, double v) {
	pImpl->setOption(optionName, v);
}

template class DPSolver<double>;
template class DPSolver<std::complex<double>>;

} // namespace tpl

} //namespace FetiLib
