//
// Created by Michel Lesoinne on 4/16/18.
//
#include <Comm.d/BaseCommunicator.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <map>
#include <cmath>
#include <vector>
#include <iostream>


namespace {

using global_subdomain_index = int;
using SubRange = std::pair<global_subdomain_index, global_subdomain_index>;

struct CompareRange {
	using is_transparent = void;

	bool operator()(const SubRange &a, const SubRange &b) const { return a < b; }

	bool operator()(global_subdomain_index a, const SubRange &b) const { return a < b.first; }

	bool operator()(const SubRange &a, global_subdomain_index b) const { return a.first < b; }
};

bool operator<(global_subdomain_index i, const SubRange &sr) { return i < sr.first; }

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

	/** \brief Obtain the process IDs that contain a given list of subdomains.
	 *
	 * @param indices The vector of subdomain indices whose owning processes are sought.
	 * @return The vector of processors containing the input subdomains, in the matching order.
	 */
	std::vector<int> processWithGlobalSubs(const std::vector<global_subdomain_index> &indices);
private:
	bool isTrivial; //!< Trivial if each process and a single subdomain.
	SubRange myRange;
	Window rangeWindow;
	std::map<SubRange, int, CompareRange> procRanges;
};

std::vector<int> ProcFinder::processWithGlobalSubs(const std::vector<global_subdomain_index> &indices) {
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

int numSubsForProc(int proc) {
	return 1+(proc%3);
}

void addSubsFor(int tg, std::vector<int> &subs, std::vector<int> &processors) {
	int fProc = 0;
	for(int i = 0; i < tg; ++i)
		fProc += numSubsForProc(i);
	int np = numSubsForProc(tg);
	for(int s = fProc; s < fProc+np; ++s) {
		subs.push_back(s);
		processors.push_back(tg);
	}
}

}



int main(int argc, char *argv[]) {
#ifdef USE_MPI
	MPI_Init(&argc,&argv);
	BaseCommunicator baseCommunicator{(MPI_Comm)MPI_COMM_WORLD};

	int numProc = baseCommunicator.commSize();
	int myProc =  baseCommunicator.rank();
	//The number of subdomains in this processor.
	int myNumSub = numSubsForProc(myProc);

	global_subdomain_index  myFirstSubc = baseCommunicator.exScan(static_cast<global_subdomain_index > (myFirstSubc));

	SubRange myRange{myFirstSubc, myFirstSubc+myNumSub};

	ProcFinder procFinder(myRange, baseCommunicator);

	auto globalNumProc = baseCommunicator.commSize();
	std::vector<int> subsToBeFound;
	std::vector<int> targetProcessors;
	for(int p = 0; p < 3; ++p) {
		int tg = (2*p+1)%globalNumProc;
		addSubsFor(tg, subsToBeFound, targetProcessors);
	}

	auto procsWithSubs = procFinder.processWithGlobalSubs(subsToBeFound);

	if(procsWithSubs != targetProcessors)
		std::cerr << "Non match" << std::endl;
	MPI_Finalize();
	return procsWithSubs == targetProcessors;
#else
	std::cout << "Test is not operational when MPI is disabled." << std::endl;
	return 0;
#endif
}
