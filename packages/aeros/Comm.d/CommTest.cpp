//
// Created by Michel Lesoinne on 2/13/18.
//
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <chrono>
#include <Utils.d/AutoTimer.h>
#include "BaseCommunicator.h"

auto now() { return std::chrono::high_resolution_clock::now(); }

using timepoint = decltype(now());
double elapsed(timepoint t0, timepoint t1) { return std::chrono::duration<double>{t1-t0}.count();}

class SubD {
public:
	SubD(std::array<int,3> subSize, std::array<int,3> nsubPerDirection, std::array<int, 3> thisIndex);

	std::map<int, std::vector<long> > getNodeSharing(const BaseCommunicator &communicator) const;
private:
	std::vector<long> nodeIndices;
	long totalNumNodes;
public:
	const std::vector<long> &getNodeIndices() const;

	long getTotalNumNodes() const;
	//!< Number of nodes in the whole model.
};

template <typename T>
T product(const std::array<T,3> &d) {
	return d[0]*d[1]*d[2];
}

template <typename T>
T productP1(const std::array<T,3> &d) {
	return (d[0]+1)*(d[1]+1)*(d[2]+1);
}

SubD::SubD(std::array<int, 3> subSize, std::array<int, 3> nsubPerDirection, std::array<int, 3> thisIndex) {
	std::array<long,3> totalSize { subSize[0]*nsubPerDirection[0]+1,
	                               subSize[1]*nsubPerDirection[1]+1,
	                               subSize[2]*nsubPerDirection[2]+1 };

	totalNumNodes = productP1(totalSize);
	nodeIndices.reserve(productP1(subSize));

	std::array<long,3> corner{ subSize[0]*thisIndex[0],
	                           subSize[1]*thisIndex[1],
	                           subSize[2]*thisIndex[2] };

	auto index = [&] (long i, long j, long k) {
		return i+(totalSize[0]*(j+totalSize[1]*k));
	};

	for(long i = 0; i <= subSize[0]; ++i)
		for(long j = 0; j <= subSize[1]; ++j)
			for(long k = 0; k <= subSize[2]; ++k)
				nodeIndices.push_back(index(i+corner[0], j+corner[1], k+corner[2]));
}

/** \brief Obtain the list of processes who are canonic owners of the nodes and for each owner the list of nodes.
 *
 * @param totalNumNodes Maximum index+1 of all the nodes in all the processors.
 * @param glNodes Array of global node numbers.
 * @param communicator Communicator over which the processors taking park communicate.
 * @return A map between process index and the array of global node indices whose canonic owner is that process.
 */
std::map<int, std::vector<long> >
getCanonicalOwners(long totalNumNodes, const std::vector<long> &glNodes, int numProc) {

	auto canonicOwner = [ nPerProc=totalNumNodes/numProc, maxRem =  totalNumNodes%numProc](long node) {
		auto baseProc = node/nPerProc;
		auto rem = node % nPerProc;
		if(rem < maxRem && rem < baseProc-1)
			return baseProc-1;
		return baseProc;
	};

	std::map<int, std::vector<long>> canonicNodePosition;
	for(auto node: glNodes) {
		auto &list = canonicNodePosition[canonicOwner(node)];
		list.push_back(node);
	}
	return canonicNodePosition;
};

std::map<int, std::vector<long> > SubD::getNodeSharing(const BaseCommunicator &communicator) const {
	int numProc = communicator.commSize();

	std::map<int, std::vector<long>> canonicNodePosition = getCanonicalOwners(totalNumNodes,
	                                                                          nodeIndices, communicator.commSize());
#ifdef USE_MPI
	long info[2] = {0,0};
	auto window = communicator.window(info, 2);
	window.open();
	long one = 1;
	for(auto &info : canonicNodePosition)
		window.accumulate(SumHandle, &one, 1, info.first, 0);
	window.close();
#endif

	return std::map<int, std::vector<long>>();
}

const std::vector<long> &SubD::getNodeIndices() const {
	return nodeIndices;
}

long SubD::getTotalNumNodes() const {
	return totalNumNodes;
}

void testGlobals(const BaseCommunicator &communicator) {
	communicator.barrier();
	double result;
	auto t1 = now();
	for(int i = 0; i < 100; ++i) {
		AutoTimer<"globalSum"_hash> timer;
		result += communicator.globalSum(i);
	}
	auto t2 = now();
	if(communicator.rank() == 0) {
		const auto &timing = AutoTimer<"globalSum"_hash>::getData();
		std::cout << "100 Global sum took: " << elapsed(t1, t2) * 1e6 << "µs. Result is " << result << "." << std::endl;
		std::cout << "Mimum time: " << timing.minTime.load()*1e-3 << "µs. Maximum: " << timing.maxTime*1e-3 << "µs."
		          << " Count:" << timing.count << std::endl;
	}
}


int main(int argc, char *argv[]) {
#ifdef USE_MPI
	int required = MPI_THREAD_MULTIPLE, provided;
	MPI_Init_thread(&argc, &argv, required, &provided);
	BaseCommunicator comm((MPI_Comm)MPI_COMM_WORLD);
	int myRank = comm.rank();
	int commSize = comm.commSize();
	if(myRank == 0)
		std::cout << "I was requiring " << required << " I got " << provided << std::endl;

	testGlobals(comm);
	int databank[4096];
	int &myData = databank[0];
	myData = 2*myRank;
	int remData;
	int operandData = myRank+3;

	int myDest = (myRank+1)%commSize;

	auto t0 = now();
	for(int i = 0; i < 16; ++i)
	{
		auto window = comm.window(&myData, 4096);

		AutoTimer<"fetchAndOp"_hash> timer;
		window.open();
		window.fetchAndOp(ProdHandle, &operandData, &remData, myDest, 0);
		window.close();
	}
	auto dt = elapsed(t0, now());
	auto &faoData = AutoTimer<"fetchAndOp"_hash>::getData();
	if(myRank == 0)
		std::cout << "Time to do the work is "
		          << AutoTimer<"fetchAndOp"_hash>::getData().minTime.load()*1e-3 << "µs minimum. "
		          << AutoTimer<"fetchAndOp"_hash>::getData().maxTime.load()*1e-3 << "µs maximum."
		          << (faoData.totalTime-faoData.maxTime)/(faoData.count-1)*1e-3 << "µs average with max rejection."
		          << std::endl;
//	for(int i = 0; i < commSize; ++i) {
//		comm.barrier();
//		if(myRank == i)
//			std::cout << i << " worked with " << myDest << " where I found " << remData << " and I was made to have " << myData << std::endl;
//	}

	if(myRank == 0)
		std::cout << "Going to try async send and receives" << std::endl;
	double sendData[4];
	double recData[2][4];
	for(int i = 0; i < 4; ++i) {
		sendData[i] = i*myRank;
	}
	int lowerProc = myRank == 0 ? commSize-1 : myRank-1;
	int higherProc = myRank+1 == commSize ? 0 : myRank+1;
	for(int i = 0; i < 10; ++i) {
		comm.barrier();
		AutoTimer<"nonBlocking"_hash> timer;
		RequestVector requests;
		auto r1 = comm.nonBlockingReceive(recData[0], 4, lowerProc, 101+i);
		auto r2 = comm.nonBlockingReceive(recData[1], 4, higherProc, 101+i);
		auto r3 = comm.nonBlockingSend(sendData, 4, lowerProc, 101+i);
		auto r4 = comm.nonBlockingSend(sendData, 4, higherProc, 101+i);
		requests.emplace_back(r1);
		requests.emplace_back(r2);
		requests.emplace_back(r3);
		requests.emplace_back(r4);
		auto results = requests.waitAll();
	}
	const auto &timerData =  AutoTimer<"nonBlocking"_hash>::getData();
	for(int i = 0; i < 4; ++i) {
		if(recData[0][i] != i*lowerProc)
			std::cerr << "Incorrect receive from lower" << std::endl;
		if(recData[1][i] != i*higherProc)
			std::cerr << "Incorrect receive from higher" << std::endl;
	}

	if(myRank == 0)
		std::cout << "Time to do the async send/receive is: " << timerData.minTime.load()*1e-3 << "µs minimum. "
		          << timerData.maxTime.load()*1e-3 << "µs maximum."
				<< (timerData.totalTime-timerData.maxTime)/(timerData.count-1)*1e-3 << "µs average with max rejection." << std::endl;
	MPI_Finalize();
#endif
}
