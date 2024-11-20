//
// Created by Michel Lesoinne on 11/3/17.
//

#include "Driver.d/Communicator.h"

void
FSCommStructure::makeSendAndRecvConnectivities()
{
	int i;
	// First step, find all the neighboring CPUs
	std::vector<int> glToLocCPU(numCPUs, -1);

	numNeighbCPUs = 0;
	for(i = 0; i < numChannels; ++i) {
		int cpuID;
		if((cpuID = channelToCPU[i]) >= 0) {
			if(glToLocCPU[cpuID] < 0) glToLocCPU[cpuID] = numNeighbCPUs++;
		}
	}
	neighbCPUs.resize(numNeighbCPUs);
	for(i = 0; i < numCPUs; ++i)
		if(glToLocCPU[i] >= 0)
			neighbCPUs[glToLocCPU[i]] = i;

	// Now count the send and receives (actually they should be the same for
	// a symmetric communication
	std::vector<size_t> sendPtr(numNeighbCPUs+1, 0);
	std::vector<size_t> recvPtr(numNeighbCPUs+1, 0);

	MapIter iter = channelMap.begin();
	while(iter != channelMap.end()) {
		Triplet key = (*iter).first;
		int cpuID = key.cpuID;
		if(cpuID >= 0) {
			if((*subToCPU)[key.from][0] == myCPU) // send pair
				sendPtr[glToLocCPU[cpuID]]++;
			else // receive
				recvPtr[glToLocCPU[cpuID]]++;
		}
		++iter;
	}

	// Make the actual pointers (shifted for easy fill in)
	for(i = 0; i < numNeighbCPUs; ++i) {
		recvPtr[i+1] += recvPtr[i];
		sendPtr[i+1] += recvPtr[i];
	}

	std::vector<int> sendTrgt(sendPtr[numNeighbCPUs]);
	std::vector<int> recvTrgt(recvPtr[numNeighbCPUs]);
	// Final fill in
	iter = channelMap.begin();
	while(iter != channelMap.end()) {
		Triplet key = (*iter).first;
		int channelID = (*iter).second;
		int cpuID = key.cpuID;
		if(cpuID >= 0) {
			if((*subToCPU)[key.from][0] == myCPU) {
				sendPtr[glToLocCPU[cpuID]]--;
				sendTrgt[sendPtr[glToLocCPU[cpuID]]] = channelID;
			}
			else {
				recvPtr[glToLocCPU[cpuID]]--;
				recvTrgt[recvPtr[glToLocCPU[cpuID]]] = channelID;
			}
		}
		++iter;
	}

	sendConnect = Connectivity(numNeighbCPUs, std::move(sendPtr), std::move(sendTrgt));
	recvConnect = Connectivity(numNeighbCPUs, std::move(recvPtr), std::move(recvTrgt));

}