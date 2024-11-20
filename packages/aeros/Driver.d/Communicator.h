#ifndef _FS_COMMUNICATOR_H_
#define _FS_COMMUNICATOR_H_
#include <Utils.d/Connectivity.h>

#ifdef USE_MPI
#ifdef RS6000_SYS
#include "/usr/lpp/ppe.poe/include/mpi.h"
#else
#ifndef MPI_NO_CPPBIND
#define MPI_NO_CPPBIND
#endif
#include <mpi.h>
#endif // USE_MPI
#endif // RS6000_SYS

#include <map>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>

#include <Comm.d/Communicator.h>

struct
FSRecInfo {
	int cpu, len;
};

/** \brief A communication object that encapsulate communcation patters. */
class FSCommunicator {
	int thisCPU;
	int numCPU;
#ifdef USE_MPI
	MPI_Comm comm;
#endif
public:
	FSCommunicator();
#ifdef USE_MPI
	FSCommunicator(Communicator *);
	MPI_Comm& getComm() { return comm; }
#endif
	int cpuNum() { return thisCPU; }
	int size() { return numCPU; }

	void exchange(int tag, int numNeighb, gsl::span<const int> cpus,
	              int *sndPtr, int *sndData, int *rcvPtr, int *rcvData);
	void exchange(int tag, int numNeighb, gsl::span<const int> cpus,
	              gsl::span<const int> sndLen, int **sndData, gsl::span<const int> rcvLen, int **rcvData);
	void exchange(int tag, int numNeighb, gsl::span<const int> cpus,
	              gsl::span<int> sndLen, double **sndData, gsl::span<int> rcvLen, double **rcvData);
	void exchange(int tag, int numNeighb, gsl::span<const int> cpus,
	              gsl::span<const int> sndLen, std::complex<double> **sndData, gsl::span<const int> rcvLen,
	              std::complex<double> **rcvData);

	void abort(int errorCode);
	void sync() {
#ifdef USE_MPI
		MPI_Barrier(comm);
#endif
	}

	// PJSA: templated functions
	template <class Type, typename X = typename std::enable_if<CommTypeTrait<Type>::isFixedSize()>::type>
	Type globalSum(Type);
	template <class Type, typename X = typename std::enable_if<CommTypeTrait<Type>::isFixedSize()>::type>
	Type globalMax(Type);
	template <class Type, typename X = typename std::enable_if<CommTypeTrait<Type>::isFixedSize()>::type>
	Type globalMin(Type);
	template <class Type>
	void globalMax(int, Type*);
	template <class Type>
	void globalMin(int, Type*);
	template <class Type>
	void globalSum(int, Type*);

	template <typename Scalar>
	void globalSum(std::vector<Scalar> &v) {
		globalSum(v.size(), v.data());
	}
	template <typename Scalar>
	void globalMin(std::vector<Scalar> &v) {
		globalMin(v.size(), v.data());
	}
	template <typename Scalar>
	void globalMax(std::vector<Scalar> &v) {
		globalMax(v.size(), v.data());
	}
	template <class Type>
	void sendTo(int cpu, int tag, Type *data, int len);
	template <class Type>
	FSRecInfo recFrom(int tag, Type *data, int len);

	template <class Type>
	void allGather(Type *send_data, int send_count, Type *recv_data,
	               int recv_count);
	template <class Type>
	void allGatherv(Type *send_data, int send_count, Type *recv_data,
	                int recv_counts[], int displacements[]);

	template <class Type>
	void gather(Type *send_data, int send_count, Type *recv_data,
	            int recv_count, int root = 0);
	template <class Type>
	void gatherv(Type *send_data, int send_count, Type *recv_data,
	             int recv_counts[], int displacements[], int root = 0);

#ifdef USE_MPI
	template <class Type>
	void globalMpiOp(int num, Type *data, MPI_Op mpi_op);
#endif

	template <class Type>
	void reduce(int num, Type*data, int root = 0);

	template <class Type>
	void broadcast(int num, Type* data, int root = 0);

private:
	Communicator communicator;
};

void initComm(int, char *[]);
void closeComm();
// extern FSCommunicator *communicator;

class Connectivity;

class Triplet {
public:
	int from, to, cpuID;
	Triplet(int _from, int _to, int _cpuID) {
		from = _from; to = _to; cpuID = _cpuID;
	}
	bool operator < (const Triplet &x) const { // required by the STL sort algorithm
		// we want to order first by cpuID (with local coms first), then by origin and then by destination
		return cpuID < x.cpuID || (cpuID  == x.cpuID &&
		                           (from < x.from || (from == x.from && to < x.to)));
	}
	bool operator == (const Triplet &x) const {
		return cpuID == x.cpuID && from == x.from && to == x.to;
	}
};

template <class T>
class FSSubRecInfo {
public:
	gsl::span<T> data;
	int len = 0;
	int leadDim = 0;
	int nvec = 0;
	bool isSend = false;

	FSSubRecInfo() = default;

	FSSubRecInfo(int _len, int _leadDim, int _nvec) : len(_len), leadDim(_leadDim), nvec(_nvec), isSend(true) {}

	FSSubRecInfo(int _len, int _leadDim, int _nvec, bool _isSend) : len(_len), leadDim(_leadDim), nvec(_nvec), isSend(_isSend) {}

	FSSubRecInfo(const FSSubRecInfo &c) = default;

};

typedef std::map<Triplet, int> MapType;
typedef MapType::value_type ValuePair;
typedef MapType::iterator MapIter;

/** \brief Data describing a communication pattern with all neighbors.
 *
 */
class FSCommStructure {
public:
	enum Mode { Share, CopyOnSend };
	enum Symmetry { Sym, NonSym };

	void makeSendAndRecvConnectivities();

	void setLen(int glFrom, int glTo, int len, int leadDim = -1, int nvec = 1) {
		setLength(glFrom, glTo, len, leadDim, nvec);
	}

//	Mode getMode() const { return mode; }
//	Symmetry getSym() const {return sym; }
//	int getMyCPU() const { return myCPU; }
//	int getNumChannels() const { return numChannels; }
//	int *getChannelToCPU() const { return channelToCPU; }
//	const MapType &getChannelMap() const { return channelMap; }
//	const Connectivity *getSubToCPU() const { return subToCPU; }
protected:
	FSCommunicator *communicator;  // PJSA
	Mode mode;
	Symmetry sym;
	int myCPU;
	int numChannels;
	MapType channelMap;
	Connectivity *subToCPU;
	int numCPUs;
	std::vector<int> reverseChannel; // corresponding reverse channel
	std::vector<int> channelToCPU; // this CPU is marked as -1
	int numNeighbCPUs;
	std::vector<int> neighbCPUs;
	Connectivity sendConnect;
	Connectivity recvConnect;
	// Buffers for cross-memory communication on a cpu per cpu basis
	std::vector<int> crossSendLen;
	std::vector<int> crossRecvLen;

	virtual void setLength(int glFrom, int glTo, int len, int leadDim, int nvec) = 0;

public:
	virtual ~FSCommStructure() {}

};
/** FSCommPattern represent a communication pattern.
 *
 * \details
 *     Communication is based on the model that a message from one
 *     subdomain to another subdomain is made of a number of vectors.
 *     Such vectors are stored in a matrix form. That matrix may have
 *     larger columns than the actual message length. This allows a subdomain
 *     to send different subparts of a given matrix to different subdomains.
 *
 *   A communication patterns contains the following information:
 *     - Length of a vector from one subdomain to its neighbors
 *     - Number of vectors in each message
 *     - Leading dimension of the matrix containing the vectors
**/
template <class T>
class FSCommPattern: public FSCommStructure {

protected:
	std::vector< FSSubRecInfo<T> > sRecInfo;
	// Buffers for local copy-communication
	T *localDBuffer;
	// Buffers for cross-memory communication on a cpu per cpu basis
	T **crossSendBuffer;
	T **crossRecvBuffer;

public:
	FSCommPattern(FSCommunicator *communicator, const Connectivity *cpuToSub,
	              int myCPU = 0,
	              Mode = Share, Symmetry = Sym);
	~FSCommPattern();
	void finalize(); // complete the internal setup
	FSSubRecInfo<T> recData(int glFrom, int glTo);
	void sendData(int glFrom, int glTo, T *);
	FSSubRecInfo<T> getSendBuffer(int glFrom, int glTo);
	void exchange();

protected:
	void setLength(int glFrom, int glTo, int len, int leadDim, int nvec) override;
};


template <class T>
FSCommPattern<T>::FSCommPattern(FSCommunicator *_communicator, const Connectivity *cpuToSub, int _myCPU, Mode _mode, Symmetry _sym)
{
	communicator = _communicator;
	myCPU = _myCPU;
	mode = _mode;
	sym = _sym;
	numChannels = 0;
	numCPUs = communicator->size();
	subToCPU = cpuToSub->alloc_reverse();
	crossSendBuffer = 0;
	crossRecvBuffer = 0;
	localDBuffer = 0;
	sendConnect = 0;
	recvConnect = 0;
}

template <class T>
FSCommPattern<T>::~FSCommPattern()
{
	if(crossSendBuffer) {
		delete [] crossSendBuffer[0]; crossSendBuffer[0] = 0;
		delete [] crossSendBuffer; crossSendBuffer = 0;
	}
	if(crossRecvBuffer) { delete [] crossRecvBuffer; crossRecvBuffer = 0; }
	if(localDBuffer) { delete [] localDBuffer; localDBuffer = 0; }
	if(subToCPU) { delete subToCPU; subToCPU = 0; }
}

template <class T>
void
FSCommPattern<T>::setLength(int glFrom, int glTo, int len, int leadDim, int nvec)
{
	// serial function

	int cpuID = (*subToCPU)[glTo][0];
	if(leadDim < 0) leadDim = len;

	if(cpuID == myCPU) {
		channelMap.insert(ValuePair(Triplet(glFrom, glTo, -1), numChannels));
		sRecInfo.insert(sRecInfo.end(), FSSubRecInfo<T>(len, leadDim, nvec));
		numChannels++;
	}
	else {
		// make send channel
		channelMap.insert(ValuePair(Triplet(glFrom, glTo, cpuID), numChannels));
		sRecInfo.insert(sRecInfo.end(), FSSubRecInfo<T>(len, leadDim, nvec));
		numChannels++;
		// make receive channel
		channelMap.insert(ValuePair(Triplet(glTo, glFrom, cpuID), numChannels));
		sRecInfo.insert(sRecInfo.end(), FSSubRecInfo<T>());  // empty for now, see ::finalize()
		numChannels++;
	}
}

template <class T>
void
FSCommPattern<T>::finalize()
{
	int i, j;
	// Step 1. build isSend, reverseChannel and channelToCPU
	reverseChannel.resize(numChannels);
	channelToCPU.resize(numChannels);
	MapIter iter = channelMap.begin();
	while(iter != channelMap.end()) {
		Triplet key = (*iter).first;
		i = (*iter).second; // channelID
		if(sRecInfo[i].isSend) reverseChannel[i] = i+1;
		else reverseChannel[i] = i-1;
		channelToCPU[i] = key.cpuID;
		++iter;
	}

	// Step 2. make send and recv connectivities
	makeSendAndRecvConnectivities();

	// Step 3. make reverse channel FSSubRecInfo's
	if(numCPUs > 1) {
		if(sym == Sym) {
			for(i = 0; i < numChannels; ++i) {
				if(channelToCPU[i] >= 0 && sRecInfo[i].isSend) {
					sRecInfo[reverseChannel[i]].len  = sRecInfo[i].len;
					sRecInfo[reverseChannel[i]].leadDim = sRecInfo[i].leadDim;
					sRecInfo[reverseChannel[i]].nvec = sRecInfo[i].nvec;
				}
			}
		}
		else {
			// send and receive the message length and number of vectors
			std::vector<int> sendMsg(2 * sendConnect.numConnect());
			std::vector<int> recvMsg(2 * recvConnect.numConnect());
			std::vector<int*> sendPtr(numNeighbCPUs);
			std::vector<int*> recvPtr(numNeighbCPUs);
			std::vector<int> sendLen(numNeighbCPUs);
			std::vector<int> recvLen(numNeighbCPUs);
			int offset = 0;
			int totLen = 0;
			for(i = 0; i < numNeighbCPUs; ++i) {
				sendPtr[i] = sendMsg.data() + offset;
				recvPtr[i] = recvMsg.data() + offset;
				sendLen[i] = 2 * sendConnect.num(i);
				recvLen[i] = 2 * sendConnect.num(i);
				for(j = 0; j < sendConnect.num(i); ++j) {
					sendPtr[i][2*j] = sRecInfo[sendConnect[i][j]].len;
					sendPtr[i][2*j+1] = sRecInfo[sendConnect[i][j]].nvec;
					totLen += sendMsg[offset] * sendMsg[offset+1];
					offset += 2;
				}
			}
			communicator->exchange(101, numNeighbCPUs, neighbCPUs, sendLen, sendPtr.data(),
			                       recvLen, recvPtr.data());
			offset = 0;
			totLen = 0;
			for(i = 0; i < numNeighbCPUs; ++i)
				for(j = 0; j < sendConnect.num(i); ++j) {
					sRecInfo[recvConnect[i][j]].len = recvMsg[offset];
					sRecInfo[recvConnect[i][j]].leadDim = recvMsg[offset];
					sRecInfo[recvConnect[i][j]].nvec = recvMsg[offset+1];
					totLen += recvMsg[offset]*recvMsg[offset+1];
					offset += 2;
				}
		}
	}

	// Step 4. allocate the cross-memory communication buffers (if necessary)
	if((mode == Share) && (numCPUs == 1)) return;
	if((numCPUs != 1) && (numNeighbCPUs > 0)) {
		int totLen = 0;
		crossSendBuffer = new T*[numNeighbCPUs];
		crossRecvBuffer  = new T*[numNeighbCPUs];
		crossSendLen.resize(numNeighbCPUs);
		crossRecvLen.resize(numNeighbCPUs);
		for(i = 0; i < numNeighbCPUs; ++i)
			crossSendLen[i] = crossRecvLen[i] = 0;
		for(i = 0; i < numNeighbCPUs; ++i) {
			for(j = 0; j < sendConnect.num(i); ++j) {
				int channel = sendConnect[i][j];
				int msgLen = sRecInfo[channel].len * sRecInfo[channel].nvec;
				crossSendLen[i] += msgLen;
				totLen += msgLen;
			}
			for(j = 0; j < recvConnect.num(i); ++j) {
				int channel = recvConnect[i][j];
				int msgLen = sRecInfo[channel].len * sRecInfo[channel].nvec;
				crossRecvLen[i] += msgLen;
				totLen += msgLen;
			}
		}
		T *crBuff = new T[totLen];
		for(i = 0; i < numNeighbCPUs; ++i) {
			crossSendBuffer[i] = crBuff;
			for(j = 0; j < sendConnect.num(i); ++j) {
				int channel = sendConnect[i][j];
				int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
				sRecInfo[channel].data = {crBuff, msgLen};
				crBuff += msgLen;
			}
			crossRecvBuffer[i] = crBuff;
			for(j = 0; j < recvConnect.num(i); ++j) {
				int channel = recvConnect[i][j];
				int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
				sRecInfo[channel].data = {crBuff, msgLen};
				crBuff += msgLen;
			}
		}
	}
	if(mode == Share) return;

	// allocate buffer for Copy On Send mode
	int len = 0;
	for(i = 0; i < numChannels; ++i)
		if((numCPUs == 1) || (channelToCPU[i] < 0))
			len += sRecInfo[i].len * sRecInfo[i].nvec;
	localDBuffer = new T[len];
	T *cBuf = localDBuffer;
	for(i = 0; i < numChannels; ++i)
		if(numCPUs == 1 || channelToCPU[i] < 0) {
			sRecInfo[i].data = {cBuf, sRecInfo[i].len * sRecInfo[i].nvec};
			cBuf += sRecInfo[i].len * sRecInfo[i].nvec;
		}
}

template <class T>
void
FSCommPattern<T>::exchange()
{
	//if(numCPUs == 1) return; // PJSA 11-7-05
	communicator->exchange(101, numNeighbCPUs, neighbCPUs, crossSendLen,
	                       crossSendBuffer, crossRecvLen, crossRecvBuffer);
}


template <class T>
FSSubRecInfo<T>
FSCommPattern<T>::recData(int glFrom, int glTo)
{
	int cpuID = (*subToCPU)[glFrom][0];
	if (cpuID == myCPU) cpuID = -1;
		int channel = channelMap.at(Triplet(glFrom, glTo, cpuID));
		FSSubRecInfo<T> ret = sRecInfo[channel];
		if (mode != Share) ret.leadDim = ret.len;
		return ret;
}

template <class T>
void
FSCommPattern<T>::sendData(int glFrom, int glTo, T *data)
{
	int cpuID = (*subToCPU)[glTo][0];
	if(cpuID == myCPU) cpuID = -1;
	int channel = channelMap.find(Triplet(glFrom, glTo, cpuID))->second;

	if((mode == Share) && ((numCPUs == 1) || (channelToCPU[channel] < 0))) {
		sRecInfo[channel].data = {data, sRecInfo[channel].len};
		return;
	}
	// For CopyOnSend, or non-local communication copy the data into the right buffer
	int i, iVec;
	FSSubRecInfo<T> chObj = sRecInfo[channel];
	for(iVec = 0; iVec < chObj.nvec; ++iVec)
		for(i = 0; i < chObj.len; ++ i)
			chObj.data[iVec*chObj.len + i] = data[iVec*chObj.leadDim + i];
}

template <class T>
FSSubRecInfo<T>
FSCommPattern<T>::getSendBuffer(int glFrom, int glTo)
{
	int cpuID = (*subToCPU)[glTo][0];
	if(cpuID == myCPU) cpuID = -1;
	int channel = channelMap.find(Triplet(glFrom, glTo, cpuID))->second;
	FSSubRecInfo<T> chObj = sRecInfo[channel];
	if(mode == Share) // Check if we have a send buffer. If not, return NULL
		chObj.data = gsl::span<T>{};
	return chObj;
}

#ifdef _TEMPLATE_FIX_
#include <Driver.d/Communicator.C>
#endif

#endif
