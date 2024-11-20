#ifndef _COMMUNICATOR_H_
#define _COMMUNICATOR_H_

#include <cstdio>
#include <Utils.d/resize_array.h>
#include <Utils.d/MyComplex.h>
#include <Utils.d/dbg_alloca.h>
#include "OpaqueHandle.h"
#include "BaseCommunicator.h"

class Connectivity;

#ifdef USE_MPI
#ifndef MPI_NO_CPPBIND
#define MPI_NO_CPPBIND
#endif
#include <mpi.h>
#include <Comm.d/MPICompatTraits.h>
#else
typedef int MPI_Comm;
typedef int MPI_Request;
typedef int MPI_Status;
typedef int MPI_Op;
typedef int MPI_Datatype;
#define MPI_COMM_WORLD 0
#define MPI_SUM 0
#endif

template <class T>
class CommTrace {
public:
	static MPI_Datatype MPIType;
};

template <>
class CommTrace<bool> {
private:
	static MPI_Datatype MPIType;
};

struct RecInfo {
	int cpu, len;
};

// Utility class to contain a message link information.
struct CPair {
	int locFrom, locTo;
	int glFrom, glTo;

};

int CPairCompare(void *a, void *b);

CommunicatorHandle getWorldComm();

class Communicator
{
public:
	explicit Communicator(MPI_Comm, FILE * = stderr);
	explicit Communicator(int _ncpu);
	explicit Communicator(const CommunicatorHandle &handle);
	Communicator(const Communicator &);
	template <class Type>
	Type globalSum(Type);
	bool globalMax(bool);
	template <class Type>
	Type globalMax(Type);
	template <class Type>
	Type globalMin(Type);
	template <class Type>
	void globalSum(int, Type*);
	template <class Type>
	void globalMax(int, Type*);
	template <class Type>
	void globalMin(int, Type*);
	template <class Type>
	void blockingSendTo(int cpu, int tag, Type *data, int len);
	template <class Type>
	void sendTo(int cpu, int tag, Type *data, int len);
	template <class Type>
	RecInfo recFrom(int tag, Type *data, int len);
	template <class Type>
	RecInfo recFrom(int cpu, int tag, Type *data, int len);
	template <class Type>
	void allGather(Type *send_data, int send_count, Type *recv_data,
	               int recv_count);
	template <class Type>
	void allGatherv(Type *send_data, int send_count, Type *recv_data,
	                int recv_counts[], int displacements[]);
	template <class Type>
	void allGather(Type *recv_data, int recv_count);
	template <class Type>
	void allGatherv(Type *recv_data, int recv_counts[], int displacements[]);
	template <class Type>
	void gather(Type *send_data, int send_count, Type *recv_data,
	            int recv_count, int root = 0);
	template <class Type>
	void gatherv(Type *send_data, int send_count, Type *recv_data,
	             int recv_counts[], int displacements[], int root = 0);

	template <class Type>
	void reduce(int num, Type *data, int root = 0, MPI_Op = MPI_SUM);
	template <class Type>
	void broadcast(int num, Type* data, int root = 0);
	void sync();
	int myID();
	int numCPUs();
	MPI_Comm* getCommunicator() { return &comm; }
	void split(int, int, Communicator**);
	int remoteSize();
	void waitForAllReq();

private:
	MPI_Comm comm;
	ResizeArray<MPI_Request> pendReq;
	ResizeArray<MPI_Status> reqStatus;
	int nPendReq;
	int glNumCPU;
	BaseCommunicator opaqueCommunicator;
};

class SysCom : public Communicator
{
	MPI_Comm *salinasCommunicator;
public:
	SysCom(int &argc, char **&argv);
	SysCom();
	~SysCom();
	SysCom(MPI_Comm *communicator);
	MPI_Comm* getCommunicator() { return salinasCommunicator; }
};

// The next routine provides tag from 100 to 200 cyclicaly
int uniqueTag();

extern SysCom *syscom;
extern Communicator *structCom;
template <class Type>
Type
Communicator::globalSum(Type data)
{
	Type buff;
	opaqueCommunicator.allReduce(&data, &buff, 1, SumHandle);
	return buff;
}

template <class Type>
Type
Communicator::globalMax(Type data)
{
	Type buff;
	opaqueCommunicator.allReduce(&data, &buff, 1, MaxHandle);
	return buff;
}

template <class Type>
Type Communicator::globalMin(Type data)
{
	Type buff;
	opaqueCommunicator.allReduce(&data, &buff, 1, MinHandle);
	return buff;
}

template <class Type>
void
Communicator::globalSum(int num, Type*data)
{
//	if(this->glNumCPU == 1)
//		return;
	Type *work;
	dbg_alloca(0);

	//int segSize = (num > 65536) ? 65536 : num;
	int segSize = (num > 4096) ? 4096 : num; // PJSA 6-19-07

	if(segSize > 5000)
		work = new Type[segSize];
	else
		work = (Type *)dbg_alloca(segSize*sizeof(Type));

	int offset;
	for(offset = 0; offset < num; offset +=segSize) {
		int msgSize = (num-offset < segSize) ? num-offset : segSize;
		opaqueCommunicator.allReduce(data+offset, work, msgSize, SumHandle);
		for(int i = 0; i < msgSize; ++i)
			data[offset+i] = work[i];
	}
	if(segSize > 5000)
		delete [] work;
}

template <class Type>
void
Communicator::globalMax(int num, Type*data)
{
	Type *work;
	dbg_alloca(0);

	//int segSize = (num > 65536) ? 65536 : num;
	int segSize = (num > 4096) ? 4096 : num; // PJSA 6-19-07

	if(segSize > 5000)
		work = new Type[segSize];
	else
		work = (Type *)dbg_alloca(segSize*sizeof(Type));

	int offset;
	for(offset = 0; offset < num; offset +=segSize) {
		int msgSize = (num-offset < segSize) ? num-offset : segSize;
		opaqueCommunicator.allReduce(data+offset, work, msgSize, MaxHandle);

		for(int i = 0; i < msgSize; ++i)
			data[offset+i] = work[i];
	}
	if(segSize > 5000)
		delete [] work;
}

template <class Type>
void
Communicator::globalMin(int num, Type*data)
{
	Type *work;
	dbg_alloca(0);

	//int segSize = (num > 65536) ? 65536 : num;
	int segSize = (num > 4096) ? 4096 : num; // PJSA 6-19-07

	if(segSize > 5000)
		work = new Type[segSize];
	else
		work = (Type *)dbg_alloca(segSize*sizeof(Type));

	int offset;
	for(offset = 0; offset < num; offset +=segSize) {
		int msgSize = (num-offset < segSize) ? num-offset : segSize;
		opaqueCommunicator.allReduce(data+offset, work, msgSize, MinHandle);
		for(int i = 0; i < msgSize; ++i)
			data[offset+i] = work[i];
	}
	if(segSize > 5000)
		delete [] work;
}

template <class Type>
void
Communicator::blockingSendTo(int cpu, int tag, Type *buffer, int len)
{
	opaqueCommunicator.blockingSend(buffer, len, cpu, tag);
}

template <class Type>
void
Communicator::sendTo(int cpu, int tag, Type *buffer, int len)
{
#ifdef USE_MPI
	int thisReq = nPendReq++;
  MPI_Request *req = pendReq+thisReq;
  MPI_Isend(buffer, len, CommTrace<Type>::MPIType,
            cpu, tag, comm, req);
#endif
}

template <class Type>
RecInfo
Communicator::recFrom(int tag, Type *buffer, int len)
{
#ifdef USE_MPI
	RecInfo rInfo;
  MPI_Status status;
  MPI_Recv(buffer, len,
           CommTrace<Type>::MPIType, MPI_ANY_SOURCE, tag, comm, &status);
  MPI_Get_count(&status, CommTrace<Type>::MPIType, &rInfo.len);
  rInfo.cpu = status.MPI_SOURCE;
  return rInfo;
#else
	return *(new RecInfo);
#endif
}

template <class Type>
RecInfo
Communicator::recFrom(int cpu, int tag, Type *buffer, int len)
{
#ifdef USE_MPI
	RecInfo rInfo;
  MPI_Status status;
  MPI_Recv(buffer, len,
           CommTrace<Type>::MPIType, cpu, tag, comm, &status);
  MPI_Get_count(&status, CommTrace<Type>::MPIType, &rInfo.len);
  rInfo.cpu = status.MPI_SOURCE;
  return rInfo;
#else
	return *(new RecInfo);
#endif
}

template <class Type>
void
Communicator::allGather(Type *send_data, int send_count,
                        Type *recv_data, int recv_count)
{
#ifdef USE_MPI
	MPI_Allgather(send_data, send_count, CommTrace<Type>::MPIType,
                recv_data, recv_count, CommTrace<Type>::MPIType, comm);
#endif
}

template <class Type>
void
Communicator::allGather(Type *recv_data, int recv_count)
{
#ifdef USE_MPI
	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                recv_data, recv_count, CommTrace<Type>::MPIType, comm);
#endif
}

template <class Type>
void
Communicator::allGatherv(Type *send_data, int send_count,
                         Type *recv_data, int recv_counts[], int displacements[])
{
#ifdef USE_MPI
	MPI_Allgatherv(send_data, send_count, CommTrace<Type>::MPIType,
                 recv_data, recv_counts, displacements,
                 CommTrace<Type>::MPIType, comm);
#endif
}

template <class Type>
void
Communicator::allGatherv(Type *recv_data, int recv_counts[], int displacements[])
{
#ifdef USE_MPI
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                 recv_data, recv_counts, displacements,
                 CommTrace<Type>::MPIType, comm);
#endif
}

template <class Type>
void
Communicator::gather(Type *send_data, int send_count, Type *recv_data, int recv_count, int root)
{
#ifdef USE_MPI
	MPI_Gather(send_data, send_count, CommTrace<Type>::MPIType,
             recv_data, recv_count, CommTrace<Type>::MPIType, root, comm);
#endif
}

template <class Type>
void
Communicator::gatherv(Type *send_data, int send_count,
                      Type *recv_data, int recv_counts[], int displacements[], int root)
{
#ifdef USE_MPI
	MPI_Gatherv(send_data, send_count, CommTrace<Type>::MPIType,
              recv_data, recv_counts, displacements,
              CommTrace<Type>::MPIType, root, comm);
#endif
}

#define _MESSAGE_SIZE 100000
template <class Type>
void
Communicator::reduce(int num, Type* data, int root, MPI_Op mpi_op)
{
#ifdef USE_MPI
	int maxSegSize = _MESSAGE_SIZE/sizeof(Type);
  int segSize = (num > maxSegSize) ? maxSegSize : num;
  Type *buffer;

  if(segSize > _MAX_ALLOCA_SIZE)
    buffer = new Type[segSize];
  else {
    dbg_alloca(0);
    buffer = (Type *)dbg_alloca(segSize*sizeof(Type));
  }

  for(int offset = 0; offset < num; offset +=segSize) {
    int count = (num-offset < segSize) ? num-offset : segSize;
    MPI_Reduce(data+offset, buffer, count, CommTrace<Type>::MPIType, /*MPI_SUM*/ mpi_op, root, comm);
    for(int i = 0; i < count; ++i) data[offset+i] = buffer[i];
  }
  if(segSize > _MAX_ALLOCA_SIZE)
    delete [] buffer;
#endif
}

template <class Type>
void
Communicator::broadcast(int num, Type* data, int root)
{
#ifdef USE_MPI
	int maxSegSize = _MESSAGE_SIZE/sizeof(Type);
  int segSize = (num > maxSegSize) ? maxSegSize : num;
  Type *buffer;

  if(segSize > _MAX_ALLOCA_SIZE)
    buffer = new Type[segSize];
  else {
    dbg_alloca(0);
    buffer = (Type *)dbg_alloca(segSize*sizeof(Type));
  }

  for(int offset = 0; offset < num; offset +=segSize) {
    int count = (num-offset < segSize) ? num-offset : segSize;
    if(myID() == root) for(int i = 0; i < count; i++) buffer[i] = data[offset+i];
    MPI_Bcast(buffer, count, CommTrace<Type>::MPIType, root, comm);
    if(myID() != root) for(int i = 0; i < count; ++i) data[offset+i] = buffer[i];
  }
  if(segSize > _MAX_ALLOCA_SIZE)
    delete [] buffer;
#endif
}


#ifdef NO_COMPLEX_MPI
// PJSA 1-7-2008 specializations of communication functions for platforms which do not support MPI_COMPLEX_DOUBLE
// implemented in Driver.d/MPIComm.C
template <>
complex<double> Communicator::globalSum(complex<double> data);

template <>
complex<double> Communicator::globalMax(complex<double> data);

template <>
complex<double> Communicator::globalMin(complex<double> data);

template <>
void Communicator::globalSum(int num, complex<double>* data);

template <>
void Communicator::sendTo(int cpu, int tag, complex<double> *buffer, int len);

template <>
RecInfo Communicator::recFrom(int tag, complex<double> *buffer, int len);

template <>
void Communicator::allGather(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_count);

template <>
void Communicator::allGatherv(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_counts[], int displacements[]);

template <>
void Communicator::allGather(complex<double> *recv_data, int recv_count);

template <>
void Communicator::allGatherv(complex<double> *recv_data, int recv_counts[], int displacements[]);

template <>
void Communicator::gather(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_count, int root);

template <>
void
Communicator::reduce(int num, complex<double> *data, int root, MPI_Op mpi_op);

template <>
void
Communicator::broadcast(int num, complex<double> *data, int root);
#endif

#endif
