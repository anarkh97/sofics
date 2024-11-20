#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <Comm.d/Communicator.h>

template <class Type, typename X>
Type
FSCommunicator::globalSum(Type data)
{
	return communicator.globalSum(data);
}

template <class Type, typename X>
Type
FSCommunicator::globalMax(Type data)
{
	return communicator.globalMax(data);
}

template <class Type, typename X>
Type FSCommunicator::globalMin(Type data)
{
	return communicator.globalMin(data);
}

template <class Type>
void
FSCommunicator::globalSum(int num, Type*data)
{
	communicator.globalSum(num, data);
}

template <class Type>
void
FSCommunicator::globalMax(int num, Type*data)
{
	communicator.globalMax(num, data);
}

template <class Type>
void
FSCommunicator::globalMin(int num, Type*data)
{
	communicator.globalMin(num, data);
}

#ifdef USE_MPI
template <class Type>
void
FSCommunicator::globalMpiOp(int num, Type *data, MPI_Op mpi_op)
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
		MPI_Allreduce(data+offset, work, msgSize,
		              CommTrace<Type>::MPIType, mpi_op, comm);
		for(int i = 0; i < msgSize; ++i)
			data[offset+i] = work[i];
	}
	if(segSize > 5000)
		delete [] work;
}
#endif

template <class Type>
void
FSCommunicator::sendTo(int cpu, int tag, Type *buffer, int len)
{
	communicator.sendTo(cpu, tag, buffer, len);
//#ifdef USE_MPI
//	MPI_Send(buffer, len, CommTrace<Type>::MPIType, cpu, tag, comm);
///*
//  int thisReq = nPendReq++;
//  MPI_Request *req = pendReq+thisReq;
//  MPI_Isend(buffer, len, CommTrace<Type>::MPIType,
//            cpu, tag, comm, req);
//*/
//#endif
}

template <class Type>
FSRecInfo
FSCommunicator::recFrom(int tag, Type *buffer, int len)
{
	FSRecInfo rInfo;
#ifdef USE_MPI
	MPI_Status status;
	MPI_Recv(buffer, len,
	         CommTrace<Type>::MPIType, MPI_ANY_SOURCE, tag, comm, &status);
	MPI_Get_count(&status, CommTrace<Type>::MPIType, &rInfo.len);
	rInfo.cpu = status.MPI_SOURCE;
#else
	rInfo.len = 0;
  rInfo.cpu = 0;
#endif
	return rInfo;
}

template <class Type>
void
FSCommunicator::allGather(Type *send_data, int send_count,
                          Type *recv_data, int recv_count)
{
#ifdef USE_MPI
	MPI_Allgather(send_data, send_count, CommTrace<Type>::MPIType,
	              recv_data, recv_count, CommTrace<Type>::MPIType, comm);
#endif
}

template <class Type>
void
FSCommunicator::allGatherv(Type *send_data, int send_count,
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
FSCommunicator::gather(Type *send_data, int send_count,
                       Type *recv_data, int recv_count, int root)
{
#ifdef USE_MPI
	MPI_Gather(send_data, send_count, CommTrace<Type>::MPIType,
	           recv_data, recv_count, CommTrace<Type>::MPIType, root, comm);
#endif
}

template <class Type>
void
FSCommunicator::gatherv(Type *send_data, int send_count,
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
FSCommunicator::reduce(int num, Type* data, int root)
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
		MPI_Reduce(data+offset, buffer, count, CommTrace<Type>::MPIType, MPI_SUM, root, comm);
		if(thisCPU == root) for(int i = 0; i < count; ++i) data[offset+i] = buffer[i];
	}
	if(segSize > _MAX_ALLOCA_SIZE)
		delete [] buffer;
#endif
}

template <class Type>
void
FSCommunicator::broadcast(int num, Type* data, int root)
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
		if(thisCPU == root) for(int i = 0; i < count; i++) buffer[i] = data[offset+i];
		MPI_Bcast(buffer, count, CommTrace<Type>::MPIType, root, comm);
		if(thisCPU != root) for(int i = 0; i < count; ++i) data[offset+i] = buffer[i];
	}
	if(segSize > _MAX_ALLOCA_SIZE)
		delete [] buffer;
#endif
}


#ifdef NO_COMPLEX_MPI
// PJSA 1-7-2008 specializations of communication functions for platforms which do not support MPI_COMPLEX_DOUBLE
// implemented in Driver.d/MPIComm.C
template <>
complex<double> FSCommunicator::globalSum(complex<double> data);

template <>
complex<double> FSCommunicator::globalMax(complex<double> data);

template <>
complex<double> FSCommunicator::globalMin(complex<double> data);

template <>
void FSCommunicator::globalSum(int num, complex<double>* data);

template <>
void FSCommunicator::globalMax(int num, complex<double>* data);

template <>
void FSCommunicator::globalMin(int num, complex<double>* data);

template <>
void FSCommunicator::sendTo(int cpu, int tag, complex<double> *buffer, int len);

template <>
FSRecInfo FSCommunicator::recFrom(int tag, complex<double> *buffer, int len);

template <>
void FSCommunicator::allGather(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_count);

template <>
void FSCommunicator::allGatherv(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_counts[], int displacements[]);

template <>
void
FSCommunicator::reduce(int num, complex<double> *data, int root);

template <>
void
FSCommunicator::broadcast(int num, complex<double> *data, int root);
#endif

