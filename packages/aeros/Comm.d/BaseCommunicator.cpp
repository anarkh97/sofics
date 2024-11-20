//
// Created by Michel Lesoinne on 2/5/18.
//

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdexcept>
#include <cstring>
#include "BaseCommunicator.h"
#include "MPICompatTraits.h"

BaseCommunicator::BaseCommunicator(const CommunicatorHandle &handle) : handle(handle) {}

class mpi_except : public std::runtime_error {
public:
	explicit mpi_except(int error_code) : std::runtime_error("MPI Error") {}
};

class not_implemented_except : public std::runtime_error {
public:
	explicit not_implemented_except(const char *reason) : std::runtime_error(reason) {}
};

	using WH = OpaqueTypedHandle<HandleType::Window>;
using CH = OpaqueTypedHandle<HandleType::Communicator>;




namespace {

void _allGather(const void *send_data, int send_count,
                void *recv_data, int recv_count, TypeHandle datatype,
                const CommunicatorHandle &comm) {
#ifdef USE_MPI
	MPI_Allgather(send_data, send_count, datatype,
	              recv_data, recv_count, datatype, comm);
#else
	size_t sz = datatype;
	memcpy(recv_data, send_data, send_count*sz);
#endif
};

WinHandle _createWindow(void *data, int size, int disp_unit, const CommunicatorHandle &comm) {
#ifdef USE_MPI
	MPI_Win winHandle;
	int err_code = MPI_Win_create(data, size, disp_unit, MPI_INFO_NULL, comm, &winHandle);
	if (err_code != MPI_SUCCESS)
		throw mpi_except(err_code);
	return WinHandle{winHandle};
#else
	static_assert(CommTypeCompatibility<WinDetails, HandleType::Window>::isCompatible, "Not compatible");
	return WinHandle{WinDetails{data,size,disp_unit}};
#endif
};

void _fence(bool openOrClose, WinHandle handle) {
#ifdef USE_MPI
	int err_code = MPI_Win_fence(openOrClose ? MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE : MPI_MODE_NOSUCCEED,
	                             handle);
	if (err_code != MPI_SUCCESS)
		throw mpi_except(err_code);
#else
	// Nothing to do
#endif
}

void _destroyWindow(WinHandle handle) {
#ifdef USE_MPI
	MPI_Win mpiHandle = handle;
	if (mpiHandle != MPI_WIN_NULL) {
		int err_code = MPI_Win_free(&mpiHandle);
		if (err_code != MPI_SUCCESS)
			throw mpi_except(err_code);
	}
#else
	// Nothing to do
#endif
}

#ifdef USE_MPI
static
void check(int err_code) {
	if (err_code != MPI_SUCCESS)
		throw mpi_except(err_code);
}
#endif

void _lock(bool isShared, int remoteRank, WinHandle handle) {
#ifdef USE_MPI
	check(
			MPI_Win_lock(isShared ? MPI_LOCK_SHARED : MPI_LOCK_EXCLUSIVE,  remoteRank, 0, handle)
	);
#else
	// Nothing to do
#endif
}

void _unlock(int remoteRank, WinHandle handle) {
#ifdef USE_MPI
	check(
			MPI_Win_unlock(remoteRank, handle)
	);
#else
	// Nothing to do
#endif
}

void _lockAll(WinHandle handle) {
#ifdef USE_MPI
	check(
			MPI_Win_lock_all(0,  handle)
	);
#else
	// Nothing to do
#endif
}

void _unlockAll(WinHandle handle) {
#ifdef USE_MPI
	check(
			MPI_Win_unlock_all(handle)
	);
#else
	// Nothing to do
#endif
}

void _flushRemote(int remoteRank, WinHandle handle) {
#ifdef USE_MPI
	check(
			MPI_Win_flush(remoteRank, handle)
	);
#else
	// Nothing to do
#endif
}

void _flushLocal(int remoteRank, WinHandle handle) {
#ifdef USE_MMPI
	check(
			MPI_Win_flush_local(remoteRank, handle)
	);
#else
	// Nothing to do
#endif
}

void _flushRemoteAll(WinHandle handle) {
#ifdef USE_MPI
	check(
			MPI_Win_flush_all(handle)
	);
#else
	// Nothing to do
#endif
}

void _flushLocalAll(WinHandle handle) {
#ifdef USE_MPI
	check(
			MPI_Win_flush_local_all(handle)
	);
#else
	// Nothing to do
#endif
}

void _fetchAndOp(WinHandle handle, OpHandle op, const void *sourceData, TypeHandle dataType,
                   void *resData, int remoteRank, int remoteOffset) {
#ifdef USE_MPI
	check(
			MPI_Fetch_and_op(const_cast<void*>(sourceData), resData, dataType, remoteRank, remoteOffset, op, handle)
	);
#else
	// Get the Window data details

	// Fetch the data from the window and put it in resData

	// Perform the operation on the window with the source data.
	throw not_implemented_except("Fetch and Op not implemented without MPI.");
#endif
}

void _accumulate(WinHandle handle, OpHandle op, const void *operand, int count, TypeHandle dataType,
                   int remoteRank, int remoteOffset) {
#ifdef USE_MPI
	check(
			MPI_Accumulate(operand, count, dataType, remoteRank, remoteOffset, count, dataType, op, handle)
	);
#else
	// Get the Window data details
	// Perform the operation on the window with the source data.
	if(count > 0)
		throw not_implemented_except("Accumulate not implemented without MPI.");
#endif
}

void _put(WinHandle handle, const void *sourceData, int count, TypeHandle datatype,
          int remoteRank, int remoteOffset) {
#ifdef USE_MPI
	check(
			MPI_Put(sourceData, count, datatype, remoteRank, remoteOffset, count, datatype, handle)
	);
#else
	// Get the Window data details
	// Perform the operation on the window with the source data.
	if(count > 0)
		throw not_implemented_except("Put not implemented without MPI.");
#endif
}

void _get(WinHandle handle, void *destData, int count, TypeHandle datatype,
          int remoteRank, int remoteOffset) {
#ifdef USE_MPI
	check(
			MPI_Get(destData, count, datatype, remoteRank, remoteOffset, count, datatype, handle)
	);
#else
	// Get the Window data details
	// Perform the operation on the window with the source data.
	if(count > 0)
		throw not_implemented_except("Get not implemented without MPI.");
#endif
}

ReqHandle _i_send(const void *buf, int count, TypeHandle datatype, RankHandle dest, int tag,
             const CommunicatorHandle &comm) {
#ifdef USE_MPI
	MPI_Request request;
	check(
		MPI_Isend(buf, count, datatype, dest.rank, tag, comm, &request)
	);
	return { OpaqueTypedHandle<HandleType::Request>{request}, datatype };
#else
	// Get the Window data details
	// Perform the operation on the window with the source data.
	if(count > 0)
		throw not_implemented_except("i_send not implemented without MPI.");
	return { OpaqueTypedHandle<HandleType::Request>{0}, datatype };
#endif
}

ReqHandle _i_receive(void *buf, int count, TypeHandle datatype, RankHandle source, int tag,
                      const CommunicatorHandle &comm) {
#ifdef USE_MPI
	MPI_Request request;
	int origin = source.rank  == RankHandle::Any.rank ? MPI_ANY_SOURCE : source.rank;
	             check(
			MPI_Irecv(buf, count, datatype, origin, tag, comm, &request)
	);
	return { OpaqueTypedHandle<HandleType::Request>{request}, datatype };
#else
	return ReqHandle{ OpaqueTypedHandle<HandleType::Request>{0}, datatype };
#endif
}

}

com_details::CommFunctions getFunctionPointers() {
	com_details::CommFunctions commFunctions{};
	commFunctions.commSize = [](const CommunicatorHandle &handle) {
#ifdef USE_MPI
		int size;
		int result = MPI_Comm_size(handle, &size);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
		return size;
#else
		return 1;
#endif
	};
	commFunctions.remoteSize = [](const CommunicatorHandle &handle) {
#ifdef USE_MPI
		int size;
		int result = MPI_Comm_remote_size(handle, &size);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
		return size;
#else
		return 0;
#endif
	};
	commFunctions.allReduce = [](const void *sendbuf, void *recvbuf, int count, TypeHandle datatype, OpHandle op,
	                      const CommunicatorHandle &comm, void(*cpData)(const void *, void *, int)) {
#ifdef USE_MPI
		int result = MPI_Allreduce(sendbuf, recvbuf, count,
		                           datatype,
		                           op,
		                           comm);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
#else
		size_t sz = datatype;
		memcpy(recvbuf, sendbuf, count*sz);
#endif
	};
	commFunctions.exScan = [](const void *sendbuf, void *recvbuf, int count, TypeHandle datatype, OpHandle op,
	                   const CommunicatorHandle &comm) {
#ifdef USE_MPI
		int result = MPI_Exscan(sendbuf, recvbuf, count,
		                        datatype,
		                        op,
		                        comm);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
#else
		size_t sz = datatype;
		memcpy(recvbuf, sendbuf, count*sz);
#endif
	};
	commFunctions.barrier = [](const CommunicatorHandle &comm) {
#ifdef USE_MPI
		int result = MPI_Barrier(comm);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
#endif // No barrier in a non mpi run.
	};
	commFunctions.rank = [](const CommunicatorHandle &handle) {
#ifdef USE_MPI
		int rank;
		int result = MPI_Comm_rank(handle, &rank);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
		return rank;
#else
		return 0;
#endif
	};
	commFunctions.blockingSend = [](const void *sendbuf, int count, TypeHandle datatype, int dest, int tag,
	                         const CommunicatorHandle &comm) {
#ifdef USE_MPI
		int result = MPI_Send(sendbuf, count, datatype, dest, tag, comm);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
#else
		if(count > 0)
			throw not_implemented_except("i_send not implemented without MPI.");
#endif
	};
	commFunctions.blockingRec = [](void *buffer, int len, TypeHandle datatype, int tag,
	                        const CommunicatorHandle &comm) {
#ifdef USE_MPI
		RecDetails rInfo;
		MPI_Status status;
		MPI_Recv(buffer, len,
		         datatype, MPI_ANY_SOURCE, tag, comm, &status);
		MPI_Get_count(&status, datatype, &rInfo.length);
		rInfo.source = status.MPI_SOURCE;
		return rInfo;
#else
		RecDetails rInfo{0,0};
		return rInfo;
#endif
	};
	commFunctions.allGather = &_allGather;
	commFunctions.createWindow = &_createWindow;
	commFunctions.destroyWindow = &_destroyWindow;
	commFunctions.fence = &_fence;
	commFunctions.lock = &_lock;
	commFunctions.unlock = &_unlock;
	commFunctions.lockAll = &_lockAll;
	commFunctions.unlockAll = &_unlockAll;
	commFunctions.flushRemote = &_flushRemote;
	commFunctions.flushLocal = &_flushLocal;
	commFunctions.flushRemoteAll = &_flushRemoteAll;
	commFunctions.flushLocalAll = &_flushLocalAll;
	commFunctions.fetchAndOp = &_fetchAndOp;
	commFunctions.accumulate = &_accumulate;
	commFunctions.get = &_get;
	commFunctions.put = &_put;
	commFunctions.i_send = &_i_send;
	commFunctions.i_receive = &_i_receive;

	return commFunctions;
}

com_details::CommFunctions BaseCommunicator::functions = getFunctionPointers();
#ifdef USE_MPI
Window BaseCommunicator::window(void *d, int nBytes, int disp_unit) const {
	MPI_Win win;
	MPI_Win_create(d, nBytes, disp_unit, MPI_INFO_NULL, handle, &win);
	return Window(WinHandle{win});
}

template<>
BaseCommunicator::BaseCommunicator(const MPI_Comm &c) : handle(c) {}

template<typename T>
std::vector<T> &RequestVector::getRequests() {
	static_assert(sizeof(requests) >= sizeof(std::vector<T>), "This system requires some more work.");
	return *reinterpret_cast<std::vector<T> *>(&requests);
}

void RequestVector::emplace_back(ReqHandle handle) {
	getRequests<MPI_Request>().emplace_back(handle.request);
	types.push_back(handle.type);
}

RequestInfo rInfoFrom(const MPI_Status &status, MPI_Datatype datatype) {
	RequestInfo requestInfo;
	MPI_Get_count(&status, datatype, &requestInfo.length);
	requestInfo.tag = status.MPI_TAG;
	requestInfo.origin = status.MPI_SOURCE;
	requestInfo.error = status.MPI_ERROR != MPI_SUCCESS;
	return requestInfo;
}

RequestInfo RequestVector::waitAny() {
	MPI_Status status;
	auto &requests = getRequests<MPI_Request>();
	int index;
	check(
			MPI_Waitany(requests.size(), requests.data(),
			            &index, &status)
	);
	if(status.MPI_ERROR != MPI_SUCCESS)
		throw mpi_except(status.MPI_ERROR);
	return rInfoFrom(status, this->types[index]);
}

std::vector<RequestInfo> RequestVector::waitAll() {
	auto &requests = getRequests<MPI_Request>();
	std::vector<MPI_Status> statuses(requests.size());
	check(
			MPI_Waitall(requests.size(), requests.data(), statuses.data())
	);
	std::vector<RequestInfo> result(requests.size());
	for(int i = 0; i < requests.size(); ++i) {
		if(statuses[i].MPI_ERROR != MPI_SUCCESS)
			throw mpi_except(statuses[i].MPI_ERROR);
		result[i] = rInfoFrom(statuses[i], types[i]);
	}
	return result;
}

RequestVector::RequestVector() {
	new (&requests) std::vector<MPI_Request>;
}

RequestVector::~RequestVector() {
	getRequests<MPI_Request>().~vector();
}

CommunicatorHandle getWorldComm() {
	return CommunicatorHandle((MPI_Comm)MPI_COMM_WORLD);
}

namespace com_details {

WinHandle Constants::nullWindow{(MPI_Win)MPI_WIN_NULL};

}

#else
CommunicatorHandle getWorldComm() {
    return CommunicatorHandle{0};
}
#endif
//#else
//
//com_details::CommFunctions BaseCommunicator::functions{
//	.commSize = [](const OpaqueHandle &handle) {
//		return 1;
//	},
//	.remoteSize = [](const OpaqueHandle &handle) {
//		return 0;
//	},
//	.allReduce = [](const void *sendbuf, void *recvbuf, int count, TypeHandle datatype, OpHandle op,
//		                CommunicatorHandle comm, void(*cpData)(const void *, void *, int)) {
//		(*cpData)(sendbuf, recvbuf, count);
//	},
//	.barrier = [](const CommunicatorHandle &comm) {
//	},
//	.rank = [](const OpaqueHandle &handle) {
//		return 0;
//	}
//};
//
//namespace com_details {
//
//template <>
//TypeHandle CommTypeTrait<int>::typeHandle() {
//	return MPI_INT;
//}
//
//template <>
//TypeHandle CommTypeTrait<double>::typeHandle() {
//	return MPI_DOUBLE;
//}
//
////
////template struct CommTypeTrait<int>;
////template struct CommTypeTrait<double>;
////template struct CommTypeTrait<std::complex<double>>;
//
//}
//#endif

