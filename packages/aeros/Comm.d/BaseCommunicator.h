//
// Created by Michel Lesoinne on 2/5/18.
//

#ifndef FEM_OPAQUECOMMUNICATOR_H
#define FEM_OPAQUECOMMUNICATOR_H


#include <type_traits>
#include <vector>
#include <algorithm>
#include "OpaqueHandle.h"

class BaseCommunicator;

struct RecDetails {
	int length;
	int source;
};

/** \brief Structure to obtain information about the underlying implementation of the standard. */
struct VersionInfo {
	int major;
	int minor;
};

namespace com_details {

struct CommFunctions {
	int (*commSize)(const CommunicatorHandle &handle) = nullptr;
	int (*remoteSize)(const CommunicatorHandle &handle){};
	void (*allReduce)(const void *sendbuf, void *recvbuf, int count,
	                  TypeHandle datatype, OpHandle op,
	                  const CommunicatorHandle &comm,
	                  void(*cpData)(const void *, void *, int)){};
	void (*exScan)(const void *sendbuf, void *recvbuf, int count,
	               TypeHandle datatype, OpHandle op,
	               const CommunicatorHandle &comm){};
	void (*barrier)(const CommunicatorHandle &){};
	int (*rank)(const CommunicatorHandle &){};
	void (*blockingSend)(const void *sendbuf, int count, TypeHandle datatype,
	                     int dest, int tag,
	                     const CommunicatorHandle &){};
	/** \brief Receive from any source
	 *
	 * @param buffer
	 * @param len
	 * @param datatype
	 * @param tag
	 * @return The information about the length received and the origin of the data.
	 */
	RecDetails (*blockingRec)(void *buffer, int len, TypeHandle datatype,
	                          int tag,
	                          const CommunicatorHandle &){};
	void (*allGather)(const void *send_data, int send_count,
	                  void *recv_data, int recv_count, TypeHandle datatype,
	                  const CommunicatorHandle &){};
	/** \brief Create a Window.
	 *
	 * @param data
	 * @param size
	 * @param disp_unit
	 * @return
	 */
	WinHandle (*createWindow)(void *data, int size, int disp_unit, const CommunicatorHandle &){};
	void (*destroyWindow)(WinHandle handle){};

	void (*fence)(bool openOrClose, WinHandle handle){};

	void (*lock)(bool isShared, int remoteRank, WinHandle handle){};
	void (*unlock)(int remoteRank, WinHandle handle){};
	void (*flushRemote)(int remoteRank, WinHandle handle){};
	void (*flushLocal)(int remoteRank, WinHandle handle){};
	void (*lockAll)(WinHandle handle){};
	void (*unlockAll)(WinHandle handle){};
	void (*flushRemoteAll)(WinHandle handle){};
	void (*flushLocalAll)(WinHandle handle){};
	void (*fetchAndOp)(WinHandle handle, OpHandle op, const void *sourceData, TypeHandle datatype,
	                   void *resData, int remoteRank, int remoteOffset){};
	void (*accumulate)(WinHandle handle, OpHandle op, const void *operand, int count, TypeHandle datatype,
	                   int remoteRank, int remoteOffset){};
	void (*put)(WinHandle handle, const void *sourceData, int count, TypeHandle datatype,
	            int remoteRank, int remoteOffset){};
	void (*get)(WinHandle handle, void *resultData, int count, TypeHandle datatype,
	            int remoteRank, int remoteOffset){};
	ReqHandle (*i_send)(const void *buf, int count, TypeHandle datatype, RankHandle dest, int tag,
	                    const CommunicatorHandle &comm){};
	ReqHandle (*i_receive)(void *buf, int count, TypeHandle datatype, RankHandle origin, int tag,
	                       const CommunicatorHandle &comm){};
};

struct Constants {
	static WinHandle nullWindow;
};

}


class Window {
public:
	class LockGuard;
	Window(const Window &) = delete;
	Window(Window &&w);
	~Window();

	/** \brief Fence operation opening the window to RMA communication. */
	void open() const;
	/** \brief Fence operation closing the window to RMA communication. */
	void close() const;

	void sharedLock(int remoteRank) const;
	void unlock(int remoteRank) const;

	void flushRemote(int remoteRank) const;
	void flushLocal(int remoteRank) const;

	void sharedLockAll() const;
	void unlockAll() const;
	void flushRemote() const;
	void flushLocal() const;

	template <typename T>
	void fetchAndOp(OpHandle op, const T *operand, T*remoteOperandResult, int remoteRank, int remoteOffset) const;

	template <typename T>
	void accumulate(OpHandle op, const T *operand, int count, int remoteRank, int remoteOffset) const;

	template <typename T>
	void put(const T* source, int count, int remoteRank, int remoteOffset) const;

	template <typename T>
	void get(T *destMemory, int count, int remoteRank, int remoteOffset) const;
private:
	/// \brief Constructor of a window.
	Window(WinHandle winHandle) : winHandle(winHandle) {}
	friend class BaseCommunicator;
	WinHandle winHandle;
};

struct RequestInfo {
	int length;
	int cancelled;
	int origin;
	int tag;
	int error;
};

///** \brief A Type erased vector of POD types. */
//class TEVector {
//public:
//	~TEVector() { delete [] d; }
//
//	size_t size() const { return numElem; }
//
//	template <typename T>
//	T &push_back(const T &t) {
//		// TODO Check we're not mixing types
//		elemSize = sizeof(T);
//		if(numElem+1 > capacity)
//			realloc(numElem+1);
//		return at(++numElem);
//	}
//
//	template <typename T>
//	T &at(size_t index) {
//		return *(data<T>()+index);
//	}
//
//	template <typename T>
//	T *data() {
//		return reinterpret_cast<T *>(d);
//	}
//private:
//	unsigned char *d = nullptr;
//	size_t capacity = 0;
//	size_t numElem = 0;
//	size_t elemSize = 0;
//
//	// TODO for non POD types, we need to have move/destructor methods.
//	void realloc(size_t newSize);
//};

class RequestVector {
public:
	RequestVector();
	~RequestVector();
	void emplace_back(ReqHandle handle);


	RequestInfo waitAny();
	std::vector<RequestInfo> waitAll();
private:
	using vec_model = std::vector<RequestInfo>;
	std::aligned_storage<sizeof(vec_model), alignof(vec_model)>::type requests;
	std::vector<TypeHandle> types;

	template <typename T>
	std::vector<T> &getRequests();
};



class Window::LockGuard {
	LockGuard(const Window &window, int remoteRank = -1) : window(window), remoteRank(remoteRank) {
		if(remoteRank >= 0)
			this->window.sharedLock(this->remoteRank);
		else
			this->window.sharedLockAll();
	}
	LockGuard(LockGuard &&);
	LockGuard(const LockGuard &) = delete;
	~LockGuard(){
		try {
			if(remoteRank >= 0)
				window.unlock(remoteRank);
			else
				window.unlockAll();
		} catch(...){}
	}
	void flushRemote() const {
		window.flushRemote(remoteRank);
	}
	void flushLocal() const {
		window.flushLocal();
	}
private:
	const Window &window;
	const int remoteRank; //!< Remote rank for a shared lock, -1 for an 'all' lock.
};

class BaseCommunicator {
public:
	template <typename T>
	BaseCommunicator(const T &t);

	BaseCommunicator(const CommunicatorHandle &handle);

	/** \brief Obtain the rank of the process in the communicator. */
	int rank() const { return (*functions.rank)(handle); }
	/** \brief Obtain the number of processes involved in the communicator. */
	int commSize() const { return (*functions.commSize)(handle); }
	/** \brief For heterogeneous communicator, obtain the number of process in the other side. */
	int remoteSize() const { return (*functions.remoteSize)(handle); }
	/** \brief Block until all the processes in the communicator have reached the barrier. */
	void barrier() const { (*functions.barrier)(handle); }

	Window window(void *d, int nBytes, int disp_unit) const;

	template <typename T>
	Window window(const T *data, int count) const {
		const void * dd = data;
		return window(const_cast<void *>(dd), count * sizeof(T), sizeof(T));
	}

	template <typename T>
	void allReduce(const T *sendbuf, T *recvbuf, int count,
	               OpHandle op) const {
		(*functions.allReduce)(sendbuf, recvbuf, count, CommTypeTrait<T>::typeHandle(), op, handle,
		                       [](const void *from, void *to, int count) {
			                       const T*f = reinterpret_cast<const T *>(from);
			                       T*t = reinterpret_cast<T *>(to);
			                       std::copy(f, f+count, t);
		                       });
	}

	template <typename T>
	void exScan(const T *sendbuf, T *recvbuf, int count,
	            OpHandle op) const {
		(*functions.exScan)(sendbuf, recvbuf, count, CommTypeTrait<T>::typeHandle(), op, handle);
	}

	template <typename T>
	void blockingSend(const T *sendbuf, int count,
	                  int dest, int tag) const {
		(*functions.blockingSend)(sendbuf, count, CommTypeTrait<T>::typeHandle(),
		                          dest, tag, handle);
	}

	template <typename T>
	RecDetails blockingRec(int tag, T *buffer, int len) const {
		return (*functions.blockingRec)(buffer, len, CommTypeTrait<T>::typeHandle(), tag, handle);
	}

	/** \brief Initiate a non-blocking send operation.
	 *
	 * @tparam T Type of data being sent
	 * @param buffer Buffer in which the data being sent can be found.
	 * @param count Number of elements to be sent.
	 * @param destination Rank of the destination process.
	 * @param tag Tag associated with the communcation.
	 * @return A handle that can be used to check if the buffer can be modified.
	 */
	template <typename T>
	ReqHandle nonBlockingSend(const T *buffer, int count, RankHandle destination, int tag) const {
		return (*functions.i_send)(buffer, count, CommTypeTrait<T>::typeHandle(), destination, tag, handle);
	}


	template <typename T>
	ReqHandle nonBlockingReceive(T *buffer, int bufferLength, RankHandle origin, int tag) const {
		return (*functions.i_receive)(buffer, bufferLength, CommTypeTrait<T>::typeHandle(), origin, tag, handle);
	}

	template <typename T>
	void allGather(const T *send_data, int send_count,
	               T *recv_data, int recv_count) {
		(*functions.allGather)(send_data, send_count, recv_data, recv_count,  CommTypeTrait<T>::typeHandle(), handle);
	}

	template <class Type>
	Type globalSum(Type) const;
	template <class Type>
	Type globalMax(Type) const;
	template <class Type>
	Type globalMin(Type) const;
	/** \brief Exclusive scan. */
	template <typename Type>
	Type exScan(Type) const;

private:
	CommunicatorHandle handle;

	friend class Window;
	static com_details::CommFunctions functions;
};

template <typename T>
T BaseCommunicator::globalSum(T t) const {
	T res;
	allReduce(&t, &res, 1, SumHandle);
	return res;
}

template <typename T>
T BaseCommunicator::globalMax(T t) const {
	T res;
	allReduce(&t, &res, 1, MaxHandle);
	return res;
}


template <typename T>
T BaseCommunicator::globalMin(T t) const {
	T res;
	allReduce(&t, &res, 1, MinHandle);
	return res;
}


template<typename Type>
Type BaseCommunicator::exScan(Type t) const {
	Type res;
	exScan(&t, &res, 1, SumHandle);
	return res;
}

inline
void Window::open() const {
	(*BaseCommunicator::functions.fence)(true, winHandle);
}

inline
void Window::close() const {
	(*BaseCommunicator::functions.fence)(false, winHandle);
}

inline Window::Window(Window &&wnd) : winHandle(wnd.winHandle) {
wnd.winHandle = com_details::Constants::nullWindow;
}

inline Window::~Window() {
	// if(winHandle != com_details::Constants::nullWindow)
	(*BaseCommunicator::functions.destroyWindow)(winHandle);
}

inline
void Window::sharedLock(int remoteRank) const {
	(*BaseCommunicator::functions.lock)(true, remoteRank, winHandle);
}

inline
void Window::unlock(int remoteRank) const {
	(*BaseCommunicator::functions.unlock)(remoteRank, winHandle);
}

inline
void Window::flushRemote(int remoteRank) const {
	(*BaseCommunicator::functions.flushRemote)(remoteRank, winHandle);
}

inline
void Window::flushLocal(int remoteRank) const {
	(*BaseCommunicator::functions.flushLocal)(remoteRank, winHandle);
}

inline
void Window::sharedLockAll() const {
	(*BaseCommunicator::functions.lockAll)(winHandle);
}

inline
void Window::unlockAll() const {
	(*BaseCommunicator::functions.unlockAll)(winHandle);
}

inline
void Window::flushRemote() const {
	(*BaseCommunicator::functions.flushRemoteAll)(winHandle);
}

inline
void Window::flushLocal() const {
	(*BaseCommunicator::functions.flushLocalAll)(winHandle);
}

template<typename T>
void Window::fetchAndOp(OpHandle op, const T *operand, T *remoteOperandResult, int remoteRank, int remoteOffset) const {
	(*BaseCommunicator::functions.fetchAndOp)(winHandle, op, operand, CommTypeTrait<T>::typeHandle(),
	                                          remoteOperandResult, remoteRank, remoteOffset);
}

template<typename T>
void Window::accumulate(OpHandle op, const T *operand, int count,
                        int remoteRank, int remoteOffset) const {
	(*BaseCommunicator::functions.accumulate)(winHandle, op, operand,
	                                          count*CommTypeTrait<T>::count(*operand),
	                                          CommTypeTrait<T>::typeHandle(),
	                                          remoteRank, remoteOffset);
}

template<typename T>
void Window::put(const T *source, int count,
                 int remoteRank, int remoteOffset) const {
	(*BaseCommunicator::functions.put)(winHandle, source,
	                                   count*CommTypeTrait<T>::count(*source),
	                                   CommTypeTrait<T>::typeHandle(),
	                                   remoteRank, remoteOffset);
}

template<typename T>
void Window::get(T *operand, int count,
                 int remoteRank, int remoteOffset) const {
	(*BaseCommunicator::functions.get)(winHandle,
	                                   operand, count*CommTypeTrait<T>::count(*operand),
	                                   CommTypeTrait<T>::typeHandle(),
	                                   remoteRank, remoteOffset);
}

#endif //FEM_OPAQUECOMMUNICATOR_H
