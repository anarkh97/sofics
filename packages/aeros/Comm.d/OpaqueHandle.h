//
// Created by Michel Lesoinne on 2/5/18.
//

#ifndef FEM_OPAQUEHANDLE_H
#define FEM_OPAQUEHANDLE_H

#include <memory>
#include <vector>
#include <type_traits>
#include <complex>

namespace com_details {
struct HandleHelper;
}

enum class HandleType {
	Communicator,
	Window,
	Type,
	Op,
	Request
};

struct RankHandle {
	static RankHandle Any;
	RankHandle(int rank) : rank(rank) {}
	int rank;
};

template <typename T, HandleType ht>
struct CommTypeCompatibility {
	static const bool isCompatible = false;
};

/** \brief Opaque handle for any implementation specific MPI communicator handle. */
class OpaqueHandle {
public:
	template <typename T>
	OpaqueHandle(const T h) {
		static_assert(sizeof(T) <= ls, "OpaqueHandle cannot store such a large object.");
		unsigned char *place = reinterpret_cast<unsigned char *>(&handle);
		T *res = new (place) T(h);
	}
private:
	static const size_t ls = sizeof(void *);
	static const size_t la = alignof(void *);
	std::aligned_storage_t<ls, la> handle;
	friend struct com_details::HandleHelper;
};

/** \brief Opaque handle for any implementation specific MPI handle.
 *
 * @tparam ht The type of handled contained.
 */
template <HandleType ht>
struct OpaqueTypedHandle {
	template <typename T, typename X = typename std::enable_if<CommTypeCompatibility<T, ht>::isCompatible, void>::type>
	explicit OpaqueTypedHandle(const T &h){
		static_assert(sizeof(T) < ls, "OpaqueHandle cannot store such a large object.");
		unsigned char *place = reinterpret_cast<unsigned char *>(&handle);
		T *res = new (place) T(h);
	}
	template <typename T, typename X = typename std::enable_if<CommTypeCompatibility<T, ht>::isCompatible, void>::type>
	operator T() const;
private:
	static const size_t ls = 32;//sizeof(void *);
	static const size_t la = alignof(void *);
	std::aligned_storage_t<ls, la> handle;
	friend struct com_details::HandleHelper;
};

template<HandleType ht>
template<typename T, typename X>
OpaqueTypedHandle<ht>::operator T() const {
	return *reinterpret_cast<const T *>(&handle);
}

using CommunicatorHandle = OpaqueTypedHandle<HandleType::Communicator>;
using TypeHandle = OpaqueTypedHandle<HandleType::Type>;
using OpHandle = OpaqueTypedHandle<HandleType::Op>;
using WinHandle = OpaqueTypedHandle<HandleType::Window>;
struct ReqHandle {
	OpaqueTypedHandle<HandleType::Request> request;
	TypeHandle type;
};

extern OpHandle MaxHandle;
extern OpHandle MinHandle;
extern OpHandle SumHandle;
extern OpHandle ProdHandle;
extern TypeHandle IntHandle;

CommunicatorHandle getWorldComm();

template <typename T>
struct CommTypeTrait {
	static TypeHandle typeHandle();
	static int count(const T &) { return 1; }
	static constexpr bool isFixedSize() { return true; }
};

template <typename T>
struct CommTypeTrait<std::pair<T, T>> {
	static TypeHandle typeHandle() {
		static_assert(CommTypeTrait<T>::isFixedSize(),
		              "Communication operations cannot handle pairs of variable size objects.");
		return CommTypeTrait<T>::typeHandle();
	}
	static int count(const std::pair<T, T> &) { return 2; }
	static constexpr bool isFixedSize() { return CommTypeTrait<T>::isFixedSize(); }
};

template <typename T>
struct CommTypeTrait<std::vector<T>> {
	static TypeHandle typeHandle() { return CommTypeTrait<T>::typeHandle(); }
	static int count(const T &v) { return v.size(); }
	static constexpr bool isFixedSize() { return false; }
};

/*
template <typename T>
struct CommTypeTrait<std::complex<T>> {
	static TypeHandle typeHandle() { return CommTypeTrait<T>::typeHandle(); }
	static int count(const T &v) { return 2; }
	static constexpr bool isFixedSize() { return true; }
};
*/

#endif //FEM_OPAQUEHANDLE_H
