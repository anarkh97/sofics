//
// Created by Michel Lesoinne on 2/5/18.
//


#include <complex>
#include "OpaqueHandle.h"
#include "MPICompatTraits.h"

#ifdef USE_MPI
#include <mpi.h>
OpHandle MaxHandle{(MPI_Op)MPI_MAX};
OpHandle MinHandle{(MPI_Op)MPI_MIN};
OpHandle SumHandle{(MPI_Op)MPI_SUM};
OpHandle ProdHandle{(MPI_Op)MPI_PROD};
TypeHandle IntHandle{(MPI_Datatype)MPI_INT};

template <>
TypeHandle CommTypeTrait<double>::typeHandle() {
	return TypeHandle{(MPI_Datatype)MPI_DOUBLE};
}

template <>
TypeHandle CommTypeTrait<char>::typeHandle() {
	return TypeHandle{(MPI_Datatype)MPI_CHAR};
}

template <>
TypeHandle CommTypeTrait<int>::typeHandle() {
	return TypeHandle{(MPI_Datatype)MPI_INT};
}

template <>
TypeHandle CommTypeTrait<long>::typeHandle() {
	return TypeHandle{(MPI_Datatype)MPI_LONG};
}

template <>
TypeHandle CommTypeTrait<size_t>::typeHandle() {
	static_assert(sizeof(size_t) == sizeof(unsigned long),
		"size_t is not compatible with unsigned long");
	return TypeHandle{(MPI_Datatype)MPI_UNSIGNED_LONG};
}

template <>
TypeHandle CommTypeTrait<std::complex<double>>::typeHandle() {
	return TypeHandle{(MPI_Datatype)MPI_DOUBLE_COMPLEX};
}

#else

OpHandle MaxHandle{1};
OpHandle MinHandle{2};
OpHandle SumHandle{3};
OpHandle ProdHandle{4};

template <>
TypeHandle CommTypeTrait<double>::typeHandle() {
    return TypeHandle{sizeof(double)};
}

template <>
TypeHandle CommTypeTrait<char>::typeHandle() {
    return TypeHandle{sizeof(char)};
}

template <>
TypeHandle CommTypeTrait<int>::typeHandle() {
    return TypeHandle{sizeof(int)};
}

template <>
TypeHandle CommTypeTrait<long>::typeHandle() {
    return TypeHandle{sizeof(long)};
}

template <>
TypeHandle CommTypeTrait<std::complex<double>>::typeHandle() {
    return TypeHandle{sizeof(std::complex<double>)};
}
#endif // USE_MPI

RankHandle RankHandle::Any{-1};
