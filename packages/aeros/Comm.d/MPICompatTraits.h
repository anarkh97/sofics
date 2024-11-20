//
// Created by Michel Lesoinne on 2/13/18.
//

#ifndef FEM_MPICOMPATTRAITS_H
#define FEM_MPICOMPATTRAITS_H

#ifdef USE_MPI

#include <mpi.h>
template <>
struct CommTypeCompatibility<MPI_Win, HandleType::Window> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<MPI_Comm, HandleType::Communicator> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<MPI_Datatype , HandleType::Type> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<MPI_Op , HandleType::Op> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<MPI_Request, HandleType::Request> {
	static const bool isCompatible = true;
};
#else
struct WinDetails {
	void *data;
	int size;
	int disp_unit;
};

template <>
struct CommTypeCompatibility<WinDetails, HandleType::Window> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<int, HandleType::Communicator> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<size_t , HandleType::Type> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<int , HandleType::Op> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<int, HandleType::Request> {
	static const bool isCompatible = true;
};
#endif


#endif //FEM_MPICOMPATTRAITS_H
