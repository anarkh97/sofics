PROJECT (FEM)

cmake_minimum_required(VERSION 3.1  FATAL_ERROR)
enable_language(Fortran)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
	cmake_policy(SET CMP0076 NEW)
endif()

option(USE_MPI "Use MPI" ON)
option(USE_CHOLMOD "Attempt to use Choldmod" ON)
option(USE_PARDISO "Attempt to use MKL Pardiso" ON)
option(TRY_ACME "Try to use ACME" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake" "/usr/local/include/eigen3/cmake/")

set(WITH_MUMPS true CACHE BOOL "Whether to use mumps or not")

# Conan approach
if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
	message("Using conan for most libraries (in build directory -- prefered method")
	include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
	set(USE_CONAN true)
	set(CONAN_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/build/conanbuildinfo.cmake) #Clion, with conanbuildinfo.cmake in build folder
	message("Using in-source conan for most libraries")
	include(${CMAKE_CURRENT_SOURCE_DIR}/build/conanbuildinfo.cmake)
	set(USE_CONAN true)
	set(CONAN_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build)
else()
	message("Using FindXXX to locate all libraries")
endif()
if(USE_CONAN)
	set(MPI_C_COMPILER ${CONAN_BUILD_DIR}/bin/mpicc)
	set(MPI_CXX_COMPILER ${CONAN_BUILD_DIR}/bin/mpic++)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix")
endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

if(USE_MARCH_NATIVE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" AND NOT WITH_CXX11 EQUAL 0)
	execute_process(
			COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
		set(CHECK_CXX_FEATURES true)
		set(CMAKE_CXX11_FLAGS "-Wno-deprecated-declarations")
	endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel" AND NOT WITH_CXX11 EQUAL 0)
	execute_process(
			COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE ICC_VERSION)
	if (ICC_VERSION VERSION_GREATER 12.1 OR ICC_VERSION VERSION_EQUAL 12.1)
		set(CHECK_CXX_FEATURES true)
		set(CMAKE_CXX11_FLAGS "-wd1478")
		add_definitions(-DHACK_INTEL_COMPILER_ITS_CPP11)
	endif ()
endif ()


SITE_NAME(HOSTNAME)
FIND_PACKAGE(LOCAL NAMES ${HOSTNAME} PATHS . QUIET)
MARK_AS_ADVANCED(HOSTNAME LOCAL_DIR)

OPTION(WITH_OPENMP "Use OpenMP" ON)
IF(WITH_OPENMP)
	FIND_PACKAGE(OpenMP)
	if( NOT (OpenMP_FLAG_DETECTED OR OpenMP_CXX_FOUND) AND APPLE )
		message( "Mac OS did not detect OpenMP. Forcing approach ")
		set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp" )
	endif()
	if(OpenMP_FLAG_DETECTED OR OpenMP_CXX_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif(OpenMP_FLAG_DETECTED OR OpenMP_CXX_FOUND)
	if(OpenMP_FLAG_DETECTED OR OpenMP_C_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	endif(OpenMP_FLAG_DETECTED OR OpenMP_C_FOUND)
	if(OpenMP_FLAG_DETECTED OR OpenMP_Fortran_FOUND)
		set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
	endif(OpenMP_FLAG_DETECTED OR OpenMP_Fortran_FOUND)
ELSE(WITH_OPENMP)
	unset(OpenMP_FLAG_DETECTED)
	unset(OpenMP_FLAG_DETECTED CACHE)
ENDIF(WITH_OPENMP)

if(NOT BLAS_FOUND)
	find_package(BLAS)
endif(NOT BLAS_FOUND)
if(NOT LAPACK_FOUND)
	FIND_PACKAGE(LAPACK)
endif(NOT LAPACK_FOUND)

SET(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(BINDIR ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(LIBDIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(EXTRALIB ${EXTRALIB} CACHE STRING "Extra link parameters")
SET(EXEC_DIR ${CMAKE_CURRENT_BINARY_DIR}/Regression.d)
SET(REG_SRC_DIR ${PROJECT_SOURCE_DIR}/Regression.d)
SET(BASELINE $EXEC_DIR/baseline)

if(USE_MPI)
	FIND_PACKAGE(MPI)
	SET(MPI_FOUND ${MPI_CXX_FOUND}) # We only care about C++
else(USE_MPI)
	set(MPI_FOUND false)
endif(USE_MPI)

if(AERO)
	MARK_AS_ADVANCED(CLEAR MPI_COMPILER)
	MARK_AS_ADVANCED(CLEAR MPI_LIBRARY)
	FIND_PACKAGE(MPI REQUIRED)
endif(AERO)

if(ASYNCHRONOUS_IO)
	add_definitions(-D_AEROS_ASYNCHRONOUS_IO)
endif(ASYNCHRONOUS_IO)

IF(UNIX)
	IF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
		add_definitions(-D__LINUX)
	ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)

	IF(${CMAKE_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
		SET(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")

		# zoltan
		IF(MPI_FOUND AND NOT WITH_ZOLTAN EQUAL 0)
			FIND_PATH(ZOLTAN_INCLUDE_PATH Zoltan_config.h
					PATHS /usr/include /usr/local/include /usr/include/trilinos)
			FIND_PATH(ZOLTAN_SRC_PATH zz/zz_const.h
					PATHS /usr/include /usr/local/include)
			FIND_LIBRARY(ZOLTAN_LIBRARY NAMES zoltan trilinos_zoltan
					PATHS /usr/lib /usr/local/lib)
			IF(ZOLTAN_INCLUDE_PATH AND ZOLTAN_SRC_PATH AND ZOLTAN_LIBRARY)
				MESSAGE("-- A library with Zoltan API found.")
				SET(ZOLTAN_FOUND TRUE)
				SET(EXTRALIBRARIES ${ZOLTAN_LIBRARY} ${EXTRALIBRARIES})
				SET(INC_DIR ${INC_DIR} ${ZOLTAN_INCLUDE_PATH}
						${ZOLTAN_SRC_PATH}/zz ${ZOLTAN_SRC_PATH}/lb ${ZOLTAN_SRC_PATH}/Utilities/shared
						${ZOLTAN_SRC_PATH}/Utilities/Timer ${ZOLTAN_SRC_PATH}/par ${ZOLTAN_SRC_PATH}/params
						${ZOLTAN_SRC_PATH}/order ${ZOLTAN_SRC_PATH}/tpls)
			ENDIF(ZOLTAN_INCLUDE_PATH AND ZOLTAN_SRC_PATH AND ZOLTAN_LIBRARY)
		ENDIF(MPI_FOUND AND NOT WITH_ZOLTAN EQUAL 0)

		# acme
		IF(TRY_ACME)
			IF(NOT MPI_FOUND)
				SET(ACME_FOUND TRUE)
				add_definitions(-DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0)
				SET(INC_DIR ${INC_DIR} "Acme.d/search" "Acme.d/enforcement")
			ELSE(NOT MPI_FOUND)
				IF(ZOLTAN_FOUND)
					SET(ACME_FOUND TRUE)
					add_definitions(-DUSE_ACME -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0)
					SET(INC_DIR ${INC_DIR} "Acme.d/search" "Acme.d/enforcement")
				ENDIF(ZOLTAN_FOUND)
			ENDIF(NOT MPI_FOUND)
		ENDIF(TRY_ACME)

		# arpack
		IF(NOT WITH_ARPACK EQUAL 0)
			FIND_LIBRARY(ARPACK_arpack_LIBRARY NAMES arpack_LINUX arpack
					PATHS /usr/lib /usr/local/lib)
			IF(NOT MPI_FOUND)
				MARK_AS_ADVANCED(FORCE ARPACK_parpack_LIBRARY)
				IF(ARPACK_arpack_LIBRARY)
					MESSAGE("-- A library with ARPACK API found.")
					SET(ARPACK_FOUND TRUE)
					add_definitions(-DUSE_ARPACK)
					SET(EXTRALIBRARIES ${ARPACK_arpack_LIBRARY} ${EXTRALIBRARIES})
				ENDIF(ARPACK_arpack_LIBRARY)
			ELSE(NOT MPI_FOUND)
				MARK_AS_ADVANCED(CLEAR ARPACK_parpack_LIBRARY)
				FIND_LIBRARY(ARPACK_parpack_LIBRARY NAMES parpack_MPI-LINUX parpack
						PATHS /usr/lib /usr/local/lib)
				IF(ARPACK_arpack_LIBRARY AND ARPACK_parpack_LIBRARY)
					MESSAGE("-- A library with ARPACK API found.")
					SET(ARPACK_FOUND TRUE)
					add_definitions(-DUSE_ARPACK)
					SET(EXTRALIBRARIES ${ARPACK_parpack_LIBRARY} ${ARPACK_arpack_LIBRARY} ${EXTRALIBRARIES})
				ENDIF(ARPACK_arpack_LIBRARY AND ARPACK_parpack_LIBRARY)
			ENDIF(NOT MPI_FOUND)
		ELSE(NOT WITH_ARPACK EQUAL 0)
			SET(ARPACK_arpack_LIBRARY "ARPACK_arpack_LIBRARY-NOTFOUND" CACHE FILEPATH "Path to a library." FORCE)
			SET(ARPACK_parpack_LIBRARY "ARPACK_parpack_LIBRARY-NOTFOUND" CACHE FILEPATH "Path to a library." FORCE)
		ENDIF(NOT WITH_ARPACK EQUAL 0)

		# spooles
		IF(NOT WITH_SPOOLES EQUAL 0)
			FIND_PATH(SPOOLES_INCLUDE_PATH FrontMtx.h
					PATHS /usr/include /usr/local/include /usr/include/spooles)
			FIND_LIBRARY(SPOOLES_spooles_LIBRARY NAMES spooles
					PATHS /usr/lib /usr/local/lib)
			FIND_LIBRARY(SPOOLES_spoolesMT_LIBRARY NAMES spoolesMT
					PATHS /usr/lib /usr/local/lib)
			if(SPOOLES_INCLUDE_PATH AND SPOOLES_spooles_LIBRARY)
				MESSAGE("-- A library with SPOOLES API found.")
				SET(SPOOLES_FOUND TRUE)
			endif(SPOOLES_INCLUDE_PATH AND SPOOLES_spooles_LIBRARY)
			IF(SPOOLES_FOUND)
				add_definitions(-DUSE_SPOOLES)
				if(SPOOLES_spoolesMT_LIBRARY)
					SET(EXTRALIBRARIES ${SPOOLES_spoolesMT_LIBRARY} ${SPOOLES_spooles_LIBRARY} ${EXTRALIBRARIES})
				else(SPOOLES_spoolesMT_LIBRARY)
					SET(EXTRALIBRARIES ${SPOOLES_spooles_LIBRARY} ${EXTRALIBRARIES})
				endif(SPOOLES_spoolesMT_LIBRARY)
				SET(INC_DIR ${INC_DIR} ${SPOOLES_INCLUDE_PATH})
			ENDIF(SPOOLES_FOUND)
		ELSE(NOT WITH_SPOOLES EQUAL 0)
			SET(SPOOLES_spooles_LIBRARY "SPOOLES_spooles_LIBRARY-NOTFOUND" CACHE FILEPATH "Path to a library." FORCE)
		ENDIF(NOT WITH_SPOOLES EQUAL 0)

		find_package(Eigen3 3.3 REQUIRED NO_MODULES)
		set(WITH_EIGEN3_EXTRA 0)
		#    include_directories(${EIGEN3_INCLUDE_DIR})
		add_definitions(-DUSE_EIGEN3 -DEIGEN_QUATERNIONBASE_PLUGIN=<Element.d/MpcElement.d/QuaternionBasePlugin.h>
				-DEIGEN_DONT_PARALLELIZE)
		set(CMAKE_MODULE_PATH "${EIGEN3_INCLUDE_DIR}/cmake" ${CMAKE_MODULE_PATH})
		IF(NOT WITH_EIGEN3_EXTRA EQUAL 0)
			find_package(Cholmod)
			if (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
				include_directories(${CHOLMOD_INCLUDES})
				add_definitions(-DEIGEN_CHOLMOD_SUPPORT)
				SET(EXTRALIBRARIES ${CHOLMOD_LIBRARIES} ${EXTRALIBRARIES})
			endif(CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
			find_package(Umfpack)
			if (UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)
				include_directories(${UMFPACK_INCLUDES})
				add_definitions(-DEIGEN_UMFPACK_SUPPORT)
				SET(EXTRALIBRARIES ${UMFPACK_LIBRARIES} ${EXTRALIBRARIES})
			endif(UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)
			find_package(SuperLU)
			if (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
				include_directories(${SUPERLU_INCLUDES})
				add_definitions(-DEIGEN_SUPERLU_SUPPORT)
				SET(EXTRALIBRARIES ${SUPERLU_LIBRARIES} ${EXTRALIBRARIES})
			endif(SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
			find_package(SPQR)
			if (SPQR_INCLUDES AND SPQR_LIBRARIES)
				include_directories(${SPQR_INCLUDES})
				add_definitions(-DEIGEN_SPQR_SUPPORT)
				set(THREADS_PREFER_PTHREAD_FLAG ON)
				find_package(Threads REQUIRED)
				SET(EXTRALIBRARIES ${SPQR_LIBRARIES} ${EXTRALIBRARIES} Threads::Threads)
			endif(SPQR_INCLUDES AND SPQR_LIBRARIES)
		ENDIF(NOT WITH_EIGEN3_EXTRA EQUAL 0)

		# sacado
		IF(NOT WITH_SACADO EQUAL 0)
			FIND_PATH(SACADO_INCLUDE_PATH NAMES Sacado.hpp
				PATHS /usr/include /usr/local/include /usr/include/trilinos)
			IF(SACADO_INCLUDE_PATH)
				set(SACADO_FOUND TRUE)
				MESSAGE("-- A library with SACADO API found. ")
				add_definitions(-DTRILINOS_NO_CONFIG_H -DUSE_SACADO)
				include_directories(${SACADO_INCLUDE_PATH})
				FIND_PATH(SACADO_MPL_INCLUDE_PATH NAMES Sacado_mpl_vector.hpp
					PATHS ${SACADO_INCLUDE_PATH} ${SACADO_INCLUDE_PATH}/mpl)
				if(SACADO_MPL_INCLUDE_PATH AND (NOT (SACADO_INCLUDE_PATH EQUAL SACADO_MPL_INCLUDE_PATH)))
					include_directories(${SACADO_MPL_INCLUDE_PATH})
				endif(SACADO_MPL_INCLUDE_PATH AND (NOT (SACADO_INCLUDE_PATH EQUAL SACADO_MPL_INCLUDE_PATH)))
				if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
					if (GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0)
						MESSAGE("Note: Sacado source code may need to be changed if Aero-S does not build "
								"(1) Sacado_trad.hpp line 657 (2) Sacado_tradvec.hpp line 468. "
								"In each case the argument to the IndepADvar assignment operator should be const.")
					endif (GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0)
				endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
			ENDIF(SACADO_INCLUDE_PATH)
		ENDIF(NOT WITH_SACADO EQUAL 0)

		# metis
		IF(NOT WITH_METIS EQUAL 0)
			FIND_PATH(METIS_INCLUDE_PATH NAMES metis.h
				PATHS /usr/include /usr/local/include $ENV{METISPATH}/include)
			FIND_LIBRARY(METIS_LIBRARY NAMES metis
					PATHS /usr/lib /usr/local/lib $ENV{METISPATH}/lib)
			IF(METIS_LIBRARY AND METIS_INCLUDE_PATH)
				set(METIS_FOUND TRUE)
				MESSAGE("-- A library with METIS API found.")
				add_definitions(-DUSE_METIS)
				include_directories(${METIS_INCLUDE_PATH})
				SET(EXTRALIBRARIES ${METIS_LIBRARY} ${EXTRALIBRARIES})
			ENDIF(METIS_LIBRARY AND METIS_INCLUDE_PATH)
		ENDIF(NOT WITH_METIS EQUAL 0)

		# scalapack
		IF(MPI_FOUND AND NOT WITH_SCALAPACK EQUAL 0)
			# blacs is required for scalapack
			FIND_LIBRARY(BLACS_blacs_LIBRARY NAMES blacs blacs-openmpi blacs_MPI-LINUX-0
					PATHS /usr/lib /usr/local/lib)
			FIND_LIBRARY(BLACS_Cinit_LIBRARY NAMES blacsCinit blacsCinit-openmpi blacsCinit_MPI-LINUX-0 blacsC
					PATHS /usr/lib /usr/local/lib)
			if(BLACS_blacs_LIBRARY AND BLACS_Cinit_LIBRARY)
				MESSAGE("-- A library with BLACS API found.")
				SET(BLACS_LIBRARIES ${BLACS_blacs_LIBRARY} ${BLACS_Cinit_LIBRARY} ${BLACS_blacs_LIBRARY})
				SET(BLACS_FOUND TRUE)
			endif(BLACS_blacs_LIBRARY AND BLACS_Cinit_LIBRARY)

			FIND_LIBRARY(SCALAPACK_LIBRARY NAMES scalapack scalapack-openmpi
					PATHS /usr/lib /usr/local/lib)
			if(SCALAPACK_LIBRARY)
				MESSAGE("-- A library with SCALAPACK API found.")
				SET(SCALAPACK_FOUND TRUE)
			endif(SCALAPACK_LIBRARY)
			if(BLACS_FOUND AND SCALAPACK_FOUND)
				add_definitions(-DUSE_SCALAPACK)
				SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${SCALAPACK_LIBRARY} ${BLACS_LIBRARIES})
				#        IF(NOT WITH_MKL EQUAL 0)
				#            MESSAGE("-- Using MKL.")
				#            add_definitions(-DUSE_MKL)
				#            include_directories($ENV{MKLPATH}/../../include)
				#        ELSE(NOT WITH_MKL EQUAL 0)
				#            include_directories("/home/avery/Codes/Scalapack/trunk/PBLAS/SRC")
				#        ENDIF(NOT WITH_MKL EQUAL 0)
			endif(BLACS_FOUND AND SCALAPACK_FOUND)
		ENDIF(MPI_FOUND AND NOT WITH_SCALAPACK EQUAL 0)

		# scotch
		IF(NOT WITH_SCOTCH EQUAL 0)
			IF(NOT SCOTCH_scotch.h_DIRS OR NOT SCOTCH_scotch_LIBRARY OR NOT SCOTCH_scotcherrexit_LIBRARY)
				find_package(Scotch)
			ELSE(NOT SCOTCH_scotch.h_DIRS OR NOT SCOTCH_scotch_LIBRARY OR NOT SCOTCH_scotcherrexit_LIBRARY)
				SET(SCOTCH_INCLUDE_DIRS ${SCOTCH_scotch.h_DIRS})
				SET(SCOTCH_LIBRARIES ${SCOTCH_scotch_LIBRARY} ${SCOTCH_scotcherrexit_LIBRARY})
			ENDIF(NOT SCOTCH_scotch.h_DIRS OR NOT SCOTCH_scotch_LIBRARY OR NOT SCOTCH_scotcherrexit_LIBRARY)
			IF(SCOTCH_LIBRARIES AND SCOTCH_INCLUDE_DIRS)
				set(SCOTCH_FOUND TRUE)
				MESSAGE("-- A library with SCOTCH API found.")
				add_definitions(-DUSE_SCOTCH)
				SET(EXTRALIBRARIES ${SCOTCH_LIBRARIES} ${EXTRALIBRARIES})
			include_directories(${SCOTCH_INCLUDE_DIRS})
			ENDIF(SCOTCH_LIBRARIES AND SCOTCH_INCLUDE_DIRS)
		ENDIF(NOT WITH_SCOTCH EQUAL 0)

		# mumps
		IF(NOT WITH_MUMPS EQUAL 0)
			FIND_PATH(MUMPS_INCLUDE_PATH dmumps_c.h
					PATHS /usr/include usr/local/include)
			if (NOT MPI_FOUND)
				UNSET(MUMPS_common_LIBRARY CACHE)
				MARK_AS_ADVANCED(CLEAR MUMPS_dmumps_seq_LIBRARY MUMPS_zmumps_seq_LIBRARY
						MUMPS_common_seq_LIBRARY MUMPS_pord_seq_LIBRARY
						MUMPS_mpiseq_LIBRARY)
				MARK_AS_ADVANCED(FORCE MUMPS_dmumps_LIBRARY MUMPS_zmumps_LIBRARY
						MUMPS_common_LIBRARY MUMPS_pord_LIBRARY)
				FIND_LIBRARY(MUMPS_dmumps_seq_LIBRARY NAMES dmumps_seq SEQ/libdmumps.a
						PATHS /usr/lib /usr/local/lib)
				FIND_LIBRARY(MUMPS_zmumps_seq_LIBRARY NAMES zmumps_seq SEQ/libzmumps.a
						PATHS /usr/lib /usr/local/lib)
				FIND_LIBRARY(MUMPS_common_seq_LIBRARY NAMES mumps_common_seq SEQ/libmumps_common.a
						PATHS /usr/lib /usr/local/lib)
				FIND_LIBRARY(MUMPS_pord_seq_LIBRARY NAMES pord_seq SEQ/libpord.a
						PATHS /usr/lib /usr/local/lib)
				FIND_LIBRARY(MUMPS_mpiseq_LIBRARY NAMES mpiseq mpiseq_seq SEQ/libmpiseq.a
						PATHS /usr/lib /usr/local/lib)
				if(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_seq_LIBRARY AND MUMPS_zmumps_seq_LIBRARY AND MUMPS_common_seq_LIBRARY AND MUMPS_pord_seq_LIBRARY AND MUMPS_mpiseq_LIBRARY)
					MESSAGE("-- A library with MUMPS API found.")
					SET(MUMPS_LIBRARIES ${MUMPS_dmumps_seq_LIBRARY} ${MUMPS_zmumps_seq_LIBRARY} ${MUMPS_common_seq_LIBRARY} ${MUMPS_pord_seq_LIBRARY} ${MUMPS_mpiseq_LIBRARY})
					SET(MUMPS_FOUND TRUE)
				endif(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_seq_LIBRARY AND MUMPS_zmumps_seq_LIBRARY AND MUMPS_common_seq_LIBRARY AND MUMPS_pord_seq_LIBRARY AND MUMPS_mpiseq_LIBRARY)
				if(MUMPS_FOUND)
					add_definitions(-DUSE_MUMPS)
					SET(INC_DIR ${INC_DIR} ${MUMPS_INCLUDE_PATH})
					SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${MUMPS_LIBRARIES})
				endif(MUMPS_FOUND)
			else(NOT MPI_FOUND)
				UNSET(MUMPS_common_seq_LIBRARY CACHE)
				MARK_AS_ADVANCED(FORCE MUMPS_dmumps_seq_LIBRARY MUMPS_zmumps_seq_LIBRARY
						MUMPS_common_seq_LIBRARY MUMPS_pord_seq_LIBRARY
						MUMPS_mpiseq_LIBRARY)
				MARK_AS_ADVANCED(CLEAR MUMPS_dmumps_LIBRARY MUMPS_zmumps_LIBRARY
						MUMPS_common_LIBRARY MUMPS_pord_LIBRARY)
				FIND_LIBRARY(MUMPS_dmumps_LIBRARY NAMES dmumps
						PATHS /usr/lib /usr/local/lib)
				FIND_LIBRARY(MUMPS_zmumps_LIBRARY NAMES zmumps
						PATHS /usr/lib /usr/local/lib)
				FIND_LIBRARY(MUMPS_common_LIBRARY NAMES mumps_common
						PATHS /usr/lib /usr/local/lib)
				FIND_LIBRARY(MUMPS_pord_LIBRARY NAMES pord
						PATHS /usr/lib /usr/local/lib)
				if(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_LIBRARY AND MUMPS_zmumps_LIBRARY AND MUMPS_common_LIBRARY AND MUMPS_pord_LIBRARY)
					SET(MUMPS_LIBRARIES ${MUMPS_dmumps_LIBRARY} ${MUMPS_zmumps_LIBRARY} ${MUMPS_common_LIBRARY} ${MUMPS_pord_LIBRARY})
					MESSAGE("-- A library with MUMPS API found.")
					SET(MUMPS_FOUND TRUE)
				endif(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_LIBRARY AND MUMPS_zmumps_LIBRARY AND MUMPS_common_LIBRARY AND MUMPS_pord_LIBRARY)

				if(BLACS_FOUND AND SCALAPACK_FOUND AND MUMPS_FOUND)
					add_definitions(-DUSE_MUMPS)
					SET(INC_DIR ${INC_DIR} ${MUMPS_INCLUDE_PATH})
					SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${MUMPS_LIBRARIES} ${SCALAPACK_LIBRARY} ${BLACS_LIBRARIES})
				endif(BLACS_FOUND AND SCALAPACK_FOUND AND MUMPS_FOUND)

			endif(NOT MPI_FOUND)
		ELSE(NOT WITH_MUMPS EQUAL 0)
			if (NOT MPI_FOUND)
				UNSET(MUMPS_common_LIBRARY CACHE)
				SET(MUMPS_common_seq_LIBRARY "MUMPS_common_seq_LIBRARY-NOTFOUND" CACHE FILEPATH "Path to a library." FORCE)
			else (NOT MPI_FOUND)
				UNSET(MUMPS_common_seq_LIBRARY CACHE)
				SET(MUMPS_common_LIBRARY "MUMPS_common_LIBRARY-NOTFOUND" CACHE FILEPATH "Path to a library." FORCE)
			endif (NOT MPI_FOUND)
		ENDIF(NOT WITH_MUMPS EQUAL 0)

		IF(NOT WITH_BOOST EQUAL 0)
			FIND_PACKAGE(Boost 1.35)
			if(Boost_FOUND)
				add_definitions(-DUSE_BOOST)
				include_directories(${Boost_INCLUDE_DIRS})
			else(Boost_FOUND)
				unset(Boost_INCLUDE_DIR CACHE)
				unset(Boost_LIBRARY_DIRS CACHE)
			endif(Boost_FOUND)
		ENDIF(NOT WITH_BOOST EQUAL 0)

	ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
ENDIF(UNIX)

IF(APPLE)
	add_definitions(-DNO_MALLOC_DOT_H)
ENDIF(APPLE)

if(MPI_FOUND)
	add_definitions(-DUSE_MPI -DDISTRIBUTED -DMPICH_IGNORE_CXX_SEEK)
	SET(EXTRALIBDYN ${EXTRALIBRARIES})
	SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${MPI_LIBRARIES})
	SET(INC_DIR ${INC_DIR} ${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

OPTION(WITH_GIT "Detect git commit ID" ON)
IF(WITH_GIT)
  FIND_PROGRAM(GIT_EXECUTABLE git)
  IF(GIT_EXECUTABLE)
	ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/GitIdentity.C
			COMMAND echo "const char *THE_VERSION = \\\""
			`${GIT_EXECUTABLE} log --pretty=format:'%h' -n 1` "\\\"" "\;" > ${CMAKE_CURRENT_BINARY_DIR}/GitIdentity.C
			DEPENDS main.C
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	ADD_DEFINITIONS(-DPRINT_CHANGESETID)
	SET(IDENTITY_FILE "${CMAKE_CURRENT_BINARY_DIR}/GitIdentity.C")
  ENDIF(GIT_EXECUTABLE)
ELSE(WITH_GIT)
  FIND_PROGRAM(MERCURIAL_EXECUTABLE hg)
  IF(MERCURIAL_EXECUTABLE)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/HgIdentity.C
			COMMAND echo "const char *THE_VERSION = \\\"" `${MERCURIAL_EXECUTABLE} identify -i` "\\\"" "\;" > ${CMAKE_CURRENT_BINARY_DIR}/HgIdentity.C
			DEPENDS main.C
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	add_definitions(-DPRINT_CHANGESETID)
	SET(IDENTITY_FILE "${CMAKE_CURRENT_BINARY_DIR}/HgIdentity.C")
  ENDIF (MERCURIAL_EXECUTABLE)
ENDIF(WITH_GIT)

add_definitions(-D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DSOWER_SURFS -DMAP_MIN_MEMORY)

IF(AEROS_NO_AD EQUAL 1)
	add_definitions(-DAEROS_NO_AD)
endif(AEROS_NO_AD EQUAL 1)

INCLUDE_DIRECTORIES( ${INC_DIR} Extermal.d/include)

add_subdirectory(Element.d)
add_subdirectory(Feti.d)
add_subdirectory(Driver.d)
add_subdirectory(Comm.d)
add_subdirectory(Corotational.d)
add_subdirectory(Dec.d)
add_subdirectory(Solvers.d)
add_subdirectory(Utils.d)
add_subdirectory(Parser.d)
add_subdirectory(Timers.d)
add_subdirectory(Threads.d)
add_subdirectory(Mortar.d)
add_subdirectory(Math.d)
add_subdirectory(Linpack.d)
add_subdirectory(Sfem.d)
add_subdirectory(Paral.d)
add_subdirectory(Problems.d)
add_subdirectory(GNU-getopt.d)
add_subdirectory(Hetero.d)
add_subdirectory(Material.d)
add_subdirectory(Rom.d)
add_subdirectory(Regression.d)

add_subdirectory(FetiLib)

SET(ALLLIBS Parser
		Driver Elem Feti solver Corotational Mortar Rom Math
		Dec Threads linpa Sfem paral Problems solver Corotational Driver Hetero Timers Comm Material Util
		FDPSolver getopt)

IF(ACME_FOUND AND (ZOLTAN_FOUND OR NOT MPI_FOUND))
	add_subdirectory(Acme.d)
	SET(ALLLIBS ${ALLLIBS} Acme)
ENDIF(ACME_FOUND AND (ZOLTAN_FOUND OR NOT MPI_FOUND))

#IF(MPI_FOUND)
#  add_subdirectory(Pita.d)
#  SET(ALLLIBS Pita ${ALLLIBS} )
#ENDIF(MPI_FOUND)

SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${EXTRALIB})
if(MPI_FOUND)
	SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${EXTRALIB_MPI})
endif(MPI_FOUND)

IF(CREATE_DSO)
	add_definitions(-DCREATE_DSO)
	ADD_LIBRARY(femlib main.C ${IDENTITY_FILE})
	# ADD_LIBRARY(femlib EXCLUDE_FROM_ALL main.C ${IDENTITY_FILE})
	target_link_libraries(femlib Parser
			Driver Elem Feti solver Corotational Mortar Rom Math
			Dec Threads linpa Sfem paral Problems Driver Hetero Timers Comm Material Util getopt Parser
			${EXTRALIBDYN} ${CMAKE_DL_LIBS} ${LAPACK_LIBRARIES})
	SET_TARGET_PROPERTIES(femlib
			PROPERTIES
			ARCHIVE_OUTPUT_DIRECTORY ${LIBDIR}
			RELEASE_OUTPUT_NAME libaeros
			DEBUG_OUTPUT_NAME libaeros_debug
			LINK_SEARCH_END_STATIC 1
			PREFIX "")
endif(CREATE_DSO)

IF(NOT CREATE_DSO)
	add_executable(femexecutable  main.C ${IDENTITY_FILE} Types.h)

	SET_TARGET_PROPERTIES(femexecutable
			PROPERTIES
			RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
			RELEASE_OUTPUT_NAME aeros
			DEBUG_OUTPUT_NAME aeros.debug
			#                       LINK_SEARCH_END_STATIC 1
			)

	target_link_libraries(femexecutable ${ALLLIBS}
			${EXTRALIBRARIES} ${CMAKE_DL_LIBS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
ENDIF(NOT CREATE_DSO)

add_executable(robcodec Rom.d/RobCodec.C)
set_target_properties(robcodec
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME rob
		DEBUG_OUTPUT_NAME rob.debug)

target_link_libraries(robcodec Rom Util)

add_executable(xpostcodec Rom.d/XPostCodec.C)
set_target_properties(xpostcodec
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME xpo
		DEBUG_OUTPUT_NAME xpo.debug)

target_link_libraries(xpostcodec Rom Util)

add_executable(xpostcodec6 Rom.d/XPostCodec6.C)
set_target_properties(xpostcodec6
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME xpo6
		DEBUG_OUTPUT_NAME xpo6.debug)

target_link_libraries(xpostcodec6 Rom Util)

add_executable(errorcodec Rom.d/ErrorCodec.C)
set_target_properties(errorcodec
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME relerr
		DEBUG_OUTPUT_NAME relerr.debug)

target_link_libraries(errorcodec Rom Util)

add_executable(impulseerrorcodec Rom.d/ImpulseErrorCodec.C)
set_target_properties(impulseerrorcodec
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME imperr
		DEBUG_OUTPUT_NAME relerr.debug)

target_link_libraries(impulseerrorcodec Rom Util)

add_executable(errorcodec6 Rom.d/ErrorCodec6.C)
set_target_properties(errorcodec6
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME relerr6
		DEBUG_OUTPUT_NAME relerr6.debug)

target_link_libraries(errorcodec6 Rom Util)

add_executable(impulseerrorcodec6 Rom.d/ImpulseErrorCodec6.C)
set_target_properties(impulseerrorcodec6
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME imperr6
		DEBUG_OUTPUT_NAME imperr6.debug)

target_link_libraries(impulseerrorcodec6 Rom Util)

add_executable(errorcodec1 Rom.d/ErrorCodec1.C)
set_target_properties(errorcodec1
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
		RELEASE_OUTPUT_NAME relerr1
		DEBUG_OUTPUT_NAME relerr1.debug)

target_link_libraries(errorcodec1 Rom Util)

#-- Add an Option to toggle the generation of the API documentation
FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND)
	IF (DEFINED ENV{MATHJAX_URL})
		SET(MJU $ENV{MATHJAX_URL})
	ELSE (DEFINED ENV{MATHJAX_URL})
		SET(MJU http://cdn.mathjax.org/mathjax/latest)
		IF (NOT DEFINED MATHJAX_URL)
			MESSAGE(STATUS "MATHJAX_URL is set to default value. Use ccmake to change it.")
		ENDIF (NOT DEFINED MATHJAX_URL)
	ENDIF (DEFINED ENV{MATHJAX_URL})
	SET(MATHJAX_URL ${MJU} CACHE STRING "URL To Mathjax. Can be set as environment.")
	configure_file(Doxyfile.in
			${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
	#-- Add a custom target to run Doxygen when ever the project is built
	add_custom_target(Docs # ALL
			COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
			SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
	# IF you do NOT want the documentation to be generated EVERY time you build the project
	# then leave out the 'ALL' keyword from the above command.
else (DOXYGEN_FOUND)
	message(WARNING
			"Doxygen not found, documentation cannot be built.")
endif (DOXYGEN_FOUND)

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
			"Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
			FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(AllTargets robcodec errorcodec errorcodec6 errorcodec1)
IF(CREATE_DSO)
	SET(AllTargets ${AllTargets} femlib)
ELSE(CREATE_DSO)
	SET(AllTargets femexecutable ${AllTargets})
ENDIF(CREATE_DSO)

INSTALL(TARGETS ${AllTargets}
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib)


MESSAGE("")
MESSAGE("=================================================")
MESSAGE("           Summary of build options")
MESSAGE("-------------------------------------------------")
IF(MPI_FOUND)
	MESSAGE("MPI                       YES")
	MESSAGE("Distributed FETI:         YES")
	MESSAGE("Aeroelastic:              YES")
	IF(MUMPS_FOUND AND SCALAPACK_FOUND AND BLACS_FOUND)
		MESSAGE("Mumps:                    YES")
	ELSE(MUMPS_FOUND AND SCALAPACK_FOUND AND BLACS_FOUND)
		MESSAGE("Mumps:                    NO ")
	ENDIF(MUMPS_FOUND AND SCALAPACK_FOUND AND BLACS_FOUND)
	IF(SCALAPACK_FOUND AND BLACS_FOUND)
		MESSAGE("Scalapack:                YES")
	ELSE(SCALAPACK_FOUND AND BLACS_FOUND)
		MESSAGE("Scalapack:                NO ")
	ENDIF(SCALAPACK_FOUND AND BLACS_FOUND)
ELSE(MPI_FOUND)
	MESSAGE("MPI                       NO ")
	MESSAGE("Distributed FETI:         NO ")
	MESSAGE("Aeroelastic:              NO ")
	IF(MUMPS_FOUND)
		MESSAGE("Mumps:                    YES")
	ELSE(MUMPS_FOUND)
		MESSAGE("Mumps:                    NO ")
	ENDIF(MUMPS_FOUND)
	MESSAGE("Scalapack:                NO ")
ENDIF(MPI_FOUND)
IF(ARPACK_FOUND)
	MESSAGE("Arpack:                   YES")
ELSE(ARPACK_FOUND)
	MESSAGE("Arpack:                   NO ")
ENDIF(ARPACK_FOUND)
IF(SPOOLES_FOUND)
	MESSAGE("Spooles:                  YES")
ELSE(SPOOLES_FOUND)
	MESSAGE("Spooles:                  NO ")
ENDIF(SPOOLES_FOUND)
IF(CHOLMOD_FOUND)
	MESSAGE("Cholmod:                  YES")
ELSE(CHOLMOD_FOUND)
	MESSAGE("Cholmod:                  NO ")
ENDIF(CHOLMOD_FOUND)
IF(PARDISO_FOUND)
	MESSAGE("Pardiso:                  YES")
ELSE(PARDISO_FOUND)
	MESSAGE("Pardiso:                  NO ")
ENDIF(PARDISO_FOUND)
IF(ACME_FOUND)
	MESSAGE("Acme:                     YES")
ELSE(ACME_FOUND)
	MESSAGE("Acme:                     NO ")
ENDIF(ACME_FOUND)
IF(METIS_FOUND)
	MESSAGE("Metis:                    YES")
ELSE(METIS_FOUND)
	MESSAGE("Metis:                    NO ")
	MESSAGE("Use METISPATH environment variable if you have METIS.")
ENDIF(METIS_FOUND)
IF(EIGEN3_FOUND)
	MESSAGE("Eigen template library:   YES")
ELSE(EIGEN3_FOUND)
	MESSAGE("Eigen template library:   NO ")
ENDIF(EIGEN3_FOUND)
if(OpenMP_FLAG_DETECTED OR OpenMP_CXX_FOUND)
	MESSAGE("OpenMP:                   YES")
else(OpenMP_FLAG_DETECTED OR OpenMP_CXX_FOUND)
	MESSAGE("OpenMP:                   NO ")
endif(OpenMP_FLAG_DETECTED OR OpenMP_CXX_FOUND)
MESSAGE("Build type:               " ${CMAKE_BUILD_TYPE})
MESSAGE("Extra libraries:          " ${EXTRALIB})
MESSAGE("=================================================")
MESSAGE("")

FIND_PATH(BASELINE REQUIRED)

MARK_AS_ADVANCED(FORCE BISON_EXECUTABLE FLEX_EXECUTABLE MERCURIAL_EXECUTABLE CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_COMPILER CMAKE_Fortran_COMPILER)
ENABLE_TESTING()
SET (COPY_FILES ${REG_SRC_DIR}/copy_files.pl ${REG_SRC_DIR} ${EXEC_DIR} )
SET (TEST_SHORT ${EXEC_DIR}/dataComp.py -l short)
SET (BUILD_SHORT ${EXEC_DIR}/buildInputs.py ALL)

SET (TEST_DURATION_SHORT $ENV{TEST_DURATION_SHORT})

ADD_CUSTOM_TARGET(test-short
		COMMAND ${COPY_FILES}
		COMMAND ${BUILD_SHORT}
		COMMAND ${TEST_SHORT}
		)

if(TEST_DURATION_SHORT)
	ADD_TEST(copyfiles ${REG_SRC_DIR}/copy_files.pl ${REG_SRC_DIR} ${EXEC_DIR} )
	ADD_TEST(buildInputs ${EXEC_DIR}/buildInputs.py ALL)
	ADD_TEST(testStatics ${EXEC_DIR}/dataComp.py -l -n statics)
	ADD_TEST(testNLstatics ${EXEC_DIR}/dataComp.py -l nlstatics)
	ADD_TEST(testEigen ${EXEC_DIR}/dataComp.py -l eigen)
	ADD_TEST(testDynamics ${EXEC_DIR}/dataComp.py -l dynamics)
	ADD_TEST(testNLDynamics ${EXEC_DIR}/dataComp.py -l nldynamics)
	ADD_TEST(testIMPE ${EXEC_DIR}/dataComp.py -l impe)
	ADD_TEST(testTempStatics ${EXEC_DIR}/dataComp.py -l tempstatics)
	ADD_TEST(testTempNLStatics ${EXEC_DIR}/dataComp.py -l tempnlstatics)
	ADD_TEST(testTempDynamics ${EXEC_DIR}/dataComp.py -l tempdynamics)
	ADD_TEST(testTempNLDynamics ${EXEC_DIR}/dataComp.py -l tempnldynamics)
	ADD_TEST(testRtest1 ${EXEC_DIR}/dataComp.py -l dsvm1)
	ADD_TEST(testRtest20 ${EXEC_DIR}/dataComp.py -l dsvm20)
	ADD_TEST(testRtest21 ${EXEC_DIR}/dataComp.py -l dsvm21)
	ADD_TEST(testRtest22 ${EXEC_DIR}/dataComp.py -l dsvm22)
	ADD_TEST(testRtest23 ${EXEC_DIR}/dataComp.py -l dsvm23)
	ADD_TEST(testRtest25 ${EXEC_DIR}/dataComp.py -l dsvm25)
	ADD_TEST(testRtest27a ${EXEC_DIR}/dataComp.py -l dsvm27a)
	ADD_TEST(testRtest27b ${EXEC_DIR}/dataComp.py -l dsvm27b)
	ADD_TEST(testRtest30 ${EXEC_DIR}/dataComp.py -l dsvm30)
	if(EIGEN3_FOUND AND HAS_CXX11_TEMPLATE_ALIAS AND SACADO_FOUND)
		ADD_TEST(testRtest34 ${EXEC_DIR}/dataComp.py -l dsvm34)
	endif(EIGEN3_FOUND AND HAS_CXX11_TEMPLATE_ALIAS AND SACADO_FOUND)
	ADD_TEST(testRtest35b ${EXEC_DIR}/dataComp.py -l dsvm35b)
	ADD_TEST(testRtest40 ${EXEC_DIR}/dataComp.py -l dsvm40)
	ADD_TEST(testvme1 ${EXEC_DIR}/dataComp.py -l vme1 )
	ADD_TEST(testvme3 ${EXEC_DIR}/dataComp.py -l vme3 )
	ADD_TEST(testvme4 ${EXEC_DIR}/dataComp.py -l vme4 )
	ADD_TEST(testvme5 ${EXEC_DIR}/dataComp.py -l vme5 )
	ADD_TEST(testvme6 ${EXEC_DIR}/dataComp.py -l vme6 )
else(TEST_DURATION_SHORT)
	ADD_TEST(copyfiles ${REG_SRC_DIR}/copy_files.pl ${REG_SRC_DIR} ${EXEC_DIR} )
	ADD_TEST(buildInputs ${EXEC_DIR}/buildInputs.py ALL)
	ADD_TEST(testStatics ${EXEC_DIR}/dataComp.py -l -n statics)
	ADD_TEST(testNLstatics ${EXEC_DIR}/dataComp.py -l nlstatics)
	ADD_TEST(testEigen ${EXEC_DIR}/dataComp.py -l eigen)
	ADD_TEST(testDynamics ${EXEC_DIR}/dataComp.py -l dynamics)
	ADD_TEST(testNLDynamics ${EXEC_DIR}/dataComp.py -l nldynamics)
	ADD_TEST(testIMPE ${EXEC_DIR}/dataComp.py -l impe)
	ADD_TEST(testFreqSweep ${EXEC_DIR}/dataComp.py -l freqsweep)
	ADD_TEST(testTempStatics ${EXEC_DIR}/dataComp.py -l tempstatics)
	ADD_TEST(testTempNLStatics ${EXEC_DIR}/dataComp.py -l tempnlstatics)
	ADD_TEST(testTempDynamics ${EXEC_DIR}/dataComp.py -l tempdynamics)
	ADD_TEST(testTempNLDynamics ${EXEC_DIR}/dataComp.py -l tempnldynamics)
	ADD_TEST(testRtest1 ${EXEC_DIR}/dataComp.py -l dsvm1)
	ADD_TEST(testRtest2 ${EXEC_DIR}/dataComp.py -l dsvm2)
	ADD_TEST(testRtest11 ${EXEC_DIR}/dataComp.py -l dsvm11)
	ADD_TEST(testRtest13 ${EXEC_DIR}/dataComp.py -l dsvm13)
	#  ADD_TEST(testRtest15 ${EXEC_DIR}/dataComp.py -l dsvm15)
	ADD_TEST(testRtest19 ${EXEC_DIR}/dataComp.py -l dsvm19)
	ADD_TEST(testRtest20 ${EXEC_DIR}/dataComp.py -l dsvm20)
	ADD_TEST(testRtest21 ${EXEC_DIR}/dataComp.py -l dsvm21)
	ADD_TEST(testRtest22 ${EXEC_DIR}/dataComp.py -l dsvm22)
	ADD_TEST(testRtest23 ${EXEC_DIR}/dataComp.py -l dsvm23)
	ADD_TEST(testRtest24 ${EXEC_DIR}/dataComp.py -l dsvm24)
	ADD_TEST(testRtest25 ${EXEC_DIR}/dataComp.py -l dsvm25)
	ADD_TEST(testRtest27a ${EXEC_DIR}/dataComp.py -l dsvm27a)
	ADD_TEST(testRtest27b ${EXEC_DIR}/dataComp.py -l dsvm27b)
	ADD_TEST(testRtest29 ${EXEC_DIR}/dataComp.py -l dsvm29)
	ADD_TEST(testRtest30 ${EXEC_DIR}/dataComp.py -l dsvm30)
	ADD_TEST(testRtest31 ${EXEC_DIR}/dataComp.py -l dsvm31)
	ADD_TEST(testRtest32 ${EXEC_DIR}/dataComp.py -l dsvm32)
	if(EIGEN3_FOUND AND HAS_CXX11_TEMPLATE_ALIAS AND SACADO_FOUND)
		ADD_TEST(testRtest34 ${EXEC_DIR}/dataComp.py -l dsvm34)
	endif(EIGEN3_FOUND AND HAS_CXX11_TEMPLATE_ALIAS AND SACADO_FOUND)
	ADD_TEST(testRtest35a ${EXEC_DIR}/dataComp.py -l dsvm35a)
	ADD_TEST(testRtest35b ${EXEC_DIR}/dataComp.py -l dsvm35b)
	ADD_TEST(testRtest37 ${EXEC_DIR}/dataComp.py -l dsvm37)
	#  ADD_TEST(testRtest38 ${EXEC_DIR}/dataComp.py -l dsvm38)
	ADD_TEST(testRtest39 ${EXEC_DIR}/dataComp.py -l dsvm39)
	ADD_TEST(testRtest40 ${EXEC_DIR}/dataComp.py -l dsvm40)
	ADD_TEST(testvme1 ${EXEC_DIR}/dataComp.py -l vme1 )
	#  ADD_TEST(testvme2 ${EXEC_DIR}/dataComp.py -l vme2 )
	ADD_TEST(testvme3 ${EXEC_DIR}/dataComp.py -l vme3 )
	ADD_TEST(testvme4 ${EXEC_DIR}/dataComp.py -l vme4 )
	ADD_TEST(testvme5 ${EXEC_DIR}/dataComp.py -l vme5 )
	ADD_TEST(testvme6 ${EXEC_DIR}/dataComp.py -l vme6 )
	ADD_TEST(testvmmech003 ${EXEC_DIR}/dataComp.py -l vmmech003 )
	ADD_TEST(testvmmech063 ${EXEC_DIR}/dataComp.py -l vmmech063 )
	ADD_TEST(testPreStress ${EXEC_DIR}/dataComp.py -l PreStressedMembrane )
	# use the -s flag to send results via email
	ADD_TEST(testPressure ${EXEC_DIR}/dataComp.py -l PlateUnderPressure )
	#  ADD_TEST(sendResults ${EXEC_DIR}/sendResults.py )
endif(TEST_DURATION_SHORT)

