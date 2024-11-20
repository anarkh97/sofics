# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH
    /home/avery/Codes/eigen
    /home/avery/Codes/MUMPS_4.10.0/include
    /home/avery/Codes/SuperLU/SuperLU_4.3/SRC
    /home/avery/Codes/trilinos/trilinos-11.4.3-Source/packages/sacado/src
    /home/avery/Codes/trilinos/trilinos-11.4.3-Obj_cmake/include
    /home/avery/Codes/trilinos/trilinos-11.4.3-Source/packages/zoltan/src)
SET(CMAKE_LIBRARY_PATH
    /home/avery/Codes/ARPACK
    /home/avery/Codes/MUMPS_4.10.0/lib
    /home/avery/Codes/SuperLU/SuperLU_4.3/lib
    /home/avery/Codes/SuperLU/SuperLU_4.3/lib
    /home/avery/Codes/trilinos/trilinos-11.4.3-Obj_cmake/lib)
SET(EXTRALIB_MPI /usr/lib/libmpif77.so
                 CACHE STRING "Extra MPI link parameters")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
SET(BLAS_blas_LIBRARY "/home/avery/Codes/eigen-build/blas/libeigen_blas.so" CACHE FILEPATH "Path to a library.")
add_definitions(-D_AEROS_ASYNCHRONOUS_IO)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix")
