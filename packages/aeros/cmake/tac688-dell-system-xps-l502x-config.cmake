# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH
    /home/tac688/eigen
    /home/tac688/trilinos-10.10.2-Source/packages/sacado/src
    /home/tac688/trilinos-10.10.2-Obj_cmake/include
    /home/tac688/trilinos-10.10.2-Source/packages/zoltan/src
#    /home/tac688/stxxl-trunk/include)
    )
SET(CMAKE_LIBRARY_PATH
#    /home/tac688/Codes/ARPACK
#    /home/tac688/stxxl-trunk/lib
    /home/tac688/trilinos-10.10.2-Obj_cmake/lib)
SET(EXTRALIB_MPI /usr/lib/libmpif77.so
                 CACHE STRING "Extra MPI link parameters")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
SET(BLAS_blas_LIBRARY "/home/tac688/eigen_build/blas/libeigen_blas.so" CACHE FILEPATH "Path to a library.")
