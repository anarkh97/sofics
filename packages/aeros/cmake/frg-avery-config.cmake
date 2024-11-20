# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /home/pavery/hg/eigen
                       /home/pavery/code/SPOOLES
                       /home/pavery/code/MUMPS/MUMPS_4.10.0/include
                       /home/pavery/code/Zoltan
                       /home/pavery/code/Zoltan/include)
SET(CMAKE_LIBRARY_PATH /home/pavery/code/ARPACK
                       /home/pavery/code/SPOOLES
                       /home/pavery/code/SPOOLES/MT/src
                       /home/pavery/code/BLACS/LIB
                       /home/pavery/code/SCALAPACK
                       /home/pavery/code/MUMPS/MUMPS_4.10.0/lib
                       /home/pavery/code/MUMPS/MUMPS_4.10.0/libseq
                       /home/pavery/code/ParMETIS3_1
                       /home/pavery/code/Zoltan/Obj_generic)
SET(BLAS_LIBRARIES /home/pavery/code/LAPACK/lapack-3.2.1/libblas.a)
SET(LAPACK_LIBRARIES /home/pavery/code/LAPACK/lapack-3.2.1/liblapack.a ${BLAS_LIBRARIES})
SET(LAPACK_FOUND true)
# add anything missed by cmake auto-detection
SET(EXTRALIB -lgfortran CACHE STRING "Extra MPI link parameters")
#SET(EXTRALIB /opt/intel/fce/10.1.015/lib/libifcore.so 
#             /opt/intel/fce/10.1.015/lib/libsvml.so
#             CACHE STRING "Extra link parameters")
#SET(EXTRALIB_MPI /usr/mpi/gcc/openmpi-1.2.6/lib64/libmpi_f77.so
#                 CACHE STRING "Extra MPI link parameters")
#
#SET(BLA_VENDOR Intel10_64lp)
