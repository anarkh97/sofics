# -----------------------------------------------------------------------------
# Find Required Dependencies
# -----------------------------------------------------------------------------

find_package(MPI REQUIRED COMPONENTS C CXX)
message(STATUS "Found MPI: ${MPI_CXX_COMPILER}")
find_package(BLAS REQUIRED)
message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
find_package(LAPACK REQUIRED)
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

# -----------------------------------------------------------------------------
# Eigen3 and Boost --- required by M2C and Aero-S
# -----------------------------------------------------------------------------
# Download Eigen if not found
find_package(Eigen3 3.4 QUIET NO_MODULE)

#if("${EIGEN3_INCLUDE_DIR}" STREQUAL "")
#  message(STATUS "Eigen3 not found, downloading ...")
#  include(FetchContent)
#  
#  FetchContent_Declare(
#      Eigen3
#      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
#      GIT_TAG 3.4.1
#      GIT_SHALLOW TRUE
#  )
#  
#  # Don't build tests/examples
#  set(BUILD_TESTING OFF)
#  set(EIGEN_BUILD_DOC OFF)
#  set(EIGEN_BUILD_TESTING OFF)
#  set(EIGEN_BUILD_PKGCONFIG OFF)
#  set(EIGEN_MPL2_ONLY ON)
#  
#  FetchContent_MakeAvailable(Eigen3)
#  message(STATUS "Eigen3 download successfull.")
#else()
#  message(STATUS "Found system Eigen3: ${EIGEN3_INCLUDE_DIR}")
#endif()

# Download Boost if not found
find_package(Boost 1.72 QUIET)

#if("${Boost_INCLUDE_DIRS}" STREQUAL "")
#  message(STATUS "Boost not found, downloading ...")
#  include(FetchContent)
#
#  #FetchContent_Declare(
#  #  Boost
#  #  GIT_REPOSITORY https://github.com/boostorg/boost.git
#  #  GIT_TAG boost-1.72.0
#  #  GIT_SHALLOW TRUE
#  #)
#
#  # TODO: Need to test with what is actually required by M2C
#  set(BOOST_INCLUDE_LIBRARIES asio regex algorithm)
#  set(BOOST_ENABLE_CMAKE ON)
#  FetchContent_Declare(
#    Boost
#    URL https://github.com/boostorg/boost/releases/download/boost-1.90.0/boost-1.90.0-b2-nodocs.7z
#    USES_TERMINAL_DOWNLOAD TRUE
#    DOWNLOAD_NO_EXTRACT FALSE
#  )
#  FetchContent_MakeAvailable(Boost)
#  message(STATUS "Boost download successfull.")
#else()
#  message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS}")
#endif()

# PETSc handled by M2C CMakeLists ... for now.
# MUMPS, METIS, SPOOLES, SCALAPACK handled by Aero-S CMakeLists ... for now.
