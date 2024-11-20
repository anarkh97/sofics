cmake_minimum_required(VERSION 3.5  FATAL_ERROR)
ADD_LIBRARY(Comm BaseCommunicator.cpp BaseCommunicator.h OpaqueHandle.cpp OpaqueHandle.h MPICompatTraits.h)
#
#SET(INC_DIR ..)
#
#INCLUDE_DIRECTORIES( ${INC_DIR} . )
TARGET_LINK_LIBRARIES(Comm ${MPI_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(Comm PUBLIC ${MPI_INCLUDE_PATH})

TARGET_SOURCES(Comm
        PRIVATE CommunicatorCore.C
        Communicator.h
        )

ADD_EXECUTABLE(comm_test CommTest.cpp)
#TARGET_SOURCES(comm_test PUBLIC )
MESSAGE("MPI Libraries: ${MPI_LIBRARIES}")
TARGET_LINK_LIBRARIES(comm_test Comm )