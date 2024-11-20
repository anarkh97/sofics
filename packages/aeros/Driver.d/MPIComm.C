#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Communicator.h>
#include <Comm.d/Communicator.h>
#include <complex>
#include <iostream>

/*
  this file implements all MPI communication functions.
  it defaults to empty functions in the case when USE_MPI is not defined
*/

void
initComm(int argc, char **argv)
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif
}

void
closeComm()
{
#ifdef USE_MPI
    MPI_Finalize();
#endif
}

FSCommunicator::FSCommunicator() : communicator(getWorldComm())
{
#ifdef USE_MPI
    comm = MPI_COMM_WORLD;
    numCPU = communicator.numCPUs();
    thisCPU = communicator.myID();
#else
    numCPU = 1;
    thisCPU = 0;
#endif
}

#ifdef USE_MPI
FSCommunicator::FSCommunicator(Communicator *c) : communicator(*c)
{
    comm = *c->getCommunicator();
    numCPU = communicator.numCPUs();
    thisCPU = communicator.myID();
}
#endif

void
FSCommunicator::exchange(int tag, int numNeighb, gsl::span<const int> cpus,
                         int *sndPtr, int *sndData, int *rcvPtr, int *rcvData) {
#ifdef USE_MPI
    int iCpu;
    MPI_Request *sndId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    MPI_Request *rcvId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    int rcvReq = 0;
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = rcvPtr[iCpu+1]-rcvPtr[iCpu];
        if(len == 0) continue;
        MPI_Irecv(rcvData+rcvPtr[iCpu], len, MPI_INT, cpus[iCpu],
                  tag, comm,  rcvId+rcvReq);
        rcvReq += 1;
    }
    // post send
    int sendReq = 0;
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = sndPtr[iCpu+1]-sndPtr[iCpu];
        if(len == 0) continue;
        // post send
        MPI_Isend(sndData+sndPtr[iCpu], len, MPI_INT, cpus[iCpu],  tag,
                  comm, sndId+sendReq);
        sendReq+=1;
    }
    MPI_Status *status = (MPI_Status *) dbg_alloca(sizeof(MPI_Status)*numNeighb);
    MPI_Waitall(rcvReq, rcvId, status);
    MPI_Waitall(sendReq, sndId, status);
#endif
}

void
FSCommunicator::exchange(int tag, int numNeighb, gsl::span<const int> cpus,
                         gsl::span<const int> sndLen, int **sndData, gsl::span<const int> rcvLen, int **rcvData) {
#ifdef USE_MPI
    int iCpu;
    MPI_Request *sndId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    MPI_Request *rcvId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    int rcvReq = 0;
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = rcvLen[iCpu];
        if(len == 0) continue;
        MPI_Irecv(rcvData[iCpu], len, MPI_INT, cpus[iCpu],
                  tag, comm,  rcvId+rcvReq);
        rcvReq += 1;
    }
    // post send
    int sendReq = 0;
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = sndLen[iCpu];
        if(len == 0) continue;
        // post send
        MPI_Isend(sndData[iCpu], len, MPI_INT, cpus[iCpu],  tag,
                  comm, sndId+sendReq);
        sendReq+=1;
    }
    MPI_Status *status = (MPI_Status *) dbg_alloca(sizeof(MPI_Status)*numNeighb);
    MPI_Waitall(rcvReq, rcvId, status);
    MPI_Waitall(sendReq, sndId, status);
#endif
}


void
FSCommunicator::exchange(int tag, int numNeighb, gsl::span<const int> cpus,
                         gsl::span<int> sndLen, double **sndData, gsl::span<int> rcvLen, double **rcvData) {
#ifdef USE_MPI
    int iCpu;
    MPI_Request *sndId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    MPI_Request *rcvId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    int rcvReq = 0;
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = rcvLen[iCpu];
        if(len == 0) continue;
        MPI_Irecv(rcvData[iCpu], len, MPI_DOUBLE, cpus[iCpu],
                  tag, comm,  rcvId+rcvReq);
        rcvReq += 1;
    }
    // post send
    int sendReq = 0;
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = sndLen[iCpu];
        if(len == 0) continue;
        // post send
        MPI_Isend(sndData[iCpu], len, MPI_DOUBLE, cpus[iCpu],  tag,
                  comm, sndId+sendReq);
        sendReq+=1;
    }
    MPI_Status *status = (MPI_Status *) dbg_alloca(sizeof(MPI_Status)*numNeighb);
    MPI_Waitall(rcvReq, rcvId, status);
    MPI_Waitall(sendReq, sndId, status);
#endif
}

void
FSCommunicator::exchange(int tag, int numNeighb, gsl::span<const int> cpus,
                         gsl::span<const int> sndLen, complex<double> **sndData, gsl::span<const int> rcvLen,
                         complex<double> **rcvData) {
#ifdef USE_MPI
    int iCpu;
    MPI_Request *sndId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    MPI_Request *rcvId = (MPI_Request *) dbg_alloca(sizeof(MPI_Request)*numNeighb);
    int rcvReq = 0;
#ifdef NO_COMPLEX_MPI
    double **tmp_rcvData = new double * [numNeighb];
  double **tmp_sndData = new double * [numNeighb];
  for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
    tmp_rcvData[iCpu] = new double[2*rcvLen[iCpu]];
    tmp_sndData[iCpu] = new double[2*sndLen[iCpu]];
    for(int j=0; j<sndLen[iCpu]; ++j) { tmp_sndData[iCpu][2*j] = sndData[iCpu][j].real(); tmp_sndData[iCpu][2*j+1] = sndData[iCpu][j].imag(); }
  }
#endif

    // post receive
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = rcvLen[iCpu];
        if(len == 0) continue;
#ifdef NO_COMPLEX_MPI
        MPI_Irecv(tmp_rcvData[iCpu], 2*len, MPI_DOUBLE, cpus[iCpu], tag, comm, rcvId+rcvReq);
#else
        MPI_Irecv(rcvData[iCpu], len, MPI_DOUBLE_COMPLEX, cpus[iCpu], tag, comm, rcvId+rcvReq);
#endif
        rcvReq += 1;
    }
    // post send
    int sendReq = 0;
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
        int len = sndLen[iCpu];
        if(len == 0) continue;
#ifdef NO_COMPLEX_MPI
        MPI_Isend(tmp_sndData[iCpu], 2*len, MPI_DOUBLE, cpus[iCpu], tag, comm, sndId+sendReq);
#else
        MPI_Isend(sndData[iCpu], len, MPI_DOUBLE_COMPLEX, cpus[iCpu], tag, comm, sndId+sendReq);
#endif
        sendReq+=1;
    }

    MPI_Status *status = (MPI_Status *) dbg_alloca(sizeof(MPI_Status)*numNeighb);
    MPI_Waitall(rcvReq, rcvId, status);
    MPI_Waitall(sendReq, sndId, status);

#ifdef NO_COMPLEX_MPI
    for(iCpu = 0; iCpu < numNeighb; ++iCpu) {
    for(int j=0; j<rcvLen[iCpu]; ++j) rcvData[iCpu][j] = complex<double>(tmp_rcvData[iCpu][2*j], tmp_rcvData[iCpu][2*j+1]);
    delete [] tmp_rcvData[iCpu];
    delete [] tmp_sndData[iCpu];
  }
  delete [] tmp_rcvData;
  delete [] tmp_sndData;
#endif

#endif
}

void
FSCommunicator::abort(int errorCode)
{
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, errorCode);
#else
    exit(errorCode);
#endif
}

#ifdef NO_COMPLEX_MPI
template <>
complex<double> FSCommunicator::globalSum(complex<double> data)
{
  double tmp_data[2] = { data.real(), data.imag() };
  globalSum(2, tmp_data);
  return complex<double>(tmp_data[0],tmp_data[1]);
}

template <>
complex<double> FSCommunicator::globalMax(complex<double> data)
{
  std::cerr << "ERROR: FSCommunicator::globalMax called with complex data\n";
  return complex<double>(0.0,0.0);
}

template <>
complex<double> FSCommunicator::globalMin(complex<double> data)
{
  std::cerr << "ERROR: FSCommunicator::globalMin called with complex data\n";
  return complex<double>(0.0,0.0);
}

template <>
void
FSCommunicator::globalSum(int num, complex<double> *data)
{
  double *tmp_data = new double[2*num];
  for(int i=0; i<num; ++i) { tmp_data[2*i] = data[i].real(); tmp_data[2*i+1] = data[i].imag(); }
  globalSum(2*num,tmp_data);
  for(int i=0; i<num; ++i) data[i] = complex<double>(tmp_data[2*i],tmp_data[2*i+1]);
  delete [] tmp_data;
}

template <>
void
FSCommunicator::globalMax(int num, complex<double> *data)
{
  std::cerr << "ERROR: FSCommunicator::globalMax called with complex data\n";
}

template <>
void
FSCommunicator::globalMin(int num, complex<double> *data)
{
  std::cerr << "ERROR: FSCommunicator::globalMin called with complex data\n";
}

template <>
void FSCommunicator::sendTo(int cpu, int tag, complex<double> *buffer, int len)
{
  std::cerr << "ERROR: FSCommunicator::sendTo called with complex data\n";
}

template <>
FSRecInfo FSCommunicator::recFrom(int tag, complex<double> *buffer, int len)
{
  std::cerr << "ERROR: FSCommunicator::recFrom called with complex data\n";
  return FSRecInfo();
}

template <>
void FSCommunicator::allGather(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_count)
{
  std::cerr << "ERROR: FSCommunicator::allGather called with complex data\n";
}

template <>
void FSCommunicator::allGatherv(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_counts[], int displacements[]) 
{
  std::cerr << "ERROR: FSCommunicator::allGatherv called with complex data\n";
}

template <>
void
FSCommunicator::reduce(int num, complex<double> *data, int root)
{ 
  std::cerr << "ERROR: FSCommunicator::reduce called with complex data\n";
}

template <>
void
FSCommunicator::broadcast(int num, complex<double> *data, int root)
{
  std::cerr << "ERROR: FSCommunicator::broadcast called with complex data\n";
}
#endif
