#ifndef ROM_BLOCKCYCLICMAP_H
#define ROM_BLOCKCYCLICMAP_H

#ifdef USE_SCALAPACK
extern "C" {
  int _FORTRAN(numroc)(const int *n, const int *nb, const int *iproc, const int *isrcproc, const int *nprocs);
  int _FORTRAN(indxl2g)(const int *indxloc, const int *nb, const int *iproc, const int *isrcproc, const int *nprocs);
}
#endif

namespace Rom {

class BlockCyclicMap {
public:
  BlockCyclicMap()
   : globalSize_(0), blockSize_(0), numCPU_(0), localNumSub_(0) {}

  BlockCyclicMap(int globalSize, int blockSize, int numCPU, int localNumSub)
   : globalSize_(globalSize), blockSize_(blockSize), numCPU_(numCPU), localNumSub_(localNumSub) {}  

  int subLen(int rank, int localSubId) const
  {
#ifdef USE_SCALAPACK
    int localCpu = rank;
    int zero = 0;
    int localSize = _FORTRAN(numroc)(&globalSize_, &blockSize_, &localCpu, &zero, &numCPU_);
    if(localNumSub_ > 1) {
      return localSize/localNumSub_ + (localSubId < localSize%localNumSub_ ? 1 : 0);
    }
    else {
      return localSize;
    }
#else
    return 0;
#endif
  }

  int localToGlobal(int rank, int localSubId, int i) const
  { 
#ifdef USE_SCALAPACK
    int localCpu = rank;
    int zero = 0;
    int localIdx;
    if(localNumSub_ > 1) {
      int localSize = _FORTRAN(numroc)(&globalSize_, &blockSize_, &localCpu, &zero, &numCPU_);
      localIdx = (localSize/localNumSub_)*localSubId + std::min(localSubId,localSize%localNumSub_) + i + 1;
    }
    else {
      localIdx = i + 1;
    }
    return _FORTRAN(indxl2g)(&localIdx, &blockSize_, &localCpu, &zero, &numCPU_) - 1; 
#else
    return 0;
#endif
  }

private:
  int globalSize_, blockSize_, numCPU_, localNumSub_; 
};

} /* end namespace Rom */

#endif /* ROM_BLOCKCYCLICALMAP_H */
