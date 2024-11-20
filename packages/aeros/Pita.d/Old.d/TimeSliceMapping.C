#include "TimeSliceMapping.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>

namespace Pita { namespace Old {

TimeSliceMapping::TimeSliceMapping(int nTS, int nCPU, int nMaxActive)
: numSlices_(nTS), numCPUs_(nCPU), numMaxActiveSlices_(nMaxActive), numCvgSlices_(0), firstInactiveSlice_(0), TStoCPU(NULL), CPUtoTS(NULL)
{
  numCvgSlices(0);
  int* arrayTS = new int[numSlices_ + 1];
  for (int i = 0; i <= numSlices_; ++i)
    arrayTS[i] = i;

  int* targetCPU = new int[numSlices_];
  
  int ratio, remain;
  if (numSlices_ < numCPUs_ * numMaxActiveSlices_)
  {
    ratio = numSlices_ / numCPUs_;
    remain = numSlices_ - ratio * numCPUs_;
  }
  else
  {
    ratio = numMaxActiveSlices_;
    remain = 0;
  }

  int currentTS = 0;
  int currentCPU = 0;
  int remainingSpotsOnCurrentCPU = ratio + (remain > 0 ? 1 : 0);
  while(currentTS < numSlices_)
  {
    if (--remainingSpotsOnCurrentCPU < 0)
    {
      currentCPU = (currentCPU + 1) % numCPUs_;
      remainingSpotsOnCurrentCPU = ratio + (remain > currentCPU ? 0 : -1);
    }
    targetCPU[currentTS++] = currentCPU; 
  }

  TStoCPU = new Connectivity(numSlices_, arrayTS, targetCPU, 1);
  CPUtoTS = TStoCPU->reverse();
  CPUtoTS->sortTargets(); 
}

TimeSliceMapping::~TimeSliceMapping()
{
  delete CPUtoTS;
  delete TStoCPU;
}

void TimeSliceMapping::numCvgSlices(int numCvgTS)
{
  numCvgSlices_ = numCvgTS;
  firstInactiveSlice_ = std::min(numSlices_, numMaxActiveSlices_ * numCPUs_ + numCvgSlices_);
}

int TimeSliceMapping::numActiveSlicesOnCPU(int CPUid) const
{
  int activeSlicesCount = 0;
  int totalSlicesOnCPU = CPUtoTS->num(CPUid);
  int *slicesOnCPU = CPUtoTS->operator[](CPUid);
  for (int i = 0; i < totalSlicesOnCPU; ++i)
  {
    if (slicesOnCPU[i] >= firstInactiveSlice_)
      break;
    if (slicesOnCPU[i] >= numCvgSlices_)
      ++activeSlicesCount;
  }
  return activeSlicesCount;
}

bool TimeSliceMapping::hasFirstInactiveSlice(int CPUid) const
{
  if (firstInactiveSlice_ >= numSlices_)
    return false;
  return (*(TStoCPU->operator[](firstInactiveSlice_)) == CPUid);
}

void TimeSliceMapping::print(std::ostream & out) const
{
  for (int i = 0; i < CPUtoTS->csize(); ++i)
  {
    out << "CPU # " << i << " ->";
    for (int j = 0; j < CPUtoTS->num(i); ++j)
    {
      out << ' ' << (*CPUtoTS)[i][j];
    }
    out << '\n';
  }
}

} /* end namespace Old */ } /* end namespace Pita */
