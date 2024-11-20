#ifndef PITA_OLD_TIMESLICEMAPPING_H
#define PITA_OLD_TIMESLICEMAPPING_H

#include <Utils.d/Connectivity.h>
#include <algorithm>
#include <iostream>

namespace Pita { namespace Old {

// Represents the mapping of time-slices on a set of available processors.
class TimeSliceMapping
{
public:
  // Constructors & destructor
  TimeSliceMapping(int nTS, int nCPU, int nMaxActive);
  ~TimeSliceMapping();

  // Mutators
  void numCvgSlices(int numCvgTS);
  int incrementNumCvgSlices() { numCvgSlices(numCvgSlices_ + 1); return numCvgSlices_; }

  // Accessors
  int numCvgSlices() const { return numCvgSlices_; } 
  int firstInactiveSlice() const { return firstInactiveSlice_; }
  int numCPUs() const { return numCPUs_; }
  int numSlices() const { return numSlices_; }
  int numMaxActiveSlices() const { return numMaxActiveSlices_; }
  int numSlicesOnCPU(int CPUid) const { return CPUtoTS->num(CPUid); } 
  int numActiveSlices() const { return firstInactiveSlice_ - numCvgSlices_; }
  int numCPU(int numTS) const { return TStoCPU->getTargetValue(numTS); } 
  int* slicesOnCPU(int CPUid) const { return CPUtoTS->operator[](CPUid); }
  bool hasFirstInactiveSlice(int CPUid) const;
  int numActiveSlicesOnCPU(int CPUid) const;

  // Display
  void print(std::ostream &) const;

private:
  // Structural variables
  int numSlices_;
  int numCPUs_;
  int numMaxActiveSlices_;
 
  // Computational status
  int numCvgSlices_;  // Also corresponds to first active slice number
  int firstInactiveSlice_;

  // Mapping slices <-> CPUs
  Connectivity *TStoCPU;
  Connectivity *CPUtoTS;
};

inline std::ostream & operator<<(std::ostream & out, const TimeSliceMapping & tsm)
{
  tsm.print(out);
  return out;
}

} /* end namespace Old */ } /* end namespace Pita */

#endif /* PITA_OLD_TIMESLICEMAPPING_H */
