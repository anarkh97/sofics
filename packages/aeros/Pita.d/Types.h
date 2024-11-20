#ifndef PITA_TYPES_H
#define PITA_TYPES_H

#include "Fwk.d/Fwk.h"

namespace Pita {

class Time;
typedef Fwk::Numeric<Time, double> Seconds;

class TimeStep;
typedef Fwk::Numeric<TimeStep, int> TimeStepCount;

class TimeSlice;
typedef Fwk::Numeric<TimeSlice, int> SliceCount;

class SliceRank : public Fwk::Interval<TimeSlice, int, SliceCount> {
public:
  explicit SliceRank(int v = 0) : Fwk::Interval<TimeSlice, int, SliceCount>(v) {}

  SliceRank next() const { return SliceRank(value() + 1); }
  SliceRank previous() const { return SliceRank(value() - 1); }
};

inline
SliceRank
operator+(const SliceRank & a, const SliceCount & d) {
  return SliceRank(a.value() + d.value());
}

inline 
SliceRank
operator-(const SliceRank & a, const SliceCount & d) {
  return SliceRank(a.value() - d.value());
}

class Cpu;
typedef Fwk::Numeric<Cpu, int> CpuCount;
class CpuTempFix;
typedef Fwk::Numeric<CpuTempFix, int> CpuRank;

inline
CpuRank
operator+(const CpuRank & a, const CpuCount & d) {
  return CpuRank(a.value() + d.value());
}

inline 
CpuRank
operator-(const CpuRank & a, const CpuCount & d) {
  return CpuRank(a.value() - d.value());
}

class IterationRank : public Fwk::Ordinal<IterationRank, int> {
public:
  explicit IterationRank(int rank = 0) : Fwk::Ordinal<IterationRank, int>(rank) {}
  IterationRank next() const { return IterationRank(value() + 1); }
};

} /* end namespace Pita */

#endif /* PITA_TYPES_H */
