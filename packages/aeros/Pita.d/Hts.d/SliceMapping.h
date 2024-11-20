#ifndef PITA_HTS_SLICEMAPPING_H
#define PITA_HTS_SLICEMAPPING_H

#include "Fwk.h"
#include "Types.h"

#include "../LoadBalancer.h"

#include <vector>

namespace Pita { namespace Hts {

class SliceMapping : public Fwk::PtrInterface<SliceMapping> {
public:
  EXPORT_PTRINTERFACE_TYPES(SliceMapping);

  class SliceIterator;

  FullSliceCount totalFullSlices() const { return FullSliceCount(totalSlices().value() / 2); }
  HalfSliceCount totalSlices() const;
  CpuCount availableCpus() const;
  HalfSliceCount maxWorkload() const;

  HalfSliceCount activeSlices() const;
  HalfSliceCount activeSlices(CpuRank cpu) const;

  HalfSliceRank firstActiveSlice() const;
  HalfSliceRank firstInactiveSlice() const;

  FullSliceCount activePrimalSlices() const;
  FullSliceCount activeDualSlices() const;

  HalfSliceCount convergedSlices() const;
  void convergedSlicesInc(HalfSliceCount increment = HalfSliceCount(1));

  CpuRank hostCpu(HalfSliceRank slice) const;
 
  SliceIterator hostedSlice(CpuRank cpu) const;

  static Ptr New(FullSliceCount totalFullSlices, CpuCount availableCpus, HalfSliceCount maxWorkload) {
    return new SliceMapping(totalFullSlices, availableCpus, maxWorkload); 
  }

protected:
  SliceMapping(FullSliceCount totalFullSlices, CpuCount availableCpus, HalfSliceCount maxWorkload);

  friend class SliceIterator;

private:
  LoadBalancer::Ptr loadBalancer_;

  DISALLOW_COPY_AND_ASSIGN(SliceMapping);
};


class SliceMapping::SliceIterator {
public:
  HalfSliceRank operator*() const { return HalfSliceRank(*impl_); }
  SliceIterator & operator++() { ++impl_; return *this; }
  SliceIterator operator++(int) { SliceIterator tmp(*this); ++(*this); return tmp; }
  operator bool() const { return impl_; }

  explicit SliceIterator(LoadBalancer::TaskIterator impl) :
    impl_(impl)
  {}

  // Default copy, assignment

private:
  LoadBalancer::TaskIterator impl_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_SLICEMAPPING_H */
