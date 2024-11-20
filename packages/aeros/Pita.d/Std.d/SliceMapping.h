#ifndef PITA_STD_SLICEMAPPING_H
#define PITA_STD_SLICEMAPPING_H

#include "Fwk.h"
#include "Types.h"

#include "../LoadBalancer.h"

namespace Pita { namespace Std {

class SliceMapping : public Fwk::PtrInterface<SliceMapping> {
public:
  EXPORT_PTRINTERFACE_TYPES(SliceMapping);

  class SliceIterator;

  SliceCount totalSlices() const;
  CpuCount availableCpus() const;
  SliceCount maxWorkload() const;

  SliceCount activeSlices() const;
  SliceCount activeSlices(CpuRank cpu) const;

  SliceRank firstActiveSlice() const;
  SliceRank firstInactiveSlice() const;

  SliceCount convergedSlices() const;
  void convergedSlicesInc(SliceCount increment = SliceCount(1));

  CpuRank hostCpu(SliceRank slice) const;
 
  SliceIterator hostedSlice(CpuRank cpu) const;

  static Ptr New(SliceCount totalSlices, CpuCount availableCpus, SliceCount maxWorkload) {
    return new SliceMapping(totalSlices, availableCpus, maxWorkload); 
  }

protected:
  SliceMapping(SliceCount totalSlices, CpuCount availableCpus, SliceCount maxWorkload);

  friend class SliceIterator;

private:
  friend OStream & operator<<(OStream &, const SliceMapping &);
  LoadBalancer::Ptr loadBalancer_;

  DISALLOW_COPY_AND_ASSIGN(SliceMapping);
};


class SliceMapping::SliceIterator {
public:
  SliceRank operator*() const { return SliceRank(*impl_); }
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

inline
OStream & operator<<(OStream & out, const SliceMapping & sm) {
  return out << *(sm.loadBalancer_);
}

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_SLICEMAPPING_H */
