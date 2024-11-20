#ifndef PITA_SEEDINITIALIZER_H
#define PITA_SEEDINITIALIZER_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"

namespace Pita {

class SeedInitializer : public Fwk::PtrInterface<SeedInitializer> {
public:
  EXPORT_PTRINTERFACE_TYPES(SeedInitializer);

  size_t vectorSize() const { return vectorSize_; }
  
  SliceRank lastSlice() const { return lastSlice_; }
  virtual DynamState initialSeed(SliceRank rank) const = 0;

protected:
  explicit SeedInitializer(size_t vectorSize, SliceRank lastSlice = SliceRank()) :
    vectorSize_(vectorSize),
    lastSlice_(lastSlice)
  {}
  
  void setSlices(SliceRank last) { lastSlice_ = last; }

private:
  size_t vectorSize_;
  SliceRank lastSlice_;

  DISALLOW_COPY_AND_ASSIGN(SeedInitializer);
};
  
} /* end namespace Pita */

#endif /* PITA_SEEDINITIALIZER_H */
