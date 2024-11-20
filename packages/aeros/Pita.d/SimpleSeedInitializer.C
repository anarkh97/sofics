#include "SimpleSeedInitializer.h"

#include <algorithm>

namespace Pita {

SimpleSeedInitializer::SimpleSeedInitializer(const DynamState & firstInitialSeed) :
  SeedInitializer(firstInitialSeed.vectorSize()),
  firstInitialSeed_(firstInitialSeed)
{}

DynamState
SimpleSeedInitializer::initialSeed(SliceRank rank) const {
  const_cast<SimpleSeedInitializer *>(this)->setSlices(std::max(rank, lastSlice()));
  return (rank == SliceRank(0)) ? firstInitialSeed_ : DynamState(vectorSize(), 0.0);
}

} /* end namespace Pita */
