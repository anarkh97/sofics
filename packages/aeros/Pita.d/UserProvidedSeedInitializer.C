#include "UserProvidedSeedInitializer.h"

#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>

namespace Pita {

inline
SliceRank
lastUserProvidedSlice(const GeoSource * gs) {
  int result = gs->getUserProvidedSeedCount() - 1;
  return SliceRank(result);
}

UserProvidedSeedInitializer::UserProvidedSeedInitializer(size_t vectorSize,
                                                         const GeoSource * geoSource,
                                                         const Domain * domain) :
  SeedInitializer(vectorSize, lastUserProvidedSlice(geoSource)),
  domain_(domain)
{}

DynamState
UserProvidedSeedInitializer::initialSeed(SliceRank seedId) const {
  if (seedId > lastSlice()) {
    return DynamState(vectorSize(), 0.0);
  }

  DynamState result(vectorSize());
  const_cast<Domain *>(domain_)->initDispVelocOnTimeSlice(result.displacement(), result.velocity(), seedId.value());
  return result;
}

} /* end namespace Pita */
