#ifndef PITA_USERPROVIDEDSEEDINITIALIZER_H
#define PITA_USERPROVIDEDSEEDINITIALIZER_H

#include "SeedInitializer.h"

class GeoSource;
class Domain;

namespace Pita {

class UserProvidedSeedInitializer : public SeedInitializer {
public:
  EXPORT_PTRINTERFACE_TYPES(UserProvidedSeedInitializer);

  virtual DynamState initialSeed(SliceRank rank) const; // overriden

  static Ptr New(size_t vectorSize, const GeoSource * geoSource, const Domain * domain) {
    return new UserProvidedSeedInitializer(vectorSize, geoSource, domain);
  }

protected:
  UserProvidedSeedInitializer(size_t vectorSize,
                              const GeoSource * geoSource,
                              const Domain * domain);

private:
  const Domain * domain_;
};
  
} /* end namespace Pita */

#endif /* PITA_USERPROVIDEDSEEDINITIALIZER_H */
