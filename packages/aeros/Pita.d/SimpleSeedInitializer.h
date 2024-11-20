#ifndef PITA_SIMPLESEEDINITIALIZER_H
#define PITA_SIMPLESEEDINITIALIZER_H

#include "SeedInitializer.h"

namespace Pita {
class SimpleSeedInitializer : public SeedInitializer {
public:
  EXPORT_PTRINTERFACE_TYPES(SimpleSeedInitializer);

  virtual DynamState initialSeed(SliceRank rank) const; // Overriden

  static Ptr New(const DynamState & firstInitialSeed) {
    return new SimpleSeedInitializer(firstInitialSeed);
  }

protected:
  explicit SimpleSeedInitializer(const DynamState & firstInitialSeed);

private:
  DynamState firstInitialSeed_;
};

} /* end namespace Pita */

#endif /* PITA_SIMPLESEEDINITIALIZER_H */
