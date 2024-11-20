#ifndef PITA_UPDATEDSEEDASSEMBLER_H
#define PITA_UPDATEDSEEDASSEMBLER_H

#include "Fwk.h"
#include "Types.h"
#include "Seed.h"
#include "SeedUpdater.h"

namespace Pita {

class UpdatedSeedAssembler : public SeedUpdater {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdatedSeedAssembler);
  typedef Fwk::GenManagerSubInterface<UpdatedSeedAssembler *, String, Fwk::GenManagerInterface<SeedUpdater *, String> > Manager;

  virtual size_t reducedBasisSize() const = 0;

  /* Sources */
  const ReducedSeed * correctionComponents() const { return correctionComponents_.ptr(); }
  
  virtual void correctionComponentsIs(const ReducedSeed * cc) { setCorrectionComponents(cc); }

protected:
  explicit UpdatedSeedAssembler(const String & name) :
    SeedUpdater(name)
  {}

  void setCorrectionComponents(const ReducedSeed * cc) { correctionComponents_ = cc; }

private:
  ReducedSeed::PtrConst correctionComponents_; 

  DISALLOW_COPY_AND_ASSIGN(UpdatedSeedAssembler);
};

} /* end namespace Pita */

#endif /* PITA_UPDATEDSEEDASSEMBLER_H */
