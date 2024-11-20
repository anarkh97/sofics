#ifndef PITA_JUMPPROJECTION_H
#define PITA_JUMPPROJECTION_H

#include "NamedTask.h"

#include "Seed.h"
#include "DynamStateBasis.h"

namespace Pita {

class JumpProjection : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpProjection);
  class Manager;

  virtual void iterationIs(IterationRank ir); // overriden
  
  // Input
  const Seed * seedJump() const { return seedJump_.ptr(); }
 
  void seedJumpIs(const Seed * sj) { setSeedJump(sj); }

  // Output
  const ReducedSeed * reducedSeedJump() const { return reducedSeedJump_.ptr(); }
  ReducedSeed * reducedSeedJump() { return reducedSeedJump_.ptr(); }

  void reducedSeedJumpIs(ReducedSeed * rsj) { setReducedSeedJump(rsj); }

  // Operator
  size_t reducedBasisSize() const { return reducedBasis()->stateCount(); }
  const DynamStateBasis * reducedBasis() const;

protected:
  JumpProjection(const String & name, const Manager * manager);

  void setSeedJump(const Seed * sj) { seedJump_ = sj; }
  void setReducedSeedJump(ReducedSeed * rsj) { reducedSeedJump_ = rsj; }
  
  friend class Manager;

private:
  Seed::PtrConst seedJump_;
  ReducedSeed::Ptr reducedSeedJump_;
 
  // Back pointer 
  const Manager * manager_;
};


class JumpProjection::Manager : public Fwk::GenManagerInterface<JumpProjection *, String>, private Fwk::GenManagerImpl<JumpProjection, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  // Overriden members
  virtual JumpProjection * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }
  virtual JumpProjection * instanceNew(const String & key) { return Impl::instanceNew(key); }
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }

  // Added members
  const DynamStateBasis * defaultReducedBasis() const { return defaultReducedBasis_.ptr(); }

  void defaultReducedBasisIs(const DynamStateBasis * drb) { defaultReducedBasis_ = drb; }

  static Ptr New(const DynamStateBasis * defaultReducedBasis) {
    return new Manager(defaultReducedBasis);
  }

protected:
  explicit Manager(const DynamStateBasis * drb);

  virtual JumpProjection * createNewInstance(const String & key);

private:
  typedef Fwk::GenManagerImpl<JumpProjection, String> Impl;
  
  DynamStateBasis::PtrConst defaultReducedBasis_;
};

inline
const DynamStateBasis *
JumpProjection::reducedBasis() const {
  return manager_->defaultReducedBasis();
}

} /* end namespace Pita */

#endif /* PITA_JUMPPROJECTION_H */
