#ifndef PITA_JUMPBUILDER_H
#define PITA_JUMPBUILDER_H

#include "Fwk.h"
#include "Seed.h"
#include "NamedTask.h"

namespace Pita {

class JumpBuilder : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpBuilder);
  typedef Fwk::GenManagerInterface<JumpBuilder *, String> Manager;
  class ManagerImpl;

  virtual void iterationIs(IterationRank ir); // overriden

  // Input
  const Seed * predictedSeed() const { return predictedSeed_.ptr(); }
  const Seed * actualSeed() const { return actualSeed_.ptr(); }

  void predictedSeedIs(const Seed * ps) { setPredictedSeed(ps); }
  void actualSeedIs(const Seed * as) { setActualSeed(as); }

  // Output
  const Seed * seedJump() const { return seedJump_.ptr(); }
  Seed * seedJump() { return seedJump_.ptr(); }

  void seedJumpIs(Seed * sj) { setSeedJump(sj); }

protected:
  explicit JumpBuilder(const String & name) :
    NamedTask(name)
  {}
 
  void setPredictedSeed(const Seed * ps) { predictedSeed_ = ps; }
  void setActualSeed(const Seed * as) { actualSeed_ = as; }
  void setSeedJump(Seed * sj) { seedJump_ = sj; }

private:
  Seed::PtrConst predictedSeed_;
  Seed::PtrConst actualSeed_;
  Seed::Ptr seedJump_;
  
  DISALLOW_COPY_AND_ASSIGN(JumpBuilder);
};

class JumpBuilder::ManagerImpl : public Fwk::GenManagerInterface<JumpBuilder *, String>, private Fwk::GenManagerImpl<JumpBuilder, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(ManagerImpl);

  // Overriden members
  virtual JumpBuilder * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }
  virtual JumpBuilder * instanceNew(const String & key) { return Impl::instanceNew(key); }
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }
  
  static Ptr New() {
    return new ManagerImpl();
  }

protected:
  ManagerImpl() {}

  virtual JumpBuilder * createNewInstance(const String & key);

private:
  typedef Fwk::GenManagerImpl<JumpBuilder, String> Impl;
}; 

} /* end namespace Pita */

#endif /* PITA_JUMPBUILDER_H */
