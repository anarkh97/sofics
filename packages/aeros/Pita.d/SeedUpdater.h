#ifndef PITA_SEEDUPDATER_H
#define PITA_SEEDUPDATER_H

#include "Seed.h"
#include "NamedTask.h"

namespace Pita {

class SeedUpdater : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(SeedUpdater);
  typedef Fwk::GenManagerInterface<SeedUpdater *, String> Manager;
  class ManagerImpl;

  virtual void iterationIs(IterationRank ir);

  virtual size_t vectorSize() const;

  /* Sources */
  const Seed * propagatedSeed() const { return propagatedSeed_.ptr(); }
  Seed * correction() const { return correction_.ptr(); } // Not guaranteed to be immutable
  
  virtual void propagatedSeedIs(const Seed * ps) { setPropagatedSeed(ps); }
  virtual void correctionIs(Seed * c) { setCorrection(c); } 

  /* Targets */
  Seed * updatedSeed() const { return updatedSeed_.ptr(); }
 
  virtual void updatedSeedIs(Seed * us) { setUpdatedSeed(us); }

protected:
  explicit SeedUpdater(const String & name) :
    NamedTask(name)
  {}

  void updateSeed();

  void setPropagatedSeed(const Seed * ps) { propagatedSeed_ = ps; }
  void setCorrection(Seed * c) { correction_ = c; }
  void setUpdatedSeed(Seed * us) { updatedSeed_ = us; }

  friend class ManagerImpl;

private:
  Seed::PtrConst propagatedSeed_;
  Seed::Ptr correction_;
  Seed::Ptr updatedSeed_;

  DISALLOW_COPY_AND_ASSIGN(SeedUpdater);
};

class SeedUpdater::ManagerImpl : public SeedUpdater::Manager, private Fwk::GenManagerImpl<SeedUpdater, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(ManagerImpl);

  virtual SeedUpdater * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }

  virtual SeedUpdater * instanceNew(const String & key) { return Impl::instanceNew(key); }
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }

  static Ptr New() { return new ManagerImpl(); }

protected:
  ManagerImpl() {}
  
  virtual SeedUpdater * createNewInstance(const String & key);

private:
  typedef Fwk::GenManagerImpl<SeedUpdater, String> Impl;
};

} /* end namespace Pita */

#endif /* PITA_SEEDUPDATER_H */
