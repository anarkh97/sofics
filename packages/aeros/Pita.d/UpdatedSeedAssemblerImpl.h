#ifndef PITA_UPDATEDSEEDASSEMBLERIMPL_H
#define PITA_UPDATEDSEEDASSEMBLERIMPL_H

#include "UpdatedSeedAssembler.h"

#include "DynamStateBasis.h"

namespace Pita {

class UpdatedSeedAssemblerImpl : public UpdatedSeedAssembler {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdatedSeedAssemblerImpl);
  class Manager;

  virtual size_t reducedBasisSize() const { return correctionBasis()->stateCount(); }
  virtual size_t vectorSize() const { return correctionBasis()->vectorSize(); }

  using UpdatedSeedAssembler::updatedSeed;
  Seed * updatedSeed() { return const_cast<Seed *>(const_cast<const UpdatedSeedAssemblerImpl *>(this)->updatedSeed()); }

  /* Added members */
  const Manager * manager() const { return manager_; }
  const DynamStateBasis * correctionBasis() const;

  virtual void iterationIs(IterationRank i);

protected:
  friend class Manager;

  explicit UpdatedSeedAssemblerImpl(const String & name, const Manager * manager);

private:
  const Manager * manager_; // Back pointer
};


class UpdatedSeedAssemblerImpl::Manager : public UpdatedSeedAssembler::Manager,
                                          private Fwk::GenManagerImpl<UpdatedSeedAssemblerImpl, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  /* Overiden members */
  virtual UpdatedSeedAssemblerImpl * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }
  virtual UpdatedSeedAssemblerImpl * instanceNew(const String & key) { return Impl::instanceNew(key); } 
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }

  /* Added members */
  const DynamStateBasis * defaultCorrectionBasis() const { return defaultCorrectionBasis_.ptr(); }
  void defaultCorrectionBasisIs(const DynamStateBasis * dcb) { defaultCorrectionBasis_ = dcb; }

  static Ptr New(const DynamStateBasis * defaultCorrectionBasis) {
    return new Manager(defaultCorrectionBasis);
  }

protected:
  explicit Manager(const DynamStateBasis * defaultCorrectionBasis);

private:
  typedef Fwk::GenManagerImpl<UpdatedSeedAssemblerImpl, String> Impl;
  
  // Overriden
  virtual UpdatedSeedAssemblerImpl * createNewInstance(const String & key);

  DynamStateBasis::PtrConst defaultCorrectionBasis_;
};

inline
const DynamStateBasis *
UpdatedSeedAssemblerImpl::correctionBasis() const {
  return manager_->defaultCorrectionBasis();
}

} /* end namespace Pita */

#endif /* PITA_UPDATEDSEEDASSEMBLERIMPL_H */
