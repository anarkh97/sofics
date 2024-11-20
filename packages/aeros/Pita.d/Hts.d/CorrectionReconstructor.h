#ifndef PITA_HTS_CORRECTIONRECONSTRUCTOR_H
#define PITA_HTS_CORRECTIONRECONSTRUCTOR_H

#include "Fwk.h"
#include "Types.h"

#include "../NamedTask.h"
#include "../Seed.h"

#include "../DynamStateReconstructor.h"

namespace Pita { namespace Hts {

class CorrectionReconstructor : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionReconstructor);
  class Manager;

  virtual void iterationIs(IterationRank ir); // overriden

  size_t vectorSize() const { return reconstructor()->vectorSize(); }
  size_t reducedBasisSize() const { return reconstructor()->reducedBasisSize(); }
  
  Seed * correction() const { return correction_.ptr(); }
  virtual void correctionIs(Seed * correction) { setCorrection(correction); }

  const ReducedSeed * correctionComponents() const { return correctionComponents_.ptr(); }
  virtual void correctionComponentsIs(const ReducedSeed * cc) { setCorrectionComponents(cc); }

protected:
  explicit CorrectionReconstructor(const String & name, DynamStateReconstructor * reconstructor) :
    NamedTask(name),
    reconstructor_(reconstructor)
  {}

  DynamStateReconstructor * reconstructor() const { return reconstructor_.ptr(); }

  void setCorrection(Seed * c) { correction_ = c; }
  void setCorrectionComponents(const ReducedSeed * cc) { correctionComponents_ = cc; }

private:
  Seed::Ptr correction_;
  ReducedSeed::PtrConst correctionComponents_;
  DynamStateReconstructor::Ptr reconstructor_;

  DISALLOW_COPY_AND_ASSIGN(CorrectionReconstructor);
};

class CorrectionReconstructor::Manager : public Fwk::PtrInterface<CorrectionReconstructor::Manager> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  CorrectionReconstructor * instance(HalfSliceRank key) const { return impl_.instance(key); }
  size_t instanceCount() const { return impl_.instanceCount(); }

  CorrectionReconstructor * instanceNew(HalfSliceRank key) { return impl_.instanceNew(key); }
  void instanceDel(HalfSliceRank key) { impl_.instanceDel(key); }

  typedef Fwk::GenManagerInterface<DynamStateReconstructor*, HalfSliceRank> OperatorManager;

  static Ptr New(OperatorManager * reconstructorManager) {
    return new Manager(reconstructorManager);
  }

protected:
  explicit Manager(OperatorManager * reconstructorManager);

  class Factory : public Fwk::InstanceFactory<CorrectionReconstructor, HalfSliceRank> {
  public:
    explicit Factory(OperatorManager * manager) :
      manager_(manager)
    {}

    CorrectionReconstructor * operator()(HalfSliceRank key);
  private:
    OperatorManager::Ptr manager_;
  }; 
private:
  typedef Fwk::FactoryManagerImpl<Factory> Impl;
  Impl impl_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_CORRECTIONRECONSTRUCTOR_H */
