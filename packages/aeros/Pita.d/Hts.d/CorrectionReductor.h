#ifndef PITA_HTS_CORRECTIONREDUCTOR_H
#define PITA_HTS_CORRECTIONREDUCTOR_H

#include "Fwk.h"
#include "Types.h"

#include "../CorrectionPropagator.h"
#include "../Seed.h"

#include "../DynamStateReductor.h"

namespace Pita { namespace Hts {

class CorrectionReductor : public CorrectionPropagator<DynamState, Vector> {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionReductor);
  class Manager;

  virtual void iterationIs(IterationRank ir); // overriden

  const DynamStateReductor * reductor() const { return reductor_.ptr(); }

protected:
  CorrectionReductor(const String & name, const DynamStateReductor * reductor);

  typedef CorrectionPropagator<DynamState, Vector> BaseClass;

private:
  DynamStateReductor::PtrConst reductor_;
};

class CorrectionReductor::Manager : public Fwk::PtrInterface<CorrectionReductor::Manager> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);
  typedef Fwk::GenManagerInterface<DynamStateReductor *, HalfSliceRank> OperatorManager;

  CorrectionReductor * instance(HalfSliceRank key) const { return impl_.instance(key); }
  size_t instanceCount() const { return impl_.instanceCount(); }

  CorrectionReductor * instanceNew(HalfSliceRank key) { return impl_.instanceNew(key); }
  void instanceDel(HalfSliceRank key) { impl_.instanceDel(key); }

  static Ptr New(OperatorManager * reductorManager) {
    return new Manager(reductorManager);
  }

protected:
  explicit Manager(OperatorManager * reductorManager);

  class Factory : public Fwk::InstanceFactory<CorrectionReductor, HalfSliceRank> {
  public:
    CorrectionReductor * operator()(HalfSliceRank key);

    Factory(OperatorManager * manager) : manager_(manager) {}

  private:
    OperatorManager::Ptr manager_;
  };

private:
  typedef Fwk::FactoryManagerImpl<Factory> Impl;
  Impl impl_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTSCORRECTIONREDUCTOR_H */
