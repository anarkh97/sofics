#ifndef PITA_HTS_NLPROPAGATOR_MANAGER
#define PITA_HTS_NLPROPAGATOR_MANAGER

#include "Fwk.h"
#include "Types.h"
#include "HalfSliceId.h"

#include "../IntegratorPropagator.h"
#include "../PostProcessingManager.h"

#include "../NlDynamTimeIntegrator.h"

#include "../ConcurrentBasisManager.h"

namespace Pita { namespace Hts {

// Create and decorate the propagators
class NlPropagatorManager : public Fwk::GenManagerInterface<DynamPropagator*, HalfSliceId>,
                            private Fwk::GenManagerImpl<GenIntegratorPropagator<NlDynamTimeIntegrator>,
                                                        HalfSliceId,
                                                        Fwk::Ptr<GenIntegratorPropagator<NlDynamTimeIntegrator> > > {
public:
  EXPORT_PTRINTERFACE_TYPES(NlPropagatorManager);

  virtual GenIntegratorPropagator<NlDynamTimeIntegrator> * instance(const HalfSliceId & id) const { return Impl::instance(id); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }

  virtual GenIntegratorPropagator<NlDynamTimeIntegrator> * instanceNew(const HalfSliceId & id) { return Impl::instanceNew(id); }
  virtual void instanceDel(const HalfSliceId & id) { Impl::instanceDel(id); }

  NlDynamOps::Ptr dynamOpsNew() { return integrator_->nlDynamOpsNew(); }

  static Ptr New(NlDynamTimeIntegrator * integrator,
                 Seconds fineTimeStep,
                 TimeStepCount halfSliceRatio,
                 Seconds initialTime) {
    return new NlPropagatorManager(integrator, fineTimeStep, halfSliceRatio, initialTime);
  }
  
  typedef Fwk::GenManager<DynamStatePlainBasis, HalfSliceId> LocalBasisManager;
  
  void postProcessingManagerIs(PostProcessing::Manager * mgr) { postProcessingMgr_ = mgr; }
  void concurrentBasisManagerIs(ConcurrentBasis::Manager * mgr) { concurrentBasisMgr_ = mgr; }
  void propagatedBasisManagerIs(LocalBasisManager * mgr) { propagatedBasisMgr_ = mgr; }

protected:
  NlPropagatorManager(NlDynamTimeIntegrator * integrator,
                      Seconds fineTimeStep,
                      TimeStepCount halfSliceRatio,
                      Seconds initialTime);

protected:
  virtual Fwk::Ptr<GenIntegratorPropagator<NlDynamTimeIntegrator> > createNewInstance(const HalfSliceId & id); // overriden

private:
  typedef Fwk::GenManagerImpl<GenIntegratorPropagator<NlDynamTimeIntegrator>, HalfSliceId, Fwk::Ptr<GenIntegratorPropagator<NlDynamTimeIntegrator> > > Impl;

  // Core members
  NlDynamTimeIntegrator::Ptr integrator_;
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  Seconds initialTime_;

  // Optional members
  PostProcessing::Manager::Ptr postProcessingMgr_;
  ConcurrentBasis::Manager::Ptr concurrentBasisMgr_;
  LocalBasisManager::Ptr propagatedBasisMgr_;

  DISALLOW_COPY_AND_ASSIGN(NlPropagatorManager);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLPROPAGATOR_MANAGER */
