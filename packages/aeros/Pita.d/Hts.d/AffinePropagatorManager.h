#ifndef PITA_HTS_AFFINEPROPAGATOR_MANAGER
#define PITA_HTS_AFFINEPROPAGATOR_MANAGER

#include "Fwk.h"
#include "Types.h"
#include "HalfSliceId.h"

#include "../AffineDynamPropagator.h"

#include "FineIntegratorManager.h"
#include "../PostProcessingManager.h"
#include "AffineBasisCollector.h"

#include "../LinearGenAlphaIntegrator.h"

namespace Pita { namespace Hts {

// DynamPropagator::Manager implementation for Pita HalfSlice problems
// Attach the necessary Reactors to the produced propagators
class AffinePropagatorManager : public Fwk::GenManagerInterface<AffineDynamPropagator*, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(AffinePropagatorManager);
 
  virtual AffineDynamPropagator * instance(const HalfSliceId & id) const;
  virtual size_t instanceCount() const;

  virtual AffineDynamPropagator * instanceNew(const HalfSliceId & id);
  virtual void instanceDel(const HalfSliceId & id);

  static Ptr New(AffineBasisCollector * collector,
                 GenFineIntegratorManager<AffineGenAlphaIntegrator> * integratorMgr,
                 PostProcessing::Manager * postProcessingMgr,
                 TimeStepCount halfSliceRatio,
                 Seconds initialTime,
                 AffineDynamPropagator::ConstantTerm defaultMode) {
    return new AffinePropagatorManager(collector, integratorMgr, postProcessingMgr, halfSliceRatio, initialTime, defaultMode);
  }

protected:
  AffinePropagatorManager(AffineBasisCollector * collector,
                          GenFineIntegratorManager<AffineGenAlphaIntegrator> * integratorMgr,
                          PostProcessing::Manager * postProcessingMgr,
                          TimeStepCount halfSliceRatio,
                          Seconds initialTime,
                          AffineDynamPropagator::ConstantTerm defaultMode);

private:
  AffineBasisCollector * collector_;
  GenFineIntegratorManager<AffineGenAlphaIntegrator>::Ptr integratorMgr_;
  PostProcessing::Manager::Ptr postProcessingMgr_;
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  Seconds initialTime_;
  AffineDynamPropagator::ConstantTerm defaultMode_;

  DISALLOW_COPY_AND_ASSIGN(AffinePropagatorManager);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_AFFINEPROPAGATOR_MANAGER */
