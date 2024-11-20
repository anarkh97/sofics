#ifndef PITA_STD_LINEARPROPAGATORMANAGER_H
#define PITA_STD_LINEARPROPAGATORMANAGER_H

#include "Fwk.h"
#include "Types.h"

#include "../LinearDynamOps.h"
#include "../AffineIntegratorPropagator.h"
#include "../AffineDynamTimeIntegrator.h"
#include "../PostProcessingManager.h"
#include "AffineBasisCollector.h"

namespace Pita { namespace Std {

class LinearPropagatorManager : public Fwk::GenManagerInterface<AffineDynamPropagator *, SliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearPropagatorManager);

  virtual AffineIntegratorPropagator * instance(const SliceRank & id) const;
  virtual size_t instanceCount() const;

  virtual AffineIntegratorPropagator * instanceNew(const SliceRank & id);
  virtual void instanceDel(const SliceRank & id);

  static Ptr New(AffineDynamTimeIntegrator * sharedIntegrator,
                 PostProcessing::Manager * postProcessingMgr,
                 AffineBasisCollector * collector,
                 TimeStepCount timeSliceRatio,
                 Seconds initialTime,
                 AffineDynamPropagator::ConstantTerm defaultMode) {
    return new LinearPropagatorManager(sharedIntegrator, postProcessingMgr, collector, timeSliceRatio, initialTime, defaultMode);
  }

protected:
  LinearPropagatorManager(AffineDynamTimeIntegrator * sharedIntegrator,
                          PostProcessing::Manager * postProcessingMgr,
                          AffineBasisCollector * collector,
                          TimeStepCount timeSliceRatio,
                          Seconds initialTime,
                          AffineDynamPropagator::ConstantTerm defaultMode);

  AffineIntegratorPropagator::Ptr createNewInstance(SliceRank rank);

private:
  AffineDynamTimeIntegrator::Ptr sharedIntegrator_;
  Seconds fineTimeStep_;
  TimeStepCount timeSliceRatio_;
  Seconds initialTime_;
  AffineDynamPropagator::ConstantTerm defaultMode_;

  PostProcessing::Manager::Ptr postProcessingMgr_;
  AffineBasisCollector * collector_;
  
  typedef std::map<SliceRank, AffineIntegratorPropagator::Ptr> InstanceMap;
  InstanceMap instance_;
};

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_LINEARPROPAGATORMANAGER_H */
