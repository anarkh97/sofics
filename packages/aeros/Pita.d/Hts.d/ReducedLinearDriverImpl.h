#ifndef PITA_HTS_REDUCEDLINEARDRIVERIMPL_H
#define PITA_HTS_REDUCEDLINEARDRIVERIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "../LinearDriverImpl.h"

#include "../DynamState.h"

#include "SliceMapping.h"

#include "../LinearGenAlphaIntegrator.h"
#include "../PostProcessingManager.h"
#include "../CorrectionPropagator.h"
#include "../SeedInitializer.h"

class SingleDomainDynamic;
class GeoSource;
class Domain;
class SolverInfo;
class Communicator;
class SDDynamPostProcessor;

namespace Pita { namespace Hts {

class ReducedLinearDriverImpl : public LinearDriverImpl {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedLinearDriverImpl);

  virtual void solve();
  
  static ReducedLinearDriverImpl::Ptr New(SingleDomainDynamic * pbDesc,
                                          GeoSource * geoSource,
                                          Domain * domain,
                                          SolverInfo * solverInfo,
                                          Communicator * baseComm) {
    return new ReducedLinearDriverImpl(pbDesc, geoSource, domain, solverInfo, baseComm);
  }

protected:
  ReducedLinearDriverImpl(SingleDomainDynamic *, GeoSource *, Domain *, SolverInfo *, Communicator *);

  void preprocess();
  void solveParallel(Communicator * timeComm, Communicator * coarseComm);
  void solveCoarse(Communicator * timeComm);
 
  PostProcessing::Manager::Ptr buildPostProcessor(CpuRank localCpu) const;
  CorrectionPropagator<DynamState>::Manager::Ptr buildCoarseCorrection(Communicator * coarseComm) const;
  LinearGenAlphaIntegrator::Ptr buildCoarseIntegrator() const; 
  DynamPropagator::Ptr buildCoarsePropagator(Communicator * coarseComm = NULL) const;
  SeedInitializer::Ptr buildSeedInitializer(Communicator * timeComm = NULL) const;

private:
  /* Space-domain */
  size_t vectorSize_;
  LinearDynamOps::Manager::Ptr dynamOpsMgr_;

  /* Time-domain */ 
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  TimeStepCount sliceRatio_;
  Seconds coarseTimeStep_;
  Seconds initialTime_;
  Seconds finalTime_;

  /* Main options */
  bool noForce_;
  bool userProvidedSeeds_; 
  bool remoteCoarse_;
  
  /* Load balancing */ 
  SliceMapping::Ptr mapping_;

  /* Other parameters */
  IterationRank lastIteration_;
  double jumpCvgRatio_;
  double projectorTolerance_;
  double coarseRhoInfinity_;
  
  /* PITA-specific output */
  bool jumpMagnOutput_;
};

} /* namespace Hts */ } /* end namespace Pita */

Pita::LinearDriver::Ptr linearReversiblePitaDriverNew(SingleDomainDynamic * pbDesc);

#endif /* PITA_HTS_REDUCEDLINEARDRIVERIMPL_H */
