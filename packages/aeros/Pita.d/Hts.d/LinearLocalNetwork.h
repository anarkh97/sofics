#ifndef PITA_HTS_LINEARLOCALNETWORK_H
#define PITA_HTS_LINEARLOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "SliceMapping.h"

#include "../Seed.h"
#include "../RemoteState.h"

#include "AffinePropagatorManager.h"
#include "ReducedCorrectionManager.h"
#include "JumpConvergenceEvaluator.h"
#include "../SeedDifferenceEvaluator.h"

#include <map>
#include <deque>
#include <cassert>

#include "LocalNetwork.h"
#include "LocalNetworkImpl.h"

namespace Pita { namespace Hts {

using namespace LocalNetworkImpl;

class LinearLocalNetwork : public LocalNetwork {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearLocalNetwork);
 
  virtual void statusIs(Status s) { /* TODO */ }

  /* Local network elements */
  TaskList halfTimeSlices() const;
  TaskList activeHalfTimeSlices() const;
  TaskList activeJumpAssemblers() const;
  TaskList activeLeftSeedSyncs() const;
  TaskList activeJumpProjectors() const;
  TaskList activeSeedAssemblers() const;

  TaskList activeFullTimeSlices() const;
  TaskList activeCorrectionSyncs() const;

  TaskList activeCoarseTimeSlices() const;
  TaskList activeFullCorrectionSyncs() const;
  
  SeedMap mainSeeds() const;
  MainSeedMap activeMainSeeds() const;

  LinearLocalNetwork(SliceMapping * mapping,
                     AffinePropagatorManager * propMgr,
                     ReducedCorrectionManager * redCorrMgr,
                     JumpConvergenceEvaluator * jumpCvgMgr,
                     RemoteState::Manager * commMgr,
                     LinSeedDifferenceEvaluator::Manager * jumpErrorMgr);

  void applyConvergenceStatus();

protected:
  void init();

  void buildForwardPropagation(HalfSliceRank sliceRank);
  void buildBackwardPropagation(HalfSliceRank sliceRank);
  void buildPrimalCorrectionNetwork(HalfSliceRank sliceRank);
  void buildDualCorrectionNetwork(HalfSliceRank sliceRank);
  
  virtual void buildCorrectionPropagator(HalfSliceRank sliceRank);
  virtual void buildCorrectionSynchronizationSend(HalfSliceRank sliceRank);
  virtual void buildCorrectionSynchronizationRecv(HalfSliceRank sliceRank);
  
  void buildPropagatedSeedSend(HalfSliceRank sliceRank);
  void buildPropagatedSeedRecv(HalfSliceRank sliceRank);
  void buildJumpEstimator(HalfSliceRank sliceRank);
  void buildJumpAssembly(HalfSliceRank sliceRank);
  void buildJumpProjection(HalfSliceRank sliceRank);
  void buildSeedUpdater(HalfSliceRank sliceRank);
  
  void buildReducedCorrectionPropagator(HalfSliceRank sliceRank);
  void buildFullCorrectionPropagator(HalfSliceRank sliceRank);
  void buildReducedCorrectionSynchronization(HalfSliceRank sliceRank);
  void buildFullCorrectionSynchronization(HalfSliceRank sliceRank);

  void eraseInactive(TaskMap & task) { task.erase(task.begin(), task.lower_bound(firstActiveSlice())); }
  TaskList getActive(const TaskMap task[2]) const;
  static TaskList getAll(const TaskMap & task);

private:
  AffinePropagatorManager::Ptr propMgr_;
  JumpBuilder::Manager::Ptr jumpBuildMgr_;
  JumpConvergenceEvaluator::Ptr jumpCvgMgr_;
  ReducedCorrectionManager::Ptr redCorrMgr_;
  
  TaskMap halfTimeSlice_[2];
  TaskMap jumpAssembler_[2];
  TaskMap jumpProjector_[2];
  TaskMap seedAssembler_[2];
  TaskMap leftSeedSync_[2];

  TaskMap correctionPropagator_[2];
  TaskMap correctionSync_[2];

  TaskMap alternateCorrectionPropagator_[2];
  TaskMap alternateCorrectionSync_[2];
  
  SeedMap mainSeed_;
  SeedMap seedCorrection_;
  ReducedSeedMap reducedSeedCorrection_;

  LinSeedDifferenceEvaluator::Manager::Ptr jumpErrorMgr_;
  
  CorrectionPropagatorBuilder<DynamState> fullCorrectionBuilder_;
  CorrectionPropagatorBuilder<Vector> reducedCorrectionBuilder_;

  DISALLOW_COPY_AND_ASSIGN(LinearLocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARLOCALNETWORK_H */
