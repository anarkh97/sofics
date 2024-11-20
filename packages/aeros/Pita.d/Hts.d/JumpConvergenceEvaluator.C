#include "JumpConvergenceEvaluator.h"

#include "../DynamStateOps.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>

namespace Pita { namespace Hts {

JumpConvergenceEvaluator::JumpConvergenceEvaluator(const String & name,
                                                   SliceMapping * mapping) :
  NamedTask(name),
  mapping_(mapping),
  localJump_()
{}

void
JumpConvergenceEvaluator::localJumpIs(HalfSliceRank rank, Seed * jump) {
  log() << "Identify Jump " << jump->name() << " at rank " << rank << "\n";
  localJump_[rank] = jump;
}

Seed *
JumpConvergenceEvaluator::localJump(HalfSliceRank rank) const {
  JumpMap::const_iterator it = localJump_.find(rank);
  return (it != localJump_.end()) ? it->second.ptr() : NULL;
}

TrivialConvergenceEvaluator::TrivialConvergenceEvaluator(SliceMapping * mapping) :
  JumpConvergenceEvaluator("Trivial Convergence", mapping)
{}

void
TrivialConvergenceEvaluator::iterationIs(IterationRank iter) {
  mapping_->convergedSlicesInc();
  log() << "First active slice is " << mapping_->firstActiveSlice() << "\n";
  setIteration(iter);
}

AccumulatedJumpConvergenceEvaluator::AccumulatedJumpConvergenceEvaluator(double targetRatio,
                                                                         const DynamOps * metric,
                                                                         SliceMapping * mapping,
                                                                         Communicator * timeComm) :
  JumpConvergenceEvaluator("Accumulated Jump Convergence", mapping),
  targetRatio_(targetRatio),
  metric_(metric),
  timeCommunicator_(timeComm),
  targetEstimate_(mapping->totalSlices().value(), 0.0),
  currentEstimate_(mapping->totalSlices().value(), 0.0),
  buffer_(mapping->totalFullSlices().value())
{}


void
AccumulatedJumpConvergenceEvaluator::iterationIs(IterationRank iter) {
  // 0) Numbering
  typedef std::vector<HalfSliceRank> SliceVector;
  SliceVector exchgToSlice(mapping_->activePrimalSlices().value()); // TODO member
  typedef std::map<HalfSliceRank, int> IndexMap;
  IndexMap sliceToExchg; // TODO member

  int cpuCount = mapping_->availableCpus().value();
  SimpleBuffer<int> recvs_counts(cpuCount); // TODO member
  std::fill_n(recvs_counts.array(), cpuCount, 0);

  HalfSliceRank firstJumpRank = mapping_->firstActiveSlice().next();
  SliceVector::iterator i_it = exchgToSlice.begin();
  for (int cpu = 0; cpu < cpuCount; ++cpu) {
    for (SliceMapping::SliceIterator s_it = mapping_->hostedSlice(CpuRank(cpu)); s_it; ++s_it) {
      HalfSliceRank rank = *s_it;
      if (rank >= mapping_->firstInactiveSlice()) {
        break;
      }
      if (rank <= mapping_->firstActiveSlice()) { 
        continue;
      }
      if ((rank.value() - firstJumpRank.value()) % 2 == 1) {
        currentEstimate_[rank.value()] = 0.0; 
      } else {
        ++recvs_counts[cpu];
        *i_it = rank;
        sliceToExchg[rank] = std::distance(exchgToSlice.begin(), i_it);
        ++i_it;
      }
    }
  }

  // 1) Compute local norms
  IterationRank jumpIter(iter.value() - 1);
  for (JumpMap::const_iterator it = localJump_.begin(); it != localJump_.end(); ++it) {
    if (jumpIter == it->second->iteration() && it->second->status() != Seed::INACTIVE) {
      buffer_[sliceToExchg[it->first]] = std::sqrt(energy(metric_.ptr(), it->second->state()));
    }
  }

  // 2) Exchange all norms (Allgatherv)
  SimpleBuffer<int> displacements(cpuCount); // TODO member
  displacements[0] = 0;
  std::partial_sum(recvs_counts.array(),
                   recvs_counts.array() + cpuCount - 1,
                   displacements.array() + 1);

  timeCommunicator_->allGatherv(buffer_.array(),
                                recvs_counts.array(),
                                displacements.array());

  // 3) Partial sums
  for (IndexMap::const_iterator it = sliceToExchg.begin(); it != sliceToExchg.end(); ++it) {
    currentEstimate_[it->first.value()] = buffer_[it->second];
  }

  std::vector<double>::iterator activeEstimateBegin = currentEstimate_.begin() + mapping_->firstActiveSlice().value();
  if (activeEstimateBegin != currentEstimate_.begin() && activeEstimateBegin != currentEstimate_.end()) {
    *activeEstimateBegin = *(activeEstimateBegin - 1);
  }
  std::partial_sum(activeEstimateBegin, currentEstimate_.end(), activeEstimateBegin);
  
  if (jumpIter == IterationRank(0)) {
    std::vector<double>::iterator last;
    last = std::transform(currentEstimate_.begin() + 1, // target is shifted back by one half-slice
                          currentEstimate_.end(),
                          targetEstimate_.begin(),
                          std::bind2nd(std::divides<double>(), targetRatio_));
    // TODO Replace HACK
    if (last != targetEstimate_.end()) {
      *last = *(last - 1) * 2.0 - *(last - 3);
    }
  }

  log() << "rank: estimate/target =\n";
  for (int i = 0; i < mapping_->totalSlices().value(); ++i) {
    log() << i << ": " << currentEstimate_[i] << "/" << targetEstimate_[i] << " ";
  }
  log() << "\n";

  // 4) Determine cvg level
  int cvgFront = mapping_->convergedSlices().value();
  int cvgFrontEnd = mapping_->firstInactiveSlice().value();
  while (cvgFront < cvgFrontEnd && currentEstimate_[cvgFront] <= targetEstimate_[cvgFront]) {
    ++cvgFront;
  }

  // 5) Apply to local jumps
  HalfSliceRank newFirstActiveSlice(cvgFront);
  
  for (JumpMap::iterator it = localJump_.begin(); it != localJump_.end(); ++it) {
    if (it->first > newFirstActiveSlice) {
      break;
    }
    if (it->second->status() == Seed::ACTIVE) {
      it->second->statusIs(Seed::CONVERGED);
    } 
  }

  // 6) Apply to mapping
  HalfSliceCount cvgIncrement = newFirstActiveSlice - mapping_->firstActiveSlice();
  mapping_->convergedSlicesInc(cvgIncrement);

  log() << "*** Convergence front at slice " << mapping_->firstActiveSlice() << "\n";

  setIteration(iter);
}

} /* end namespace Hts */ } /* end namespace Pita */
