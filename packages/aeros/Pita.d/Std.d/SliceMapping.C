#include "SliceMapping.h"

#include <algorithm>

namespace Pita { namespace Std {

SliceMapping::SliceMapping(SliceCount totalSlices, CpuCount availableCpus, SliceCount maxWorkload) :
  loadBalancer_(LoadBalancer::New(totalSlices.value(), availableCpus.value(), maxWorkload.value()))
{}

SliceCount
SliceMapping::totalSlices() const {
  return SliceCount(loadBalancer_->totalTasks());
}

CpuCount
SliceMapping::availableCpus() const {
  return CpuCount(loadBalancer_->availableWorkers());
}


SliceCount
SliceMapping::maxWorkload() const {
  return SliceCount(loadBalancer_->maxWorkload());
}

SliceCount
SliceMapping::activeSlices() const {
  return SliceCount(loadBalancer_->currentGlobalWorkload());
}

SliceCount
SliceMapping::activeSlices(CpuRank cpu) const {
  return SliceCount(loadBalancer_->currentWorkload(cpu.value()));
}

SliceRank
SliceMapping::firstActiveSlice() const {
  return SliceRank(loadBalancer_->firstCurrentTask());
}

SliceRank
SliceMapping::firstInactiveSlice() const {
  return SliceRank(loadBalancer_->firstWaitingTask());
}

SliceCount
SliceMapping::convergedSlices() const {
  return SliceCount(loadBalancer_->completedTasks());
}

void
SliceMapping::SliceMapping::convergedSlicesInc(SliceCount increment) {
  loadBalancer_->completedTasksInc(increment.value());
}

CpuRank
SliceMapping::hostCpu(SliceRank slice) const {
  return CpuRank(loadBalancer_->worker(slice.value()));
}

SliceMapping::SliceIterator
SliceMapping::hostedSlice(CpuRank cpu) const {
  return SliceIterator(loadBalancer_->tasks(cpu.value()));
}


} /* end namespace Std */ } /* end namespace Pita */
