#include "LinearProjectionNetwork.h"

#include "../DynamStateBasisWrapper.h"

#include "../DynamStateOps.h"
#include <Comm.d/Communicator.h>

#include <stdexcept>

#include <Timers.d/GetTime.h>

namespace Pita { namespace Hts {

LinearProjectionNetwork::LinearProjectionNetwork(size_t vSize,
                                                         Communicator * timeComm,
                                                         const SliceMapping * mapping,
                                                         const DynamOps * metric,
                                                         RankDeficientSolver * solver) :
  vectorSize_(vSize),
  timeCommunicator_(timeComm),
  mapping_(mapping),
  metric_(metric),
  gBuffer_(),
  mBuffer_(),
  mpiParameters_(),
  localBasis_(),
  metricBasis_(DynamStatePlainBasis::New(vectorSize_)),
  finalBasis_(DynamStatePlainBasis::New(vectorSize_)),
  originalMetricBasis_(DynamStatePlainBasis::New(vectorSize_)),
  originalFinalBasis_(DynamStatePlainBasis::New(vectorSize_)),
  normalMatrix_(),
  transmissionMatrix_(),
  reprojectionMatrix_(),
  solver_(solver),
  collector_(AffineBasisCollector::New()),
  globalExchangeNumbering_()
{}

void
LinearProjectionNetwork::prepareProjection() {
  GlobalExchangeNumbering::Ptr numbering = new GlobalExchangeNumbering(mapping_.ptr());
  globalExchangeNumbering_.push_back(numbering);
}

void
LinearProjectionNetwork::buildProjection() {
#ifndef NDEBUG
  double tic = getTime(), toc;
#endif /* NDEBUG*/

  const int previousMatrixSize = normalMatrix_.dim();

  // Build buffer
  const size_t stateSize = 2 * vectorSize_;
  const size_t targetBufferSize = stateSize * globalExchangeNumbering_.back()->stateCount();
  if (gBuffer_.size() < targetBufferSize) {
    gBuffer_.sizeIs(targetBufferSize);
  }

#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Build buffer: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/
  
  // Collect data from time-slices
  // 1) Final states
  AffineBasisCollector::CollectedState cs = collector_->firstForwardFinalState();
  while (cs.second.vectorSize() != 0) {
    const int inBufferRank = globalExchangeNumbering_.back()->globalIndex(HalfSliceId(cs.first, FORWARD));
    if (inBufferRank >= 0) {
      double * targetBuffer = gBuffer_.array() + (inBufferRank * stateSize);
      bufferStateCopy(cs.second, targetBuffer); 
      const int accumulatedIndex = inBufferRank + (2 * previousMatrixSize);
      localBasis_.insert(std::make_pair(accumulatedIndex, cs.second));
    }
    collector_->firstForwardFinalStateDel();
    cs = collector_->firstForwardFinalState();
  }

  // 2) Initial states
  cs = collector_->firstBackwardFinalState();
  while (cs.second.vectorSize() != 0) {
    const int inBufferRank = globalExchangeNumbering_.back()->globalIndex(HalfSliceId(cs.first, BACKWARD));
    if (inBufferRank >= 0) {
      double * targetBuffer = gBuffer_.array() + (inBufferRank * stateSize);
      if (metric_) {
        // Pre-multiply initial local states by the metric first
        mult(metric_.ptr(), cs.second, targetBuffer);
      } else {
        log() << "Warning, no metric found !\n";
        bufferStateCopy(cs.second, targetBuffer); 
      }
      const int accumulatedIndex = inBufferRank + (2 * previousMatrixSize);
      localBasis_.insert(std::make_pair(accumulatedIndex, cs.second));
    }
    collector_->firstBackwardFinalStateDel();
    cs = collector_->firstBackwardFinalState();
  }

#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Collect data from local time-slices: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/
  
  // Setup parameters for global MPI communication
  const int numCpus = mapping_->availableCpus().value(); 
  mpiParameters_.sizeIs(2 * numCpus);
  int * recv_counts = mpiParameters_.array();
  int * displacements = mpiParameters_.array() + numCpus;

  for (int i = 0; i < numCpus; ++i) {
    recv_counts[i] = globalExchangeNumbering_.back()->stateCount(CpuRank(i)) * stateSize;
  }

  displacements[0] = 0;
  for (int i = 1; i < numCpus; ++i) {
    displacements[i] = displacements[i-1] + recv_counts[i-1];
  }

  timeCommunicator_->allGatherv(gBuffer_.array(), recv_counts, displacements);
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Allgather: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  // Add new states to projection bases
  DynamStateBasisWrapper::Ptr receivedBasis = DynamStateBasisWrapper::New(vectorSize_, globalExchangeNumbering_.back()->stateCount(), gBuffer_.array());

  for (GlobalExchangeNumbering::IteratorConst it = globalExchangeNumbering_.back()->globalIndex(); it; ++it) {
    std::pair<Direction, int> p = *it;
    DynamStatePlainBasis * targetBasis = (p.first == BACKWARD) ?
                                        originalMetricBasis_.ptr() : // Initial state (premultiplied)
                                        originalFinalBasis_.ptr();   // Final state
    targetBasis->lastStateIs(receivedBasis->state(p.second));
  }
 
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Add new states: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/
 
  // Assemble normal and reprojection matrices in parallel (local rows)
  const int matrixSizeIncrement = globalExchangeNumbering_.back()->stateCount(BACKWARD);
  const int newMatrixSize = previousMatrixSize + matrixSizeIncrement;
 
  log() << "*** Matrix size: previous = " << previousMatrixSize << ", increment = " << matrixSizeIncrement << ", newSize = " << newMatrixSize << "\n";
  
  const size_t matrixBufferSize = (2 * newMatrixSize) * newMatrixSize; // For normalMatrix and reprojectionMatrix data
  if (mBuffer_.size() < matrixBufferSize) {
    mBuffer_.sizeIs(matrixBufferSize);
  }
  for (int i = 0; i < numCpus; ++i) {
    recv_counts[i] = 0;
    for (NumberingList::const_iterator it = globalExchangeNumbering_.begin(); it != globalExchangeNumbering_.end(); ++it) {
      recv_counts[i] += (*it)->stateCount(CpuRank(i));
    }
    recv_counts[i] *= newMatrixSize;
  }

  displacements[0] = 0;
  for (int i = 1; i < numCpus; ++i) {
    displacements[i] = displacements[i-1] + recv_counts[i-1];
  }
  
  // TODO Better efficiency: Do not recompute every coefs
  const int myCpu = timeCommunicator_->myID();
  double * rowBuffer = mBuffer_.array() + displacements[myCpu];
  for (LocalBasis::const_iterator it = localBasis_.begin();
      it != localBasis_.end(); ++it) {
    for (int i = 0; i < newMatrixSize; ++i) {
      rowBuffer[i] = originalMetricBasis_->state(i) * it->second;
    }
    rowBuffer += newMatrixSize;
  }
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Compute normal and reprojection matrices: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  // Exchange normal/reprojection matrix data  
  timeCommunicator_->allGatherv(mBuffer_.array(), recv_counts, displacements);
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Exchange normal matrix: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  // Assemble updated normal & reprojection matrices
  normalMatrix_.reSize(newMatrixSize);
  transmissionMatrix_.reSize(newMatrixSize);

  // Loop on iterations
  int originRowIndex = 0; // Target matrix base row index for current iteration 
  for (NumberingList::const_iterator it = globalExchangeNumbering_.begin(); it != globalExchangeNumbering_.end(); ++it) {
    GlobalExchangeNumbering::PtrConst numbering = *it;
    GlobalExchangeNumbering::IteratorConst jt_i = numbering->globalHalfIndex(BACKWARD);
    GlobalExchangeNumbering::IteratorConst jt_f = numbering->globalHalfIndex(FORWARD);
   
    // Loop on cpus
    for (int cpu = 0; cpu < numCpus; ++cpu) {
      const double * originBufferBegin = mBuffer_.array() + displacements[cpu]; // Position in AllgatherBuffer
      
      const int initialStateCountInIter = numbering->stateCount(CpuRank(cpu), BACKWARD);
      for (int s = 0; s < initialStateCountInIter; ++s) {
        int rowIndex = (*jt_i).second + originRowIndex; // Row in target normal matrix
        std::copy(originBufferBegin, originBufferBegin + newMatrixSize, normalMatrix_[rowIndex]);
        originBufferBegin += newMatrixSize;
        ++jt_i;
      }

      const int finalStateCountInIter = numbering->stateCount(CpuRank(cpu), FORWARD);
      for (int s = 0; s < finalStateCountInIter; ++s) {
        const int rowIndex = (*jt_f).second + originRowIndex; // Row in target reprojection matrix
        std::copy(originBufferBegin, originBufferBegin + newMatrixSize, transmissionMatrix_[rowIndex]);
        originBufferBegin += newMatrixSize;
        ++jt_f;
      }
      // Update buffer cpu base displacement for next iteration
      displacements[cpu] += newMatrixSize * (initialStateCountInIter + finalStateCountInIter);
    } 
    
    originRowIndex += numbering->stateCount(BACKWARD); // Udpate target matrix base row for next iteration
  }

#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Assemble normal matrix: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  solver_->transposedMatrixIs(normalMatrix_);

  reprojectionMatrix_.reSize(solver_->factorRank());
  metricBasis_->stateBasisDel();
  finalBasis_->stateBasisDel();

  for (int compactIndex = 0; compactIndex < solver_->factorRank(); ++compactIndex) {
    const int originalIndex = solver_->factorPermutation(compactIndex);

    for (int j = 0; j < solver_->factorRank(); ++j) {
      reprojectionMatrix_[compactIndex][j] = transmissionMatrix_[originalIndex][solver_->factorPermutation(j)];
    }

    metricBasis_->lastStateIs(originalMetricBasis_->state(originalIndex));
    finalBasis_->lastStateIs(originalFinalBasis_->state(originalIndex));
  }
  
  solver_->orderingIs(RankDeficientSolver::COMPACT);
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Factor: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  log() << "*** Projector rank = " << solver_->factorRank() << "\n";
  
}
  
// LinearProjectionNetwork::GlobalExchangeNumbering

LinearProjectionNetwork::GlobalExchangeNumbering::GlobalExchangeNumbering(const SliceMapping * mapping)
{
  this->initialize(mapping);
}

void
LinearProjectionNetwork::GlobalExchangeNumbering::initialize(const SliceMapping * mapping) {
  // Count only full timeslices
  const int fullSliceCount = std::max(0, (mapping->activeSlices().value() - 1) / 2);
  const int numCpus = mapping->availableCpus().value();

  // Avoid reallocations
  stateCount_.reserve(numCpus);
  initialStateCount_.reserve(numCpus);
  finalStateCount_.reserve(numCpus);
 
  stateId_.reserve(fullSliceCount * 2);
  initialStateId_.reserve(fullSliceCount);
  finalStateId_.reserve(fullSliceCount);

  // Determine location of the HalfTimeSlices updated since last correction
  typedef std::set<HalfSliceRank> HalfStateSet;
  typedef std::pair<HalfStateSet, HalfStateSet> FullStateSet; // pair<initialStates, finalStates>
  typedef std::vector<FullStateSet> TempMapping; 
  TempMapping currentStateMapping(mapping->availableCpus().value());
  
  {
    HalfSliceRank r(HalfSliceRank(1) + mapping->convergedSlices());
    for (int fts = 0; fts < fullSliceCount; ++fts) {
      // Initial states <=> Backward
      {
        const CpuRank cpu = mapping->hostCpu(r);
        currentStateMapping[cpu.value()].first.insert(r);
        r = r + HalfSliceCount(1); 
      }

      // Initial states <=> Forward
      {
        const CpuRank cpu = mapping->hostCpu(r);
        currentStateMapping[cpu.value()].second.insert(r);
        r = r + HalfSliceCount(1);
      }
    }
  }

  // Fill the member data structures
  size_t currentFullGlobalIndex = 0;
  size_t currentInitialGlobalIndex = 0;
  size_t currentFinalGlobalIndex = 0;

  for (TempMapping::const_iterator it = currentStateMapping.begin();
       it != currentStateMapping.end();
       ++it) {
    
    initialStateCount_.push_back(it->first.size());
    finalStateCount_.push_back(it->second.size());
    stateCount_.push_back(it->first.size() + it->second.size());

    // Initial states
    for (HalfStateSet::const_iterator jt = it->first.begin(); jt != it->first.end(); ++jt) {
      const HalfSliceId stateId(*jt, BACKWARD);
      
      initialGlobalIndex_.insert(std::make_pair(stateId, currentInitialGlobalIndex++));
      globalIndex_.insert(std::make_pair(stateId, currentFullGlobalIndex++));
      
      initialStateId_.push_back(stateId);
      stateId_.push_back(stateId);
    }
    
    // Final states
    for (HalfStateSet::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
      const HalfSliceId stateId(*jt, FORWARD);
      
      finalGlobalIndex_.insert(std::make_pair(stateId, currentFinalGlobalIndex++));
      globalIndex_.insert(std::make_pair(stateId, currentFullGlobalIndex++));
      
      finalStateId_.push_back(stateId);
      stateId_.push_back(stateId);
    }

  }

}

size_t
LinearProjectionNetwork::GlobalExchangeNumbering::stateCount() const {
  return stateId_.size();
}

size_t
LinearProjectionNetwork::GlobalExchangeNumbering::stateCount(Direction d) const {
  switch (d) {
    case NO_DIRECTION:
      return 0;
    case FORWARD:
      return finalStateId_.size();
    case BACKWARD:
      return initialStateId_.size();
  }
  throw Fwk::InternalException("in LinearProjectionNetwork::GlobalExchangeNumbering::stateCount");
}

size_t
LinearProjectionNetwork::GlobalExchangeNumbering::stateCount(CpuRank c) const {
  try {
    return stateCount_.at(c.value());
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return 0; 
}

size_t
LinearProjectionNetwork::GlobalExchangeNumbering::stateCount(CpuRank c, Direction d) const {
  try {
    switch (d) {
      case NO_DIRECTION:
        return 0;
      case FORWARD:
        return finalStateCount_.at(c.value());
      case BACKWARD:
        return initialStateCount_.at(c.value());
      default:
        throw Fwk::InternalException("in LinearProjectionNetwork::GlobalExchangeNumbering::stateCount");
    }
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return 0;
}

int
LinearProjectionNetwork::GlobalExchangeNumbering::globalIndex(const HalfSliceId & id) const {
  IndexMap::const_iterator it = globalIndex_.find(id);
  return (it != globalIndex_.end()) ? static_cast<int>(it->second) : -1; 
}

LinearProjectionNetwork::GlobalExchangeNumbering::IteratorConst
LinearProjectionNetwork::GlobalExchangeNumbering::globalIndex() const {
  return IteratorConst(this->globalIndex_.begin(), this->globalIndex_.end());
}

int
LinearProjectionNetwork::GlobalExchangeNumbering::globalHalfIndex(const HalfSliceId & id) const {
  int result = -1;
  
  IndexMap::const_iterator it;
  switch (id.direction()) {
    case NO_DIRECTION:
      break;
    case FORWARD:
      it = finalGlobalIndex_.find(id);
      if (it != finalGlobalIndex_.end())
        result = static_cast<int>(it->second);
      break;
    case BACKWARD:
      it = initialGlobalIndex_.find(id);
      if (it != initialGlobalIndex_.end())
        result = static_cast<int>(it->second);
  }

  return result;
}

LinearProjectionNetwork::GlobalExchangeNumbering::IteratorConst
LinearProjectionNetwork::GlobalExchangeNumbering::globalHalfIndex(Direction d) const {
  switch (d) {
    case NO_DIRECTION:
      break;
    case FORWARD:
      return IteratorConst(this->finalGlobalIndex_.begin(), this->finalGlobalIndex_.end());
    case BACKWARD:
      return IteratorConst(this->initialGlobalIndex_.begin(), this->initialGlobalIndex_.end());
  }
  return IteratorConst();
}

HalfSliceId
LinearProjectionNetwork::GlobalExchangeNumbering::stateId(int gfi) const {
  try {
    return stateId_.at(gfi);
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return HalfSliceId();
}

HalfSliceId
LinearProjectionNetwork::GlobalExchangeNumbering::stateId(int ghi, Direction d) const {
  try {
    switch (d) {
      case NO_DIRECTION:
        break;
      case FORWARD:
        return finalStateId_.at(ghi);
      case BACKWARD:
        return initialStateId_.at(ghi);
    }
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return HalfSliceId();
}

} /* end namespace Hts */ } /* end namespace Pita */
