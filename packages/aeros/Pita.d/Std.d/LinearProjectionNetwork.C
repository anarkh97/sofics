#include "LinearProjectionNetwork.h"

#include "../DynamStateOps.h"
#include "../DynamStateBasisWrapper.h"

#include <Comm.d/Communicator.h>

#include <map>
#include <vector>
#include <numeric>

namespace Pita { namespace Std {

inline
OStream & operator<<(OStream & out, LinearProjectionNetwork::Kind k) {
  char c = (k == LinearProjectionNetwork::INITIAL) ? 'i' : 'f';
  return out << c;
}

inline
OStream & operator<<(OStream & out, const LinearProjectionNetwork::IterStateId & isi) {
  return out << isi.rank() << isi.type();
}

inline
OStream & operator<<(OStream & out, const LinearProjectionNetwork::StateId & si) {
  return out << "{" << si.iteration() << "," << si.slice() << si.type() << "}";
}

// Auxiliary class for in-iteration global indexing

class LinearProjectionNetwork::StateExchgNumbering : public Fwk::PtrInterface<StateExchgNumbering> {
public:
  EXPORT_PTRINTERFACE_TYPES(StateExchgNumbering);

  CpuCount cpuCount() const;

  // State count
  int stateCount() const;
  int stateCount(CpuRank c) const;
  int stateCount(CpuRank c, Kind k) const;
  
private: 
  typedef std::map<IterStateId, int> IndexMap;

public:
  class IndexIterator {
  public:
    int operator*() const { return it_->second; }
    IndexIterator & operator++() { ++it_; return *this; }
    IndexIterator operator++(int) { IndexIterator tmp(*this); this->operator++(); return tmp; }
    operator bool() const { return it_ != endIt_; }
   
    IndexIterator() {
      it_ = endIt_;
    }

  protected:
    typedef IndexMap::const_iterator ItImpl;

    explicit IndexIterator(ItImpl beginIt, ItImpl endIt) : 
      it_(beginIt),
      endIt_(endIt)
    {}

    friend class StateExchgNumbering;

  private:
    ItImpl it_;
    ItImpl endIt_;
  };

  // Numbering for all states (both initial AND final)
  int index(IterStateId id) const;
  IndexIterator index() const; // Indices sorted by increasing IterStateId = (SliceRank, Kind)
  IterStateId stateId(int idx) const;
  
  explicit StateExchgNumbering(const SliceMapping * m) {
    initialize(m);
  }

private:
  void initialize(const SliceMapping * m);

  // CpuRank to local State count
  std::vector<int> stateCount_;
 
  // IterStateId to index 
  IndexMap index_;
  
  // Index to IterStateId
  std::vector<IterStateId> stateId_;

  DISALLOW_COPY_AND_ASSIGN(StateExchgNumbering);
};

void
LinearProjectionNetwork::StateExchgNumbering::initialize(const SliceMapping * mapping) {
  int cpuCount = mapping->availableCpus().value();
  int sliceCount = mapping->activeSlices().value();

  // Avoid reallocations
  stateCount_.reserve(cpuCount);
 
  stateId_.reserve(2 * sliceCount);

  // Determine location of the TimeSlices updated since last correction
  typedef std::set<SliceRank> SliceSet;
  typedef std::vector<SliceSet> TempMapping; 
  TempMapping currentStateMapping(mapping->availableCpus().value());

  for (SliceRank sr = mapping->firstActiveSlice();
       sr < mapping->firstInactiveSlice();
       sr = sr.next()) {
    CpuRank cpu = mapping->hostCpu(sr);
    currentStateMapping[cpu.value()].insert(sr);
  }

  // Fill the member data structures
  int currIndex = 0;

  for (TempMapping::const_iterator it = currentStateMapping.begin();
       it != currentStateMapping.end();
       ++it) {
    
    stateCount_.push_back(2 * it->size());

    // Initial states
    for (SliceSet::const_iterator jt = it->begin(); jt != it->end(); ++jt) {
      IterStateId i_stateId(INITIAL, *jt);
      index_.insert(std::make_pair(i_stateId, currIndex++));
      stateId_.push_back(i_stateId);
    }
      
    // Final states
    for (SliceSet::const_iterator jt = it->begin(); jt != it->end(); ++jt) {
      IterStateId f_stateId(FINAL, *jt);
      index_.insert(std::make_pair(f_stateId, currIndex++));
      stateId_.push_back(f_stateId);
    }
  }
}

CpuCount
LinearProjectionNetwork::StateExchgNumbering::cpuCount() const {
  return CpuCount(stateCount_.size());
}

int
LinearProjectionNetwork::StateExchgNumbering::stateCount() const {
  return stateId_.size();
}

int
LinearProjectionNetwork::StateExchgNumbering::stateCount(CpuRank c) const {
  return stateCount_[c.value()];
}

int
LinearProjectionNetwork::StateExchgNumbering::stateCount(CpuRank c, Kind k) const {
  return stateCount(c) / 2; 
}

int
LinearProjectionNetwork::StateExchgNumbering::index(IterStateId id) const {
  IndexMap::const_iterator it = index_.find(id);
  return (it != index_.end()) ? static_cast<int>(it->second) : -1; 
}


LinearProjectionNetwork::StateExchgNumbering::IndexIterator
LinearProjectionNetwork::StateExchgNumbering::index() const {
  return IndexIterator(index_.begin(), index_.end());
}

LinearProjectionNetwork::IterStateId
LinearProjectionNetwork::StateExchgNumbering::stateId(int index) const {
  return stateId_[index];
}


// Auxiliary class for all-iteration exchange numbering

class LinearProjectionNetwork::MatrixExchgNumbering : public Fwk::PtrInterface<MatrixExchgNumbering> {
public:
  EXPORT_PTRINTERFACE_TYPES(MatrixExchgNumbering);

  int iterationCount() const;
  const StateExchgNumbering * iterNumbering(IterationRank iter) const;
  const StateExchgNumbering * lastIterNumbering() const;
  void lastIterNumberingIs(const StateExchgNumbering * iterNumbering);

  int stateCount() const;
  int stateCount(CpuRank cpu) const;
  int stateCount(CpuRank cpu, Kind k) const;

private:
  typedef std::map<StateId, int> IndexMap;

public:
  class IndexIterator {
  public:
    int operator*() const { return it_->second; }
    IndexIterator & operator++() { ++it_; return *this; }
    IndexIterator operator++(int) { IndexIterator tmp(*this); this->operator++(); return tmp; }
    operator bool() const { return it_ != endIt_; }
   
    IndexIterator() {
      it_ = endIt_;
    }

  protected:
    typedef IndexMap::const_iterator ItImpl;

    explicit IndexIterator(ItImpl beginIt, ItImpl endIt) : 
      it_(beginIt),
      endIt_(endIt)
    {}

    friend class MatrixExchgNumbering;

  private:
    ItImpl it_;
    ItImpl endIt_;
  };

  int index(StateId id) const;  
  IndexIterator index() const; // Indices sorted by increasing StateId = (Iteration, StateId) 
  StateId stateId(int index) const;

  explicit MatrixExchgNumbering(CpuCount cpus);

private:
  void update();
  void appendState(StateId id);

  typedef std::vector<StateExchgNumbering::PtrConst> ExchgHistory;
  ExchgHistory history_;

  // CpuRank to local State count
  std::vector<int> stateCount_;
  
  // StateId to index
  IndexMap index_; 
  
  // Index to StateId
  std::vector<StateId> stateId_;
};

LinearProjectionNetwork::MatrixExchgNumbering::MatrixExchgNumbering(CpuCount cpus) :
  stateCount_(cpus.value(), 0)
{}

int
LinearProjectionNetwork::MatrixExchgNumbering::iterationCount() const {
  return history_.size();
}

const LinearProjectionNetwork::StateExchgNumbering *
LinearProjectionNetwork::MatrixExchgNumbering::iterNumbering(IterationRank iter) const {
  return history_[iter.value()].ptr();
}

const LinearProjectionNetwork::StateExchgNumbering *
LinearProjectionNetwork::MatrixExchgNumbering::lastIterNumbering() const {
  return !history_.empty() ? history_.back().ptr() : NULL;
}

LinearProjectionNetwork::MatrixExchgNumbering::IndexIterator
LinearProjectionNetwork::MatrixExchgNumbering::index() const {
  return IndexIterator(index_.begin(), index_.end());
}

void
LinearProjectionNetwork::MatrixExchgNumbering::lastIterNumberingIs(const LinearProjectionNetwork::StateExchgNumbering * iterNumbering) {
  IterationRank iteration(history_.size());
  int cpuCount = iterNumbering->cpuCount().value();
  
  std::vector<StateId> prevStateId;
  prevStateId.swap(stateId_);
  stateId_.reserve(prevStateId.size() + iterNumbering->stateCount());

  int prevStateBegin = 0;
  int addedStateBegin = 0;
  for (int cpu = 0; cpu < cpuCount; ++cpu) {
    // Previous states
    int prevStateCount = stateCount_[cpu];
    int prevStateEnd = prevStateBegin + prevStateCount;
    for (int ps = prevStateBegin; ps < prevStateEnd; ++ps) {
      StateId id = prevStateId[ps];
      appendState(id);
    }
    prevStateBegin = prevStateEnd;

    // Current states
    int addedStateCount = iterNumbering->stateCount(CpuRank(cpu));
    int addedStateEnd = addedStateBegin + addedStateCount;
    for (int as = addedStateBegin; as < addedStateEnd; ++as) {
      StateId id(iteration, iterNumbering->stateId(as));
      appendState(id);
    }
    addedStateBegin = addedStateEnd;

    stateCount_[cpu] += addedStateCount;
  }
  
  history_.push_back(iterNumbering);
}


void
LinearProjectionNetwork::MatrixExchgNumbering::appendState(LinearProjectionNetwork::StateId id) {
  int newIndex = stateId_.size();
  stateId_.push_back(id);
  index_[id] = newIndex; 
}

int
LinearProjectionNetwork::MatrixExchgNumbering::index(LinearProjectionNetwork::StateId id) const {
  IndexMap::const_iterator it = index_.find(id);
  return (it != index_.end()) ? static_cast<int>(it->second) : -1; 
}

LinearProjectionNetwork::StateId
LinearProjectionNetwork::MatrixExchgNumbering::stateId(int index) const {
  return stateId_[index];
}

int
LinearProjectionNetwork::MatrixExchgNumbering::stateCount() const {
  return stateId_.size();
}

int
LinearProjectionNetwork::MatrixExchgNumbering::stateCount(CpuRank c) const {
  return stateCount_[c.value()];
}

int
LinearProjectionNetwork::MatrixExchgNumbering::stateCount(CpuRank c, Kind k) const {
  return stateCount(c) / 2; 
}

OStream &
operator<<(OStream & out, const LinearProjectionNetwork & n) {
  n.print(out);
  return out;
}

void
LinearProjectionNetwork::print(OStream & out) const {
  out << "Mapping:\n" << *mapping_;
  
  const StateExchgNumbering * numbering = numbering_->lastIterNumbering();
  out << "State Numbering:\n";
  out << "State count: Total = " << numbering->stateCount() << " / ";
  for (CpuRank c(0); c < CpuRank(0) + mapping_->availableCpus(); c = c + CpuCount(1)) {
    out << numbering->stateCount(c) << " ";
  }
  out << "\n";

  out << "Exchange indices (I and F):\n";
  for (SliceRank sr = mapping_->firstActiveSlice(); sr < mapping_->firstInactiveSlice(); sr = sr.next()) {
    out << IterStateId(INITIAL, sr) << ":" << numbering->index(IterStateId(INITIAL, sr)) << " "
        << IterStateId(FINAL, sr)   << ":" << numbering->index(IterStateId(FINAL, sr))   << " ";
  }
  out << "\n";

  out << "Matrix Numbering:\n";
  out << "State count: Total = " << numbering_->stateCount() << " / ";
  for (CpuRank c(0); c < CpuRank(0) + mapping_->availableCpus(); c = c + CpuCount(1)) {
    out << numbering_->stateCount(c) << " ";
  }
  out << "\n";
  out << "Exchange indices (I and F):\n";
  for (int index = 0; index < numbering_->stateCount(); ++index) {
    out << numbering_->stateId(index) << ":" << numbering_->index(numbering_->stateId(index)) << " ";
  }
  out << "\n";
}

// Constructor

LinearProjectionNetwork::LinearProjectionNetwork(
    const SliceMapping * mapping,
    Communicator * timeComm,
    const DynamOps * metric,
    size_t vecSize,
    RankDeficientSolver * solver) :
  vectorSize_(vecSize),
  metric_(metric),
  mapping_(mapping),
  numbering_(new MatrixExchgNumbering(mapping->availableCpus())),
  timeCommunicator_(timeComm),
  projectionBasis_(DynamStatePlainBasis::New(vecSize)),
  propagatedBasis_(DynamStatePlainBasis::New(vecSize)),
  normalMatrixSolver_(solver),
  reprojectionMatrix_(),
  originalProjectionBasis_(DynamStatePlainBasis::New(vecSize)),
  originalPropagatedBasis_(DynamStatePlainBasis::New(vecSize)),
  normalMatrix_(),
  transmissionMatrix_(),
  collector_(AffineBasisCollector::New()),
  localState_()
{}

// Implementation

void
LinearProjectionNetwork::buildProjection() {
  // Prelude 
  StateExchgNumbering::Ptr iterNumbering = new StateExchgNumbering(mapping_.ptr());
  numbering_->lastIterNumberingIs(iterNumbering.ptr());
  IterationRank iteration(numbering_->iterationCount() - 1);

  int previousMatrixSize = normalMatrix_.dim();
  int previousStateCount = 2 * previousMatrixSize;
  int cpuCount = mapping_->availableCpus().value();

  // 1) Collect new local states
#ifndef NDEBUG
  log() << "^^ Collecting new local states\n";
#endif /* NDEBUG */

  size_t stateSize = 2 * vectorSize_;
  size_t targetBufferSize = stateSize * iterNumbering->stateCount();
  SimpleBuffer<double> stateBuffer(targetBufferSize);

  AffineBasisCollector::CollectedState cs;

  //    a) Initial states to be premultiplied
  cs = collector_->firstInitialState();
  while (cs.state.vectorSize() != 0) {
    IterStateId iterState(INITIAL, cs.sliceId);
    int inBufferRank = iterNumbering->index(iterState);
    if(inBufferRank >= 0) {
      double * inBufferAddrBegin = stateBuffer.array() + (inBufferRank * stateSize);
      if (metric_) {
        mult(metric_.ptr(), cs.state, inBufferAddrBegin);
      } else {
        bufferStateCopy(cs.state, inBufferAddrBegin);
      }
      StateId id(iteration, iterState);
      localState_.insert(std::make_pair(id, cs.state));
    }
    collector_->firstInitialStateDel();
    cs = collector_->firstInitialState();
  }
  
  //    b) Final states as is
  cs = collector_->firstFinalState();
  while (cs.state.vectorSize() != 0) {
    IterStateId iterState(FINAL, cs.sliceId);
    int inBufferRank = iterNumbering->index(iterState);
    if(inBufferRank >= 0) {
      double * inBufferAddrBegin = stateBuffer.array() + (inBufferRank * stateSize);
      bufferStateCopy(cs.state, inBufferAddrBegin);
      StateId id(iteration, iterState);
      localState_.insert(std::make_pair(id, cs.state));
    }
    collector_->firstFinalStateDel();
    cs = collector_->firstFinalState();
  }

  // 2) Consolidate new states
#ifndef NDEBUG
  log() << "^^ Consolidating new local states\n";
#endif /* NDEBUG */
  SimpleBuffer<int> recvs_counts(cpuCount);
  SimpleBuffer<int> displacements(cpuCount);
 
  for (int cpu = 0; cpu < cpuCount; ++cpu) {
    recvs_counts[cpu] = iterNumbering->stateCount(CpuRank(cpu)) * stateSize;
  }
  displacements[0] = 0;
  std::partial_sum(recvs_counts.array(),
                   recvs_counts.array() + cpuCount - 1,
                   displacements.array() + 1);

  timeCommunicator_->allGatherv(stateBuffer.array(),
                                recvs_counts.array(),
                                displacements.array());

  DynamStateBasisWrapper::Ptr receivedStates = DynamStateBasisWrapper::New(
      vectorSize_,
      iterNumbering->stateCount(),
      stateBuffer.array());

  // Get in-iteration exchange indices by increasing SliceRank / Kind
  for (StateExchgNumbering::IndexIterator it = iterNumbering->index(); it; ++it) {
    int index = *it;
    IterStateId id = iterNumbering->stateId(index);
    DynamStatePlainBasis * targetBasis = (id.type() == INITIAL) ?
                                         originalProjectionBasis_.ptr() :
                                         originalPropagatedBasis_.ptr();
    targetBasis->lastStateIs(receivedStates->state(index));
  }

  // 3) Assemble [new] local ROWS of the TRANSPOSED reduced operators
#ifndef NDEBUG
  log() << "^^ Assembling local rows\n";
#endif /* NDEBUG */
  int matrixSizeIncrement = iterNumbering->stateCount() / 2;
  int newMatrixSize = previousMatrixSize + matrixSizeIncrement;

  log() << "*** Matrix size: previous = " << previousMatrixSize
        << ", increment = " << matrixSizeIncrement
        << ", newSize = " << newMatrixSize << "\n";

  size_t matrixBufferSize = 2 * (newMatrixSize * newMatrixSize); // Buffer covers 2 matrices
  SimpleBuffer<double> matrixBuffer(matrixBufferSize);
  
  // Fill buffer with dot products:
  // -> Row ordering = allgather-exchange (cpu / kind / slice)
  // -> Col ordering = rank-deficient operators (iteration / slice)
  // TODO: Efficiency: Do not recompute already known coefs
  for (LocalStateMap::const_iterator row_it = localState_.begin();
       row_it != localState_.end();
       ++row_it) {
    int row_index = numbering_->index(row_it->first);
    double * inBufferAddrBegin = matrixBuffer.array() + (row_index * newMatrixSize);
    for (DynamStatePlainBasis::IteratorConst col_it = originalProjectionBasis_->state(); col_it; ++col_it) {
      *inBufferAddrBegin++ = (*col_it) * row_it->second;
    }
  }

  // 4) Consolidate [new] rows of reduced operators
#ifndef NDEBUG
  log() << "^^ Consolidating reduced operators\n";
#endif /* NDEBUG */

  for (int cpu = 0; cpu < cpuCount; ++cpu) {
    recvs_counts[cpu] = numbering_->stateCount(CpuRank(cpu)) * newMatrixSize; // (# rows) * (size of row)
  }
  displacements[0] = 0;
  std::partial_sum(recvs_counts.array(),
                   recvs_counts.array() + cpuCount - 1,
                   displacements.array() + 1);

  timeCommunicator_->allGatherv(matrixBuffer.array(),
                                recvs_counts.array(), 
                                displacements.array());

  // 5) Assemble rank-deficient operators
#ifndef NDEBUG
  log() << "^^ Assembling reduced operators\n";
#endif /* NDEBUG */
  normalMatrix_.reSize(newMatrixSize);
  transmissionMatrix_.reSize(newMatrixSize);

  int initialStateIndex = 0;
  int finalStateIndex = 0;

  const double * originBuffer = matrixBuffer.array();
  for (MatrixExchgNumbering::IndexIterator idx_it = numbering_->index(); idx_it; ++idx_it) {
    // Source
    int index = *idx_it;
    const double * originBufferRowBegin = originBuffer + index * (newMatrixSize);
    const double * originBufferRowEnd = originBufferRowBegin + newMatrixSize;
    
    // Target
    StateId id = numbering_->stateId(index);
    Kind kind = id.type();
    double * targetBufferRowBegin;
    if (kind == INITIAL) {
      targetBufferRowBegin = normalMatrix_[initialStateIndex++];
    } else {
      targetBufferRowBegin = transmissionMatrix_[finalStateIndex++];
    }
    
    std::copy(originBufferRowBegin, originBufferRowEnd, targetBufferRowBegin); 
  }

  // 6) Factor rank-deficient normal matrix
#ifndef NDEBUG
  log() << "^^ Factorizing rank-deficient matrix\n";
#endif /* NDEBUG */
  normalMatrixSolver_->transposedMatrixIs(normalMatrix_);

  // 7) Assemble full-rank compact operators and bases
#ifndef NDEBUG
  log() << "^^ Assembling full-rank operators\n";
#endif /* NDEBUG */
  int numericalRank = normalMatrixSolver_->factorRank();
  log() << "*** Rank = " << numericalRank << "\n";
  reprojectionMatrix_.reSize(numericalRank);
  projectionBasis_->stateBasisDel();
  propagatedBasis_->stateBasisDel();

  for (int compactIndex = 0; compactIndex < numericalRank; ++compactIndex) {
    int originalIndex = normalMatrixSolver_->factorPermutation(compactIndex);

    for (int j = 0; j < numericalRank; ++j) {
      int p = normalMatrixSolver_->factorPermutation(j);
      reprojectionMatrix_[compactIndex][j] = transmissionMatrix_[originalIndex][p];
    }

    projectionBasis_->lastStateIs(originalProjectionBasis_->state(originalIndex));
    propagatedBasis_->lastStateIs(originalPropagatedBasis_->state(originalIndex));
  }

  normalMatrixSolver_->orderingIs(RankDeficientSolver::COMPACT);
  
  // Perform notification
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onProjectionOperators);
}

class LinearProjectionNetwork::Task : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(Task);

  virtual void iterationIs(IterationRank iter);

  Task(LinearProjectionNetwork * parent) :
    NamedTask("Reduced Projection Assembly"),
    parent_(parent)
  {}

private:
    LinearProjectionNetwork::Ptr parent_;
};

void
LinearProjectionNetwork::Task::iterationIs(IterationRank iter) {
  assert(iter > iteration());
  parent_->buildProjection();
  setIteration(iter);
}

NamedTask::Ptr 
LinearProjectionNetwork::projectionTaskNew() {
  return new Task(this);
}

LinearProjectionNetwork::~LinearProjectionNetwork() {}

} /* end namespace Std */ } /* end namespace Pita */
