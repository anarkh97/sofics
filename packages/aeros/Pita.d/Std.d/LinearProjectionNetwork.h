#ifndef PITA_STD_LINEARPROJECTIONNETWORK_H
#define PITA_STD_LINEARPROJECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "../DynamOps.h"

class Communicator;
#include "SliceMapping.h"

#include "../DynamStatePlainBasis.h"
#include "../RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

#include "AffineBasisCollector.h"

#include "../NamedTask.h"

#include <list>

namespace Pita { namespace Std {

class LinearProjectionNetwork : public Fwk::PtrInterface<LinearProjectionNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearProjectionNetwork);

  // Data collection
  AffineBasisCollector * collector() const { return collector_.ptr(); }

  // Projection operators
  size_t reducedBasisSize() const { return projectionBasis_->stateCount(); }
  const DynamOps * metric() const { return metric_.ptr(); }
  const DynamStateBasis * projectionBasis() const { return projectionBasis_.ptr(); }
  const DynamStateBasis * propagatedBasis() const { return propagatedBasis_.ptr(); }
  const RankDeficientSolver * normalMatrixSolver() const { return normalMatrixSolver_.ptr(); }
  const FullSquareMatrix * reprojectionMatrix() const { return &reprojectionMatrix_; }

  // Execution
  NamedTask::Ptr projectionTaskNew();

  // Notification
  class NotifieeConst : public Fwk::BaseMultiNotifiee<const LinearProjectionNetwork, NotifieeConst> {
  public:
    EXPORT_PTRINTERFACE_TYPES(NotifieeConst);

    virtual void onProjectionOperators() {}

  protected:
    explicit NotifieeConst(const LinearProjectionNetwork * notifier = NULL) :
      Fwk::BaseMultiNotifiee<const LinearProjectionNetwork, NotifieeConst>(notifier)
    {}
  };

  void lastNotifieeIs(NotifieeConst * notifiee) const { notifierDelegate().lastNotifieeIs(notifiee); }
  void notifieeDel(NotifieeConst * notifiee) const { notifierDelegate().notifieeDel(notifiee); }

  // Debug
  friend OStream & operator<<(OStream &, const LinearProjectionNetwork &);

  static Ptr New(const SliceMapping * mapping, Communicator * timeComm, const DynamOps * metric, size_t vecSize, RankDeficientSolver * solver) {
    return new LinearProjectionNetwork(mapping, timeComm, metric, vecSize, solver);
  }

protected:
  LinearProjectionNetwork(const SliceMapping * mapping, Communicator * timeComm, const DynamOps * metric, size_t vecSize, RankDeficientSolver * solver);
  ~LinearProjectionNetwork();
 
  // Numbering for Allgather 
  class StateExchgNumbering;
  class MatrixExchgNumbering;

  // Notification
  GenNotifierDelegate<NotifieeConst> & notifierDelegate() const { return const_cast<LinearProjectionNetwork *>(this)->notifierDelegate_; }
  
  // Execution
  class Task;
  friend class Task;
  void buildProjection();

private:
  // Problem characteristics 
  size_t vectorSize_;
  DynamOps::PtrConst metric_;

  // Network & Communication
  SliceMapping::PtrConst mapping_;
  Fwk::Ptr<MatrixExchgNumbering> numbering_;
  Communicator * timeCommunicator_;

  // Projection operators
  DynamStatePlainBasis::Ptr projectionBasis_; 
  DynamStatePlainBasis::Ptr propagatedBasis_;
  RankDeficientSolver::Ptr normalMatrixSolver_;
  FullSquareMatrix reprojectionMatrix_;

  // Rank-deficient operators
  DynamStatePlainBasis::Ptr originalProjectionBasis_; 
  DynamStatePlainBasis::Ptr originalPropagatedBasis_; 
  FullSquareMatrix normalMatrix_;
  FullSquareMatrix transmissionMatrix_;

  // Local data collection
  AffineBasisCollector::Ptr collector_;
  
  enum Kind {
    INITIAL = 0,
    FINAL = 1
  };
  friend OStream & operator<<(OStream &, LinearProjectionNetwork::Kind);

  typedef GenId<Kind> IterStateId;
  friend OStream & operator<<(OStream &, const LinearProjectionNetwork::IterStateId &);
  
  class StateId {
  public:
    IterationRank iteration() const { return iteration_; }
    SliceRank slice() const { return inIterId_.rank(); }
    Kind type() const { return inIterId_.type(); }

    StateId(IterationRank iter, IterStateId inIterId) :
      iteration_(iter), inIterId_(inIterId)
    {}

    bool operator==(const StateId & other) const {
      return iteration() == other.iteration() && inIterId_ == other.inIterId_;
    }

    bool operator<(const StateId & other) const {
      return iteration() == other.iteration() ? inIterId_ < other.inIterId_ : iteration() < other.iteration();
    }

  private:
    IterationRank iteration_;
    IterStateId inIterId_;
  };
  friend OStream & operator<<(OStream &, const LinearProjectionNetwork::StateId &);
  
  typedef std::map<StateId, DynamState> LocalStateMap;
  LocalStateMap localState_; // Accumulated index to local state

  void print(OStream &) const;
  
  // Notifier implementation
  GenNotifierDelegate<NotifieeConst> notifierDelegate_;

  DISALLOW_COPY_AND_ASSIGN(LinearProjectionNetwork);
};
  
OStream &
operator<<(OStream &, const LinearProjectionNetwork &);

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_LINEARPROJECTIONNETWORK_H */
