#ifndef PITA_HTS_NLPROJECTIONNETWORK_H
#define PITA_HTS_NLPROJECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"
#include "HalfSliceId.h"

#include "../NlDynamOps.h"
#include "../DynamStatePlainBasis.h"
#include "../NearSymmetricSolver.h"
#include "../BasisTransform.h"
#include "CorrectionReductor.h"
#include "CorrectionReconstructor.h"

#include "../ConcurrentBasisManager.h"
#include "GlobalStateSharing.h"

namespace Pita { namespace Hts {

// Root Containers / Managers

class MidBasisManager : public Fwk::GenManager<const DynamStateBasis, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(MidBasisManager);

  const DynamStateBasis * commonBasis() { return commonBasis_.ptr(); }
  void commonBasisIs(const DynamStateBasis * commonBasis) { commonBasis_ = commonBasis; } // TODO Update the instances !

  explicit MidBasisManager(const DynamStateBasis * commonBasis) :
    commonBasis_(commonBasis)
  {}

protected:
  virtual DynamStateBasis::PtrConst createNewInstance(const HalfSliceRank & key) {
    return commonBasis_;
  }

private:
  DynamStateBasis::PtrConst commonBasis_;
};


// Internal Repositories 

class SimpleBasisManager : public Fwk::GenManager<DynamStatePlainBasis, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(SimpleBasisManager);

  size_t vectorSize() const { return vectorSize_; }
  
  static Ptr New(size_t vectorSize) {
    return new SimpleBasisManager(vectorSize);
  }

protected:
  explicit SimpleBasisManager(size_t vectorSize) : vectorSize_(vectorSize) {}

  virtual DynamStatePlainBasis::Ptr createNewInstance(const HalfSliceId & key) {
    return DynamStatePlainBasis::New(vectorSize());
  }

private:
  size_t vectorSize_;
};

template <typename SolverType>
class SimpleSolverManager : public Fwk::GenManager<SolverType, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(SimpleSolverManager);

  double tolerance() const { return tolerance_; }

  static Ptr New(double tol) {
    return new SimpleSolverManager(tol);
  }

protected:
  explicit SimpleSolverManager(double tol) : tolerance_(tol) {}

  virtual typename SolverType::Ptr createNewInstance(const HalfSliceId & key) {
    return SolverType::New(tolerance());
  }

private:
  double tolerance_;
};

// External Projection Managers

class ReductorManager : public Fwk::GenManagerInterface<DynamStateReductor *, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(ReductorManager);
  
  virtual DynamStateReductor * instance(const HalfSliceRank & key) const;
  virtual size_t instanceCount() const;

  virtual DynamStateReductor * instanceNew(const HalfSliceRank & key);
  virtual void instanceDel(const HalfSliceRank & key);

  typedef SimpleBasisManager BasisManager;
  typedef SimpleSolverManager<NearSymmetricSolver> SolverManager;

  static Ptr New(BasisManager * basisMgr, SolverManager * solverMgr) {
    return new ReductorManager(basisMgr, solverMgr);
  }

protected:
  ReductorManager(BasisManager * basisMgr, SolverManager * solverMgr) :
    basisMgr_(basisMgr),
    solverMgr_(solverMgr)
  {}
  
  static HalfSliceId keyToKey(HalfSliceRank slice) {
    return HalfSliceId(slice, BACKWARD);
  }


private:
  BasisManager::Ptr basisMgr_;
  SolverManager::Ptr solverMgr_;

  DISALLOW_COPY_AND_ASSIGN(ReductorManager);
};

class ReconstructorManager : public Fwk::GenManagerInterface<DynamStateReconstructor *, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(ReconstructorManager);
  
  virtual DynamStateReconstructor * instance(const HalfSliceRank & key) const;
  virtual size_t instanceCount() const;

  virtual DynamStateReconstructor * instanceNew(const HalfSliceRank & key);
  virtual void instanceDel(const HalfSliceRank & key);

  typedef SimpleBasisManager BasisManager;

  static Ptr New(BasisManager * basisMgr) {
    return new ReconstructorManager(basisMgr);
  }

protected:
  ReconstructorManager(BasisManager * basisMgr) :
    basisMgr_(basisMgr)
  {}

  static HalfSliceId keyToKey(HalfSliceRank slice) {
    return HalfSliceId(slice.previous(), FORWARD);
  }

private:
  BasisManager::Ptr basisMgr_;

  DISALLOW_COPY_AND_ASSIGN(ReconstructorManager);
};

// Tasks

class BasisCondensation : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(BasisCondensation);

  virtual void iterationIs(IterationRank ir); // overriden

  BasisCondensation(const String & name, const DynamStateBasis * origin,
                    DynamStatePlainBasis * target, BasisTransform * operation);

private:
  DynamStateBasis::PtrConst origin_;
  DynamStatePlainBasis::Ptr target_;
  BasisTransform::Ptr operation_;
};

class BasisCondensationManager : public Fwk::GenManager<BasisCondensation, HalfSliceId, BasisCondensation*> {
public:
  EXPORT_PTRINTERFACE_TYPES(BasisCondensationManager);

  BasisCondensationManager(size_t vectorSize, double tolerance, SimpleBasisManager * endBasisMgr, MidBasisManager * midBasisMgr);

protected:
  virtual BasisCondensation * createNewInstance(const HalfSliceId & sliceId); // overriden

private:
  size_t vectorSize_;
  double tolerance_;

  SimpleBasisManager::Ptr endBasisMgr_;
  MidBasisManager::Ptr midBasisMgr_; 
};


class ProjectionBuilding : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(ProjectionBuilding);
 
  virtual void iterationIs(IterationRank ir); // Overriden

  // Solver info
  const NlDynamOps * dynamOps() const { return dynamOps_.ptr(); }
  size_t vectorSize() const { return refDisp_.size(); }
  const Vector & referenceDisplacement() const { return refDisp_; }
  Seconds referenceTime() const { return refTime_; }
  double tolerance() const { return solver()->tolerance(); }

  // Input
  const DynamStateBasis * projectionBasis() const { return projectionBasis_.ptr(); }

  // Output
  const DynamStatePlainBasis * dualBasis() const { return dualBasis_.ptr(); }
  const RankDeficientSolver * solver() const { return solver_.ptr(); }

  void referenceDisplacementIs(const GenVector<double> & disp, Seconds time = Seconds(0.0));
  void toleranceIs(double tol);

  ProjectionBuilding(const String & taskName, NlDynamOps * dynOps,
                     const DynamStateBasis * projectionBasis,
                     DynamStatePlainBasis * dualBasis, RankDeficientSolver * solver);

private:
  NlDynamOps::Ptr dynamOps_;
  Vector refDisp_;
  Seconds refTime_;

  DynamStateBasis::PtrConst projectionBasis_;
  DynamStatePlainBasis::Ptr dualBasis_;
  RankDeficientSolver::Ptr solver_;
};

class ProjectionBuildingFactory : public Fwk::PtrInterface<ProjectionBuildingFactory> {
public:
  EXPORT_PTRINTERFACE_TYPES(ProjectionBuildingFactory);

  ProjectionBuilding * instanceNew(HalfSliceRank key);

  ProjectionBuildingFactory(NlDynamOps * dynamOps,
                            SimpleBasisManager * endBasisMgr,
                            SimpleBasisManager * projBasisMgr,
                            SimpleSolverManager<NearSymmetricSolver> * solverMgr);
private:
  NlDynamOps::Ptr dynamOps_;
  SimpleBasisManager::Ptr endBasisMgr_;
  SimpleBasisManager::Ptr projBasisMgr_;
  SimpleSolverManager<NearSymmetricSolver>::Ptr solverMgr_;
};


class NlProjectionNetwork : public Fwk::PtrInterface<NlProjectionNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(NlProjectionNetwork);

  // Data sharing (global only for the moment, TODO other)
  DynamStatePlainBasis * sharedProjectionBasis() { return sharedProjectionBasis_.ptr(); }
  GlobalStateSharing * sharing() { return sharing_.ptr(); }

  // Concurrent propagation
  ConcurrentBasis::Manager * concurrentMgr() { return concurrentMgr_.ptr(); }
  SimpleBasisManager * endBasisMgr() { return endBasisMgr_.ptr(); } // TODO Remove ?

  // Projection-based correction tasks
  BasisCondensationManager * condensMgr() { return condensMgr_.ptr(); }
  ProjectionBuildingFactory * projBuildMgr() { return projBuildMgr_.ptr(); }
  CorrectionReductor::Manager * corrRedMgr() { return corrRedMgr_.ptr(); }
  CorrectionReconstructor::Manager * corrReconMgr() { return corrReconMgr_.ptr(); }

  NlProjectionNetwork(GlobalStateSharing * sharing,
                      NlDynamOps * dynamOps,
                      size_t vectorSize, double projectionTolerance);

private:
  size_t vectorSize_;
  GlobalStateSharing::Ptr sharing_;
  
  DynamStatePlainBasis::Ptr sharedProjectionBasis_;
  MidBasisManager::Ptr midBasisMgr_;
  
  ConcurrentBasis::Manager::Ptr concurrentMgr_;
  SimpleBasisManager::Ptr endBasisMgr_;
  
  SimpleBasisManager::Ptr projBasisMgr_;
  SimpleSolverManager<NearSymmetricSolver>::Ptr solverMgr_;

  BasisCondensationManager::Ptr condensMgr_;
  ProjectionBuildingFactory::Ptr projBuildMgr_;
  CorrectionReductor::Manager::Ptr corrRedMgr_;
  CorrectionReconstructor::Manager::Ptr corrReconMgr_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLPROJECTIONNETWORK_H */
