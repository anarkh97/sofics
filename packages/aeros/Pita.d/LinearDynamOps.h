#ifndef PITA_LINEARDYNAMOPS_H
#define PITA_LINEARDYNAMOPS_H

#include "DynamOps.h"
#include "Types.h"
#include <Driver.d/Dynam.h>
#include <map>

class SingleDomainDynamic;;

namespace Pita {
 
class LinearDynamOps : public DynamOps {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearDynamOps);

  typedef double Scalar;
  typedef GenDynamMat<Scalar> DynOpsType;
  typedef ::GenSparseMatrix<Scalar> MatrixType;
  typedef GenSolver<Scalar> SolverType;

  virtual const MatrixType * massMatrix()        const;
  virtual const MatrixType * stiffnessMatrix()   const;
  virtual const MatrixType * dampingMatrix()     const;
  virtual const SolverType * dynamicMassSolver() const;
  
  MatrixType * massMatrix();
  MatrixType * stiffnessMatrix();
  MatrixType * dampingMatrix();
  SolverType * dynamicMassSolver();

  // Handle to internal implementation, owned by the instance
  const DynOpsType * dynamMat() const { return allOps_; }

  class Manager;
  friend class Manager;

protected:
  explicit LinearDynamOps(DynOpsType *);
  virtual ~LinearDynamOps();
  
private:
  DynOpsType * allOps_;

  DISALLOW_COPY_AND_ASSIGN(LinearDynamOps);
};


class GeneralizedAlphaParameter {
public:
  double alpham() const { return alpham_; }
  double alphaf() const { return alphaf_; }
  double beta()   const { return beta_;   }
  double gamma()  const { return gamma_;  }

  // rhoInfinity = dissipation coefficient of forward time-integration when frequency -> infinity
  double rhoInfinity() const { return rhoInfinity_; }
  void rhoInfinityIs(double rhoInf); // 0.0 <= rhoInf <= 1.0
  
  Seconds timeStepSize() const { return timeStepSize_; }
  void timeStepSizeIs(Seconds dt);

  explicit GeneralizedAlphaParameter(Seconds stepSize, double rhoInf = 0.0);
  // Default GeneralizedAlphaParameter(const GeneralizedAlphaParameter &);
  // Default GeneralizedAlphaParameter & operator=(const GeneralizedAlphaParameter &);

  bool operator==(const GeneralizedAlphaParameter &) const;
  bool operator<(const GeneralizedAlphaParameter &) const; // (Arbitrary) total ordering defined for std::map

private:
  double alpham_, alphaf_, beta_, gamma_;
  double rhoInfinity_;
  Seconds timeStepSize_;
};


class LinearDynamOps::Manager : public PtrInterface<LinearDynamOps::Manager> {
public:
  typedef Fwk::Ptr<LinearDynamOps::Manager> Ptr;
  typedef Fwk::Ptr<const LinearDynamOps::Manager> PtrConst;

  typedef double Scalar;
  typedef SingleDomainDynamic ProblemDescriptor;

  LinearDynamOps::Ptr dynOpsNew(const GeneralizedAlphaParameter & param);
  LinearDynamOps::Ptr dynOps(const GeneralizedAlphaParameter & param) const;
  void dynOpsDel(const GeneralizedAlphaParameter & param);
  size_t dynOpsCount() const;

  const ProblemDescriptor * probDesc() const { return probDesc_; }
  ProblemDescriptor * probDesc() { return probDesc_; }

  static Manager::Ptr New(ProblemDescriptor * pbDesc) {
    return new Manager(pbDesc);
  }

protected:
  explicit Manager(ProblemDescriptor * pbDesc);

private:
  ProblemDescriptor * probDesc_;

  typedef std::map<GeneralizedAlphaParameter, LinearDynamOps::Ptr> DynOpsMap;
  DynOpsMap dynOpsMap_;

  // Optimization possible when no damping
  bool noPhysicalDamping_;
};

} // end namespace Pita

#endif /* PITA_LINEARDYNAMOPS_H */
