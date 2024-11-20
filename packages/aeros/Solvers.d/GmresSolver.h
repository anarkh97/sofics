#ifndef __GMRESSOLVER_H__
#define __GMRESSOLVER_H__

#include <Solvers.d/Solver.h>
#include <Solvers.d/MultiDomainSolver.h>

template <class Scalar> class GmresOrthoSet;
class FSCommunicator;

template<class Scalar, class AnyVector, class AnyOperator, class LeftPreconditioner, class RightPreconditioner>
class GmresSolver : public GenSolver<Scalar>, public MultiDomainSolver<Scalar>
{
  private: 
    int maxit;
    double tol;
    AnyOperator *op;
    void (AnyOperator::*matvec)(AnyVector &, AnyVector &);
    LeftPreconditioner *leftprec;
    void (LeftPreconditioner::*applyLeft)(const AnyVector &, AnyVector &);
    RightPreconditioner *rightprec;
    void (RightPreconditioner::*applyRight)(const AnyVector &, AnyVector &);
    mutable GmresOrthoSet<Scalar> *oSetGMRES;
    int rank;
    /// !< Iteration counter. TODO return the number of interations instead!
    mutable int m_info[1];
  protected:
    using MultiDomainSolver<Scalar>::com;

  public:
    // Constructor
    GmresSolver(int maxit, double tol, AnyOperator *_op, void (AnyOperator::*_matvec)(AnyVector &, AnyVector &),
                LeftPreconditioner *_leftprec, void (LeftPreconditioner::*_applyLeft)(const AnyVector &, AnyVector &),
                RightPreconditioner *_rightprec, void (RightPreconditioner::*_applyRight)(const AnyVector &, AnyVector &),
                FSCommunicator* _com = NULL);

    // Destructor
    ~GmresSolver();

    int maxortho, printNumber, verbose;

    // Linear solution function
    void solve(const AnyVector &b, AnyVector &x) override { x=b; reSolve(x); }
    void reSolve(AnyVector &x) override;
    void reset() { oSetGMRES->reset(); }
    int info(int i) { return m_info[i]; }
    long size() const override { return 0; }
    int neqs() const override { return op->neqs(); }
    void factor() override { }

    Timings& getTimers() override { return GenSolver<Scalar>::getTimers(); }
    double getSolutionTime() const override { return ((GenSolver<Scalar> *)this)->getSolutionTime(); }
    void solve(const Scalar *rhs, Scalar *solution) override {
      std::cerr << "GmresSolver::solve(Scalar *rhs, Scalar *solution) is not implemented\n"; }
};

#ifdef _TEMPLATE_FIX_
  #include <Solvers.d/GmresSolver.C>
#endif

#endif
