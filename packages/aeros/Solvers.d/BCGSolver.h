#ifndef _BCGSOLVER_H_
#define _BCGSOLVER_H_

#include <Solvers.d/Solver.h>

template <class Scalar, class AnyVector, class AnyOperator, class AnyPreconditioner>
class GenBCGSolver : public GenSolver<Scalar> {
 protected:
   AnyOperator *A;
   AnyPreconditioner *P;
   double tolerance;
   int maxiter;
   mutable double solveTime;
 public:
   int printNumber, verbose;

   GenBCGSolver(int _maxit, double _tol, AnyOperator *_A, AnyPreconditioner *__P = 0)
     { maxiter = _maxit; tolerance = _tol; A = _A; P = __P; solveTime = 0;
       printNumber = 1; verbose = 1; }
   ~GenBCGSolver() {};
   int neqs() const override { return A->neqs(); }
   void solve(const AnyVector &, AnyVector &) override;
   void reSolve(AnyVector &rhs) override { AnyVector rhs_copy(rhs); solve(rhs_copy, rhs); }
   long size() const override { return 0; }
   void factor() override {}
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/BCGSolver.C>
#endif

#endif
