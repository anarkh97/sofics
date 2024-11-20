#ifndef _BASEPCG_H_
#define _BASEPCG_H_

#include <Threads.d/Paral.h>

template <class AnyVector>
class NullPrec {
   public:
     void apply(AnyVector &r, AnyVector &pr) { pr = r; }
};
 
// NOTE: SP2 does not support this default template arguement currently!

template <class Scalar, class AnyVector, class AnyOperator, 
	  class AnyProjector, class AnyPreconditioner >
class BasePCG {
  protected:
    AnyOperator       *A;
    AnyProjector      *proj;
    AnyPreconditioner *prec;

    AnyVector         &res1;
    AnyVector         &res2;
    AnyVector         &p;
    AnyVector         &z1;
    AnyVector         &ap;
    AnyVector         &z2;
    double tolpcg;
    int maxitr;
    int verbose;

    mutable Scalar finalNorm;
    mutable int numIterations;

  public:
    BasePCG(int _maxitr, double _tolpcg, int _verbose,
          AnyOperator *AA, AnyPreconditioner *_prec = 0, AnyProjector *_proj = 0)
          : res1(*new AnyVector(AA->dim())), 
            res2(*new AnyVector(AA->dim())),
            p(*new AnyVector(AA->dim())),  z1(*new AnyVector(AA->dim())),
            ap(*new AnyVector(AA->dim())), z2(*new AnyVector(AA->dim()))
         { A = AA; proj = _proj; prec = _prec; maxitr = _maxitr; 
           tolpcg = _tolpcg; verbose = _verbose; }
    BasePCG() { };
    ~BasePCG() { if(A) delete A; if(proj) delete proj; if(prec) delete prec;};
    AnyOperator * getOperator() { return A; }
    int neqs() const { return A->neqs(); }
    int doSolve(const AnyVector& rhs, AnyVector& sol);
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/BasePCG.C>
#endif

#endif
