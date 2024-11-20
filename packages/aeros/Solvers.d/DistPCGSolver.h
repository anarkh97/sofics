#ifndef _DIST_PCG_SOLVER_H
#define _DIST_PCG_SOLVER_H

#include <Solvers.d/PCGSolver.h>

template<class Scalar, class AnyVector, class AnyOperator>
class GenDistPCGSolver : public PCGSolver<Scalar, AnyVector, AnyOperator> 
{
  int numUncon;
  int firstRow;
  int numRows;
  double *nlines;
 public:
   GenDistPCGSolver(AnyOperator *K, int precno, int maxiter,
                    double tolerance, int fRow, int nRow);
   virtual ~GenDistPCGSolver();

   void factor();
   void reSolve(Scalar *rhs);

}; 

#ifdef _TEMPLATE_FIX_
  #include <Solvers.d/DistPCGSolver.C>
#endif

#endif

