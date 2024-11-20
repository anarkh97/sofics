#include <cstdio>
#include <Math.d/Vector.h>
#include <Timers.d/GetTime.h>

/*
  CR is for real symmetric indefinite linear systems
*/

template<class Scalar, class AnyVector, class AnyOperator, class AnyPreconditioner>
void
GenCRSolver<Scalar, AnyVector, AnyOperator, AnyPreconditioner>::solve(const AnyVector &rhs, AnyVector &sol)
{
 solveTime -= getTime();
 AnyVector r(rhs);

 // ... INITIALIZE SOL TO ZERO
 sol.zero();

 double r0r0 = r.squareNorm();
 if(verbose && printNumber > 0)
   fprintf(stderr," ... Iteration #%d\tTwo norm = %1.7e\t  Rel. residual  %1.7e \n",0,r0r0,1.0);
 if(r0r0 == 0.0) {
   solveTime += getTime();
   return;
 }

 // ... FIRST ITERATION
 AnyVector q ( r );
 AnyVector Aq( q );

 int niter = 1;

 if(P) P->apply(r, q);

 A->mult(q, Aq);

 AnyVector p  ( q   );
 AnyVector Ap ( Aq  );
 AnyVector PAp( Ap  );

 if(P) P->apply(PAp, PAp);

 Scalar qAq = q*Aq, qAq_prev;

 while( 1 )
 {
   if(P) P->apply(Ap, PAp); else PAp = Ap;

   Scalar alpha = (qAq)/(PAp*Ap);

   // sol = sol + alpha*p
   sol.linAdd(alpha, p);

   // r = r - alpha*Ap 
   r.linAdd(-alpha, Ap);

   // CHECK CONVERGENCE
   double r2 = r.squareNorm();

   if(verbose && printNumber > 0 && (niter%printNumber == 0))
     fprintf(stderr," ... Iteration #%d\tTwo norm = %1.7e\t  Rel. residual  %1.7e \n",niter,r2,sqrt(r2/r0r0));
   if( r2 <= r0r0*tolerance*tolerance || niter >= maxiter ) {
      if(verbose) fprintf(stderr,"\n ... Total # Iterations = %d",niter);
      if(verbose) fprintf(stderr,"\n ...     Final Two norm = %1.7e\n",r2);
      solveTime += getTime();
      return;
   }

   // q = q - alpha*PAp
   q.linAdd(-alpha, PAp);

   A->mult(q, Aq);

   qAq_prev = qAq;
   qAq = q*Aq;

   Scalar beta = qAq/qAq_prev;

   // p = q + beta*p
   p.linC(q, beta, p);

   // Ap = Aq + beta*Ap
   Ap.linC(Aq, beta, Ap);

   ++niter;
 }
}

