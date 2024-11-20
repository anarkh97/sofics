#include <Timers.d/GetTime.h>
#include <Utils.d/MyComplex.h>
#include <Utils.d/print_debug.h>
#include <Sfem.d/Sfem.h>
#include <cstdio>
#include <iostream>
extern Sfem* sfem;

template <class Scalar, class AnyVector, class AnyOperator, 
	  class AnyProjector, class AnyPreconditioner >
int
BasePCG<Scalar,AnyVector,AnyOperator,AnyProjector,AnyPreconditioner>
::doSolve(const AnyVector& f, AnyVector& sol)
{
 double t1 = getTime();
 res1.zero();
 res2.zero();
 // ... INITIALIZE SOL TO ZERO
 sol.zero(); 

 if(proj) {
   proj->newSystem();
   proj->project(&f, &sol);
   int hasInit = 1;
   if(hasInit)
     A->mult( sol, res1 ); // MATRIX MULTIPLY: res1 = A*sol
   // ... INITIAL RESIDUAL (compute res1 = f - res1) 
   res1.linC( f, -1.0, res1 );
 } else {
   res1 = f;
 }

 Scalar r0tr0 = res1*res1;
 if(verbose)
   std::cerr << " ... Iteration #  0\t Two norm = " << sqrt(r0tr0) << "\t Rel. residual = 1.0" << std::endl;
 if(r0tr0 == 0.0) {
   return 0;
 }

 // ... FIRST ITERATION
 prec->apply(res1,z1);
 
 // ... APPLY KRYLOV CORRECTION
 if(proj) proj->project(&res1, &z1); 
 p = z1;

 // ... MATRIX-VECTOR MULTIPLY: ap = k*p  )
 A->mult( p, ap ); // MATRIX MULTIPLY: ap = k*p
 Scalar ptap = p*ap;
 if(proj) proj->addDirection( p, ap, ptap);
 Scalar z1tr1 = z1*res1;
 Scalar alpha = z1tr1/ptap;

 // ... UPDATE SOLUTION AND RESIDUALS
 sol.linAdd(alpha,p);
 res2 = res1;
 z2   = z1;
 res1.linAdd(-alpha,ap);

 // OTHER ITERATIONS 
 int niter = 1;
 while( 1 ) {
   // FIRST CHECK CONVERGENCE
   Scalar r1tr1 = res1*res1;

   if(verbose)
     std::cerr << " ... Iteration #  " << niter << "\t Two norm = " << sqrt(r1tr1) << "\t Rel. residual = "  << ScalarTypes::norm(sqrt(r1tr1/r0tr0)) << std::endl;
   if( ScalarTypes::norm(r1tr1) <= ScalarTypes::norm(r0tr0*tolpcg*tolpcg) || niter >= maxitr ) {
      Scalar twonr = sqrt(r1tr1);
      if(verbose) {
        std::cerr << " ...     Final Two norm = " << twonr << std::endl;
        std::cerr << " ...     Final residual = " << r1tr1 << std::endl;
        if(niter >= maxitr && ScalarTypes::norm(r1tr1) >= ScalarTypes::norm(r0tr0*tolpcg*tolpcg))
          std::cerr << " ... Achieved a rel. residual of " << ScalarTypes::norm(sqrt(r1tr1/r0tr0)) 
               << " in " << niter << " iter." << std::endl;
      }
      finalNorm     = r1tr1;
      numIterations = niter;
      return ( ScalarTypes::norm(r1tr1) <= ScalarTypes::norm(r0tr0*tolpcg*tolpcg) ) ? 0 : -1;
   }

   prec->apply(res1,z1);

   // ... APPLY KRYLOV CORRECTION
   if(proj) proj->project(&res1, &z1);
   z1tr1        = z1*res1;
   Scalar beta  = z1tr1/(z2*res2);
   p.linC( z1, beta, p );
   if(proj) proj->ortho(&p, niter);
   A->mult( p, ap ); // SPARSE MATRIX MULTIPLY: ap = k*p
   ptap  = p*ap;
   if(proj) proj->addDirection( p, ap, ptap);
   alpha = z1tr1/ptap;

   // UPDATE SOLUTION AND RESIDUALS THEN INCREMENT COUNTER.
   sol.linAdd(alpha,p);
   res2 = res1;
   z2   = z1;
   res1.linAdd(-alpha,ap);
   ++niter;
 }
}

