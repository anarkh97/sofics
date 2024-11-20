#ifndef _GEOM_NL_SOLVER_H_
#define _GEOM_NL_SOLVER_H_

#include <Driver.d/Domain.h>
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;

class GeomNLSolver {
   const SolverInfo  &sinfo;		// Contains Solver Information
   Vector            force;		// external force
   Vector            residual; 		// external force - internal force
   Vector            pVec;
   Vector            arcLenResid;
   Vector	     elementInternalForce;
   Solver           *solver;
   int               numele;		// number of elements
   CoordSet         &nodes;
   Corotator       **allCorot;
   FullSquareMatrix *kelArray;
   DofSetArray      *c_dsa;
   int             **dofsArray;
   Connectivity     *allDOFs;
   int               maxNumDof;
   int               maxit;
   double	     forceNorm;

   double	     lmin;
   double	     lmax;
 public:
   // Constructors
   GeomNLSolver(int nele, Domain &domain, CoordSet &cs, int numdofs, 
             SolverInfo &_sinfo, DofSetArray *_c_dsa, Vector &f, 
  	     Corotator **_allCorot);

   GeomNLSolver(int _numele, Connectivity* _allDOFs,     int _maxNumDof,
           CoordSet &_nodes,            int numdofs, SolverInfo &_sinfo,
FullSquareMatrix *_kelArray,    DofSetArray *_c_dsa,          Vector &f,
      Corotator **_allCorot, Solver *_solver);

   void extendedNewton(GeomState &u,     Vector &dU, double &lambda,
               double deltaLambda, double &deltaS, double w,  int &numExtIter);

   void newton(GeomState &u, const double lambda=1.0);

   int **getDofArray() { return dofsArray;}

   void reBuild(int iteration);

   void getStiffAndForce(GeomState &u, Vector &residual);

   int  checkConvergence(int iteration, double residualNorm, double firstNorm);

};

#endif
