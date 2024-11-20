#ifndef _NONLINEAR_INFO_
#define _NONLINEAR_INFO_

#include <limits>

// 
// Nonlinear information used in solution algorithms: static, dynamic problems
// with Newton-Raphson iteration or Arclength method
//

struct LinesearchInfo {
   int  type;
   double minstep;
   double maxstep;
   bool verbose;
   int maxit;
   double c1; // constant used to define sufficient decrease condition
   double c2; // constant used to define curvature condition
   double tau; // constant used for backtracking linesearch

   LinesearchInfo() { setDefaults(); }

   void setDefaults() {
     type       = 0; // no linesearch
     maxit      = 10;
     minstep    = std::numeric_limits<double>::epsilon();
     maxstep    = std::numeric_limits<double>::infinity();
     verbose    = false;
     c1         = 1e-4;
     c2         = 0.1;
     tau        = 0.8;
   }
};

struct NonlinearInfo {

   int updateK;         // number of Newton iterations between K updates in time/load step
   int stepUpdateK;     // number of time/load steps between K updates
   int kryflg;          // Krylov correction flag
   int initflg;         // initialization flag
   int reorthoflg;      // full reorthogonalization flag
   int maxiter;         // maximum # of Newton iterations
   int maxVec;          // maximum # of vectors to store
   int fitAlgShell;     // fit algorithm for shell elements
   int fitAlgBeam;      // fit algorithm for beam elements

   double tolRes;       // Newton iteration relative force residual tolerance
   double tolInc;       // Newton iteration relative displacement increment tolerance
   double absTolRes;    // Newton iteration absolute force residual tolerance
   double absTolInc;    // Newton iteration absolute displacement increment tolerance
   double dlambda;      // load step increment
   double maxLambda;    // maximum load step value to attain
   double lfactor;      // scaling factor to determine step size in the load
   int    extMin;       // lower iteration limit scaling in arclength 
   int    extMax;       // maximum iteration limit scaling in arclength

   bool unsymmetric;
   LinesearchInfo linesearch;
   bool failsafe;
   double failsafe_tol;

   int linearelastic;   // 1: blanket small displacement assumption, excepting external follower forces.
                        // 2: blanket small displacement assumption, including external follower forces.
                        //    joints with freeplay and other inequality constraints with penalty enforcment
                        //    also get special treatment in both case 1 and case 2.

   NonlinearInfo() { setDefaults(); }

   void setDefaults() {
                     updateK     = 1; kryflg     =   0; initflg =   0;
                     reorthoflg  = 0; maxiter    = 100; maxVec  =   1;
                     fitAlgShell = 2; fitAlgBeam =   2; dlambda = 1.0;
                     tolRes = 1.0E-6; tolInc     = std::numeric_limits<double>::infinity();
                     absTolRes = 0.0; absTolInc  = std::numeric_limits<double>::infinity();
                     maxLambda = 1.0; lfactor    = 1.0; extMin  =   4;
                     extMax      = 6; unsymmetric = false;
                     failsafe = false; failsafe_tol = std::numeric_limits<double>::epsilon();
                     stepUpdateK = 1; linearelastic = 0; }

};

#endif
