#ifndef _NL_STATICPROBTYPE_H_
#define _NL_STATICPROBTYPE_H_

#include <Driver.d/StateUpdater.h>

template < class OpSolver, 
           class VecType, 
           class PostProcessor, 
           class ProblemDescriptor, 
           class GeomType,
           class StateUpdate = IncrUpdater<ProblemDescriptor, VecType, GeomType> >
class NLStaticSolver {
     ProblemDescriptor *probDesc;
     GeomType *geomState;
     typename StateUpdate::RefState *refState;
     typename StateUpdate::StateIncr *stateIncr;

   public:
     // Constructor
     NLStaticSolver(ProblemDescriptor *PrbD) 
       { probDesc = PrbD; geomState = NULL; }
     // Destructor
     ~NLStaticSolver();

     void solve();
     void arclength();
     static int newton(VecType& force, VecType& residual, VecType &glResid,
                       VecType& elementInternalForce, ProblemDescriptor *probDesc,
                       typename StateUpdate::RefState *refState,
                       GeomType* geomState, typename StateUpdate::StateIncr *stateIncr,
                       int& numIter, double& residualNorm, double lambda=1.0, int step=1);

     void extendedNewton(GeomType &u, GeomType &u0, VecType &dU, double &lambda, 
                         double deltaLambda,
                         double &deltaS, double w, int &numExtIter,
                         VecType& force, VecType& residual, VecType& glRes,
                         VecType& arcLenResid,
                         double forceNorm, VecType& elementInternalForce, VecType& pVec, int step=1);

     void predictorStep(GeomType &u, GeomType &u0, VecType &dU, double &lambda, double &deltaLambda,
                        double &deltaS, double &deltaS0, double w,
                        VecType& force, VecType& residual, VecType &totRes, VecType& elementInternalForce,
                        VecType& duds, int step=1);

     GeomType *getGeomState() { return geomState; }
};

#ifdef _TEMPLATE_FIX_
#include <Driver.d/NLStaticProbType.C>
#endif

#endif
