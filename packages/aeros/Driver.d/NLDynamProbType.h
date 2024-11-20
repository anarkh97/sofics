#ifndef _NL_DYNAMPROBTYPE_H_
#define _NL_DYNAMPROBTYPE_H_

#include <cstdio>
#include <cstdlib>
#include <Driver.d/StateUpdater.h>

template < class OpSolver, 
           class VecType, 
           class PostProcessor, 
           class ProblemDescriptor, 
           class GeomType,
           class StateUpdate = IncrUpdater<ProblemDescriptor,VecType, GeomType> >
class NLDynamSolver {
     ProblemDescriptor *probDesc;
     typename StateUpdate::RefState *refState;
     typename StateUpdate::StateIncr *stateIncr;

     double beta, gamma, alphaf, alpham;
   public:

     // Constructor
     NLDynamSolver(ProblemDescriptor *PrbD) :
       probDesc(PrbD), refState(0), stateIncr(0)
     {} 

     ~NLDynamSolver() {
       delete refState;
       delete stateIncr;
     }

     void solve();
};

#ifdef _TEMPLATE_FIX_
#include <Driver.d/NLDynamProbType.C>
#endif

#endif
