#ifndef _CONDPROBTYPE_H_
#define _CONDPROBTYPE_H_

template < class ConditionOps, class PostProcessor, class ProblemDescriptor >
class CondSolver {
     ProblemDescriptor *probDesc;
   public:
     CondSolver(ProblemDescriptor *PrbD) 
       { probDesc = PrbD; }
     void solve();
};

#ifdef _TEMPLATE_FIX_
#include <Driver.d/CondProbType.C>
#endif


#endif
