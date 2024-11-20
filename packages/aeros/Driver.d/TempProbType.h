#ifndef _TEMPPROBTYPE_H_
#define _TEMPPROBTYPE_H_

template <class VecType>
class TempState {
   VecType &d_n, &v_n, &v_p;
 public:
   TempState(VecType &d, VecType &v, VecType &vp) :
     d_n(d), v_n(v), v_p(vp) { }
   VecType &getDisp()  { return d_n; }
   VecType &getVeloc() { return v_n; }
   VecType &getPrevVeloc() { return v_p; }
};


template < 
     class DynOps, 
     class VecType, 
     class PostProcessor, 
     class ProblemDescriptor>
class TempSolver {

     ProblemDescriptor *probDesc;
     PostProcessor     *postProcessor;
     double minVel, maxVel, qsbeta;
     //int    stepNum;
     int    steadyFlag,steadyMin,steadyMax;
     double steadyTol;


   public:
     TempSolver(ProblemDescriptor *PrbD) { probDesc = PrbD; }
     void solve();
     void NewmarkTempLoop(TempState<VecType>&, double, double);
     void quasistaticLoop(TempState<VecType>&, double, double);
};

#ifdef _TEMPLATE_FIX_
#include <Driver.d/TempProbType.C>
#endif

#endif
