#ifndef _LINE_SEARCH_H_
#define _LINE_SEARCH_H_

#include <Utils.d/DistHelper.h>
#include <limits>

template<class ProbDescr, class VecType, class GeomType, template<typename P, typename V, typename G> class Updater>
class StaticResidualEvaluator {

public:
  typedef typename Updater<ProbDescr,VecType,GeomType>::StateIncr StateIncr;
  typedef typename Updater<ProbDescr,VecType,GeomType>::RefState RefState;

private:
  ProbDescr *pbd;
  RefState *refState;
  StateIncr *du;
  VecType &elementInternalForce;
  VecType &gRes;
  double lambda;
  VecType &force;

public:
  StaticResidualEvaluator(ProbDescr *_pbd, RefState *_refState, StateIncr *_du, VecType &_elementInternalForce,
                          VecType &_gRes, double _lambda, VecType &_force)
   : pbd(_pbd), refState(_refState), du(_du), elementInternalForce(_elementInternalForce), gRes(_gRes),
     lambda(_lambda), force(_force) {}

  void operator() (GeomType *geomState, VecType &inc, VecType &residual) {
    Updater<ProbDescr,VecType,GeomType>::updateIncr(du, inc);
    residual.linC(force, lambda);
    Updater<ProbDescr,VecType,GeomType>::integrate(1, pbd, refState, geomState, du, residual, elementInternalForce,
                                                   gRes, lambda, true);
  }
};

template<class ProbDescr, class VecType, class GeomType, template<typename P, typename V, typename G> class Updater>
class DynamicResidualEvaluator {

public:
  typedef typename Updater<ProbDescr,VecType,GeomType>::StateIncr StateIncr;
  typedef typename Updater<ProbDescr,VecType,GeomType>::RefState RefState;

private:
  ProbDescr *pbd;
  RefState *refState;
  StateIncr *du;
  VecType &elementInternalForce;
  VecType &gRes;
  double midtime;
  VecType &ext_force;
  VecType &vel_n;
  VecType &accel;
  StateIncr &inc_displac;
  double delta;
  bool zeroRot;
  VecType &toto;

public:
  DynamicResidualEvaluator(ProbDescr *_pbd, RefState *_refState, StateIncr *_du, VecType &_elementInternalForce,
                           VecType &_gRes, double _midtime, VecType &_ext_force, VecType &_vel_n, VecType &_accel,
                           StateIncr &_inc_displac, double _delta, bool _zeroRot, VecType &_toto)
   : pbd(_pbd), refState(_refState), du(_du), elementInternalForce(_elementInternalForce), gRes(_gRes),
     midtime(_midtime), ext_force(_ext_force), vel_n(_vel_n), accel(_accel), inc_displac(_inc_displac), delta(_delta),
     zeroRot(_zeroRot), toto(_toto) {}

  void operator() (GeomType *geomState, VecType &inc, VecType &residual) {
    Updater<ProbDescr,VecType,GeomType>::updateIncr(du, inc);
    toto = ext_force;
    Updater<ProbDescr,VecType,GeomType>::integrate(1, pbd, refState, geomState, du, toto, elementInternalForce, gRes,
                                                   vel_n, accel, midtime, true);
    Updater<ProbDescr,VecType,GeomType>::get_inc_displacement(pbd, geomState, inc_displac, refState, zeroRot);
    Updater<ProbDescr,VecType,GeomType>::formRHScorrector(pbd, inc_displac, vel_n, accel, toto, residual, geomState,
                                                          delta);
  }
};

template<class ProbDescr, class VecType, class GeomType, class ResidualEvaluatorType>
void linesearchImpl(ProbDescr *pbd, GeomType *geomState, VecType &residual, VecType &p, ResidualEvaluatorType &resEval)
{
  // XXX need to consider case where p is not a descent direction
  union { double s0; double r0; };
  double sL, sU, sJ, sJm1, alphaL, alphaU, alphaJ, alphaJm1;
  int &type = pbd->linesearch().type;
  int &maxit = pbd->linesearch().maxit;
  double &c1 = pbd->linesearch().c1;
  double &c2 = pbd->linesearch().c2;
  double &maxEta = pbd->linesearch().maxstep;
  double &minEta = pbd->linesearch().minstep;
  double &tau = pbd->linesearch().tau;
  bool &verbose = pbd->linesearch().verbose;
  if(type == 1) {
    r0 = residual.sqNorm();
    if(verbose) filePrint(stderr, "Backtracking Line Search, r0 = %-g\n", r0);
  }
  else {
    s0 = p*residual;
    if(verbose) filePrint(stderr, "Bisection Line Search, s0 = %-g\n", s0);
  }
  VecType tmp(pbd->solVecInfo());
  double alpha = 1;
  for(int i = 0; i < maxit; ++i) {
    GeomType *tmpState = new GeomType(*geomState);
    tmp.linC(p, alpha);
    resEval(tmpState, tmp, residual);
    delete tmpState;

    if(type == 1) {
      // type 1 is a backtracking linesearch using merit function 0.5*(L2 norm of residual)^2
      // and first Wolfe condition as the stopping criterion (Amijo rule, sufficient decrease condition)
      double r = residual.sqNorm();
      if(verbose) filePrint(stderr, "i = %3d, r = %#8.3e, r/r0 = %#8.3e, alpha = %#8.3g, (1-2*c1*alpha)*r0 = %#8.3e\n",
                            i, r, r/r0, alpha, (1-2*c1*alpha)*r0);
      if(r > (1-2*c1*alpha)*r0) { // here assuming that the Jacobian matrix was exact.
        alpha = alpha*tau;
      }
      else {
        break;
      }
    }
    else {
      // types 2, 3, 4 and 5 use root-finding algorithms such as bisection to find the step length which 
      // inexactly minimizes the objective function for conservative systems (i.e. when the stiffness matrix
      // is symmetric and the residual is the gradient of a scalar-valued potential), using the 2nd strong
      // Wolf condition as a stopping criterion which forces the step length to lie in at least a broad
      // neighborhood of a local minimizer or stationary point.
      // Note: even for non-conservative systems these methods may still have a meaningful interpretation, but in
      // this case the line search is not a minimization problem, rather it inexactly solves for a step length
      // such that the updated residual is orthogonal to the previous search direction.

      double s = p*residual;

      // set some variables
      if(i == 0) {
        alphaU   = 1.0;
        alphaL   = 0.0;
        sU     = s;
        sL     = s0;
        alphaJ   = 1.0;
        alphaJm1 = 0.0;
        sJ     = s;
        sJm1   = s0;
        if(pbd->linesearch().type == 2 || pbd->linesearch().type == 3) {
          // first search for a bracket to a solution, i.e. we want sU * sL < 0.0
          int count = 0;
          while ((sU * sL > 0.0) && (alphaU < maxEta)) {
          
            count++;
            alphaU *= 2.0;
            tmp.linC(p, alphaU);
      
            // new residual
            tmpState = new GeomType(*geomState);
            resEval(tmpState, tmp, residual);
            delete tmpState;

            sU = p*residual;
          
            if(verbose) filePrint(stderr, "bracketing: count = %3d, alphaU = % #8.3g, alphaL = % #8.3g, sU = % #8.3e, sL = % #8.3e\n",
                                  count, alphaU, alphaL, sU, sL);
          }

          // return if no bracket for a solution found, resetting to initial values
          if (sU * sL > 0.0) {
            filePrint(stderr, " *** WARNING: Bisection Line Search could not bracket solution\n");
            alpha = 1;
            break;
          }
        }
      }
      else {
        if(s*sU < 0) { alphaL = alpha; sL = s; }
        else { alphaU = alpha; sU = s; }
        alphaJm1 = alphaJ; sJm1 = sJ;
        alphaJ = alpha; sJ = s;
      }
      if(verbose) filePrint(stderr, "i = %3d, s = % #8.3e, |s/s0| = %#8.3e, alpha = % #8.3g, sL = % #8.3e, sU = % #8.3e\n",
                            i, s, fabs(s/s0), alpha, sL, sU);
      if(fabs(s) > std::max(c2*fabs(s0), std::numeric_limits<double>::epsilon())) {
        switch(type) {
          case 2 : // Bisection Line Search
            alpha = (alphaL + alphaU)/2;
            break;
          case 3 : // RegulaFalsi Line Search
            alpha = alphaU - sU*(alphaL-alphaU)/(sL-sU);
            break;
          case 4 : // Secant Line Search
            alpha = alphaJ - sJ*(alphaJm1-alphaJ)/(sJm1-sJ);
            break;
          default : 
            filePrint(stderr, " *** ERROR: invalid linesearch type\n");
            exit(-1);
        }
      }
      else {
        break;
      }
    }
  }
  if(verbose) filePrint(stderr,"alpha = %-g\n", alpha);
  if(alpha != 1.0) p *= alpha;
}

#endif
