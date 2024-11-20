#ifdef USE_EIGEN3
#include <Element.d/Joint.d/DotConstraintType1a.h>
#include <Element.d/Joint.d/ElementaryFunction.h>

#if (__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)
#include <cmath>
using std::remainder;
#else
double remainder(double x, double y) { return ((x > 0 && y > 0) || (x < 0 && y < 0)) ? x - int((x+0.5*y)/y)*y : x - int((x-0.5*y)/y)*y; }
#endif

DotConstraintType1a::DotConstraintType1a(int* _nn, int _axis1, int _axis2)
 : DotType1ConstraintElement(_nn, _axis1, _axis2, 0)
{
  axis1_copy = _axis1;
  t_reparam = -1;
  offset = M_PI/2;
}

void
DotConstraintType1a::buildFrame(CoordSet& cs)
{
  double theta;
  if(prop) {
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    theta = f(0);
  }
  else {
    theta = 0;
  }
  d0 = std::cos(theta-offset);
  DotType1ConstraintElement::buildFrame(cs);
}

void 
DotConstraintType1a::update(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
  double theta = f(t);
  double TOL = M_PI/4;
  if(std::fabs(remainder(theta-offset,M_PI)) < TOL && t > t_reparam) { // reparameterize to stay away from singularity
    if(axis1 != axis2) {
      axis1 = axis2;
      offset = 0;
    }
    else {
      axis1 = axis1_copy;
      offset = M_PI/2;
    }
    t_reparam = t;
  }
  d0 = std::cos(theta-offset);

  DotType1ConstraintElement::update(refState, gState, cs, t);
}

double
DotConstraintType1a::getVelocityConstraintRhs(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  double vel_rhs = 0;
  // g(t)   = -cos(f(t)-c) + d
  // g'(t)  = f'(t)*(-sin(c-f(t)))
  if(prop) {
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    double theta = f(t);
    double dfdt = f.firstDerivative(t);
    vel_rhs -= dfdt*(-std::sin(offset-theta));
  }
  return vel_rhs;
}

double
DotConstraintType1a::getAccelerationConstraintRhs(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  double acc_rhs = MpcElement::getAccelerationConstraintRhs(refState, gState, cs, t);
  // g(t)   = -cos(f(t)-c) + d
  // g''(t) = f'(t)^2*cos(c-f(t)) - f''(t)*sin(c-f(t))
  if(prop) {
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    double theta  = f(t);
    double dfdt   = f.firstDerivative(t);
    double d2fdt2 = f.secondDerivative(t);
    acc_rhs -= std::pow(dfdt,2)*std::cos(offset-theta) - d2fdt2*std::sin(offset-theta);
  }
  return acc_rhs;
}

void
DotConstraintType1a::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                          GeomState *gs, int cflg, double t)
{
  ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
  d0 = f(t);

  DotType1ConstraintElement::computePressureForce(cs, elPressureForce, gs, cflg, t);
}
#endif
