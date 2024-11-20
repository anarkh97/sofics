#ifdef USE_EIGEN3
#include <Element.d/Joint.d/DotConstraintType2a.h>
#include <Element.d/Joint.d/ElementaryFunction.h>

DotConstraintType2a::DotConstraintType2a(int* _nn, int _axis)
 : DotType2ConstraintElement(_nn, _axis) {}

void
DotConstraintType2a::buildFrame(CoordSet& cs)
{
  if(prop) {
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    d0 = f(0);
  }
}

void 
DotConstraintType2a::update(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
  d0 = f(t);

  DotType2ConstraintElement::update(refState, gState, cs, t);
}

double
DotConstraintType2a::getVelocityConstraintRhs(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  double vel_rhs = 0;
  if(prop) {
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    vel_rhs -= f.firstDerivative(t);
  }
  return vel_rhs;
}

double
DotConstraintType2a::getAccelerationConstraintRhs(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  double acc_rhs = MpcElement::getAccelerationConstraintRhs(refState, gState, cs, t);
  if(prop) {
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    acc_rhs -= f.secondDerivative(t);
  }
  return acc_rhs;
}

void
DotConstraintType2a::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                          GeomState *gs, int cflg, double t)
{
  ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
  d0 = f(t);

  DotType2ConstraintElement::computePressureForce(cs, elPressureForce, gs, cflg, t);
}
#endif
