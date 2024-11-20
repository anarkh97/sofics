#ifdef USE_EIGEN3
#include <Element.d/Joint.d/NonlinearTranslationalSpring2.h>

NonlinearTranslationalSpring2::NonlinearTranslationalSpring2(int* _nn, int _axis)
 : DotType3ConstraintElement(_nn, _axis)
{}

void
NonlinearTranslationalSpring2::setProp(StructProp *p, bool _myProp)
{
  StructProp *prop = (_myProp) ? p : new StructProp(*p);
  prop->penalty = prop->k1;
  prop->lagrangeMult = false;
  DotType3ConstraintElement::setProp(prop, true);
}

void 
NonlinearTranslationalSpring2::buildFrame(CoordSet& cs)
{
  DotType3ConstraintElement::buildFrame(cs);
  sp0 = -rhs.r_value;
  rhs.r_value = 0;
}

void 
NonlinearTranslationalSpring2::update(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  DotType3ConstraintElement::update(refState, gState, cs, t);
  rhs.r_value += sp0;
}
#endif
