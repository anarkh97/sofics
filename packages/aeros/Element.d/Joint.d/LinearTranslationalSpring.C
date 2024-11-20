#ifdef USE_EIGEN3
#include <Element.d/Joint.d/LinearTranslationalSpring.h>

LinearTranslationalSpring::LinearTranslationalSpring(int* nn, int type)
 : ConstantDistanceConstraint(nn, type)
{
}

void
LinearTranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  StructProp *prop = (_myProp) ? p : new StructProp(*p);
  prop->penalty = prop->k1;
  if(type == 1) f0 += p->freeplay[0].ul;
  prop->lagrangeMult = false;
  ConstantDistanceConstraint::setProp(prop, true);
}
#endif
