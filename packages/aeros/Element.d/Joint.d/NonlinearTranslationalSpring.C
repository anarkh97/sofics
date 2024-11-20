#ifdef USE_EIGEN3
#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>

NonlinearTranslationalSpring::NonlinearTranslationalSpring(int* _nn, int _axis, int _propIndex, int _type, int _ieqtype)
 : DotType2ConstraintElement(_nn, _axis, _type, _ieqtype)
{
  propIndex = _propIndex;
}

void
NonlinearTranslationalSpring::setProp(StructProp *p, bool _myProp)
{
  StructProp *prop = (_myProp) ? p : new StructProp(*p);

  const double k[3] = { p->k1, p->k2, p->k3 };
  const int &i = propIndex;
  if(type == 1 && ieqtype == 1) {
    d0 += p->freeplay[i].ul;
    prop->penalty = (p->freeplay[i].uz-p->freeplay[i].dz)*k[i];
  }
  else if(type == 1 && ieqtype == 2) {
    d0 += p->freeplay[i].ll;
    prop->penalty = (p->freeplay[i].lz-p->freeplay[i].dz)*k[i];
  }
  else if(type == 0) {
    prop->penalty = p->freeplay[i].dz*k[i];
  }
  prop->lagrangeMult = false;

  DotType2ConstraintElement::setProp(prop, true);
}

void 
NonlinearTranslationalSpring::buildFrame(CoordSet& cs)
{
  DotType2ConstraintElement::buildFrame(cs);
  sp0 = -rhs.r_value;
  rhs.r_value = 0;
}

void 
NonlinearTranslationalSpring::update(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  DotType2ConstraintElement::update(refState, gState, cs, t);
  rhs.r_value += sp0;
}
#endif
