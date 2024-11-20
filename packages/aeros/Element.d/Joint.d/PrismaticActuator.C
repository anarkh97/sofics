#ifdef USE_EIGEN3
#include <Element.d/Joint.d/PrismaticActuator.h>
#include <Element.d/Joint.d/PrismaticJointSpringCombo.h>
#include <Element.d/Force.d/FollowerForceElement.h>

PrismaticActuator::PrismaticActuator(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];

  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new PrismaticJointSpringCombo(nnloc);
  subElems[1] = new FollowerForceElement(&nnloc[0]);
  subElems[2] = new FollowerForceElement(&nnloc[1]);
}

void
PrismaticActuator::setProp(StructProp *p, bool myProp)
{
  StructProp *p1 = new StructProp(*p);
  p1->amplitude *= -1.0;
  p1->offset *= -1.0;

  SuperElement::setProp(p, myProp);
  subElems[1]->setProp(p1, true);
}

int 
PrismaticActuator::getTopNumber() const
{ 
  return 106; 
}
#endif
