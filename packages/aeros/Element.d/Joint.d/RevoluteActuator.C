#ifdef USE_EIGEN3
#include <Element.d/Joint.d/RevoluteActuator.h>
#include <Element.d/Joint.d/RevoluteJointSpringCombo.h>
#include <Element.d/Force.d/FollowerMomentElement.h>

RevoluteActuator::RevoluteActuator(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];

  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new RevoluteJointSpringCombo(nnloc);
  subElems[1] = new FollowerMomentElement(&nnloc[0]);
  subElems[2] = new FollowerMomentElement(&nnloc[1]);
}

void
RevoluteActuator::setProp(StructProp *p, bool myProp)
{
  StructProp *p1 = new StructProp(*p);
  p1->amplitude *= -1.0;
  p1->offset *= -1.0;

  SuperElement::setProp(p, myProp);
  subElems[1]->setProp(p1, true);
}

int 
RevoluteActuator::getTopNumber() const
{ 
  return 106; 
}
#endif
