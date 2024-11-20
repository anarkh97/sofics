#ifdef USE_EIGEN3
#include <Element.d/Joint.d/RevoluteJointSpringCombo.h>
#include <Element.d/Joint.d/RevoluteJoint.h>
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>

RevoluteJointSpringCombo::RevoluteJointSpringCombo(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new RevoluteJoint(nnloc);              // ↓ propIndex
  subElems[1] = new NonlinearTorsionalSpring(nnloc, 2, 1, 0);
}

int 
RevoluteJointSpringCombo::getTopNumber() const
{ 
  return 106; 
}
#endif
