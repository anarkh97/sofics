#ifdef USE_EIGEN3
#include <Element.d/Joint.d/SphericalJointSpringCombo.h>
#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>

SphericalJointSpringCombo::SphericalJointSpringCombo(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 4;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new SphericalJoint(nnloc);             // ↓ propIndex
  subElems[1] = new NonlinearTorsionalSpring(nnloc, 2, 1, 0);
  subElems[2] = new NonlinearTorsionalSpring(nnloc, 2, 0, 1);
  subElems[3] = new NonlinearTorsionalSpring(nnloc, 1, 0, 2);
}

int 
SphericalJointSpringCombo::getTopNumber() const
{ 
  return 106; 
}
#endif
