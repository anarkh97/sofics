#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidSpring.h>
#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/TranslationalJoint.h>

RigidSpring::RigidSpring(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };
  subElems[0] = new SphericalJoint(indices);
  subElems[1] = new TranslationalJoint(indices);
}

int
RigidSpring::getTopNumber() const
{
  return 101;
}
#endif
