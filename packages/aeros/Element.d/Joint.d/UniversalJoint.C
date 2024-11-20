#ifdef USE_EIGEN3
#include <Element.d/Joint.d/UniversalJoint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/CommonPointConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.h>

UniversalJoint::UniversalJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new CommonPointConstraint(nnloc);
  subElems[1] = new RotationBlockerConstraint(nnloc, 2, 1);
}

int 
UniversalJoint::getTopNumber() const
{ 
  return 106; 
}
#endif
