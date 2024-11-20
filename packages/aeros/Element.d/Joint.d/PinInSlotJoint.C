#ifdef USE_EIGEN3
#include <Element.d/Joint.d/PinInSlotJoint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/StraightLinePointFollowerConstraint.h>

PinInSlotJoint::PinInSlotJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new ParallelAxesConstraint(nnloc, 1);
  subElems[1] = new StraightLinePointFollowerConstraint(nnloc);
}

int 
PinInSlotJoint::getTopNumber() const
{ 
  return 106; 
}
#endif
