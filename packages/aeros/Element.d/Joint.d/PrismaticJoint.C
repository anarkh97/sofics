#ifdef USE_EIGEN3
#include <Element.d/Joint.d/PrismaticJoint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/StraightLinePointFollowerConstraint.h>

PrismaticJoint::PrismaticJoint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new ParallelAxesConstraint(nnloc);
  subElems[1] = new RotationBlockerConstraint(nnloc, 2, 1);
  subElems[2] = new StraightLinePointFollowerConstraint(nnloc);
}

int 
PrismaticJoint::getTopNumber() const
{ 
  return 106; 
}
#endif
