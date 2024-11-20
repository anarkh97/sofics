#ifdef USE_EIGEN3
#include <Element.d/Joint.d/RevoluteDriver.h>
#include <Element.d/Joint.d/BuildingBlocks.d/CommonPointConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.h>
#include <Element.d/Joint.d/DotConstraintType1a.h>

RevoluteDriver::RevoluteDriver(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new CommonPointConstraint(nnloc);
  subElems[1] = new ParallelAxesConstraint(nnloc);
  subElems[2] = new DotConstraintType1a(nnloc, 2, 1);
}

int 
RevoluteDriver::getTopNumber() const
{ 
  return 106; 
}
#endif
