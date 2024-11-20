#ifdef USE_EIGEN3
#include <Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.h>

ParallelAxesConstraint::ParallelAxesConstraint(int* _nn, int axis)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  switch(axis) {
    case 0:
      subElems[0] = new RotationBlockerConstraint(nnloc, 1, 0);
      subElems[1] = new RotationBlockerConstraint(nnloc, 2, 0);
      break;
    case 1: 
      subElems[0] = new RotationBlockerConstraint(nnloc, 0, 1);
      subElems[1] = new RotationBlockerConstraint(nnloc, 2, 1);
      break;
    case 2:
      subElems[0] = new RotationBlockerConstraint(nnloc, 0, 2);
      subElems[1] = new RotationBlockerConstraint(nnloc, 1, 2);
      break;
    default:
      std::cerr << "*** ERROR: invalid argument value in ParallelAxesConstraint, axis = " << axis << std::endl;
      exit(-1);
  }
}

int 
ParallelAxesConstraint::getTopNumber() const
{ 
  return 106; 
}
#endif
