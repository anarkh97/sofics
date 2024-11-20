#ifdef USE_EIGEN3
#include <Element.d/Joint.d/BuildingBlocks.d/StraightLinePointFollowerConstraint.h>
#include <Element.d/MpcElement.d/DotType2ConstraintElement.h>
#include <iostream>

StraightLinePointFollowerConstraint::StraightLinePointFollowerConstraint(int* _nn, int axis)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  switch(axis) {
    case 0 :
      subElems[0] = new DotType2ConstraintElement(nnloc, 1);
      subElems[1] = new DotType2ConstraintElement(nnloc, 2);
      break;
    case 1 :
      subElems[0] = new DotType2ConstraintElement(nnloc, 0);
      subElems[1] = new DotType2ConstraintElement(nnloc, 2);
      break;
    case 2 :
      subElems[0] = new DotType2ConstraintElement(nnloc, 0);
      subElems[1] = new DotType2ConstraintElement(nnloc, 1);
      break;
    default :
      std::cerr << "*** ERROR: invalid argument in StraightLinePointFollowerConstraint, axis = " << axis << std::endl;
      exit(-1);
  }
}

int 
StraightLinePointFollowerConstraint::getTopNumber() const
{ 
  return 106; 
}
#endif
