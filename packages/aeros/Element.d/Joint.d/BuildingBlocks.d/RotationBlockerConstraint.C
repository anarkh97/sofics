#ifdef USE_EIGEN3
#include <Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.h>

RotationBlockerConstraint::RotationBlockerConstraint(int* _nn, int _axis1, int _axis2)
 : DotType1ConstraintElement(_nn, _axis1, _axis2)
{
}

void 
RotationBlockerConstraint::buildFrame(CoordSet& cs)
{
  if(DotType1ConstraintElement::C0) {
    Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType1ConstraintElement::C0[0][0]);
    d0 = C0.row(axis1).dot(C0.row(axis2));
  }

  DotType1ConstraintElement::buildFrame(cs);
}

int 
RotationBlockerConstraint::getTopNumber() const
{ 
  return 106; 
}
#endif
