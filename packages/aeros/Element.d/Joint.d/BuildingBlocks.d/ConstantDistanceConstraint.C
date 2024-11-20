#ifdef USE_EIGEN3
#include <Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.h>

ConstantDistanceConstraint::ConstantDistanceConstraint(int* _nn, int _type)
 : DistanceConstraintElement(_nn, 0.0, _type)
{
}

void 
ConstantDistanceConstraint::buildFrame(CoordSet& cs)
{
  Eigen::Vector3d a;
  a << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z;

  f0 = a.norm();

  DistanceConstraintElement::buildFrame(cs);
}

int 
ConstantDistanceConstraint::getTopNumber() const
{ 
  return 106; 
}
#endif
