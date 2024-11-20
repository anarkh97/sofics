#ifdef USE_EIGEN3
#include <Element.d/Joint.d/SphericalJoint.h>

SphericalJoint::SphericalJoint(int* _nn)
 : CommonPointConstraint(_nn)
{
}

int 
SphericalJoint::getTopNumber() const
{ 
  return 106; 
}
#endif
