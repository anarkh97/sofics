#ifndef _LINEARTRANSLATIONALSPRING_H_
#define _LINEARTRANSLATIONALSPRING_H_

#include <Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.h>

// this element is a translational spring for small displacements and rotations, implemented using
// the penalized constraint method, in which the penalty parameter is the value of the spring stiffness.
// The constraint which is penalized in the distance l from node A to node B minus
// the constant value of same distance in the undeformed configuration l0.
// i.e f(x) = l - l0

class LinearTranslationalSpring : public ConstantDistanceConstraint
{
  public:
    LinearTranslationalSpring(int*, int=0);

	int getElementType() const override { return 200; }
    void setProp(StructProp *p, bool _myProp) override;

    bool isSpring() const override { return true; }
};

class LinearTranslationalSpring203 : public LinearTranslationalSpring
{
public:
	using LinearTranslationalSpring::LinearTranslationalSpring;

	int getElementType() const override { return 203; }

};

#endif
