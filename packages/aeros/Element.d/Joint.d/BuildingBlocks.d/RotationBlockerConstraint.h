#ifndef _ROTATIONBLOCKERCONSTRAINT_H_
#define _ROTATIONBLOCKERCONSTRAINT_H_

#include <Element.d/MpcElement.d/DotType1ConstraintElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, section 5.2.1.4
// Constraint Building Block type BB4, also known as a rotation blocker constraint
// one constrained rotational DOF

class RotationBlockerConstraint : public DotType1ConstraintElement
{
public:
	RotationBlockerConstraint(int*, int, int);

	int getElementType() const override { return 113; }
	Category getCategory() const override { return Category::Structural; }
	void buildFrame(CoordSet&) override;
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
};

#endif
