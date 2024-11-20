#ifndef _COMMONPOINTCONSTRAINT_H_
#define _COMMONPOINTCONSTRAINT_H_

#include <Element.d/SuperElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, sections 5.2.1.1 and 5.2.2.1
// Constraint Building Block type BB1, also known as a common point constraint
// constrains three translational dofs

class CommonPointConstraint : public SuperElement
{
public:
	explicit CommonPointConstraint(int*);

	int getElementType() const override { return 114; }
	Category getCategory() const override { return Category::Structural; }
	void buildFrame(CoordSet& cs) override;
	int getTopNumber() const override;
};

#endif
