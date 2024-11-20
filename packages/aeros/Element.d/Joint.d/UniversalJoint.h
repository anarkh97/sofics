#ifndef _UNIVERSALJOINT_H_
#define _UNIVERSALJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.4
// constrains three translational and one rotational dof

class UniversalJoint : public SuperElement
{
public:
	UniversalJoint(int*);

	int getElementType() const override { return 122; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
