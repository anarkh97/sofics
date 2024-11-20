#ifndef _CYLINDRICALJOINT_H_
#define _CYLINDRICALJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.6
// constrains two translational and two rotational dofs

class CylindricalJoint : public SuperElement
{
public:
	CylindricalJoint(int*);

	int getElementType() const override { return 124; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
