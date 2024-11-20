#ifndef _REVOLUTEJOINT_H_
#define _REVOLUTEJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.5
// constrains three translational and two rotational dofs

class RevoluteJoint : public SuperElement
{
public:
	RevoluteJoint(int*);

	int getElementType() const override { return 123; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
