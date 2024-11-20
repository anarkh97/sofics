#ifndef _PRISMATICJOINT_H_
#define _PRISMATICJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.7
// constrains three rotational and two translational dofs

class PrismaticJoint : public SuperElement
{
public:
	PrismaticJoint(int*);

	int getElementType() const override { return 125; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
