#ifndef _SPHERICALJOINTSPRINGCOMBO_H_
#define _SPHERICALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class SphericalJointSpringCombo : public SuperElement
{
public:
	SphericalJointSpringCombo(int*);

	int getElementType() const override { return 220; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
