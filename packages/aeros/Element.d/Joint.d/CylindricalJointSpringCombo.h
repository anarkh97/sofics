#ifndef _CYLINDRICALJOINTSPRINGCOMBO_H_
#define _CYLINDRICALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class CylindricalJointSpringCombo : public SuperElement
{
public:
	CylindricalJointSpringCombo(int*);

	int getElementType() const override { return 224; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
