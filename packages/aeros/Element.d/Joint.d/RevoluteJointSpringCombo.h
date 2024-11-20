#ifndef _REVOLUTEJOINTSPRINGCOMBO_H_
#define _REVOLUTEJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class RevoluteJointSpringCombo : public SuperElement
{
public:
	RevoluteJointSpringCombo(int*);

	int getElementType() const override { return 223; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
