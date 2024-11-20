#ifndef _PININSLOTJOINTSPRINGCOMBO_H_
#define _PININSLOTJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class PinInSlotJointSpringCombo : public SuperElement
{
public:
	PinInSlotJointSpringCombo(int*);

	int getElementType() const override { return 227; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
