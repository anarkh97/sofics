#ifndef _TRANSLATIONALJOINTSPRINGCOMBO_H_
#define _TRANSLATIONALJOINTSPRINGCOMBO_H_

#include <Element.d/SuperElement.h>

class TranslationalJointSpringCombo : public SuperElement
{
public:
	TranslationalJointSpringCombo(int*);

	int getElementType() const override { return 221; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
