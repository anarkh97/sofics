#ifndef _REVOLUTEACTUATOR_H_
#define _REVOLUTEACTUATOR_H_

#include <Element.d/SuperElement.h>

class RevoluteActuator : public SuperElement
{
public:
	explicit RevoluteActuator(int*);

	int getElementType() const override { return 226; }
	Category getCategory() const override { return Category::Structural; }
	void setProp(StructProp *p, bool myProp) override;
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
