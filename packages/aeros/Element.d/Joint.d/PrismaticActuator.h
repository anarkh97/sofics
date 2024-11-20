#ifndef _PRISMATICACTUATOR_H_
#define _PRISMATICACTUATOR_H_

#include <Element.d/SuperElement.h>

class PrismaticActuator : public SuperElement
{
public:
	explicit PrismaticActuator(int*);

	int getElementType() const override { return 234; }
	Category getCategory() const override { return Category::Structural; }
	void setProp(StructProp *p, bool myProp) override;
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
