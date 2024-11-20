#ifndef _PRISMATICDRIVER_H_
#define _PRISMATICDRIVER_H_

#include <Element.d/SuperElement.h>

class PrismaticDriver : public SuperElement
{
public:
	PrismaticDriver(int*);

	int getElementType() const override { return 134; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
