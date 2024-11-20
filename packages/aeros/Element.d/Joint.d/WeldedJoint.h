#ifndef _WELDEDJOINT_H_
#define _WELDEDJOINT_H_

#include <Element.d/SuperElement.h>

class WeldedJoint : public SuperElement
{
	EFrame *elemframe;
	bool myframe;
public:
	explicit WeldedJoint(int*);
	~WeldedJoint() override;

	int getElementType() const override { return 119; }
	void setFrame(EFrame *_elemframe) override;
	void buildFrame(CoordSet&) override;
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
