#ifndef _RIGIDSOLID_H_
#define _RIGIDSOLID_H_

#include <Element.d/SuperElement.h>

class RigidSolid : public SuperElement
{
public:
	explicit RigidSolid(int, int*);

	int getElementType() const override { return 71; }
	void buildFrame(CoordSet& cs) override;
	int getTopNumber() const override;
	int numTopNodes() const override;
	bool isRigidElement() const override { return true; }
	bool isSafe() const override;
	Category getCategory() const override { return Element::Structural; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif

