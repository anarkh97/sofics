#ifndef _TRANSLATIONALJOINT_H_
#define _TRANSLATIONALJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.3
// constrains three rotational dofs

class TranslationalJoint : public SuperElement
{
public:
	TranslationalJoint(int*);

	int getElementType() const override { return 121; }
	Category getCategory() const override { return Category::Structural; }
	int getTopNumber() const override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;
};

#endif
