#ifndef _ROTNSPRLINK_H_
#define _ROTNSPRLINK_H_

#include <Element.d/Element.h>

class RotnSprlink : public Element
{
	int nn[2];

public:

	explicit RotnSprlink(int*);

	int getElementType() const override { return 22; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double* kel, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double* mel, int cmflg) const override;

	void markDofs(DofSetArray&) const override;
	int* dofs(DofSetArray&, int*) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int*) const override;
	Corotator* getCorotator(CoordSet&, double*, int, int) override;

	int getTopNumber() const override;
	bool isSafe() const override { return false; }
	bool isSpring() const override { return true; }
	PrioInfo examine(int sub, MultiFront*) override;

};
#endif
