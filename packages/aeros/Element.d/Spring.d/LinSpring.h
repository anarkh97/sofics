#ifndef _LINSPRING_H_
#define _LINSPRING_H_

#include <Element.d/Element.h>

class LinSpring : public Element {

	int nn[1];
public:

	explicit LinSpring(int*);

	int getElementType() const override { return 12; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *kel, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	bool isSafe() const override {return true;}
	bool isSpring() const override {return true;}
	int getTopNumber() const override {return 111;}
	PrioInfo examine(int sub, MultiFront *) override;
	Corotator *getCorotator(CoordSet &, double*, int, int) override;
};
#endif
