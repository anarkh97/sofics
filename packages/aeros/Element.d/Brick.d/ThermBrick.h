#ifndef _THERMBRICK_H_
#define _THERMBRICK_H_

#include <Element.d/Element.h>

class ThermBrick: public Element {

	int nn[8];
public:
	explicit ThermBrick(int*);

	int getElementType() const override { return 51; }
	Category getCategory() const override { return Category::Thermal; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;


	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override;

	Corotator *	getCorotator(CoordSet &cs, double* kel, int, int) override;

};
#endif

