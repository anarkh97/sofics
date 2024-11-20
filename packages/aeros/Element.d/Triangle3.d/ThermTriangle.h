#ifndef _THERMTRIANGLE_H_
#define _THERMTRIANGLE_H_

#include <Element.d/Element.h>

class ThermTriangle: public Element {

	int nn[3];
public:
	explicit ThermTriangle(int*);

	int getElementType() const override { return 53; }
	Category getCategory() const override { return Category::Thermal; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet& cs, double *d, int flg) const override;
	FullSquareMatrix  massMatrix(const CoordSet& cs, double *mel, int cmflg) const override;
	double  getMass(const CoordSet& cs) const override;


	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;
	int getTopNumber() const override;

	void computeTemp(CoordSet&cs, State &state, double gp[2], double*res) override;
	void getFlFlux(double gp[2], double *flF, double *resF) override;
	PrioInfo examine(int sub, MultiFront *) override;
	Corotator * getCorotator(CoordSet &, double*, int, int) override { return 0; }
};
#endif

