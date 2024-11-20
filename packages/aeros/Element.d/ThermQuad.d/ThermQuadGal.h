#ifndef _THERMQUADGAL_H_
#define _THERMQUADGAL_H_

#include <Element.d/Element.h>

class ThermQuadGal: public Element {

	int nn[4];
public:
	explicit ThermQuadGal(int*);

	int getElementType() const override { return 10; }
	Category getCategory() const override { return Category::Thermal; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flag) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg) const override;
	double           getMass(const CoordSet&) const override;

	void             markDofs(DofSetArray &) const override;
	int*             dofs(DofSetArray &, int *p) const override;
	int              numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override;
	void computeTemp(CoordSet&cs, State &state, double gp[2], double*res) override;
	void getFlFlux(double gp[2], double *flF, double *resF) override;
	void computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp,
	                       int hflInd) override;
	Corotator * getCorotator(CoordSet &, double*, int, int) override { return nullptr; }
};
#endif

