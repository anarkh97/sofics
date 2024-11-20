#ifndef _THERM3DQUAD_H_
#define _THERM3DQUAD_H_

#include <Element.d/Element.h>

class Therm3DQuad: public Element {

	int nn[4];
public:
	explicit Therm3DQuad(int*);
	int getElementType() const override { return 3; }
	Category getCategory() const override { return Category::Thermal; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs,double *d, int cmflg) const override;
	double  getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p=0) const override;
	int  numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override;
	void computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp,
	                       int hflInd) override;
	Corotator * getCorotator(CoordSet &, double*, int, int) override { return 0; }
};

#endif
