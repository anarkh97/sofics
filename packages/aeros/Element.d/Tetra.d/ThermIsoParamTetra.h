#ifndef _THERMISOPARAMTETRA_H_
#define _THERMISOPARAMTETRA_H_

#include <Element.d/Element.h>

class ThermIsoParamTetra: public Element {
public:
	int order;
	int *nn;
	ThermIsoParamTetra(const ThermIsoParamTetra& e);

public:
	ThermIsoParamTetra(int,int*);

	int getElementType() const override { return 50; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	double  getMass(const CoordSet& cs) const override;

	Category getCategory() const override { return Category::Thermal; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return (order*(order+1)*(order+2))/6; }
	int numNodes() const override;
	int* nodes(int * = 0) const override;

	PrioInfo examine(int sub, MultiFront *mf) override;
	int getTopNumber() const override;

	Corotator * getCorotator(CoordSet &, double*, int, int) override { return 0; }
};
#endif
