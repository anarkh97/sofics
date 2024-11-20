#ifndef _SLOSHQUADGAL_H_
#define _SLOSHQUADGAL_H_

#include <Element.d/Element.h>
#include <cstdio>

class SloshQuadGal: public Element {

	int nn[4];
public:
	SloshQuadGal(int*);

	int getElementType() const override { return 301; }
	Category getCategory() const override { return Category::Fluid; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flag=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	double           getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override {
		fprintf(stderr,"SloshQuadGal.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C");
		return *(new PrioInfo);
	};

	void computeSloshDisp(Vector& heatflux, CoordSet &cs, Vector& elTemp,
						  int hflInd) override;
	void computeSloshDispAll(Vector& heatflux, CoordSet &cs, Vector& elTemp) override;
};

#endif
