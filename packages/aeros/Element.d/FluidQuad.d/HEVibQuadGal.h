#ifndef _HEVIBQUADGAL_H_
#define _HEVIBQUADGAL_H_

#include <Element.d/Element.h>
#include <cstdio>

class HEVibQuadGal: public Element {

	int nn[4];
public:
	explicit HEVibQuadGal(int*);

	int getElementType() const override { return 321; }
	Category getCategory() const override { return Category::Fluid; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flag) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg) const override;
	double           getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() const override;

	bool isHEVFluidElement() override { return true; }

	PrioInfo examine(int sub, MultiFront *) override {
		fprintf(stderr,"HEVibQuad.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C"); return *(new PrioInfo);
	};
};
#endif

