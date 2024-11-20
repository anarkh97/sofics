// ---------------------------------------------------------------------
// HB - 05-24-05
// ---------------------------------------------------------------------
// 32 nodes brick element
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#ifndef _HELMBRICK32_H_
#define _HELMBRICK32_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmBrick32: public HelmElement, public Element {

	int nn[32];
public:
	explicit HelmBrick32(int*);

	int getElementType() const override { return 93; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix  acousticm(CoordSet&, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;

	void addFaces(PolygonSet *pset) override;

	int getTopNumber() const override;
	int numTopNodes() const override;

	PrioInfo examine(int sub, MultiFront *mf) override;
};
#endif

