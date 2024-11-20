#ifndef _HELMPENTA_H_
#define _HELMPENTA_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmPenta: public HelmElement, public Element {

	int nn[6];
public:
	HelmPenta(int*);

	int getElementType() const override { return 90; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg = 1) const override;
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

	PrioInfo examine(int sub, MultiFront *mf) override;
};
#endif

