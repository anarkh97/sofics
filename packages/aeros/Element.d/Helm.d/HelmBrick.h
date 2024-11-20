#ifndef _HELMBRICK_H_
#define _HELMBRICK_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmBrick: public HelmElement, public Element {
private:
	int nn[8];
public:
	HelmBrick(int*);

	int getElementType() const override { return 45; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg = 1) const override;
	FullSquareMatrix acousticm(CoordSet&, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int * = 0) const override;

	void addFaces(PolygonSet *pset) override;

	int getTopNumber() const override;

	PrioInfo examine(int sub, MultiFront *mf) override;
	int nDecFaces() const override { return 6;}
	int getDecFace(int iFace, int *fn) override;
};
#endif

