#ifndef _HELMBRICKGLS_H_
#define _HELMBRICKGLS_H_

#include <cmath>
#include <Element.d/Helm.d/HelmElement.h>

class HelmBrickGLS: public HelmElement, public Element {

	int nn[8];
	mutable double coef; // TODO Get rid of this variable!
public:
	HelmBrickGLS(int*);

	int getElementType() const override { return 44; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix  acousticm(CoordSet&, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const override;
	double getMass(const CoordSet& cs) const override;


	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;

	void addFaces(PolygonSet *pset) override;

	int getTopNumber() const override;

	double helmCoef() override { return coef; }

	PrioInfo examine(int sub, MultiFront *mf) override;
	int nDecFaces() const override { return 6;}
	int getDecFace(int iFace, int *fn) override;
};
#endif

