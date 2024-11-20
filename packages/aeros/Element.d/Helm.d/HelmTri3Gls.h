#ifndef _HELMTRI3GLS_H_
#define _HELMTRI3GLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmTri3Gls: public HelmElement, public Element {

	int nn[3];
	mutable double coef; // TODO Get rid of this variable.
public:
	explicit HelmTri3Gls(int*);

	int getElementType() const override { return 36; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet& cs, double *d, int flg) const override;
	FullSquareMatrix  acousticm(CoordSet& cs, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg) const override;

	double getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;

	void addFaces(PolygonSet *pset) override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override;

	virtual double getHelmCoef() { return coef; }

};
#endif

