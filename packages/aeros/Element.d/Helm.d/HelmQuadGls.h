#ifndef _HELMQUADGLS_H_
#define _HELMQUADGLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuadGls: public HelmElement, public Element {

	int nn[4];
	mutable double coef; // TODO Get rid of this variable.
public:
	explicit HelmQuadGls(int*);

	int getElementType() const override { return 31; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet& cs, double *d, int flg=1) const override;
	FullSquareMatrix  acousticm(CoordSet& cs, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet& cs,double *d, int cmflg) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int *nodes(int *) const override;
	int getTopNumber() const override;

	PrioInfo examine(int sub, MultiFront *) override;

	void addFaces(PolygonSet *pset) override;

	double helmCoef() override { return coef; }
};
#endif

