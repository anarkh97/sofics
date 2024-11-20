#ifndef _HELMQUADGAL_H_
#define _HELMQUADGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuadGal: public HelmElement, public Element {

	int nn[4];
public:
	explicit HelmQuadGal(int*);

	int getElementType() const override { return 30; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg = 1) const override;
	FullSquareMatrix acousticm(CoordSet&, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	void getHelmForce(CoordSet&, ComplexVector &, ComplexVector &) override;


	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;

	int getTopNumber() const override;

	PrioInfo examine(int sub, MultiFront *) override;

	void addFaces(PolygonSet *pset) override;

};
#endif

