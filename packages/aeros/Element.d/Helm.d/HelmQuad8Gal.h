#ifndef _HELMQUAD8GAL_H_
#define _HELMQUAD8GAL_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmQuad8Gal: public HelmElement, public Element {

	int nn[8];
public:
	explicit HelmQuad8Gal(int*);

	int getElementType() const override { return 32; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix acousticm(CoordSet&, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg) const override;
	double getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() const override { return(132); }
	void addFaces(PolygonSet *pset) override;

};
#endif

