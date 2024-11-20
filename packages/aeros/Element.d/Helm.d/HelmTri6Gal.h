#ifndef _HELMTRI6GAL_H_
#define _HELMTRI6GAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix22[2][2];

class HelmTri6Gal: public HelmElement, public Element {

	int nn[6];
public:
	HelmTri6Gal(int*);

	int getElementType() const override { return 38; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix acousticm(CoordSet&, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;

	double           getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int             numNodes() const override;
	int * nodes(int *) const override;

	void addFaces(PolygonSet *pset) override;

	void computedxdxi(CoordSet &cs, int nint, double (*derivatives)[6][2],
					  Matrix22 *dxdxi, double *det);
	void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
						double kappa, double *waveDir) override;
	int getTopNumber() const override {return 138;}
};
#endif

