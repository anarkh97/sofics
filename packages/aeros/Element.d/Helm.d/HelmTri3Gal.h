#ifndef _HELMTRI3GAL_H_
#define _HELMTRI3GAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix22[2][2];

class HelmTri3Gal: public HelmElement, public Element {

	int nn[3];
public:
	explicit HelmTri3Gal(int*);

	int getElementType() const override { return 35; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix  acousticm(CoordSet&, double *d) override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;

	void getHelmForce(CoordSet& cs, ComplexVector &vc, ComplexVector &force) override;

	void computedxdxi(CoordSet &cs, int nint, double (*derivatives)[3][2],
					  Matrix22 *dxdxi, double *det);
	void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
						double kappa, double *waveDir) override;

	double           getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;

	void addFaces(PolygonSet *pset) override;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() const override;

};

#endif

