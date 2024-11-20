#ifndef _TETRAHELMGAL_H_
#define _TETRAHELMGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix33[3][3];

class TetraHelmGal: public HelmElement, public Element {

	int nn[4];
	static double tetra4_weights[4];
	static double tetra4_values[4][4];
	static double tetra4_derivatives[4][4][3];
public:
	TetraHelmGal(int*);

	int getElementType() const override { return 40; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const override;
	FullSquareMatrix acousticm(CoordSet&, double *kel) override;
	void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
	                    double kappa, double *waveDir) override;
	FullSquareMatrix massMatrix(const CoordSet&,double *mel,int cmflg=1) const override;
	double getMass(const CoordSet& cs) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() const override { return 4;}
	int getDecFace(int iFace, int *fn) override;

private:

	void addFaces(PolygonSet *pset) override;

	void computedxdxi(const CoordSet &cs, int nint,
	                  double (*derivatives)[4][3], Matrix33 *dxdxi, double *det) const;
};
#endif

