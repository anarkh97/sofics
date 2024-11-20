#ifndef _TETRA10HELMGAL_H_
#define _TETRA10HELMGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix33[3][3];

class Tetra10HelmGal: public HelmElement, public Element {

	int nn[10];
	static double tetra10_weights[15];
	static double tetra10_values[15][10];
	static double tetra10_derivatives[15][10][3];
	static double tetra10_vertex_derivatives[10][10][3];

public:
	explicit Tetra10HelmGal(int*);

	int getElementType() const override { return 42; }
	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const override;
	FullSquareMatrix acousticm(CoordSet&, double *kel) override;
	void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*) override;
	void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
	                    double kappa, double *waveDir) override;
	FullSquareMatrix massMatrix(const CoordSet&,double *mel,int cmflg) const override;
	double getMass(const CoordSet& cs) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int             numNodes() const override;
	int * nodes(int *) const override;

	void addFaces(PolygonSet *pset) override;
private:
	void computedxdxi(const CoordSet &cs, int nint, double (*derivatives)[10][3], Matrix33 *dxdxi, double *det) const;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() const override {return (142);}
};
#endif
