#ifndef _HELMLAGQUADGAL_H_
#define _HELMLAGQUADGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

class HelmLagQuadGal: public HelmElement, public Element {

	int order;
	int *nn;
	void shapeFunctions(double xi, double eta, double *N) const;
	HelmLagQuadGal(const HelmLagQuadGal& e);

public:
	HelmLagQuadGal(int,int*);

	int getElementType() const override { return 43; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg) const override;
	FullSquareMatrix acousticm(CoordSet&, double *d) override;
	void wErrors(CoordSet&,
	             double *l2e, double *h1e, double *l2, double *h1,
	             ComplexD *u, double kappa, double *waveDir) override;
	double           getMass(const CoordSet&) const override;
	void edgeShapeFunctions(int n1, int n2, int *ng,
	                        double **gw, double **N) override;

	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return order*order; }
	int numNodes() const override { return order*order; }
	int * nodes(int *) const override;
	void addFaces(PolygonSet *pset) override;
	int getTopNumber() const override {return 163;}

	PrioInfo examine(int sub, MultiFront *mf) override;
	double weight() const override;
	double trueWeight() const override;

};

class HelmLagQuadGal63 : public HelmLagQuadGal {
public:
	using HelmLagQuadGal::HelmLagQuadGal;

	int getElementType() const override { return 63; }

};
#endif

