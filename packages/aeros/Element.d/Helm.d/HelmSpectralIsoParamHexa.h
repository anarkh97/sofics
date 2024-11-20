#ifndef _HELMSPECTRALISOPARAMHEXA_H_
#define _HELMSPECTRALISOPARAMHEXA_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmSpectralIsoParamHexa: public HelmElement, public Element {

	int order;
	int *nn;
	HelmSpectralIsoParamHexa(const HelmSpectralIsoParamHexa& e);

public:
	HelmSpectralIsoParamHexa(int,int*);
	~HelmSpectralIsoParamHexa() { delete [] nn; }

	int getElementType() const override { return 105; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	double  getMass(const CoordSet& cs) const override;

	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int getTopNumber() const override {return 195;}
	int numTopNodes() const override {return order*order*order;}
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return order*order*order; }
	int numNodes() const override;
	int* nodes(int * = 0) const override;
	void addFaces(PolygonSet *pset) override;

	PrioInfo examine(int sub, MultiFront *mf) override;
	double weight() const override;
	double trueWeight() const override;

};
#endif
