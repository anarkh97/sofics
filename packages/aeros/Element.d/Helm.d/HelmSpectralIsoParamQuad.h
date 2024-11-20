#ifndef _HELMSPECTRALISOPARAMQUAD_H_
#define _HELMSPECTRALISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmSpectralIsoParamQuad: public HelmElement, public Element {
private:
	int order;
	int *nn;
	HelmSpectralIsoParamQuad(const HelmSpectralIsoParamQuad& e);

public:
	HelmSpectralIsoParamQuad(int,int*);

	int getElementType() const override { return 108; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	double  getMass(const CoordSet& cs) const override;

	Category getCategory() const override { return Category::Acoustic; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int getTopNumber() const override {return 163;}
	int numTopNodes() const override {return order*order;}
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return order*order; }
	int numNodes() const override;
	int* nodes(int * = 0) const override;

	PrioInfo examine(int sub, MultiFront *mf) override;
	double weight() const override;
	double trueWeight() const override;

};
#endif
