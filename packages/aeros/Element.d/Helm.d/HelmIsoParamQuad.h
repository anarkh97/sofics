#ifndef _HELMISOPARAMQUAD_H_
#define _HELMISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamQuad: public HelmElement, public Element {

	int order;
	int *nn;
	HelmIsoParamQuad(const HelmIsoParamQuad& e);

public:
	HelmIsoParamQuad(int,int*);

	int getElementType() const override { return 98; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	FullSquareMatrixC stiffness(const CoordSet&, complex<double> *d) const override;
	FullSquareMatrixC massMatrix(const CoordSet&, complex<double> *d) const override;
	double  getMass(const CoordSet& cs) const override;

	Category getCategory() const override { return Category::Acoustic; }
	int getTopNumber() const override;
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
//	int getTopNumber() const override {return 195;}
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
