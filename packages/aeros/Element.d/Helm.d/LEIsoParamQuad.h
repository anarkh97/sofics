#ifndef _LEISOPARAMQUAD_H_
#define _LEISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamQuad: public Element {

	int order;
	int *nn;
	LEIsoParamQuad(const LEIsoParamQuad& e);

public:
	LEIsoParamQuad(int,int*);

	int getElementType() const override { return 100; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	double  getMass(const CoordSet& cs) const override;
	double getMassThicknessSensitivity(CoordSet&) override;

	int getTopNumber() const override;
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
//	int getTopNumber() const override {return 195;}
	int numTopNodes() const override {return order*order;}
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return 2*order*order; }
	int numNodes() const override;
	int* nodes(int * = 0) const override;

//        PrioInfo examine(int sub, MultiFront *mf) override;
	double weight() const override;
	double trueWeight() const override;
};
#endif
