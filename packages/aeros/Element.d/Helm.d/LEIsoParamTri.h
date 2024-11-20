#ifndef _LEISOPARAMTRI_H_
#define _LEISOPARAMTRI_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamTri: public Element {

	int order;
	int *nn;
	LEIsoParamTri(const LEIsoParamTri& e);

public:
	LEIsoParamTri(int,int*);

	int getElementType() const override { return 101; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg) const override;
	double  getMass(const CoordSet& cs) const override;
	double getMassThicknessSensitivity(CoordSet&) override;

	int getTopNumber() const override;
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
//	int getTopNumber() const override {return 195;}
	int numTopNodes() const override {return (order*(order+1))/2;}
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return (order*(order+1)); }
	int numNodes() const override;
	int* nodes(int *) const override;

//        PrioInfo examine(int sub, MultiFront *mf) override;
	double weight() const override;
	double trueWeight() const override;

};
#endif
