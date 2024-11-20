#ifndef _LEISOPARAMHEXA_H_
#define _LEISOPARAMHEXA_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamHexa: public Element {

	int order;
	int *nn;
	LEIsoParamHexa(const LEIsoParamHexa& e);

public:
	LEIsoParamHexa(int,int*);

	int getElementType() const override { return 102; }
	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
	void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega) override;
	double  getMass(const CoordSet& cs) const override;
	void   getGravityForce(CoordSet&,double *gravity,Vector &force,
						   int gravflg, GeomState *gs) override;


	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	void markDofs(DofSetArray &) const override;
	int getTopNumber() const override {return 195;}
	int numTopNodes() const override {return order*order*order;}
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override { return 3*order*order*order; }
	int numNodes() const override;
	int* nodes(int * = 0) const override;

	PrioInfo examine(int sub, MultiFront *mf) override;
	double weight() const override;
	double trueWeight() const override;
	int nDecFaces() const override { return 6;}
	int getDecFace(int iFace, int *fn) override;

};
#endif
