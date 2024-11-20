#ifndef _SPECTRALISOPARAMQUADSOMMER_H_
#define _SPECTRALISOPARAMQUADSOMMER_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class SpectralIsoParamQuadSommer : public SommerElement {

	int *nn;
	int order;
public:
	SpectralIsoParamQuadSommer(int, int *, Element *_el = 0, int eType = 10);

	int getElementType() const override { return 14; }
	int numNodes() const override { return order * order; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override { return order * order; }

	int numSolidDofs() const override { return 3 * order * order; }

	int numWetDofs() const override { return 4 * order * order; }

	int dim() const override { return 3; }

	int *dofs(DofSetArray &, int *p = 0) const override;

	virtual SpectralIsoParamQuadSommer *clone() override;

	int *wetDofs(DofSetArray &, int *p = 0) const override;

	int *solidDofs(DofSetArray &, int *p = 0) const override;

	void flipNormal() override;

	void neumVector(CoordSet &, ComplexVector &,
	                double, double, double, double);

	void neumVectorDeriv(CoordSet &cs, ComplexVector &cv, double k,
	                     double dx, double dy, double dz, int n);

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	//GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &cs,
	//                                                double *d);
	virtual void ffp(CoordSet &cs, int numFFP, double *dirFFP,
	                 complex<double> *sol, complex<double> *ffpv);

	FullSquareMatrixC sommer2Matrix(CoordSet &, complex<double> *);

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	ComplexD ffpCoef(double k) const override { return {0.25 / M_PI, 0.0}; }

	void getNormal(const CoordSet &, double [3]) const override;

	void markDofs(DofSetArray &) const override;
};

#endif
