#ifndef _ISOPARAMTRISOMMER_H_
#define _ISOPARAMTRISOMMER_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class IsoParamTriSommer : public SommerElement {

	int *nn;
	int order;
public:
	IsoParamTriSommer(int, int *, Element *_el = nullptr, int etype = 11);

	int getElementType() const override { return 11; }
	int numNodes() const override { return (order * (order + 1)) / 2; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override { return (order * (order + 1)) / 2; }

	int numSolidDofs() const override { return 3 * (order * (order + 1)) / 2; }

	int numWetDofs() const override { return 2 * (order * (order + 1)); }

	int dim() const override { return 3; }

	int *dofs(DofSetArray &, int *p) const override;

	int *wetDofs(DofSetArray &, int *p) const override;

	int *solidDofs(DofSetArray &, int *p) const override;

	IsoParamTriSommer *clone() override;

	int nFaceCorners() const override { return 3; }

	int *faceCorners() const override {
		auto *fc = new int[3];
		fc[0] = nn[0];
		fc[1] = nn[order - 1];
		fc[2] = nn[(order * (order + 1)) / 2 - 1];
		return fc;
	}

	void flipNormal() override;

	void neumVector(CoordSet &, ComplexVector &,
	                double, double, double, double, int pflag) override;

	void neumVectorDeriv(CoordSet &cs, ComplexVector &cv, double k,
	                     double dx, double dy, double dz, int n,
	                     int pflag) override;

	void wetInterfaceVector(CoordSet &, ComplexVector &,
	                        double, double, double, double, int, int) override;

	void wetInterfaceVector(CoordSet &, ComplexVector &,
	                        complex<double> (*)[3], complex<double> *) override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &cs, double *d) override;

	void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) override;

	void ffp(CoordSet &cs, int numFFP, double *dirFFP,
	                 complex<double> *sol, complex<double> *ffpv, bool direction) override;

	FullSquareMatrixC sommer2Matrix(CoordSet &, complex<double> *);

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	ComplexD ffpCoef(double k) const override { return {0.25 / M_PI, 0.0}; }

	void getNormal(const CoordSet &, double [3]) const override;

	void markDofs(DofSetArray &) const override;
};

#endif
