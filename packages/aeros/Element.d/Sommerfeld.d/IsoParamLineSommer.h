#ifndef _ISOPARAMLINESOMMER_H_
#define _ISOPARAMLINESOMMER_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class IsoParamLineSommer : public SommerElement {

	int *nn;
	int order;
public:
	IsoParamLineSommer(int, int *, Element *_el = 0);

	int getElementType() const override { return 12; }
	int numNodes() const override { return order; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override { return order; }

	int numWetDofs() const override { return 3 * order; }

	int dim() const override { return 2; }

	int *dofs(DofSetArray &, int *p = 0) const override;

	virtual IsoParamLineSommer *clone() override;

	virtual int nFaceCorners() const override { return 2; }

	virtual int *faceCorners() const override {
		int *fc = new int[2];
		fc[0] = nn[0];
		fc[1] = nn[order - 1];
		return fc;
	}

	int *wetDofs(DofSetArray &, int *p = 0) const override;

	void neumVector(CoordSet &, ComplexVector &,
	                double, double, double, double, int pflag = 0) override;

	void neumVectorDeriv(CoordSet &, ComplexVector &,
	                     double, double, double, double, int, int pflag = 0) override;

	void wetInterfaceVector(CoordSet &, ComplexVector &,
	                        double, double, double, double, int, int) override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &cs,
	                                                double *d) override;

	void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) override;

	FullSquareMatrixC sommer2Matrix(CoordSet &, complex<double> *);

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	ComplexD ffpCoef(double k) const override { return exp(ComplexD(0.0, M_PI / 4.0)) / sqrt(8.0 * M_PI * k); }

	void getNormal(const CoordSet &, double [3]) const override;

	void markDofs(DofSetArray &) const override;

};

#endif
