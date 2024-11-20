#ifndef _LINESOMMERBC_H_
#define _LINESOMMERBC_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

#ifdef WINDOWS
#include <cmath>
#endif

class LineSommerBC : public SommerElement {

	int nn[2];
public:
	LineSommerBC(int, int, Element *_el = 0, int etype = 1);

	int getElementType() const override { return 1; }
	int numNodes() const override { return 2; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int *dofs(DofSetArray &, int *p) const override;

	int numDofs() const override { return 2; }

	int dim() const override { return 2; }

	LineSommerBC *clone() override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	void sommerVector(CoordSet &, ComplexVector &, ComplexVector &) override;

	void btVector(CoordSet &, ComplexVector &, ComplexVector &) override;

	void ffpDir(int, ComplexD *, CoordSet &, ComplexD *, ComplexD *,
	            double, double(*)[3], double *) override;

	ComplexD ffpCoef(double k) const override {
		return exp(complex<double>(0.0, M_PI / 4.0)) / sqrt(8.0 * M_PI * k);
	}

	void getNormal(const CoordSet &, double [3]) const override;

	double getSize(CoordSet &) override;

	void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) override;
};

#endif
