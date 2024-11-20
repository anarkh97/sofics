#ifndef _LAGLINESOMMER_H_
#define _LAGLINESOMMER_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class LagLineSommer : public SommerElement {

	int *nn;
	int order;
public:
	LagLineSommer(int, int *, Element *_el = 0);

	int getElementType() const override { return 9; }
	int numNodes() const override { return order; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override { return order; }

	int dim() const override { return 2; }

	int *dofs(DofSetArray &, int *p = 0) const override;

	virtual LagLineSommer *clone() override;

	void neumVector(CoordSet &, ComplexVector &,
	                double, double, double, double);

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	ComplexD ffpCoef(double k) const override { return exp(ComplexD(0.0, M_PI / 4.0)) / sqrt(8.0 * M_PI * k); }

	void getNormal(const CoordSet &, double [3]) const override;
};

#endif

