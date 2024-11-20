#ifndef _TRIANGLE6SOMMERBC_H_
#define _TRIANGLE6SOMMERBC_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

typedef double Matrix22[2][2];

class Triangle6SommerBC : public SommerElement {

	int nn[6];
public:
	Triangle6SommerBC(int, int, int, int, int, int, Element *_el = 0, int etype = 6);

	int getElementType() const override { return 6; }
	int numNodes() const override { return 6; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override;

	int dim() const override { return 3; }

	int *dofs(DofSetArray &, int *p) const override;

	void flipNormal() override;

	Triangle6SommerBC *clone() override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	void sommerMatrixEllipsoid(CoordSet &cs, double kappa, double H[3], double K[3], ComplexD *d) override;

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	ComplexD ffpCoef(double k) const override { return {0.25 / M_PI, 0.0}; }

	void getNormal(const CoordSet &, double[3]) const override;

	void BT2(CoordSet &cs, double *e, double *f, double *g,
	         double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d) override;

	void sphereBT2(CoordSet &cs, double r, double k, ComplexD *d) override;

	void ellipsoidBT2(CoordSet &cs, double a, double b, double k, ComplexD *d) override;

	void markDofs(DofSetArray &) const override;
//        FullSquareMatrix  stiffness(const CoordSet&, double *d, int flg = 1) const override;
//        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
//        int* nodes(int * = 0) const override;

//        bool isSommerElement() { return true; }
private:
	static double tri6_derivatives[12][6][2];
	static double tri6_values[12][6];
	static double tri6_weights[12];
	static double tri6_coord[12][2];
	static double tri3_values[12][3];

	void computedxdxi(double x[6], double y[6], Matrix22 *dxdxi, double *det);

	void get_basis(int, int, int, double (*)[3], double *, double *);

	double getArea(CoordSet &, int *);

	void getLocalCoordinates(CoordSet &, double xx[6], double yy[6], double zz[6]) const;

	void getLocalCoordinatesNew(CoordSet &, double T[3][3], double xi1[3], double xi2[3], double xi3[3]);


};

#endif
