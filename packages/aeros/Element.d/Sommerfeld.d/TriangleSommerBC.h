#ifndef _TRIANGLESOMMERBC_H_
#define _TRIANGLESOMMERBC_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class TriangleSommerBC : public SommerElement {

	int nn[3];
public:
	TriangleSommerBC(int, int, int, Element *_el = nullptr, int eType = 3);

	int getElementType() const override { return 3; }
	int numNodes() const override { return 3; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override;

	int dim() const override { return 3; }

	int *dofs(DofSetArray &, int *p = 0) const override;

	void markDofs(DofSetArray &) const override;

	TriangleSommerBC *clone() override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	FullSquareMatrix refinedSommerMatrix(CoordSet &, double *) override;

	//FullSquareMatrix surfStiffMatrix(CoordSet&, double *);
	FullSquareMatrix HSommerMatrix(const CoordSet &cs, double *d) const override;
	//FullSquareMatrix HKSommerMatrix(CoordSet&, double *);

	ComplexD ffpCoef(double k) const override { return {0.25 / M_PI, 0.0}; }

	void getNormal(const CoordSet &, double[3]) const override;

	double getSize(CoordSet &) override;

	void BT2(CoordSet &cs, double *e, double *f, double *g,
	         double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d) override;

	void BT2n(CoordSet &cs, double *e, double *f, double *g,
	          double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d, int n) override;

	void sphereBT2(CoordSet &cs, double r, double k, ComplexD *d) override;

	void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) override;

private:
	void get_basis(int, int, int, double (*)[3], double *, double *);

	double getArea(CoordSet &, int *);

	void getLocalCoordinates(CoordSet &, double xx[3], double yy[3], double zz[3]) const;

	void SurfaceRefinement(int nNo, double *x, double *y, double *z, double *xx, double *yy, double *zz) const;

	void GaussCoordinates(int Ngp, double *uPg, double *vPg, double *weight) const;

};

#endif

