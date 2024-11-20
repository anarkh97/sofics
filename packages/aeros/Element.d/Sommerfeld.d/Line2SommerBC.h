#ifndef _LINE2SOMMERBC_H_
#define _LINE2SOMMERBC_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class Line2SommerBC : public SommerElement {

	int nn[3];
public:
	Line2SommerBC(int, int, int, Element *_el = 0, int etype = 2);

	int getElementType() const override { return 2; }
	int numNodes() const override { return 3; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override { return 3; }

	int dim() const override { return 2; }

	int *dofs(DofSetArray &, int *p = 0) const override;

	virtual Line2SommerBC *clone() override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	ComplexD ffpCoef(double k) const override { return exp(ComplexD(0.0, M_PI / 4.0)) / sqrt(8.0 * M_PI * k); }

	void getNormal(const CoordSet &, double [3]) const override;

	void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) override;
};

#endif

