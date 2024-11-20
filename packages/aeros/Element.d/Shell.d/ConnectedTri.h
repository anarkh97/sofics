#ifndef _CONNECTED_TRI_H_
#define _CONNECTED_TRI_H_

#include	<Element.d/Element.h>

class ConnectedTri : public Element {

	int nn[4];
public:
	ConnectedTri(int*);

	int getElementType() const override { return 80; }
	Category getCategory() const override { return Category::Structural; }
	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;

	void computeDisp(CoordSet&, State &, const InterpPoint &,
	                     double*, GeomState *gs) override;
	void getFlLoad(CoordSet &, const InterpPoint &,
	                       double *flF, double *resF, GeomState *gs) override;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() const override;

};

#endif
