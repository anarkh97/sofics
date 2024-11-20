#ifndef _TRIANGLERADIATION_H_
#define _TRIANGLERADIATION_H_

#include <Element.d/Element.h>

class TriangleRadiation: public virtual Element {

	int nn[3];
public:
	explicit TriangleRadiation(int*);
	~TriangleRadiation() override;

	int getElementType() const override { return 57; }
	Category getCategory() const override { return Category::Thermal; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs,double *d, int cmflg) const override;
	double getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	Corotator* getCorotator(CoordSet &, double*, int, int) override;

	int numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override;

	bool isRadiationElement() override { return true; }

	void computeTemp(CoordSet&, State &, double[2], double*) override;
	void getFlFlux(double[2], double *, double *) override;
};

#endif
