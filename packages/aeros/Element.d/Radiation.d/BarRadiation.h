#ifndef _BARRADIATION_H_
#define _BARRADIATION_H_

#include <Element.d/Element.h>

class BarRadiation: public virtual Element {

	int nn[2];
public:
	explicit BarRadiation(int*);
	~BarRadiation() override;

	int getElementType() const override { return 56; }
	Category getCategory() const override { return Category::Thermal; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	Corotator* getCorotator(CoordSet &, double*, int, int) override;

	int numNodes() const override;
	int * nodes(int *) const override;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() const override;

	bool isRadiationElement() override { return true; }

	void computeTemp(CoordSet&, State &, double[2], double*) override;
	void getFlFlux(double[2], double *, double *) override;
};

#endif
