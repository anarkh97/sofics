#ifndef _BARSLOSHFS_H_
#define _BARSLOSHFS_H_

#include <Element.d/Element.h>
#include <cstdio>

class BarSloshFS: public Element {
private:
	int nn[2];
public:

	explicit BarSloshFS(int*);

	int getElementType() const override { return 302; }
	Category getCategory() const override { return Category::Fluid; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	PrioInfo examine(int sub, MultiFront *) override {
		fprintf(stderr,"BarSloshFS.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C.");
		return *(new PrioInfo);
	};
	int getTopNumber() const override;
};

#endif
