#ifndef _HEVIBTETRA_H_
#define _HEVIBTETRA_H_

#include <Element.d/Element.h>
#include <cstdio>

class HEVibTetra: public Element {

	int nn[4];
public:
	HEVibTetra(int*);

	int getElementType() const override { return 331; }
	Category getCategory() const override { return Category::Fluid; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&,double *mel, int cmflg=1) const override;
	double getMass(const CoordSet& cs) const override;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() const override;

	bool isHEVFluidElement() override { return true; }

	PrioInfo examine(int sub, MultiFront *) override {
		fprintf(stderr,"HEVibTetra.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C\n"); return *(new PrioInfo);
	};
	//void getThermalForce(CoordSet &, Vector &, Vector &force, int, GeomState *geomState=0) { force.zero(); }
private:

	double volume(double dOmega) const {return dOmega/6.0;}
	double computeMetrics(const CoordSet &, double gN[4][3]) const;
	void buildTetraStiff(double m[4][4], double gN[4][3], double dOmega) const;

};
#endif

