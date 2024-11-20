#ifndef _TRIANGLE3_H_
#define _TRIANGLE3_H_

#include <Element.d/Element.h>
#include <Element.d/Triangle3.d/Triangle3ElementTemplate.hpp>

class Triangle3: public Element,
                 public Triangle3ElementTemplate<double>
{

	int nn[3];
public:
	explicit Triangle3(int*);
	int getElementType() const override { return 4; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg = 1) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const override;
	double           getMass(const CoordSet&) const override;
	double getMassThicknessSensitivity(CoordSet&) override;

	void getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                 GeomState *gs) override;
	void getGravityForceThicknessSensitivity(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                         GeomState *gs) override;

	void getVonMises (Vector &stress, Vector &weight,
	                              CoordSet &cs, Vector &elDisp,
	                              int strInd, int surface,
	                              double *ndTemps,
	                              double ylayer, double zlayer, int avgnum) override;

	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
	                                        CoordSet &cs, Vector &elDisp, int strInd, int surface,
	                                        double *ndTemps, int avgnum, double ylayer, double zlayer) override;

	void getAllStress(FullM &stress, Vector &weight,
	                  CoordSet &cs, Vector &elDisp,
	                  int strInd, int surface,
	                  double *ndTemps) override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int  numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;

	int getTopNumber() const override;

	// Routines for the decomposer
	PrioInfo examine(int sub, MultiFront *) override;

};
#endif

