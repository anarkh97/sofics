#ifndef _SHEARPANEL_H_
#define _SHEARPANEL_H_

#include <Element.d/Element.h>
#include <Element.d/Shear.d/ShearPanelTemplate.hpp>

class ShearPanel: public Element,
                  public ShearPanelTemplate<double> {

	int nn[4];
public:
	explicit ShearPanel(int*);

	int getElementType() const override { return 18; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double  getMass(const CoordSet& cs) const override;
	double getMassThicknessSensitivity(CoordSet&) override;

	void getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                     GeomState *geomState) override;
	void getGravityForceThicknessSensitivity(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                         GeomState *geomState) override;

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
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() const override;
	bool hasRot() const override {return true;}

	int getMassType() const override { return 0; } // lumped only
};
#endif

