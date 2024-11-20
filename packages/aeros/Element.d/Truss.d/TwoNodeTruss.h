#ifndef _TWONODETRUSS_H_
#define _TWONODETRUSS_H_

#include <Element.d/Element.h>
#include <Element.d/Truss.d/TrussElementTemplate.hpp>

class TwoNodeTruss : public virtual Element,
                     public TrussElementTemplate<double>
{
	int nn[2];
	double preload;
public:
	explicit TwoNodeTruss(int*);
	int getElementType() const override { return 1; }
	Element *clone() override;
	Category getCategory() const override { return Category::Structural; }
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg) const override;
	void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs) override;
	int getMassType() const override { return 2; } // both consistent and lumped
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet&) const override;
	void getMassNodalCoordinateSensitivity(CoordSet &cs, Vector &dMassdx);
	void getLengthNodalCoordinateSensitivity(CoordSet &cs, Vector &dLengthdx);
	void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration) override;
	void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg, GeomState *gs) override;
	void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
	                                               GenFullM<double> &dGfdx, int gravflg, GeomState *geomState) override;
	void getIntrnForce(Vector &elForce, CoordSet& cs,
	                   double *elDisp, int forceIndex, double *ndTemps) override;
	void getVonMises(Vector& stress, Vector& weight,CoordSet &cs,
	                 Vector& elDisp, int strInd,int surface,
	                 double *ndTemps, double ylayer, double zlayer, int avgnum) override;
	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp, CoordSet &cs,
	                                        Vector &elDisp, int strInd, int surface,
	                                        double *ndTemps, int avgnum, double ylayer, double zlayer) override;
	void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
	                                           double *ndTemps, int avgnum, double ylayer, double zlayer) override;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;
	int numNodes() const override;
	int* nodes(int *) const override;
	Corotator* getCorotator(CoordSet &cs, double *kel,int,int) override;
	int getTopNumber() const override;
	void getThermalForce(CoordSet &cs, Vector &ndTemps,
	                     Vector &ThermalForce, int glflag,
	                     GeomState *gs) override;
	void setPreLoad(std::vector<double> &load) override;
	std::vector<double> getPreLoad() override;
	bool isSafe() const override { return false; }
	bool isStart() override {return false; }
#ifndef SALINAS
	PrioInfo examine(int sub, MultiFront *) override;
#endif
private:
	void buildBarFrame(CoordSet&, double xg[3][3], double xl[3][3]);
};
#endif
