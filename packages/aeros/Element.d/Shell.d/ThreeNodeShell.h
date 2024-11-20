#ifndef _THREENODESHELL_H_
#define _THREENODESHELL_H_

#include <Element.d/Element.h>
#include <Driver.d/MultiFront.h>

class GeomState;
class Shell3Corotator;

class ThreeNodeShell : public Element
{
protected:
	int nn[3];
	double w;
	Shell3Corotator *corot;
	PressureBCond *pbc;
public:
	explicit ThreeNodeShell(int*, double _w=3);

	int getElementType() const override { return 8; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;
	void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg,
	                     GeomState *gs) override;
	void getVonMises(Vector& stress, Vector& weight, CoordSet& cs,
	                 Vector& elDisp, int strInd, int surface,
	                 double *ndTemps, double ylayer, double zlayer, int avgnum) override;
	void getAllStress(FullM& stress, Vector& weight, CoordSet& cs,
	                  Vector& elDisp, int strInd, int surface,
	                  double *ndTemps) override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override ;
	int* nodes(int *) const override;
	Corotator *getCorotator(CoordSet &, double *, int , int) override;

	void computeDisp(CoordSet&, State &, const InterpPoint &,
	                 double*, GeomState *gs) override;
	void printInfo(CoordSet&, State &, double[2]);
	void getFlLoad(CoordSet &, const InterpPoint &,
	               double *flF, double *resF, GeomState *gs) override;

	int getTopNumber() const override;
	int nDecFaces() const override { return 1;}
	int getDecFace(int iFace, int *fn) override { for(int i; i<3; i++) fn[i] = nn[i]; return 3; }

	int getFace(int iFace, int *fn) override { return getDecFace(iFace,fn); }

	void setPressure(PressureBCond *_pbc) override { pbc = _pbc; }
	PressureBCond* getPressure() override { return pbc; }
	void computePressureForce(CoordSet&, Vector& elPressureForce,
	                          GeomState *gs, int cflg, double t) override;

	void getThermalForce(CoordSet& cs, Vector& ndTemps,Vector &elThermalForce,
	                     int glfag, GeomState *gs) override;

	bool isShell() const override { return true; }

	int getMassType() const override; // lumped only

	// DEC
	bool hasRot() const override {return true;}
	PrioInfo examine(int sub, MultiFront *mf) override;

#ifdef USE_EIGEN3
	// NEW STRUCTOPT 
	double getMassThicknessSensitivity(CoordSet& cs) override;
	void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet&, double *gravityAcceleration) override;
	void getGravityForceThicknessSensitivity(CoordSet&, double *gravity, Vector&, int gravflg,
	                                         GeomState *gs) override;
	void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
	                                               GenFullM<double> &dGfdx, int gravflg, GeomState *gs) override;
	void getStiffnessThicknessSensitivity(CoordSet& cs, FullSquareMatrix &dStiffdThick, int flg) override;
	void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs) override;
	void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
	                                           Vector &elDisp, int strInd, int surface,
	                                           double *ndTemps, int avgnum, double ylayer, double zlayer) override;
	void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp,
	                                     int strInd, int surface, double *ndTemps, int avgnum,
	                                     double ylayer, double zlayer) override;
	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *, CoordSet &cs,
	                                        Vector &elDisp, int strInd, int surface,
	                                        double *ndTemps, int avgnum, double ylayer, double zlayer) override;
#endif
};
#endif

