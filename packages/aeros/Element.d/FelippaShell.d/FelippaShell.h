#ifndef _FELIPPASHELL_H_
#define _FELIPPASHELL_H_

#ifdef USE_EIGEN3
#include <Element.d/Element.h>
#include <Corotational.d/Shell3Corotator.h>

template <typename doublereal> class ShellMaterial;

class FelippaShell : public Element,
					 public Shell3Corotator
{
	static int sflg, tflg; // see comments in ShellElementTemplate.hpp
	int      nn[3];
	int      type;
	double  *cFrame;
	PressureBCond *pbc;
	ShellMaterial<double> *gpmat; // gauss points' material
	ShellMaterial<double> *nmat;  // nodes material

public:
	explicit FelippaShell(int*);
	~FelippaShell() override;

	int getElementType() const override { return 15; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;

	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg) const override;

	void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg,
						 GeomState *gs) override;

	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
					 Vector &elDisp, int strInd, int surface = 0,
					 double *ndTemps = 0, double ylayer = 0, double zlayer = 0,
					 int avgnum = 0) override;

	void getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
					  Vector &elDisp, int strInd, int surface = 0,
					  double *ndTemps = 0) override;

	void setProp(StructProp *p, bool myProp) override;
	void setCompositeData(int _type, int nlays, double *lData,
						  double *coefs, double *frame) override;
	double * setCompositeData2(int _type, int nlays, double *lData,
							   double *coefs, CoordSet &cs, double theta) override;
	void getCFrame(CoordSet &cs, double cFrame[3][3]) const override;
	void setMaterial(NLMaterial *) override;
	int numStates() override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;
	int getTopNumber() const override;

	int numNodes() const override;
	int* nodes(int * ) const override;
	double  getMass(const CoordSet& cs) const override;

	void computeDisp(CoordSet&, State &, const InterpPoint &,
					 double*, GeomState *gs = 0) override;
	void getFlLoad(CoordSet &, const InterpPoint &, double *flF,
				   double *resF, GeomState *gs) override;

	void setPressure(PressureBCond *_pbc) override { pbc = _pbc; }
	PressureBCond* getPressure() override { return pbc; }
	void computePressureForce(CoordSet&, Vector& elPressureForce,
							  GeomState *gs = 0, int cflg = 0, double t = 0) override;
	void getThermalForce(CoordSet&, Vector&, Vector &, int, GeomState *) override;

	// Nonlinear
	Corotator* getCorotator(CoordSet&, double*, int, int) override;
	void getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
						  FullSquareMatrix &elK, double *f, double dt, double t) override;
	void getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
						  FullSquareMatrix &elK, double *f, double dt, double t) override;
	void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0, double dt = 0) override;
	bool checkElementDeletion(GeomState &) override;
	void initStates(double *) override;
	double getDissipatedEnergy(GeomState &, CoordSet &) override;
	void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
							 int &nlflag) override;
	void getNLVonMises(Vector& stress, Vector& weight, GeomState &curState,
					   GeomState *refState, CoordSet& c0, int strIndex, int surface = 0,
					   double ylayer = 0, double zlayer = 0, int avgnum = 0, int measure = -1) override;
	void getNLAllStress(FullM &stress, Vector &weight, GeomState &curState,
						GeomState *refState, CoordSet &c0, int strInd, int surface = 0,
						int measure = -1) override;

	// Routines for the decomposer
	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() const override { return 1; }
	int getDecFace(int iFace, int *fn) override { for(int i=0; i<3; i++) fn[i] = nn[i]; return 3; }

	// Miscellaneous
	int getFace(int iFace, int *fn) override { return getDecFace(iFace,fn); }
	bool hasRot() const override { return true; }
	int getMassType() const override { return 0; } // lumped only

	// NEW STRUCTOPT
	double getMassThicknessSensitivity(CoordSet& cs) override;
	void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet&, double *gravityAcceleration) override;
	void getGravityForceThicknessSensitivity(CoordSet&, double *gravity, Vector&, int gravflg,
											 GeomState *gs = 0) override;
	void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
												   GenFullM<double> &dGfdx, int gravflg, GeomState *gs = 0) override;
	void getStiffnessThicknessSensitivity(CoordSet& cs, FullSquareMatrix &dStiffdThick, int flg = 1) override;
	void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs) override;
	void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
											   Vector &elDisp, int strInd, int surface,
											   double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0) override;
	void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp,
										 int strInd, int surface, double *ndTemps = 0, int avgnum = 1,
										 double ylayer = 0, double zlayer = 0) override;
	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *, CoordSet &cs,
											Vector &elDisp, int strInd, int surface,
											double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0) override;
	void getInternalForceThicknessSensitivity(GeomState *refState, GeomState &geomState, CoordSet &cs, Vector &dFintdThick,
											  double dt, double t) override;
private:
	void getVonMisesImpl(Vector &stress, Vector &weight, CoordSet &cs,
						 Vector &elDisp, int strInd, int surface,
						 double *ndTemps, double ylayer, double zlayer,
						 int flag, double *staten = 0, double *statenp = 0);
	void getAllStressImpl(FullM &stress, Vector &weight, CoordSet &cs,
						  Vector &elDisp, int strInd, int surface,
						  double *ndTemps, int flag, double *staten = 0,
						  double *statenp = 0);
};
#endif
#endif

