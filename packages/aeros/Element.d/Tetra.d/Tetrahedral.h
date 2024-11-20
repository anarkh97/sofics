#ifndef _TETRAHEDRAL_H_
#define _TETRAHEDRAL_H_

#include <Element.d/Element.h>
#include <Element.d/Tetra.d/TetraElementTemplate.hpp>

class Tetrahedral: public Element,
				   public TetraElementTemplate<double>
{
	int nn[4];
	double *cCoefs;
	double *cFrame;
	NLMaterial *mat;
	void computeDjDx(double x[4], double y[4], double z[4], double J, double djdx[12]);

public:
	explicit Tetrahedral(int*);
	~Tetrahedral() override;

	int getElementType() const override { return 23; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const override;
	void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs) override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration) override;
	double getMass(const CoordSet& cs) const override;
	void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega) override;

	void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs) override;
	void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
												   GenFullM<double> &dGfdx, int gravflg, GeomState *geomState) override;
	void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs) override;

	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
					 int surface, double *ndTemps, double ylayer, double zlayer, int avgnum) override;

	void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight,
											   CoordSet &cs, Vector &elDisp, int strInd,
											   int surface, double* ndTemps,
											   int avgnum, double ylayer, double zlayer) override;

	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *,
											CoordSet &cs, Vector &elDisp, int strInd, int surface,
											double *ndTemps, int avgnum, double ylayer, double zlayer) override;

	void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
					  int surface, double *ndTemps) override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;

	int getTopNumber() const override;

	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() const override { return 4; }
	int getDecFace(int iFace, int *fn) override;

	int getFace(int iFace, int *fn) override { return getDecFace(iFace,fn); }

	void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame) override
	{ cCoefs = coefs; cFrame = frame; }

	double* setCompositeData2(int, int, double*, double*, CoordSet&, double) override
	{ fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
				  "              for Tetrahedral el.\n"); return (double *) 0;
	}
	void getCFrame(CoordSet &cs, double cFrame[3][3]) const override;

	void getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
						  int surface, double *ndTemps, double ylayer, double zlayer, int avgnum);

	void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
						   Vector &elDisp, int strInd, int surface, double *ndTemps);

	void setMaterial(NLMaterial *) override;
	int numStates() override;
	void initStates(double *st) override;
	Corotator *getCorotator(CoordSet &cs, double *kel, int, int) override;
};

#endif
