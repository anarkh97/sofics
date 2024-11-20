#ifndef _EIGHTNODEBRICK_H_
#define _EIGHTNODEBRICK_H_

#include <Element.d/Element.h>

#include <complex>
using std::complex;

class EightNodeBrick: virtual public Element
{
protected:
	int nn[8];
	double *cCoefs;
	double *cFrame;
	NLMaterial *mat;

public:
	explicit EightNodeBrick(int*);
	~EightNodeBrick() override;

	int getElementType() const override { return 17; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega) override;

	double getMass(const CoordSet& cs) const override;
	double weight(CoordSet&, double *, int);

	void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs) override;
	void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs) override;

	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
					 int surface, double *ndTemps, double ylayer, double zlayer, int avgnum) override;

	void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
					  int surface, double *ndTemps) override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;

	int getTopNumber() const override;

	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() const override { return 6; }
	int getDecFace(int iFace, int *fn) override;

	int getFace(int iFace, int *fn) override;

	void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame) override
	{ cCoefs = coefs; cFrame = frame; }

	double* setCompositeData2(int, int, double*, double*, CoordSet&, double) override
	{ fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
				  "              for Hexahedral el.\n"); return (double *) 0;
	}
	void getCFrame(CoordSet &cs, double cFrame[3][3]) const override;

	void getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
						  int surface, double *ndTemps, double ylayer, double zlayer, int avgnum);

	void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
						   Vector &elDisp, int strInd, int surface, double *ndTemps);

	void setMaterial(NLMaterial *) override;
	int numStates() override;
	void initStates(double *) override;
	Corotator *getCorotator(CoordSet &cs, double *kel, int, int) override;
};

#endif