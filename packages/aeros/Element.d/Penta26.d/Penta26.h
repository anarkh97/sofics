#ifndef _PENTA26_H_
#define _PENTA26_H_

#include <Element.d/Element.h>

class Penta26: public Element
{
	int nn[26];
	double *cCoefs;
	double *cFrame;
	NLMaterial *mat;

public:
	explicit Penta26(int*);
	~Penta26() override;

	int getElementType() const override { return 92; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet& cs) const override;

	void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs) override;
	void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs=0) override;

	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
					 int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0) override;

	void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
					  int surface=0, double *ndTemps=0) override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;

	int getTopNumber() const override;
	int numTopNodes() const override;

	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() const override { return 5; }
	int getDecFace(int iFace, int *fn) override;

	int getFace(int iFace, int *fn) override;

	void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame) override
	{ cCoefs = coefs; cFrame = frame; }

	double* setCompositeData2(int, int, double*, double*, CoordSet&, double) override
	{ fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
				"              for Penta26 el.\n"); return (double *) 0;
	}
	void getCFrame(CoordSet &cs, double cFrame[3][3]) const override;

	void setMaterial(NLMaterial *) override;
	int numStates() override;
	void initStates(double *st) override;
	Corotator *getCorotator(CoordSet &cs, double *kel, int, int) override;
};

#endif
