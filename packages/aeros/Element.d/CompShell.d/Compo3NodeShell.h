#ifndef _COMPO3NODESHELL_H_
#define _COMPO3NODESHELL_H_

#include <Element.d/Element.h>

class GeomState;

class Compo3NodeShell : public Element {
private:
	int      nn[3];
	int      type;
	int numLayers;
	int     (*idlay)[5];
	double  *layData;
	double  *cCoefs;
	double  *cFrame;
	PressureBCond *pbc;
	static bool Wzero_density, Wthermal_force;

public:
	explicit Compo3NodeShell(int*);

	int getElementType() const override { return 20; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;

	void renum(const int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs,double *mel,int cmflg) const override;

	void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
						 GeomState *gs) override;

	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
					 Vector &elDisp, int strInd, int surface,
					 double *ndTemps, double ylayer, double zlayer, int avgnum) override;

	void getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
					  Vector &elDisp, int strInd, int surface,
					  double *ndTemps) override;

	void setCompositeData(int _type, int nlays, double *lData,
						  double *coefs, double *frame) override;
	double * setCompositeData2(int _type, int nlays, double *lData,
							   double *coefs, CoordSet &cs, double theta) override;
	void getCFrame(CoordSet &cs, double cFrame[3][3]) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;
	int getTopNumber() const override;

	int numNodes() const override;
	int* nodes(int *) const override;
	Corotator *getCorotator(CoordSet &, double *,int,int) override;
	double  getMass(const CoordSet& cs) const override;
	double getMassThicknessSensitivity(CoordSet &) override;

	void computeDisp(CoordSet&, State &, const InterpPoint &,
					 double*, GeomState *gs) override;
	void getFlLoad(CoordSet &, const InterpPoint &, double *flF,
				   double *resF, GeomState *gs) override;

	void setPressure(PressureBCond *_pbc) override { pbc = _pbc; }
	PressureBCond* getPressure() override { return pbc; }
	void computePressureForce(CoordSet&, Vector& elPressureForce,
							  GeomState *gs, int cflg, double t) override;
	void getThermalForce(CoordSet&, Vector&, Vector&, int glflag,
						 GeomState *gs) override;

	int getCompositeLayer() override { return numLayers;  }

	double * getMidPoint(CoordSet &) override;
	double * getCompositeData(int nl) override { return layData+(nl*9); }
	double * getCompositeFrame() override { return cFrame;  }

	// Routines for the decomposer
	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() const override { return 1; }
	int getDecFace(int iFace, int *fn) override { for(int i=0; i<3; i++) fn[i] = nn[i]; return 3; }

	int getFace(int iFace, int *fn) override { return getDecFace(iFace,fn); }

	bool hasRot() const override { return true; }

	int getMassType() const override { return 0; } // lumped only

};
#endif

