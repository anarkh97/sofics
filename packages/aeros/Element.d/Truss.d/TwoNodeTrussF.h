#ifndef _TWONODETRUSSF_H_
#define _TWONODETRUSSF_H_

#include <Element.d/Element.h>
class BarFCorotator;

class TwoNodeTrussF : public virtual Element {

	int nn[2];
	double preload;
	BarFCorotator *myCorot;
public:
	explicit TwoNodeTrussF(int*);

	int getElementType() const override { return 111; }
	Category getCategory() const override { return Category::Structural; }
	Element *clone() override;
	void renum(const int *) override;
	void renum(EleRenumMap&) override;
	FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double  getMass(const CoordSet& cs) const override;
	void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg,
						 GeomState *gs) override;
	void getIntrnForce(Vector &elForce, CoordSet& cs,
					   double *elDisp, int forceIndex, double *ndTemps=0) override;
	void getVonMises(Vector& stress, Vector& weight,CoordSet &cs,
					 Vector& elDisp, int strInd,int surface=0,
					 double *ndTemps=0,double ylayer=0.0, double zlayer=0.0, int avgnum=0) override;
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

private:

	void buildBarFrame(CoordSet&, double xg[3][3], double xl[3][3]);

#ifndef SALINAS
	PrioInfo examine(int sub, MultiFront *) override;
#endif
};
#endif
