#ifndef _MATNL_COROTATOR_H_
#define _MATNL_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class MatNLElement;
class NLMaterial;

class MatNLCorotator : public Corotator {

	MatNLElement *ele;
	bool own;

public:
	MatNLCorotator(MatNLElement *, bool own = true);
	~MatNLCorotator() override;

	void getStiffAndForce(GeomState &curState, CoordSet &cs,
						  FullSquareMatrix &elk, double *f, double dt, double t) override
	{ getStiffAndForce((GeomState *) nullptr, curState, cs, elk, f, dt, t); }

	void getInternalForce(GeomState &curState, CoordSet &cs,
						  FullSquareMatrix &elk, double *f, double dt, double t) override 
	{ getInternalForce((GeomState *) nullptr, curState, cs, elk, f, dt, t); }

	void getStiffAndForce(GeomState *refState, GeomState &curState, CoordSet &cs,
						  FullSquareMatrix &elk, double *f, double dt, double t) override;

	void getInternalForce(GeomState *refState, GeomState &curState, CoordSet &cs,
						  FullSquareMatrix &elk, double *f, double dt, double t) override;

	void extractDeformations(GeomState &geomState, CoordSet &cs,
									 double *vld, int &nlflag) override;

	void getNLVonMises(Vector& stress, Vector& weight, GeomState &,
					   GeomState *, CoordSet &, int strIndex, int surface,
					   double ylayer, double zlayer, int avgnum, int measure) override;

	void getNLAllStress(FullM &stress, Vector &weight, GeomState &curState,
						GeomState *refState, CoordSet &c0, int strInd, int surface,
						int measure) override;

	void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0, double dt) override;

	bool checkElementDeletion(GeomState &curState) override;

	double getElementEnergy(GeomState &, CoordSet &) override;
	double getDissipatedEnergy(GeomState &geomState, CoordSet &cs) override;

	int getNumGaussPoints() const;
};

#endif
