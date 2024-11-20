#ifndef _SUPERCOROTATOR_H_
#define _SUPERCOROTATOR_H_

#include <Corotational.d/Corotator.h>
#include <Element.d/SuperElement.h>

class SuperCorotator : public Corotator
{
	int nSubElems;
	Corotator **subElemCorotators;
	SuperElement *superElem;
	FullM *origK;
	double **sub_vld;
	double **sub_vlr;
	double **sub_dvld;

public:
	explicit SuperCorotator(SuperElement *_superElem);
	~SuperCorotator() override;

	void setSubCorotator(int i, Corotator *subCorotator)
	{ subElemCorotators[i] = subCorotator; }
	double* getPreviouslyExtractedSubDeformations(int i) { return (sub_vld) ? sub_vld[i] : 0; }
	double* getPreviouslyExtractedSubDeformationsSensitivities(int i) { return (sub_dvld) ? sub_dvld[i] : 0; }
	double* getPreviouslyExtractedSubRigidBodyMotion(int i) { return (sub_vlr) ? sub_vlr[i] : 0; }

	void getStiffAndForce(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t) override;
	void getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t) override;
	void getDExternalForceDu(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f) override;
	void getInternalForce(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t) override;
	void getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f, double dt, double t) override;
	void getExternalForce(GeomState &geomState, CoordSet &cs,  double *f) override;

	void formGeometricStiffness(GeomState &geomState, CoordSet &cs, FullSquareMatrix &k, double *f) override;
	double* getOriginalStiffness() override;
	void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld, int &nlflag) override;
	void getNLVonMises(Vector& stress, Vector& weight, GeomState &curState, GeomState *refState, CoordSet& c0, int strIndex,
					   int surface = 0, double ylayer = 0, double zlayer = 0, int avgnum = 0, int measure = -1) override;
	void getNLAllStress(FullM &stress, Vector &weight, GeomState &curState, GeomState *refState, CoordSet &c0, int strInd,
						int surface = 0, int measure = -1) override;
	double getElementEnergy(GeomState &geomState, CoordSet &cs) override;
	double getDissipatedEnergy(GeomState &geomState, CoordSet &cs) override;
	void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs, double *vlr) override;
	void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0, double dt = 0) override;
	bool checkElementDeletion(GeomState &curState) override;

	void getResidualCorrection(GeomState &gs, double *r) override;
	void initMultipliers(GeomState& c1) override;
	void updateMultipliers(GeomState& c1) override;
	double getError(GeomState& c1) override;

	bool useDefaultInertialStiffAndForce() override;
	void getInertialStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0,
								  FullSquareMatrix &elK, double *f, double dt, double t,
								  double beta, double gamma, double alphaf, double alpham) override;


	void getInternalForceThicknessSensitivity(GeomState *refState, GeomState &geomState, CoordSet &cs, Vector &dFintdThick,
											  double dt, double t) override;
	void getInternalForceNodalCoordinateSensitivity(GeomState *refState, GeomState &geomState, CoordSet &cs, Vector *&dFintdx,
													double dt, double t) override;
	void extractDeformationsDisplacementSensitivity(GeomState &geomState, CoordSet &cs, double *dvld) override;
};

#endif
