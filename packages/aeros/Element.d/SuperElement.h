#ifndef _SUPERELEMENT_H_
#define _SUPERELEMENT_H_

#include <Element.d/Element.h>
#include <map>
class SuperCorotator;

class SuperElement : public Element 
{
private:
	Elemset *eset;
	DofSetArray *dsa;
	SuperCorotator *superCorotator;
	int nInternalNodes;
	double **sub_extf;

protected:
	CoordSet *css;
	Element **subElems;
	int nSubElems;
	int **subElemDofs;
	int **subElemNodes;
	int nnodes; // not including internal nodes
	int ndofs;
	int *nn; // all the node numbers
	bool localFlag;

	FullSquareMatrix stiffness(const CoordSet& cs, double *k, int flg=1) const override;
public:
	explicit SuperElement(bool = false);

	~SuperElement() override;
	Category getCategory() const override;

	double * getPreviouslyComputedSubExternalForce(int i) { return (sub_extf) ? sub_extf[i] : 0; }

	int getNumSubElems() { return nSubElems; }
	int getSubElemNumDofs(int i) { return subElems[i]->numDofs(); }
	int getSubElemNumNodes(int i) { return subElems[i]->numNodes(); }
	int* getSubElemDofs(int i) { return subElemDofs[i]; }
	int* getSubElemNodes(int i) { return subElemNodes[i]; }

	void setPressure(PressureBCond *) override;
	PressureBCond* getPressure() override;

	void renum(const int *table) override;
	void renum(EleRenumMap&) override;
	void setGlNum(int gn, int sn = 0) override;

	void setProp(StructProp *p, bool _myProp) override;
	void setPreLoad(std::vector<double> &load) override;
	std::vector<double> getPreLoad() override;
	void setFrame(EFrame *frame) override;
	void buildFrame(CoordSet &cs) override;
	void setOffset(double *o) override;
	void setCompositeData(int _type, int nlays, double *lData,
	                      double *coefs, double *frame) override;
	double * setCompositeData2(int _type, int nlays, double *lData,
	                               double *coefs, CoordSet &cs, double theta) override;
	void setMaterial(NLMaterial *) override;

	FullSquareMatrix massMatrix(const CoordSet& cs, double *m, int cmflg=1) const override;
	void getStiffnessThicknessSensitivity(CoordSet& cs, FullSquareMatrix &dStiffdThick, int flg=1) override;
	void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs) override;

	double getMass(const CoordSet&) const override;
	double getMassThicknessSensitivity(CoordSet&) override;
	double weight(CoordSet&, double *) override;
	double getWeightThicknessSensitivity(CoordSet&, double *) override;
	void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration) override;
	void getGravityForce(CoordSet &cs, double *gravity, Vector &force,
	                     int gravflg, GeomState *gs=0) override;
	void getGravityForceThicknessSensitivity(CoordSet &cs, double *gravity, Vector &forceSen,
	                                         int gravflg, GeomState *gs=0) override;
	void getGravityForceNodalCoordinateSensitivity(CoordSet &cs, double *gravity, GenFullM<double> &,
	                                               int gravflg, GeomState *gs=0) override;
	void getThermalForce(CoordSet &cs, Vector &ndT, Vector &force,
	                     int glflag, GeomState *gs=0) override;
	void getIntrnForce(Vector &elForce, CoordSet &cs,
	                   double *elDisp, int Index, double *ndTemps) override;
	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
	                 Vector &elDisp, int strInd, int surface=0,
	                 double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0) override;
	void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
	                                     int surface, double *, int avgnum, double ylayer, double zlayer) override;
	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
	                                        CoordSet &cs, Vector &elDisp, int strInd, int surface,
	                                        double *ndTemps, int avgnum, double ylayer, double zlayer) override;
	void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp,
	                                           int strInd, int surface, double * = 0, int avgnum=1, double ylayer=0,
	                                           double zlayer=0) override;
	void getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
	                  Vector &elDisp, int strInd, int surface=0,
	                  double *ndTemps=0) override;
	void computeHeatFluxes(Vector &heatflux, CoordSet &cs, Vector &elTemp, int hflInd) override;
	void trussHeatFluxes(double &trussflux, CoordSet &cs, Vector &elTemp, int hflInd) override;
	void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0) override;
	void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *, double *res, GeomState *gs=0) override;
	void computeTemp(CoordSet &cs, State &state, double[2], double *res) override;
	void getFlFlux(double[2], double *flF, double *res) override;

	void markDofs(DofSetArray &dsa) const override;
	int* dofs(DofSetArray &dsa, int *p=0) const override;
	int numDofs() const override;
	int numNodes() const override;
	int* nodes(int *p=0) const override;

	Corotator *getCorotator(CoordSet &, double *, int = 2, int = 2) override;
	void computePressureForce(CoordSet& cs, Vector& elPressureForce,
	                          GeomState *gs=0, int cflg = 0, double t = 0) override;

	double* getMidPoint(CoordSet &cs) override;
	double* getCompositeData(int nl) override;
	double* getCompositeFrame() override;
	int getCompositeLayer() override;
	int dim() const override ;
	void addFaces(PolygonSet *pset) override;
	int numInternalNodes() const override;
	void setInternalNodes(int *in) override;
	bool isSafe() const override;
	bool isRotMidSideNode(int iNode) override;
	bool isMpcElement() override;
	//bool isRigidMpcElement(const DofSet & = DofSet::nullDofset, bool forAllNodes=false);
	bool isConstraintElement() override;
	bool isFreeplayElement() const override;

	int getMassType() const override;
	int getNumMPCs() override;
	LMPCons** getMPCs() override;
	void makeAllDOFs();

	int numStates() override;
	void setStateOffset(int) override;
	void initStates(double *) override;
};

#endif
