#ifndef _SUB_DOMAIN_H_
#define _SUB_DOMAIN_H_

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>
#include <Feti.d/DistrVector.h>
#include <Feti.d/FetiSub.h>
#include <Corotational.d/Corotator.h>
#include <Math.d/DistVector.h>
#include <Utils.d/MyComplex.h>
#include <Math.d/FsiSparse.h>
#include <Driver.d/SComm.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/GlobalToLocalMap.h>
#include <Utils.d/MathUtils.h>
#include <Rom.d/VecBasis.h>
#include <Rom.d/DistrVecBasis.h>
#include <vector>
#include <list>

extern GeoSource *geoSource;

template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class Scalar> class GenDistrDomain;
typedef GenDistrDomain<double> DistrDomain;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenFetiSolver;
typedef GenFetiSolver<double> FetiSolver;
template <class Scalar> class GenFetiDPSolver;
typedef GenFetiDPSolver<double> FetiDPSolver;
template <class Scalar> class GenFetiOp;
typedef GenFetiOp<double> FetiOp;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class IntFullM;
class GeomState;
class DistrGeomState;
template <class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
template <class Scalar> class GenSkyMatrix;
typedef GenSkyMatrix<double> SkyMatrix;
template <class Scalar> class GenDomainGroupTask;
typedef GenDomainGroupTask<double> DomainGroupTask;
template <class Scalar> class GenAssembledFullM;
typedef GenAssembledFullM<double> AssembledFullM;
template <class Scalar> class GenSparseSet;
typedef GenSparseSet<double> SparseSet;
class FetiSubCornerHandler;
template <class Type> class FSCommPattern;
class DistrComplexVector;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenMpcSparse;
//template <class Scalar, class GenVecType> class GenVecBasis;
//typedef GenVecBasis<double, DistrVector> DistrVecBasis;

class BaseSub : virtual public Domain , virtual public FetiBaseSub
{
protected:

	GlobalToLocalMap glToLocalElem;
	int *glElems = nullptr;
	int glNumNodes;
	double *bcx = nullptr;
	DComplex *bcxC = nullptr; // FETI-H
	double *vcx = nullptr, *acx = nullptr;
	int *locToGlSensorMap = nullptr;
	int *locToGlActuatorMap = nullptr;
	int *locToGlUserDispMap = nullptr;
	int *locToGlUserForceMap = nullptr;
	int *crnPerNeighb = nullptr;
#ifdef DISTRIBUTED
	// for distributed output of single nodes
	int numNodalOutput = 0;
	int *outputNodes = nullptr;
	int *outIndex = nullptr;
#endif
	int globalNMax;  // highest global node number of all the nodes in this subdomain
	int globalEMax;  // highest global element number of all the elements in this subdomain


public:
	BaseSub();
	BaseSub(Domain &dom, int sn, Connectivity &con, Connectivity &nds, int gn);
	BaseSub(Domain &dom, int sn, int nNodes, int *nds,
	        int nElems, int *elems, int gn);
	BaseSub(Domain &dom, int sn, CoordSet* nodes, Elemset* elems, int *glNodeNums,
	        int *glElemNums, int gn); // PJSA: for new sower
	virtual ~BaseSub();

	const FetiInfo &getFetiInfo() const override { return solInfo().getFetiInfo(); }
	int *localToGlobalFSI = nullptr;  // local to global FSI numbering

	int crnDofLen() const  { return crnDofSize; }
	IntFullM* getC(int &crnDofSize, FSCommPattern<int> *sPat);
	void showExchangeData();
	void countCornerDofs(int *cWeight);
	void applyAuxData();
	void distributeBCs(int *);
	void setControlData(ControlLawInfo *_claw, int *, int *, int *, int *);
#ifdef DISTRIBUTED
	void setOutputNodes(int, int *, int *);
	int *getOutputNodes()       { return outputNodes; }
	int *getOutIndex()          { return outIndex; }
	int getNumNodalOutput() const { return numNodalOutput; }
#endif
	int getBC(BCond *, int, int *, BCond *&);
	int *getGlElems() const    { return glElems; }
	int *getGlMPCs()  const     { return localToGlobalMPC; }
	int glToPackElem(int e) const { return (geoSource->glToPackElem(e) > globalEMax) ? -1 : glToLocalElem[geoSource->glToPackElem(e)]; }
	int *getSensorDataMap() const { return locToGlSensorMap; }
	int *getUserDispDataMap() const { return locToGlUserDispMap; }
	int countElemNodes();
	int globalNumNodes();
	int numNodes() const        { return numnodes; }
	const Connectivity *getNodeToNode() const override { return nodeToNode.get(); }
	int findProfileSize();
	int renumberBC(int *);
	void makeGlobalToLocalNodeMap();
	void makeGlobalToLocalElemMap();
	int globalToLocal(int i)    { return (i < 0 || i > globalNMax) ? -1 : glToLocalNode[i]; }
	int localToGlobal(int i)    { return glNums[i]; }
	int globalToLocalElem(int i) { return (i < 0 || i > globalEMax) ? -1 : glToLocalElem[i]; }
	int localToGlobalElem(int i) { return glElems[i]; }
	int getGlobalNMax()         { return globalNMax; }
	int getNumUncon() const override { return numUncon(); }
	void putNumMPC(int *ptr) { ptr[subNumber] = numMPC; }
	void putLocalToGlobalMPC(int *ptr, int *tg) { for(int i=0; i<numMPC; ++i) tg[ptr[subNumber]+i] = localToGlobalMPC[i]; }
	void putNumMPC_primal(int *ptr) { ptr[subNumber] = numMPC_primal; }
	void putLocalToGlobalMPC_primal(int *ptr, int *tg) { for(int i=0; i<numMPC_primal; ++i) tg[ptr[subNumber]+i] = localToGlobalMPC_primal[i]; }
	void putNumFSI(int *ptr) { ptr[subNumber] = numFSI; }
	void putLocalToGlobalFSI(int *ptr, int *tg) { for(int i=0; i<numFSI; ++i) tg[ptr[subNumber]+i] = localToGlobalFSI[i]; }

	void makeMpcInterface(Connectivity *subToMpc, const Connectivity &lmpcToSub,
	                      Connectivity *subToSub_mpc);
	void makeFsiInterface(const Connectivity *subToFsi, const Connectivity &fsiToSub,
	                      const Connectivity *subToSub_fsi);

	bool checkForColinearCrossPoints(int numCornerPoints, int *localCornerPoints);
	void addCornerPoints(int *glCornerList);

	const DofSetArray * getDsa() const override { return dsa; }
	const ConstrainedDSA * get_c_dsa() const override { return c_dsa; }
	double getShiftVal() const override { return geoSource->shiftVal(); }

	void addFsiElements();
	void addSingleFsi(LMPCons *localFsi);
	void addNodeXYZ(double *centroid, double* nNodes);

	void setCorners(gsl::span<const lc_node_idx> crnList) override;

	const bool* getInternalMasterFlag();

	void setDofPlusCommSize(FSCommStructure *) const;

	// for timing file
	double getSharedDofCount();
	int getTotalDofCount();
	FetiSubCornerHandler *getCornerHandler();

	void initHelm(Domain &dom);

	// DPH functions
	int isFluid(int i=0);
	void setWaveNumbers(double *waveNumbers);
	void computeWaveNumbers() override;
	void averageMatProps() override;

	void setDirichletBC(std::list<BCond *> *_list);
	void setNeumanBC(std::list<BCond *> *_list);
	void setInitialDisplacement(std::list<BCond *> *_list);
	void setInitialDisplacement6(std::list<BCond *> *_list);
	void setInitialVelocity(std::list<BCond *> *_list);
	void setSensor(std::list<BCond *> *_list);
	void setActuator(std::list<BCond *> *_list);
	void setUsdd(std::list<BCond *> *_list);
	void setUsdf(std::list<BCond *> *_list);
	void setClaw(char* _fileName, char* _routineName) {
		claw = new ControlLawInfo;
		claw->fileName = _fileName;
		claw->routineName = _routineName;
	}
	void setComplexDirichletBC(std::list<ComplexBCond *> *_list);
	void setComplexNeumanBC(std::list<ComplexBCond *> *_list);
	void setDnb(std::list<SommerElement *> *_list);
	void setScat(std::list<SommerElement *> *_list);
	void setArb(std::list<SommerElement *> *_list);
	void setWet(std::list<SommerElement *> *_list);
	// coupled_dph
protected:
	std::vector<bool> wetInterfaceMark;
	std::vector<bool> wetInterfaceFluidMark;
	std::vector<bool> wetInterfaceStructureMark;
	std::vector<bool> wetInterfaceCornerMark;

	int *wDofToNode = nullptr; //HB
	Connectivity *drySharedNodes = nullptr;
	int numFsiNeighb;
	int *fsiNeighb = nullptr;

public:
	/// \brief Non-owning pointer to the global(?) node to subdomain connectivity.
	Connectivity *nodeToSub;
	SparseConnectivityType2 *nodeToSub_sparse;
	void setnodeToSubConnectivity(Connectivity *nTsubConn) { nodeToSub = nTsubConn; }
	void setnodeToSubConnectivity(SparseConnectivityType2 *nTsubConn) { nodeToSub_sparse = nTsubConn; }
	void markWetInterface(int nWI, int *wiNum, bool soweredInput);
	bool onWetInterface(int iNode) const override { return wetInterfaceMark[iNode]; }
	bool onWetInterfaceFluid(int iNode) { return wetInterfaceFluidMark[iNode]; }
	bool onWetInterfaceStructure(int iNode) { return wetInterfaceStructureMark[iNode]; }
	bool isWetInterfaceCorner(int iNode) const override { return wetInterfaceCornerMark[iNode]; }
	void setWetInterface(int nWI, int *wiNum, bool soweredInput);
	void makeDSA();
	void makeCDSA();
	void makeCCDSA();

#ifdef HB_COUPLED_PRECOND
	Connectivity* precNodeToNode;
#endif
};


template<class Scalar>
class GenSubDomain : public BaseSub , public FetiSub<Scalar>
{
private:
	int *cornerWeight;
	void initialize();

protected:
	/// \brief Send to neighbors the list of DOFs on the shared nodes
	void sendDOFList(FSCommPattern<int> *pat) const;

public:

	std::shared_ptr<GenSparseMatrix<Scalar>>   MPCsparse;
	Corotator           	    **corotators;

	Scalar *bcx_scalar;

public:
	GenSubDomain(int, int);
	GenSubDomain(Domain &, int sn, Connectivity &con, Connectivity &nds, int gn);
	GenSubDomain(Domain &, int sn, int nNodes, int *nds,
	             int nElems, int *elems, int gn);
	GenSubDomain(Domain &, int sn, CoordSet* nodes, Elemset* elems, int *glNodeNums,
	             int *glElemNums, int gn); // PJSA: for new sower
	~GenSubDomain();

	long getMemoryK() { Scalar s; return memK*long(sizeof(s))/long(1024); /* PJSA return memory usage in KB */ }
	long getMemoryPrec() { Scalar s; return memPrec*long(sizeof(s))/long(1024); /* PJSA return memory usage in KB */ }
	void extractControlData(Scalar *, Scalar *, Scalar *,
	                        Scalar *, Scalar *, Scalar *);
	void addUserForce(Scalar *, Scalar *);
	void addCtrl(Scalar *, Scalar *);
	Scalar *getBcx()  { if(!bcx_scalar) makeBcx_scalar(); return bcx_scalar; }
	double *getVcx()  { return vcx; }
	double *getAcx()  { return acx; }
	void setUserDefBC(double *, double *, double *, bool nlflag);
	void reBuildKbb(FullSquareMatrix *kel);
	void addDMass(int glNum, int dof, double m);

	void extractAndSendInterf(const Scalar *subvec, FSCommPattern<Scalar> *pat) const;
	void assembleInterf(Scalar *subvec, FSCommPattern<Scalar> *pat) const;
	void assembleInterfInvert(Scalar *subvec, FSCommPattern<Scalar> *pat) const;
	/** \brief Renumber the element nodes to local numbers. */
	void renumberElements();
	void renumberElementsGlobal();
	void renumberSharedNodes();
	void renumberDirichlet();
	/** \brief Renumber nodes in boundary condition data to use local numbers. */
	void renumberBCsEtc();
	void renumberControlLaw();
	void renumberMPCs();
	void sendNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat);
	void collectNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat);

	void getSRMult(const Scalar *lvec, const Scalar *lbvec, int nRBM, const double *locRBMs, Scalar *alpha) const;
	void sendDeltaF(const Scalar *deltaF, FSCommPattern<Scalar> *vPat);
	double collectAndDotDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat);
	void multFi(GenSolver<Scalar> *s, Scalar *, Scalar *);
	void assembleLocalComplexEls(GenSparseMatrix<Scalar> *Kas, GenSolver<Scalar> *smat = 0);
	void mergePrimalError(Scalar* error, Scalar* primal);
	void mergeStress(Scalar *stress, Scalar *weight,
	                 Scalar *globStress, Scalar *globWeight, int glNumNodes);
	void mergeElemStress(Scalar *loc, Scalar *glob, const Connectivity *glElemToNode);
	void mergeDisp(Scalar (*xyz)[11], GeomState* locGS, Scalar (*xyz_loc)[11] = NULL);
	void mergeAllDisp(Scalar (*xyz)[11], Scalar *d, Scalar (*xyz_loc)[11] = NULL);
	void mergeAllVeloc(Scalar (*xyz)[11], Scalar *v, Scalar (*xyz_loc)[11] = NULL);
	void mergeAllAccel(Scalar (*xyz)[11], Scalar *a, Scalar (*xyz_loc)[11] = NULL);
	void forceContinuity(Scalar *locdisp, Scalar (*xyz)[11]);
	void mergeDistributedNLDisp(Scalar (*xyz)[11], GeomState* u, Scalar (*xyz_loc)[11] = NULL);
	void mergeForces(Scalar (*mergedF)[6], Scalar *subF);
	void mergeReactions(Scalar (*mergedF)[11], Scalar *subF);
	void mergeDistributedForces(Scalar (*mergedF)[6], Scalar *subF);
	void mergeDistributedReactions(Scalar (*mergedF)[11], Scalar *subF);
	void mergeElemProps(double* props, double* weights, int propType);
	template<class Scalar1> void dispatchNodalData(FSCommPattern<Scalar> *pat, NewVec::DistVec<Scalar1> *);
	template<class Scalar1> void addNodalData(FSCommPattern<Scalar> *pat, NewVec::DistVec<Scalar1> *);
	void dispatchInterfaceGeomState(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void collectInterfaceGeomState(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void dispatchInterfaceGeomStateDynam(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void collectInterfaceGeomStateDynam(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void dispatchInterfaceNodalInertiaTensors(FSCommPattern<double> *pat);
	void collectInterfaceNodalInertiaTensors(FSCommPattern<double> *pat);
	void dispatchGeomStateData(FSCommPattern<double> *, GeomState *);
	void collectGeomStateData(FSCommPattern<double> *, GeomState *);
	void computeElementForce(int, Scalar *u, int Findex, Scalar *force);
	void computeElementForce(GeomState *gs, Corotator **allCorot,
                                 int, int Findex, Scalar *force);
	void computeStressStrain(int, Scalar *u, int Findex,
	                         Scalar *stress, Scalar *weight = 0);
	void computeStressStrain(GeomState *gs, Corotator **allCorot,
	                         int, int Findex, Scalar *glStress, Scalar *glWeight = 0, GeomState *refState = NULL);
	void updatePrescribedDisp(GeomState *geomState, Scalar deltaLambda);
	void updatePrescribedDisp(GeomState *geomState);
	Scalar displacementNorm(Scalar *displacement);
	void firstAssemble(GenSparseMatrix<Scalar> *K);
	void makeKccDofsExp2(int nsub, GenSubDomain<Scalar> **sd, int augOffset,
	                     Connectivity *subToEdge);
	void deleteKcc();
	void multQt(int glMPCnum, const Scalar *x, Scalar *result) const;

	friend class GenDistrDomain<Scalar>;
	friend class GenDecDomain<Scalar>;
	friend class GenFetiOp<Scalar>;
	friend class GenFetiSolver<Scalar>;

	void precondGrbm();
	void setMpcSparseMatrix();
	void constructKrc();
	void initSrc();
	// MPC and contact functions
	void extractMPCs(int glNumMPC, ResizeArray<LMPCons *> &lmpc);
	void extractMPCs_primal(int glNumMPC, ResizeArray<LMPCons *> &lmpc);
	void printLMPC();
	void applySplitting();
	void applyDmassSplitting();
	void applyForceSplitting();
	void applyMpcSplitting();
	void constraintProduct(int num_vect, const double* R[], Scalar** V, int trans);
	void addConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &f);
	void addCConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &fc, double s);
	void dualConstraintProjection(std::vector<std::map<int,double> > &W, Rom::DistrVecBasis &CtW, Eigen::Matrix<double,Eigen::Dynamic,1> &WtRhs,
	                              int startCol, int blockCols);
	void locateMpcDofs();
	void deleteMPCs();

	void bmpcQualify(std::vector<LMPCons *> *bmpcs, int *pstatus, int *nstatus);

	void printMpcStatus();
	void computeContactPressure(Scalar *globStress, Scalar *globWeight);
	void computeLocalContactPressure(Scalar *stress, Scalar *weight);

	void getConstraintMultipliers(std::map<std::pair<int,int>,double> &mu, std::vector<double> &lambda);
	void getLocalMultipliers(std::vector<double> &lambda);
	void setMpcRhs(Scalar *interfvec, double t, int flag);
	void updateMpcRhs(Scalar *interfvec);

	const CoordSet &getNodeSet() const override { return getNodes(); }

protected:
	double *mpcForces;

	Scalar Bcx(int i);
	void makeBcx_scalar();
#ifdef HB_COUPLED_PRECOND
	Scalar* kSumWI;
#endif

public:
	void addSommer(SommerElement *ele); // XDEBUG

	// frequency sweep
	GenSparseMatrix<Scalar> *M;
	GenCuCSparse<Scalar> *Muc;
	GenSparseMatrix<Scalar> *C;
	GenCuCSparse<Scalar> *Cuc;
	GenSparseMatrix<Scalar> **C_deriv;
	GenSparseMatrix<Scalar> **Cuc_deriv;
	int numC_deriv;
	GenSparseMatrix<Scalar> **K_deriv;
	GenSparseMatrix<Scalar> **Kuc_deriv;
	int numK_deriv;
	int num_K_arubber;
	GenSparseMatrix<Scalar> **K_arubber_l;
	GenSparseMatrix<Scalar> **K_arubber_m;
	GenSparseMatrix<Scalar> **Kuc_arubber_l;
	GenSparseMatrix<Scalar> **Kuc_arubber_m;

private:
	// frequency sweep
	void makeFreqSweepLoad(Scalar *load, int iRHS, double omega);
	GenVector<Scalar> **a, **b;  // pade P, Q coefs
	int ia, ib;
	GenVector<Scalar> *P, *Q; // pade P(x), Q(x)
	bool rebuildPade;

public:
	void multM(Scalar *localrhs, GenStackVector<Scalar> **u, int k);
	void multMCoupled1(Scalar *localrhs, GenStackVector<Scalar> **u, int k,
	                   FSCommPattern<Scalar> *wiPat);
	void multMCoupled2(Scalar *localrhs, FSCommPattern<Scalar> *wiPat);
	void multWCAWE(Scalar *localrhs, GenStackVector<Scalar> **u,
	               Scalar *pU, Scalar *pb, int maxRHS, int iRHS);

	void pade(GenStackVector<Scalar> *sol,  GenStackVector<Scalar> **u, double *h, double x);
	void setRebuildPade(bool _rebuildPade) { rebuildPade = _rebuildPade; }

	int *l2g;
	void makeLocalToGlobalDofMap();
	void multAddLT(const Scalar *localvec, Scalar *globalvec);
	void multAddLinv(const Scalar *localvec, Scalar *globalvec);
	void multLTinv(const Scalar *globalvec, Scalar *localvec);
	void multL(const Scalar *globalvec, Scalar *localvec);

};

typedef GenSubDomain<double> SubDomain;


#endif
