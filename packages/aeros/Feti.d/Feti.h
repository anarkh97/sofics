#ifndef _FETI_H_
#define _FETI_H_
#include <iostream>

#include <Feti.d/DistrVector.h>
#include <Threads.d/Paral.h>
#include <Timers.d/Timing.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Communicator.h>
#include <Utils.d/MyComplex.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/SComm.h>
#include <Utils.d/dbg_alloca.h>
#include <Feti.d/DistrVectorSet.h>
#include <Solvers.d/ParallelSolver.h>
#include <Feti.d/CCtSolver.d/CCtSolver.h>
#include <Feti.d/FetiBaseClass.h>

template <class Scalar> class GenFetiOp;
typedef GenFetiOp<double> FetiOp;
template <class Scalar> class GenFetiOpControler;
typedef GenFetiOpControler<double> FetiOpControler;
template <class Scalar> class GenFetiSolver;
typedef GenFetiSolver<double> FetiSolver;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class IntFullM;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class OrthoSet;
template <class Scalar> class GenGMRESOrthoSet;
template <class Scalar> class GenGCROrthoSet;
template <class Scalar> class GenCGOrthoSet;
template <class Scalar> class GenFetiWorkSpace;
typedef GenFetiWorkSpace<double> FetiWorkSpace;
template <class Scalar> class GenSymFullMatrix;
typedef GenSymFullMatrix<double> SymFullMatrix;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class FetiInfo;
class DistrGeomState;
template <class Scalar> class GenBigMatrix;
typedef GenBigMatrix<double> BigMatrix;
template <class Scalar> class GenSkyMatrix;
typedef GenSkyMatrix<double> SkyMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
class Rbm;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenBlockSky;
typedef GenBlockSky<double> BlockSky;
template <class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
template <class Type> class FSCommPattern;
template <class Scalar> class CCtSolver;

template <typename Scalar>
class FetiSub;

// The FetiSolver class is too big. Future versions need it to be leaner
template <class Scalar>
class GenFetiSolver  : public GenParallelSolver<Scalar>
{
protected:
	std::vector<FetiSub<Scalar> *> subdomains;
	GenSubDomain<Scalar> **sd = nullptr;
	int nsub = 0, glNumSub = 0;
	FetiInfo *fetiInfo = nullptr;
	int verboseFlag;
	FSCommunicator *fetiCom = nullptr;
	int myCPU, numCPUs;
	Connectivity *subToSub = nullptr, *mpcToSub = nullptr, *mpcToSub_primal = nullptr;
	Connectivity *edgeToSub = nullptr, *subToEdge = nullptr;
	Connectivity *coarseConnect = nullptr;  // first level coarse prob. connectivity
	compStruct *renum = nullptr;
	compStruct renumber;
	SimpleNumberer *eqNums = nullptr;
	SimpleNumberer *PFcNums = nullptr;
	int gOffset;
	int mOffset;
	GenSparseMatrix<Scalar> *singleCoarse = nullptr;
	GenSolver<Scalar> *singleCoarseSolver = nullptr;
	double epsilon2;
	int maxiter;
	mutable GenCGOrthoSet<Scalar> *oSetCG = nullptr; // Workspace
	mutable GenGMRESOrthoSet<Scalar> *oSetGMRES = nullptr;
	mutable GenGCROrthoSet<Scalar> *oSetGCR = nullptr;
	int numP;
	int numrbms = 0, halfSize = 0; // number of  rbms, half interface size
	int glNumRBM = 0;
	DistrInfo internalDI, interface;
	FSCommPattern<Scalar> *vPat = nullptr;
	FSCommPattern<Scalar> *rbmPat = nullptr;
	FSCommPattern<int> *sPat = nullptr;
	int *glSubToLoc = nullptr;
	int crns = 0; // sum of corner lambdas of all subdomains
	TaskDescr **fetiTasks = nullptr;
	GenFetiOp<Scalar> **fetiOps = nullptr;
	mutable GenFetiOpControler<Scalar> *opControl = nullptr;
	GenSolver<Scalar> *GtGsolver = nullptr;
	GenBigMatrix<Scalar> *PCtFPC = nullptr;
	GenFetiWorkSpace<Scalar> *wksp = nullptr;
	int QGisLocal;  // Whether or not QG is local
	int isDynamic;  // Whether or not we are in dynamics
	int isFeti2;    // Whether or not we are using FETI2
	mutable int numSystems = 0; // Nonlinear additions
	mutable Timings times;
	GenSymFullMatrix<Scalar> *GtQGs = nullptr;
	GenFullM<Scalar> *GtFCs = nullptr;
	Scalar *gtqglocal = nullptr;
	int numNodes = 0;
	Connectivity *cpuToSub = nullptr;
	int glNumMpc = 0, glNumMpc_primal = 0;

	void makeGtG();
	void makeDistGtG(int *glSubToLocal);
	void assembleDistGtQGs(int i, int *);
	bool isLowestLocalNeighbor(int subI, int subJ) const;
	void addNonLocalGtQG(int subI, int subJ);
	void addNonLocalGContrib(int subI, int subJ);
	void addNonLocalCContrib(int subI, int subJ);
	void getNonLocalGtQMult(int myNum, int neighbN, Scalar *va,
	                        GenDistrVector<Scalar> *dv) const;
	void getNonLocalFCtMult(int myNum, int neighbN, Scalar *va,
	                        GenDistrVector<Scalar> *dv) const;
	void getNonLocalSubAlphaGtQ(int subI, int subJ, Scalar *va,
	                            GenDistrVector<Scalar> *dv) const;
	void getNonLocalGtQMult(int subI, int subJ);
	void getGtQMult(int iSub, Scalar *, GenDistrVector<Scalar> *) const;
	void getFCMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void reBuildGtG();
	void makePCtFPC();
	void reBuildPCtFPC();
	void makelocalFcoarse();
	int  collectIntGlobalSum();
	void subRgcMult(int i, int nThreads, GenVector<Scalar>* alpha, GenVector<Scalar>* result) const;
	void subRgcTransMult(int i, int nThreads, GenVector<Scalar>* alpha, GenVector<Scalar>* result) const;
	void fSend(int i,  GenDistrVector<Scalar> &) const;
	void fScale(int i, GenDistrVector<Scalar> &) const;
	void fSplit(int iSub, GenDistrVector<Scalar> &force) const;
	void fSendCoupled(int iSub, GenDistrVector<Scalar> &force, GenDistrVector<Scalar> &fw) const;
	void fScaleCoupled(int iSub, GenDistrVector<Scalar> &force, GenDistrVector<Scalar> &fw) const;

	// Routines to create the single coarse problem
	void makeGandFG();
	void makeSingleCoarse();
	void singleCoarseAssembly();
	void singleCoarseAssembleG(int isub);
	void singleCoarseAssembleMPCs(int iSub);
	Connectivity getCoarseToSubConnect() const;

	Connectivity * makeSingleConnect(const Connectivity *coarseConnect,
	                                 const Connectivity *coarseToSub,
	                                 const Connectivity *subToCoarse, int gOffset);
	void singleCoarseSolve(const GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &u) const;
	void setAllOffsets(int iSub, int gOffset);
	void getSRMult(int iSub, GenDistrVector<Scalar> *r, GenDistrVector<Scalar> *lambda,
	               Scalar *alpha) const;
	void getSGtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *alpha) const;
	void getSCtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void getSQtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void getFGMult(int iSub, GenDistrVector<Scalar> *r, Scalar *alpha) const;
	void addG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void addSG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void addGs(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w, GenVector<Scalar> &) const;
	void singlePr(GenDistrVector<Scalar> &y, GenDistrVector<Scalar> &p, GenVector<Scalar> &) const;
	void addMpcRhs(int iMPC, Scalar *singleC) const;
	void preprocessMPCs(int iSub);
	// corner preprocessing
	void preProcessCorners();

	void addC(int iSub, GenDistrVector<Scalar> *lambda, Scalar *sv) const;
	void getQtKpBMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void getNeighbFGs(int iSub);

public:
	GenFetiSolver(int nsub, GenSubDomain<Scalar> **, Connectivity *,
	              FetiInfo *finfo, FSCommunicator *fetiCom, int *glToLoc,
	              Connectivity *mpcToSub, Connectivity *cpuToSub,
	              GenSolver<Scalar> **sysMatrices = 0, GenSparseMatrix<Scalar> **sysMat = 0,
	              Rbm **_rbms = 0, int verboseFlag = 0);
	GenFetiSolver(int _nsub, GenSubDomain<Scalar> **subs, int _numThreads, int _verboseFlag);
	virtual ~GenFetiSolver();

	void sendDeltaF(int iSub, GenDistrVector<Scalar>& deltaF) const;
	void normDeltaF(int iSub, double * subDots, GenDistrVector<Scalar>* deltaF) const;
	void factorMatrices(int isub);
	void sendScale(int isub);
	void collectScale(int isub);
	void getErrorEstimator(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &es);
	void interfaceDiff(int iSub, GenDistrVector<Scalar> &v) const;
	void rebuildInterface(int iSub, GenDistrVector<Scalar> &v) const;
	void multKbb(int iSub, const GenDistrVector<Scalar>& v, GenDistrVector<Scalar> &interfvec,
	             GenDistrVector<Scalar>& deltaU, GenDistrVector<Scalar> &deltaF, bool &errorFlag) const;
	void localSolve(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4) const;
	void localSolve2(int iSub, GenDistrVector<Scalar> *v1, GenDistrVector<Scalar> *v2,
	                 GenVector<Scalar> *beta, GenDistrVector<Scalar> *v3,
	                 GenDistrVector<Scalar> *v4) const;
	void interfSend(int iSub, GenDistrVector<Scalar> &dv1) const;
	//void interfDiff(int iSub, GenDistrVector<Scalar> &dv1);
	void interfDiffAndDot(int iSub, GenDistrVector<Scalar> &dv1, GenDistrVector<Scalar> &dv2) const;
	void getRMult(int iSub, GenDistrVector<Scalar> *localvec, Scalar *alpha) const;
	void getGtMult(int iSub, const GenDistrVector<Scalar> *localvec, Scalar *alpha) const;
	void addRP(int iSub, GenDistrVector<Scalar> * localvec, Scalar *alpha) const;
	void addRS(int iSub, GenDistrVector<Scalar> * localvec, Scalar *alpha) const;
	void solve(const GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &x) override;
	void distributeForce(GenDistrVector<Scalar> &force) const;
	void distributeForce(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw) const;
	void reSolve(GenDistrVector<Scalar> &) override;
	Timings& getTimers() override { return times; }
	void makeRbmPat();
	void makeSingleIntPat();

	// NonLinear Functions
	virtual void reBuild(FullSquareMatrix **kel, DistrGeomState& gs, int iter=0,
	                     int step = 1);
	void reBuildMatrices(int isub, FullSquareMatrix **kel);
	void reBuildErrorEstimator(GenFullSquareMatrix<Scalar> **kel);
	void subdomainReBuild(int isub, FullSquareMatrix **kel, DistrGeomState *gs);
	int nlPreCondition(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Pr) const;
	double getSolutionTime() const override { return times.solve; }
	void reSendInterfaceRBM(int iSub);
	void reGetNeighbQGs(int iSub);
	void reMakeLocalFcoarse();
	void reGetNeighbFC(int iSub);
	void reComputeFiBC(int iSub);
	void makeCompatible(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &deltaU) const;
	Scalar localSolveAndJump(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
	                         GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) const;
	Scalar localSolveAndJump(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
	                         GenVector<Scalar> &, GenDistrVector<Scalar> &,
	                         GenDistrVector<Scalar> &) const;
	void project(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &, int isDirect = 1) const;
	void tProject(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &, int isDirect = 1) const;
	void computeDynamL0(GenDistrVector<Scalar> &, GenVector<Scalar> &,
	                    GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) const;
	void computeL0(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &) const;
	void computeL0(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
	               GenVector<Scalar> &, GenDistrVector<Scalar> &) const;
	void addR(GenDistrVector<Scalar> &, GenVector<Scalar> &) const;
	void addRSingle(GenDistrVector<Scalar> &, GenVector<Scalar> &) const;
	void outputPrimalResidual(int iter, GenDistrVector<Scalar> &deltaF) const;

	// preconditioner returns the error estimator
	double preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag = true) const;

	void resetOrthoSet();
	void orthoAddCG(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp) const;
	void orthogonalize(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) const;

	// GMRES functions
	void initGMRES(GenDistrVector<Scalar> &p);
	double orthoAddGMRES(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
	void GMRESSolution(GenDistrVector<Scalar> &p);

	// GCR functions
	int predictGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0);
	void orthogonalizeGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fr,
	                      GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
	void orthoAddGCR(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar FpFp);

	/** \brief Create a guess for lambda0.
	 *
	 * @param r RHS for which a guess is desired.
	 * @param lambda0 The best guess lambda.
	 * @return Whether a prediction was made.
	 */
	bool predict(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0) const;

	DistrInfo &interfInfo() { return interface; }
	DistrInfo &localInfo()  { return internalDI; }
	GenDistrVector<Scalar> & getLambda() { return wksp->ret_lambda(); }
	int neq() { return internalDI.len; }

	void updateFeti2lambda(GenDistrVector<Scalar> &, GenDistrVector<Scalar>&,
	                       GenDistrVector<Scalar> &,
	                       GenVector<Scalar>&, GenVector<Scalar>&) const;
	void updateFeti2y(GenDistrVector<Scalar> &, GenDistrVector<Scalar>&,
	                  GenVector<Scalar>&, GenVector<Scalar>&) const;
	void computeRgcTransMult(GenVector<Scalar>&, GenVector<Scalar>&) const;
	void computeRgcMult(GenVector<Scalar>& beta, GenVector<Scalar>& result) const;
	void updateFeti2Vector(GenDistrVector<Scalar> &, GenVector<Scalar>&, GenVector<Scalar>&) const;

	void addAllFcoarse(GenFullM<Scalar> &);
	void finishRgc(int, int);

	void getCtMult(GenDistrVector<Scalar> &w, GenVector<Scalar> &gamma) const;

	void scatterHalfInterface(int iSub, Scalar *v1, GenDistrVector<Scalar> *v2) const;
	void gatherHalfInterface(int iSub, const GenDistrVector<Scalar> *v1,
	                         const GenDistrVector<Scalar> *v2,
	                         Scalar *v3, Scalar *v4) const;

	void setAndStoreInfo(int iter, double finalPrimal2, double finalDual2 ) const;
	void findProfileSizes(int iSub, int *subSizes);

	// For eigen problem
	int numRBM() override;
	void getRBMs(GenDistrVectorSet<Scalar> &) override;
	void getRBMs(Scalar *) override;
	void Ksolve(int iSub, GenStackDistVector<Scalar> &R);

	int halfOffset(int iSub) const { return fetiOps[iSub]->halfOffset; }
	int numNeighbor(int iSub) const;
	Scalar *interfaceBuffer(int iSub) { return fetiOps[iSub]->interfBuff; }
	virtual void clean_up();
	double getFNormSq(GenDistrVector<Scalar> &f) override;

	virtual void getLocalMpcForces(int iSub, double *mpcLambda) { };  // only implemented for DP

protected:
	FSCommPattern<Scalar> *wiPat = nullptr;
};


template<class Scalar>
class GenFetiDPSolver : public FetiBaseClass<Scalar>
{
	using FetiBaseClass<Scalar>::fetiInfo;
	using FetiBaseClass<Scalar>::verboseFlag;

	DistrInfo internalR, internalC, internalWI;
	DistrInfo *coarseInfo = nullptr;
	GenSolver<Scalar>         *KccSolver = nullptr;
	GenParallelSolver<Scalar> *KccParallelSolver = nullptr;
	Scalar *kccrbms = nullptr;
	GenSparseMatrix<Scalar> *KccSparse = nullptr;
	int glNumCorners = 0;
	Connectivity *cornerToSub = nullptr;
	DofSetArray *cornerEqs = nullptr;
	int mpcOffset = 0, augOffset = 0; // mpc equation offset for coarse grid
	bool rbmFlag;
	bool geometricRbms;
	enum StepType { CG, PROPORTIONING, EXPANSION };
	bool proportional;

public:
	GenFetiDPSolver(int nsub, int glNumSub, std::vector<FetiSub<Scalar> *> subdomains, const Connectivity *subToSub,
	                FetiInfo *finfo, FSCommunicator *fetiCom, int *glToLoc, const Connectivity *mpcToSub,
	                const Connectivity *mpcToSub_primal,
	                Connectivity *mpcToMpc, const Connectivity *mpcToCpu, const Connectivity *cpuToSub,
	                const Connectivity *bodyToSub = 0, std::vector<std::unique_ptr<GenSolver<Scalar>>> sysMatrices = {},
	                GenSparseMatrix<Scalar> **sysMat = 0, Rbm **rbms = 0, bool rbmFlag = 0,
	                bool geometricRbms = true, int verboseFlag = 0);
	virtual ~GenFetiDPSolver();

	void makeFc(int iSub, const GenDistrVector<Scalar> &fr, const GenDistrVector<Scalar> &lambda) const;
	void makeFcB(int iSub, GenDistrVector<Scalar> &bf) const;
	void KrrReSolve(int iSub, GenDistrVector<Scalar> &ur);
	void makeKcc();
	void makeMultiLevelDP(std::unique_ptr<const Connectivity> subToCorner);
	void makeMultiLevelDPNew(const Connectivity &subToCorner);
	void assembleFcStar(GenVector<Scalar> &FcStar) const;
	void mergeSolution(GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, GenDistrVector<Scalar> &u,
	                   GenDistrVector<Scalar> &lambda) const;
	void mergeUr(int iSub, GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, GenDistrVector<Scalar> &u,
	             GenDistrVector<Scalar> &lambda) const;
	/** \brief Compute \f$ f_r = K_{rc} u_c \f$ for subdomain iSub */
	void multKrc(int iSub, GenDistrVector<Scalar> &fr, const GenVector<Scalar> &uc) const;
	void extractFc(int iSub, const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fc) const;
	void getFc(const GenDistrVector<Scalar> &f, GenVector<Scalar> &fc) const;
	void makeEdgeConnectivity();
	void countEdges(int iSub, int *edges) const;
	void numberEdges(int iSub, int *eP, int *ep2, int *edges, FSCommPattern<int> *sPat);
	void receiveNeighbEdgeNums(int iSub, int *eP, int *edges, FSCommPattern<int> *sPat);
	void factorLocalMatrices(int isub);
	/** \brief Extract the components of f into corner, remainder and ?, modifies f.
	 *
	 * @param[inout] f
	 * @param[out] fr Remainder nodes components of the force vector.
	 * @param fc Corner nodes components of the force vector.
	 * @param fw
	 * @return \f$ \| f \|^2 if it is not zero, 1.0 otherwise.
	 */
	double extractForceVectors(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr,
	                           GenVector<Scalar> &fc, GenDistrVector<Scalar> &fw) const;
	void printSummary(int iter) const;
	double getFNormSq(GenDistrVector<Scalar> &f);
	void getRBMs(Scalar *);
	void getRBMs(GenDistrVectorSet<Scalar> &);
	void getGlobalRBM(int iSub, int &iRBM, GenDistrVector<Scalar> &R);
	int numRBM();
	void clean_up();
	void subdomainSolve(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                    GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                    GenVector<Scalar> &v5) const;
	void subdomainSolveCoupled1(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                            GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                            GenVector<Scalar> &v5) const;
	void subdomainSolveCoupled2(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                            GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                            GenVector<Scalar> &v5, GenDistrVector<Scalar> &fw) const;
	void subdomainSolveCoupled2b(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                            GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                            GenVector<Scalar> &v5) const;
	void localSolveAndJump(GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda,
	                       GenDistrVector<Scalar> &ur, GenVector<Scalar> &fc,
	                       GenVector<Scalar> &uc, GenDistrVector<Scalar> &r,
	                       GenDistrVector<Scalar> &fw) const;
	Scalar localSolveAndJump(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &dur,
	                         GenVector<Scalar> &duc, GenDistrVector<Scalar> &Fp) const;
	void solve(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	void solveCG(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	void solveGMRES(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	void solveGCR(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	double preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag = true) const;

	void subtractMpcRhs(int iSub, GenDistrVector<Scalar> &dv1) const;
	bool updateActiveSet(GenDistrVector<Scalar> &v, int flag, double tol = 0.0);
	void subUpdateActiveSet(int iSub, GenDistrVector<Scalar> &v, double tol, int flag, bool *statusChange);
	void subRecvMpcStatus(int iSub, FSCommPattern<int> *mpcPat, int flag, bool *statusChange);

	void projectActiveIneq(const GenDistrVector<Scalar> &x, GenDistrVector<Scalar> &y) const;
	void split(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &v_f, GenDistrVector<Scalar> &v_c) const;
	void update(Scalar nu, GenDistrVector<Scalar> &lambda,
	            GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fp,
	            GenDistrVector<Scalar> &ur, GenDistrVector<Scalar> &dur,
	            GenVector<Scalar> &uc, GenVector<Scalar> &duc) const;
	void saveStep() const;
	void restoreStep() const;

private:
	bool globalFlagCtc;
	mutable double t7;
	int *ngrbmGr = nullptr;
	int nGroups = 0, nGroups1 = 0;
	int *groups = nullptr;
	const Connectivity *groupToSub = nullptr, *bodyToSub = nullptr, *subToGroup = nullptr;
	GenSolver<Scalar> *GtGtilda = nullptr;
	GenSparseMatrix<Scalar> *GtGsparse = nullptr;
	const Connectivity *subToBody = nullptr;
	/// Statistic variables.
	mutable int nSubIterDual = 0, nSubIterPrimal = 0, nMatVecProd = 0, nRebuildGtG = 0, nRebuildCCt = 0;
	mutable int nLinesearch = 0, nLinesearchIter = 0, nStatChDual = 0, nStatChPrimal = 0;
	mutable bool dualStatusChange = false, primalStatusChange = false, stepLengthChange = false;
	int ngrbms = 0;
	Connectivity coarseConnectGtG;
	SimpleNumberer *eqNumsGtG = nullptr;
	Connectivity *mpcToMpc = nullptr;
	CCtSolver<Scalar> *CCtsolver = nullptr;
	bool mpcPrecon;  // mpc preconditioner flag, true = use generalized preconditioner with CCt
	// false = use scaling method (diagonal CCt) or no preconditioning
	const Connectivity *mpcToCpu = nullptr;
	int numSubsWithMpcs = 0;
	int *mpcSubMap = nullptr;
	void singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow, int &rank, double tol, FullM *V = 0);
	FSCommPattern<int> *mpcPat = nullptr;
	mutable double lastError;

	// Contact functions
public:
	void makeGtG();
	void deleteGtG() { std::cerr << "deleteGtG is not implemented\n"; };
	void trMultC(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f);
	void multC(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu);
private:
	void makeE(GenDistrVector<Scalar> &f) const; // Modifies the workspace.
	void assembleE(int iGroup, GenVector<Scalar> &e, GenDistrVector<Scalar> &f) const;
	void assembleGtG(int iGroup);
	void rebuildGtGtilda();
	void computeL0(GenDistrVector<Scalar> &lambda0, GenDistrVector<Scalar> &f) const;
	void normalizeC();
	void subTrMultC(int iSub, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f);
	void subMultC(int iSub, GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu) const;
	void project(GenDistrVector<Scalar> &z, GenDistrVector<Scalar> &y, int eflag = 0) const;
	double tProject(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w) const;
	void multG(const GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, double beta) const;
	void trMultG(const GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, double beta) const;
	void subTrMultG(int iGroup, const GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha) const;
	void addRalpha(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &alpha) const;
	void computeProjectedDisplacement(GenDistrVector<Scalar> &u) const;
	void addRstar_gT(int iGroup, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta) const;
	void subtractRstar_g(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta) const;
	bool checkStoppingCriteria(int iter, double error, double fnorm) const;

	// MPC & WI functions
public:
	void buildCCt();
	void rebuildCCt();
	void deleteCCt() { if(CCtsolver) delete CCtsolver; CCtsolver = 0; }
	Connectivity * getBlockToMpc();
	void cctSolveMpc(GenDistrVector<Scalar> &v) const;
	void getLocalMpcForces(int iSub, double *mpcLambda);
private:
	void addMpcRHS(int iMPC, Scalar *singleC) const;
	void wetInterfaceComms();  // coupled_dph
	void computeLocalWaveNumbers();
public:
	void reconstruct();
	void refactor();
	void reconstructMPCs(Connectivity *_mpcToSub, Connectivity *_mpcToMpc, Connectivity *_mpcToCpu);
	void zeroG();
	void deleteG();

	/** \brief Build the sub to corner connectivity into the member variable and return corner to subdomain. */
	Connectivity makeCornerToSub();
};

typedef GenFetiSolver<double> FetiSolver;
typedef GenFetiDPSolver<double> FetiDPSolver;
typedef GenFetiWorkSpace<double> FetiWorkSpace;

// HB: MPC stuff -> for creating "superblocks"
struct BlockPair {

	int Id;
	double cost;
	double bandwidth;

	BlockPair() { Id = 0; cost = 0.0; bandwidth = 0.0; }
	BlockPair(int _Id, double _cost, double _bandwidth)
	{ Id = _Id; cost = _cost; bandwidth = _bandwidth; }
	BlockPair(int _Id, double _cost)
	{ Id = _Id; cost = _cost; bandwidth = 0.0; }
	BlockPair(int _Id)
	{ Id = _Id; cost = 0.0; bandwidth = 0.0; }
	BlockPair(double _cost)
	{ Id = 0;  cost = _cost; bandwidth = 0.0; }
	bool operator <(const BlockPair &bp) const {
		return (cost < bp.cost) || ( (cost == bp.cost) && (Id < bp.Id) ) ;
	}
	void print() {
		filePrint(stderr," --- block ---\n");
		filePrint(stderr,"  # Id       : %d\n",Id);
		filePrint(stderr,"  # cost     : %e\n",cost);
		filePrint(stderr,"  # bandwidth: %e\n",bandwidth);
	}
};

#ifdef _TEMPLATE_FIX_
#ifdef DISTRIBUTED
#include <Dist.d/DistFeti.C>
#endif
#endif

#endif
