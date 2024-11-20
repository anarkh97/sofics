//
// Created by Michel Lesoinne on 11/3/17.
//

#ifndef FEM_FETUSUB_H
#define FEM_FETUSUB_H
#include <vector>
#include <gsl/span>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Math.d/matrix.h>
#include <Utils.d/GlobalToLocalMap.h>
#include <Driver.d/SComm.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/MpcSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/FsiSparse.h>
#include <Driver.d/Domain.h>
#include <Feti.d/FetiSub.h>

class FSCommStructure;
template <typename Scalar>
class FSCommPattern;
template <typename Scalar>
class GenVector;
template <typename Scalar>
class GenFullM;
template <typename Scalar>
class GenSolver;
template <typename Scalar>
class GenSparseMatrix;
template <typename Scalar>
class GenAssembledFullM;
template <typename Scalar>
class GenCuCSparse;
template <typename Scalar>
class SubLMPCons;
template <typename Scalar>
class GenSparseSet;

class Connectivity;
class CoordSet;
class DofSet;
class FetiInfo;

//#if (__GNUC__ < 7)
using vec_const_int = std::vector<int>;
//#else
//using vec_const_int = std::vector<const int>;
//#endif

/** \brief Pure Interface of what the notion of Subdomain provides for FETI solver. */
class FetiBaseSub {
public:
	virtual ~FetiBaseSub() = default;
	/// \brief Obtain the solver settings. TODO Get rid of this. Why should the subdomain data know all the solver details?
	virtual const FetiInfo &getFetiInfo() const = 0;
	/** \brief Obtain the size of the interface of this subdomain. */
	int interfLen() const;
	/** \brief Obtain the size of the half interface for which this subdomain is the master. */
	int halfInterfLen() const;
	/** \brief Obtain the index of this subdomain in the global system. */
	gl_sub_idx subNum() const { return subNumber; }
	/** \brief Set the communication size in the pattern.
	 * TODO abstract FSCommPattern to isolate the scalar type.
	 *
	 * @param numRBM The number of RBMs for which the setup should be done.
	 * @param pattern The pattern to adjust.
	 */
	void setRbmCommSize(int numRBM, FSCommStructure *pattern) const;
	/** \brief Set the communication pattern from this subdomain to its neighbors for a RHS or solution vector. */
	void setDofCommSize(FSCommStructure *) const;

	void setCommSize(FSCommStructure *pat, int size) const;

	void setMpcCommSize(FSCommStructure *mpcPat) const;

	void setMpcNeighbCommSize(FSCommPattern<int> *pt, int size) const;

	// Pure virtual defined in FetiSub<Scalar>.
	virtual void setGCommSize(FSCommStructure *pat) const = 0;

	lc_sub_idx localSubNum() const { return localSubNumber; }
	int localLen() const { return (cc_dsa) ? cc_dsa->size() : get_c_dsa()->size(); }
	int localRLen() const { return cc_dsa->size(); }

	/** \brief Obtain the global node numbers. */
	auto &getGlNodes() const { return glNums; }
	/** \brief Obtain the number of unconstrained degrees of freedom (includes R and C DOFs). */
	virtual int getNumUncon() const = 0;

	/** \brief Make the subdomain determine its master flag, i.e. whether it is the primary holder of a variable. */
	void computeMasterFlag(const Connectivity &mpcToSub);

	void mergeInterfaces();

	virtual void setCorners(gsl::span<const lc_node_idx> cNum);

	/// \brief Obtain the number of MPC constraints.
	int numMPCs() const  { return numMPC; }
	int numMPCs_primal() const  { return numMPC_primal; }
	/// \brief Obtain the number of coarse dofs. This method computes a cached value.
	int numCoarseDofs();

	int numCornerDofs()	const { return numCRNdof; }
	/// \brief Obtain the number of corner nodes.
	int numCorners() const  { return numCRN; }

	const std::vector<int> &getLocalCornerNodes() const { return cornerNodes; };

	int numWetInterfaceDofs() const { return numWIdof; }

	/// \brief Send to neighbors the list of DOFs on the shared nodes
	void sendDOFList(FSCommPattern<int> *pat) const;
	void gatherDOFListPlus(FSCommPattern<int> *pat);

	virtual const CoordSet& getNodeSet() const = 0;

	bool isEdgeNeighbor(int neighb) const { return scomm->isEdgeNeighb[neighb]; }

	virtual double getShiftVal() const = 0;

	void computeInternalMasterFlag();
	const std::vector<bool> &getMasterFlag() const { return masterFlag; }

	virtual const Connectivity *getNodeToNode() const = 0;

	// TODO Remove this from here. Only temporary to help refactor of CornerMaker into CornerSelector
	virtual bool onWetInterface(int iNode) const { return false; }
	virtual bool isWetInterfaceCorner(int iNode) const { return false; }

	const std::vector<int> &getCornerNodes() const { return glCornerNodes; }
	std::vector<int> &getCornerNodes() { return glCornerNodes; }
	void markCornerDofs(gsl::span<int> glCornerDofs) const;
	void makeKccDofs(DofSetArray *cornerEqs, int augOffset, Connectivity *subToEdge, int mpcOffset = 0);
	/** \brief Create the vector cornerEqNums with the corner degrees of freedom for 2 level FETI-DP.
	 *
	 * @param sd Vector of upper-level (coarse) subdomains.
	 * @param augOffset Offset of the augmentation
	 * @param subToEdge Connectivity between lower level subs and the edges between them.
	 */
	void makeKccDofsMultiLevel(gsl::span<FetiBaseSub *> sd, int augOffset, Connectivity *subToEdge);

	int numEdgeDofs(int i) const { return edgeDofSize[i]; }

	const std::vector<int> &getWeights() const { return weight; }
	std::vector<int> &getWeights() { return weight; }
	int dofWeight(int i) const { return weight[i]; }



	gsl::span<const int> getCornerEqNums() const { return cornerEqNums; }
	int getGroup() const { return group; }

	SComm *getSComm() { return scomm; }
	void setSComm(SComm *sc);
	/** \brief Obtain the number of neighbor subdomains. */
	int numNeighbors() const { return scomm->numNeighb;}
	int numEdgeNeighbors() const { return scomm->numEdgeNeighb; }

	bool* getMpcMaster() const { return mpcMaster; }
	// Multiple Point Constraint (MPC) functions
	int getNumMpc() const       { return numMPC; }

	/// \brief Obtain the DofSetArray constructed with all the DOFs of the subdomain.
	virtual const DofSetArray * getDsa() const = 0;
	/// \brief Obtain the DofSetArray after constrained DOFs have been removed.
	virtual const ConstrainedDSA * get_c_dsa() const = 0;
	/// \brief Obtain the DofSetArray after constrained and corner DOFs have been removed.
	const ConstrainedDSA * getCCDSA() const;

	/// \copydoc
	int getLocalMPCIndex(int globalMpcIndex) const;
	/// \copydoc
	int getGlobalMPCIndex(int localMpcIndex) const;
	void makeLocalMpcToGlobalMpc(const Connectivity *mpcToMpc);

	void addMPCsToGlobalZstar(FullM &globalZstar, int startRow, int startCol, int numCol) const;
	void addSPCsToGlobalZstar(FullM &globalZstar, int &zRow, int zColOffset) const;

	void setWIoneCommSize(FSCommStructure *pat) const;
	void setWICommSize(FSCommStructure *wiPat);
	void setWImapCommSize(FSCommPattern<int> *pat);
	bool isWetInterfaceDof(int d) const { return (wetInterfaceMap.size() > 0) ? (wetInterfaceMap[d] > -1) : false; }
	void GramSchmidt(double *Q, bool *isUsed, DofSet desired, int nQPerNeighb, bool isPrimalAugmentation);

	void setLocalMpcToBlock(const Connectivity *mpcToBlock, const Connectivity *blockToMpc);
	void findEdgeNeighbors();
	void zeroEdgeDofSize();

	virtual void computeWaveNumbers() = 0;
	void sendWaveNumbers(FSCommPattern<double> *kPat);
	void collectWaveNumbers(FSCommPattern<double> *kPat);
	virtual void averageMatProps() = 0;
	void sendMatProps(FSCommPattern<double> *matPat);
	void collectMatProps(FSCommPattern<double> *matPat);

public:
	std::vector<DofSet> cornerDofs;
	std::vector<DofSet> edgeDofs; //<! \brief "Or" of all the DOFs for each edge with a neighbor.

protected:
	gl_sub_idx subNumber;
	lc_sub_idx localSubNumber; // relevant when running in distributed

	/// \brief Vector of global indices for local nodes.
	vec_const_int glNums;

	SComm *scomm = nullptr;
	bool isCoupled = false; // TODO Ensure this is set or derived from some other info.
	int boundLen = 0;
	int internalLen = 0;
	int totalInterfSize;
	std::vector<int> allBoundDofs;
	std::vector<int> boundMap;
	std::vector<int> internalMap;

	std::vector<int> glBoundMap;
	std::vector<int> glInternalMap;

	/// \brief Corner nodes in local numbering.
	std::vector<int> cornerNodes;
	std::vector<bool> isCornerNode;   //<! \brief True for node which is a corner node; false otherwise.
	std::vector<int> glCornerNodes; //!< \brief Corner nodes in global numbering.
	int numCRN = 0;
	int numCRNdof = 0;
	int crnDofSize = 0;
	std::vector<std::vector<DofSet>> boundaryDOFs;
	std::vector<int> edgeDofSizeTmp;   // XXXX
	std::vector<int> edgeDofSize;      //<! \brief Number of edge DOF per neighbor.
	std::vector<int> cornerEqNums; //<! \brief unique equation numbers for subdomain corner dofs
	std::unique_ptr<ConstrainedDSA> cc_dsa; //!< DOFSet for non-corner free DOFs.
	std::vector<int> ccToC; //!< Mapping from cc_dsa to c_dsa. All indices are >= 0
	std::vector<int> cToCC; //!< Mapping from c_dsa to cc_dsa. Indices for corner DOFs are < 0.
	bool isMixedSub = false;
	bool isThermalSub = false;
	bool isUndefinedSub = false;
	bool isFluidSub = false;

	std::vector<int> weight; ///!< \brief DOF weights (i.e. number of subd sharing that dof).
	std::vector<int> weightPlus; ///!< \brief DOF weights (i.e. number of subd sharing that dof) including corner DOFs.

	// MPC related data
	Connectivity *localMpcToMpc = nullptr;
	Connectivity *localMpcToGlobalMpc = nullptr;
	bool *faceIsSafe = nullptr;
	int *localToGroupMPC = nullptr;
	std::vector<int> boundDofFlag;  // boundDofFlag[i] = 0 -> perfect interface dof  (not contact or mpc)
	// boundDofFlag[i] = 1 -> node-to-node contact interface dof
	// boundDofFlag[i] = 2 -> mpc dof, only used for rixen method, domain->mpcflag = 1
	// note: boundDofFlag[i] > 2 can be used for future extensions, eg mortar contact
	std::vector<bool> masterFlag; // masterFlag[i] = true if this sub is the "master" of allBoundDofs[i]
	bool *internalMasterFlag = nullptr;
	int masterFlagCount = 0;
	bool *mpcMaster = nullptr;  // mpcMaster[i] = true if this subd contains masterdof for mpc i
	Connectivity *mpcToDof = nullptr;
	Connectivity *localMpcToBlock = nullptr;
	Connectivity *blockToLocalMpc = nullptr;
	std::unique_ptr<Connectivity> blockToBlockMpc;
	Connectivity *localMpcToBlockMpc = nullptr;
	Connectivity *mpcToBoundDof = nullptr;
	double *localLambda = nullptr;  // used for contact pressure output
	std::vector<int> invBoundMap;
	int *mpclast = nullptr;
	int nCDofs = -1; // TODO Clean this up

	mutable int *mpcStatus;
	mutable bool *mpcStatus1, *mpcStatus2;

	std::vector<int> makeBMaps(const DofSetArray *dofsetarray=0);
	std::vector<int> makeIMaps(const DofSetArray *dofsetarray=0);
public:
	void setGroup(const Connectivity *subToGroup) { this->group = (*subToGroup)[subNum()][0]; }
	void setNumGroupRBM(int *ngrbmGr);
	void getNumGroupRBM(int *ngrbmGr);
	void makeLocalToGroupMPC(const Connectivity &groupToMPC);

	void setNodeCommSize(FSCommStructure *, int d = 1) const;

	GlobalToLocalMap &getGlobalToLocalNode() { return glToLocalNode; }

	int group = 0;
	/// @name MPC Data
	///@g{
	// Multiple Point Constraint (MPC) Data
	int numMPC = 0;             // number of local Multi-Point Constraints
	int *localToGlobalMPC = nullptr;  // local to global MPC numbering
	GlobalToLocalMap globalToLocalMPC; // alternative data structure for global to local MPC numbering
	// not a pointer so don't have to de-reference before using [] operator

	int numMPC_primal = 0;
	int *localToGlobalMPC_primal = nullptr;
	GlobalToLocalMap globalToLocalMPC_primal;
	///@}

	std::vector<int> cornerMap;

	void sendNumNeighbGrbm(FSCommPattern<int> *pat);
	void recvNumNeighbGrbm(FSCommPattern<int> *pat);

	void sendNumWIdof(FSCommPattern<int> *sPat) const;
	void recvNumWIdof(FSCommPattern<int> *sPat);
	void sendWImap(FSCommPattern<int> *pat);
	void recvWImap(FSCommPattern<int> *pat);

	void sendNeighbGrbmInfo(FSCommPattern<int> *pat);
	void receiveNeighbGrbmInfo(FSCommPattern<int> *pat);

	bool isWetInterfaceNode(int n) const { return (wetInterfaceNodeMap.size() > 0) ? (wetInterfaceNodeMap[n] > -1) : false; }

	// variables and routines for parallel GRBM algorithm and floating bodies projection
	// and MPCs (rixen method)
protected:
	int nGrbm = 0;
	int *neighbNumGRBMs = nullptr;

	int numGroupRBM = 0, groupRBMoffset = 0;
	int *neighbNumGroupGrbm = nullptr;
	int *neighbGroupGrbmOffset = nullptr;
	int numGlobalRBMs = 0;
	int *dualToBoundary = nullptr;

	std::unique_ptr<Rbm> rigidBodyModesG;

	int numWIdof = 0;  // number of dofs on the wet interface (both fluid and structure)
	int numWInodes = 0;  // number of nodes on the wet interface (both fluid and structure)
	std::vector<int> wetInterfaceMap;  // dof map
	std::vector<int> wetInterfaceNodeMap;
	std::vector<int> wetInterfaceNodes;
	std::vector<int> numNeighbWIdof;
	std::vector<int> wiInternalMap;
	std::vector<DofSet> wetInterfaceDofs;
	std::vector<bool> wiMaster;

	GlobalToLocalMap glToLocalWImap;

	GlobalToLocalMap *neighbGlToLocalWImap = nullptr;

	GlobalToLocalMap glToLocalNode; //!< This seems to be for coarse problem only.

	/// \brief store indices for possible rebuild (multiple LHS freq sweep)
	int edgeQindex[2] = {-1, -1};

	double prev_cscale_factor;

	/// \brief wave numbers for FETI-DPH for this subdomain only in SALINAS.
	double k_f = 0.0, k_p = 0.0, k_s = 0.0, k_s2 = 0.0;
	double *neighbK_p = nullptr, *neighbK_s = nullptr, *neighbK_s2 = nullptr, *neighbK_f = nullptr;  // neighbors' wave numbers
	double Ymod = 0.0, Prat = 0.0, Dens = 0.0, Thih = 0.0, Sspe = 0.0;  // Young's modulus, Poisson ration, density, thickness, speed of sound
	double *neighbYmod = nullptr, *neighbPrat = nullptr, *neighbDens = nullptr, *neighbThih = nullptr, *neighbSspe = nullptr;  // neighbor's values

	long memK = 0;       //!< memory necessary to store K(s)
	long memPrec = 0;    //!< memory necessary to store Preconditioner
};

template <typename Scalar>
class _AVMatrix : public Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> {
public:
	using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

	_AVMatrix() = default;

	using MatrixType::operator=;
	const Scalar *operator[](int i) const { return &(*this)(0, i); }
	Scalar *operator[](int i) { return &(*this)(0, i); }
};

/** \brief Pure Interface of what a the notion of Subdomain provides for FETI solver. */
template <typename Scalar>
class FetiSub : virtual public FetiBaseSub {
public:
	FetiSub() = default;

	void gatherDOFList(FSCommPattern<int> *pat);

	double densProjCoefficient(int dof) { return 1.0; } // Defined as virtual in SubDomain, but never defined otherwise.
	void multMFi(GenSolver<Scalar> *s, Scalar *, Scalar *, int numRHS) const;
	void getQtKQ(GenSolver<Scalar> *s);
	void getQtKQ(int iMPC, Scalar *QtKQ);
	Scalar getMpcRhs(int iMPC) const;
	Scalar getMpcRhs_primal(int iMPC) const;
	// TODO Figure out how to make this const.
	void sendDiag(GenSparseMatrix<Scalar> *s, FSCommPattern<Scalar> *vPat);
	void factorKii();
	void sendInterf(const Scalar *interfvec, FSCommPattern<Scalar> *vPat) const;
	void scatterHalfInterf(const Scalar *s, Scalar *loc) const;
	void getHalfInterf(const Scalar *s, Scalar *t) const;
	void getHalfInterf(const Scalar *s, Scalar *t, const Scalar *ss, Scalar *tt) const;
	void interfaceJump(Scalar *interfvec, FSCommPattern<Scalar> *vPat) const;
	void rebuildInterf(Scalar *v, FSCommPattern<Scalar> *vPat) const;
	void fSend(const Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw = 0) const;
	void splitInterf(Scalar *subvec) const;
	void sendDeltaF(const Scalar *deltaF, FSCommPattern<Scalar> *vPat);
	double collectAndDotDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat);
	/** \brief Basic FETI operation.
	 * \details Computes localvec = K-1 (localvec -B interfvec)
	 * then    interfvec = B^T localvec and sends local data to neighbors
	 *
	 * @param s
	 * @param localvec
	 * @param interfvec
	 */
	void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const;
	void fetiBaseOp(Scalar *uc,GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const;
	void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec, Scalar *beta) const;
	void fetiBaseOpCoupled2(const Scalar *uc, const Scalar *localvec, Scalar *interfvec,
	                        FSCommPattern<Scalar> *wiPat, const Scalar *fw = nullptr) const;
	void fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, const Scalar *interfvec,
	                                FSCommPattern<Scalar> *wiPat) const;
	void multQt(int glMPCnum, const Scalar *V, int numV, Scalar *QtV) const ;
	void multQtKBt(int glNumMPC, const Scalar *G, Scalar *QtKBtG, Scalar alpha=1.0, Scalar beta=1.0) const;
	/// \brief Retrieve the number of RBMs.
	int numRBM() const { return nGrbm; }
	const Scalar *getQtKpBt() const { return QtKpBt.data(); };
	void split(const Scalar *v, Scalar *v_f, Scalar *v_c) const;

	/// \brief Generate the B matrices for the problem
	void makeBs();

	void makeLocalMpcToDof(); //HB: create the LocalMpcToDof connectivity for a given DofSetArray
	void makeLocalMpcToMpc();
	void updateActiveSet(Scalar *v, double tol, int flag, bool &statusChange);

	// (G^T*G) matrix assembly
	void assembleGtGsolver(GenSparseMatrix<Scalar> *GtGsolver);
	void getLocalMpcForces(double *mpcLambda, DofSetArray *cornerEqs,
	                       int mpcOffset, GenVector<Scalar> &uc);
	void useKrrNullspace();
	// R matrix construction and access
	void makeLocalRstar(FullM **Qtranspose); // this is used by decomposed domain GRBM algorithm
	// R matrix-vector multiplication
	void addRalpha(Scalar *u, GenVector<Scalar> &alpha) const;  // u += R_g*alpha

	int zColDim() { return rigidBodyModesG->Zstar->numCol(); }
	int zRowDim() { return rigidBodyModesG->Zstar->numRow(); }

	void deleteLocalRBMs() { rigidBodyModesG.reset(nullptr); }
	void setBodyRBMoffset(int _boff) { bodyRBMoffset = _boff; }
	void assembleE(GenVector<Scalar> &e, Scalar *f) const; // e = R^T*f
	// G matrix construction and destruction
	void makeG();
	void makeTrbmG(Scalar *rbms, int nrbms, int glNumCDofs);

	void setGCommSize(FSCommStructure *pat) const override;
	void sendG(FSCommPattern<Scalar> *rbmPat);
	void receiveG(FSCommPattern<Scalar> *rbmPat);
	void zeroG();
	void deleteG();

	void getFr(const Scalar *f, Scalar *fr) const;
	void getFc(const Scalar *f, Scalar *fc) const;
	void getFw(const Scalar *f, Scalar *fw) const;
	// G matrix-vector multiplication
	void multG(const GenVector<Scalar> &x, Scalar *y, Scalar alpha) const;  // y = alpha*G*x
	void trMultG(const Scalar *x, GenVector<Scalar> &y, Scalar alpha) const; // y = alpha*G^T*x


	// R_g matrix construction and access
	void buildGlobalRBMs(GenFullM<Scalar> &Xmatrix, const Connectivity *cornerToSub); // use null space of (G^T*P_H*G) ... trbm method !!!
	void getGlobalRBM(int iRBM, Scalar *Rvec) const;
	// R_g matrix-vector multiplication
	void subtractRstar_g(Scalar *u, GenVector<Scalar> &beta) const; // u -= R_g*beta
	void addRstar_gT(Scalar *u, GenVector<Scalar> &beta) const; // u += R_g*beta
	// (R_g^T*R_g) matrix assembly
	void assembleRtR(GenFullM<Scalar> &RtRu);

	void makeKbbMpc();
	void makeKbb(const DofSetArray *dofsetarray = 0);
	void rebuildKbb();

	// new B operators
	void multBr(const Scalar *localvec, Scalar *interfvec, const Scalar *uc = 0, const Scalar *uw = 0) const;
	void multAddBrT(const Scalar *interfvec, Scalar *localvec, Scalar *uw = nullptr) const;

	void multAddCT(const Scalar *interfvec, Scalar *localvec) const;
	void multC(const Scalar *localvec, Scalar *interfvec) const;

	/** \brief Compute \f$ f_r = K_{rc} u_c \f$. */
	void multKrc(Scalar *fr, const Scalar *uc) const;

	/** \brief Form the \f$ K_{cc}^\star \f$ contribution for the subdomain. */
	void formKccStar();
	void multKbb(const Scalar *u, Scalar *Pu, Scalar *delta_u = 0, Scalar * delta_f= 0, bool errorFlag = true);
	void multKbbMpc(const Scalar *u, Scalar *Pu, Scalar *deltaU, Scalar *deltaF, bool errorFlag = true);
	void multKbbCoupled(const Scalar *u, Scalar *Pu, Scalar *deltaF, bool errorFlag = true);
	void multDiagKbb(const Scalar *u, Scalar *Pu) const;

	void collectScaling(FSCommPattern<Scalar> *vPat);
	void fScale(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw = 0);
	void initMpcScaling();
	void initScaling();
	void weightEdgeGs();
	void makeQ();

	void applyBtransposeAndScaling(const Scalar *u, Scalar *v, Scalar *deltaU = 0, Scalar *localw = 0) const;
	void applyScalingAndB(const Scalar *res, Scalar *Pu, Scalar *localw = 0) const;
	void setMpcDiagCommSize(FSCommStructure *mpcDiagPat) const;
	void sendMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);
	void collectMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);
	void sendMpcScaling(FSCommPattern<Scalar> *mpcPat);
	void collectMpcScaling(FSCommPattern<Scalar> *mpcPat);
	void assembleMpcIntoKcc();

	void makeEdgeVectorsPlus(bool isFluidSub = false, bool isThermalSub = false,
	                         bool isUndefinedSub = false);
	void makeAverageEdgeVectors();
	void extractInterfRBMs(int numRBM, Scalar *locRBMs, Scalar *locInterfRBMs);
	void sendInterfRBMs(int numRBM, Scalar *locInterfRBMs, FSCommPattern<Scalar> *rbmPat);
	void recvInterfRBMs(int iNeighb, int numNeighbRBM, Scalar *neighbInterfRBMs, FSCommPattern<Scalar> *rbmPat);
	void sendInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
	void receiveInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
	// templated R and G functions
	// note #1: we use feti to solve global domain problem: min 1/2 u_g^T*K_g*u_g - u_g^T*f_g subj. to C_g*u_g <= g
	//          by solving an equivalent decomposed domain problem: min 1/2 u^T*K*u - u^T*f subj to B*u = 0, C*u <= g
	// the columns of R_g span the left null space of [ K_g & C_gtilda^T // C_gtilda & 0 ] ... C_gtilda = gtilda are the active constraints
	// the columns of R span the left null space of K
	// G = [B^T C^T]^T*R
	// e = R^T*f
	// note #2: the null space of a matrix must be templated (i.e. it is real if the matrix is real or complex if the matrix is complex)
	// note #3: the geometric rigid body modes (GRBMs) or heat zero energy modes (HZEMs) can be used to construct the null space SOMETIMES, NOT ALWAYS !!!!
	// note #4: the GRBMs are not always computed correctly when there are mechanisms
	// note #5: the GRBMs/HZEMs are always real

	void addTrbmRalpha(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *alpha, Scalar *ur) const; // u += R_g*alpha
	void assembleTrbmE(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *e, Scalar *fr) const; // e = R^T*f

	void projectActiveIneq(Scalar *v) const;
	void normalizeCstep1(Scalar *cnorm);
	void normalizeCstep2(Scalar *cnorm);
	void recvMpcStatus(FSCommPattern<int> *mpcPat, int flag, bool &statusChange);

	void mergeUr(Scalar *ur, Scalar *uc, Scalar *u, Scalar *lambda);

	void multfc(const VectorView<const Scalar> &fr, /*Scalar *fc,*/ const VectorView<const Scalar> &lambda) const;
	void multFcB(Scalar *bf);

	void sendMpcStatus(FSCommPattern<int> *mpcPat, int flag);
	void subtractMpcRhs(Scalar *interfvec);
	void insertBlockMpcResidual(Scalar *subv, GenVector<Scalar> **mpcv, const Connectivity *mpcToBlock,
	                            SimpleNumberer **blockMpcEqNums);
	void extractBlockMpcResidual(int block, Scalar *subv, GenVector<Scalar> *mpcv,
	                             SimpleNumberer *blockMpcEqNums);
	void sendMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec);
	void combineMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec);
	void setLocalLambda(Scalar *_localLambda);
	double getMpcError() const;
	void applyMpcSplitting();

	void initMpcStatus();
	void saveMpcStatus();
	void restoreMpcStatus();
	void saveMpcStatus1() const; // const is a lie but we have mutable.
	void saveMpcStatus2();
	void cleanMpcData();

	void constructKcc();
	void constructKcw();
	void scaleAndSplitKww();
	void reScaleAndReSplitKww();
	void precondGrbm();

	void makeZstarAndR(double *centroid);  // makes Zstar and R

	const std::vector<Scalar> &getfc() const { return fcstar; }

	/// \brief Solver for the remainder DOFs.
	std::unique_ptr<GenSolver<Scalar>> Krr;
	/// \brief Sparse view alias of the Krr solver. \details Typically used to fill the matrix before calling factor.
	GenSparseMatrix<Scalar>   *KrrSparse = nullptr;
	/// \brief Solver the internal DOFs used in the Dirichlet Preconditionner. Alias to the sparse matrix.
	GenSolver<Scalar>         *KiiSolver = nullptr;
	std::unique_ptr<GenSparseMatrix<Scalar>> KiiSparse;
	std::unique_ptr<GenCuCSparse<Scalar> >     Kib;
	std::unique_ptr<GenAssembledFullM<Scalar>> Kcc;
	std::shared_ptr<GenCuCSparse<Scalar>>      Krc;
	std::shared_ptr<GenCuCSparse<Scalar>>      Grc;
	_AVMatrix<Scalar> Ave;
	_AVMatrix<Scalar> Eve;

	void assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver); //HB: add the subdomain contributions to global CCt
	void assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver, SimpleNumberer *mpcEqNums);
	void computeSubContributionToGlobalCCt(SimpleNumberer *mpcEqNums); //HB: only compute the subdomain contribution to global CCt
	void constructLocalCCtsolver();
	void deleteLocalCCtsolver() { localCCtsolver.reset(); }
	void solveLocalCCt(Scalar *subv);
	void assembleBlockCCtsolver(int iBlock, GenSolver<Scalar> *CCtsolver, SimpleNumberer *blockMpcEqNums);
	void assembleLocalCCtsolver();
	void setCCtCommSize(FSCommPattern<Scalar> *cctPat);
	void sendNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, const Connectivity *mpcToSub);
	void recNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, const Connectivity *mpcToSub);
	void factorLocalCCtsolver();
	void zeroLocalCCtsolver();

	void extractMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums);
	void insertMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums);

public:
	std::vector<Scalar> rbms;
	std::vector<Scalar> interfaceRBMs;
	// MPC related data
	mutable std::vector<std::unique_ptr<SubLMPCons<Scalar>>> mpc; // multiple point constraints
	std::vector<std::unique_ptr<SubLMPCons<Scalar>>> mpc_primal;
	/// \brief Augmentation matrix.
	std::unique_ptr<GenSparseSet<Scalar>> Src;
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> BKrrKrc;
	std::unique_ptr<GenDBSparseMatrix<Scalar>> Kbb;    //!< Boundary to boundary stiffness matrix.

	GenVector<Scalar> diagCCt;

	std::unique_ptr<GenFullM<Scalar>> qtkq;
	std::vector<Scalar> QtKpBt;
protected:
	std::unique_ptr<GenSolver<Scalar>> localCCtsolver;
	GenSparseMatrix<Scalar> *localCCtsparse; // Alias to localCCtsolver.
	int lengthCCtData = 0;
	std::vector<int> CCtrow, CCtcol;
	std::vector<Scalar> CCtval;
protected:
	GenVector<Scalar> scaling;
	mutable GenVector<Scalar> deltaFmpc;
	GenVector<Scalar> deltaFwi;
	GenVector<Scalar> wweight;
	GenVector<Scalar> kweight; //!< stiffness weights (i.e. sum of Kii for all subd sharing that dof)

	// templated RBMs
	GenFullM<Scalar> Rstar;
	GenFullM<Scalar> Rstar_g;
	std::unique_ptr<GenFullM<Scalar>> sharedRstar_g;
	std::unique_ptr<GenFullM<Scalar>> tmpRstar_g;
	std::vector<std::unique_ptr<GenFullM<Scalar>>> G;
	std::vector<std::unique_ptr<GenFullM<Scalar>>> neighbG;

	int bodyRBMoffset = 0;
	std::unique_ptr<Rbm> rigidBodyModes;

	mutable std::vector<Scalar> fcstar; // TODO Move this out!

	std::unique_ptr<GenFsiSparse<Scalar>> neighbKww;
	mutable std::vector<Scalar> localw_copy;
	// coupled_dph
	std::unique_ptr<GenDBSparseMatrix<Scalar>> Kww;
	std::unique_ptr<GenCuCSparse<Scalar>>      Kcw;
	std::unique_ptr<GenMpcSparse<Scalar>> Kcw_mpc;
	std::unique_ptr<GenCuCSparse<Scalar>> Krw;
	mutable std::vector<Scalar> localw;

	Eigen::SparseMatrix<double> B, Bw;
	Eigen::SparseMatrix<Scalar> Bm, Bc;

	double tolsvd = 1e-10; // TODO get the value from the input.
};

#endif //FEM_FETUSUB_H
