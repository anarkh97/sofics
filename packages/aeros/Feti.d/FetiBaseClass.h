//
// Created by Michel Lesoinne on 1/18/18.
//

#ifndef FEM_FETIBASECLASS_H
#define FEM_FETIBASECLASS_H

#include <iostream>
#include <vector>

#include <Feti.d/DistrVector.h>
#include <Timers.d/Timing.h>
#include <Utils.d/dofset.h>
#include <Solvers.d/ParallelSolver.h>

template <class Scalar> class GenFetiOp;
template <class Scalar> class GenFetiOpControler;
template <class Scalar> class GenGMRESOrthoSet;
template <class Scalar> class GenGCROrthoSet;
template <class Scalar> class GenCGOrthoSet;
template <class Scalar> class GenFetiWorkSpace;
template <class Scalar> class FetiSub;
class FetiInfo;

template <class Type> class FSCommPattern;

template <typename Scalar>
class FetiSub;


// Class FetiWorkSpace is used to allocate and
// organize Distributed Vectors for FETI.
template<class Scalar>
class GenFetiWorkSpace
{
	GenDistrVector<Scalar> *r;
	GenDistrVector<Scalar> *lambda;
	GenDistrVector<Scalar> *w;
	GenDistrVector<Scalar> *y;
	GenDistrVector<Scalar> *z;
	GenDistrVector<Scalar> *p;
	GenDistrVector<Scalar> *pr;
	GenDistrVector<Scalar> *Fp;
	GenDistrVector<Scalar> *du;
	GenDistrVector<Scalar> *uzero;
	GenDistrVector<Scalar> *fr;
	GenDistrVector<Scalar> *fr2;
	GenDistrVector<Scalar> *ur;
	GenDistrVector<Scalar> *fw;
	GenDistrVector<Scalar> *wrk1;
	GenDistrVector<Scalar> *wrk2;
	GenDistrVector<Scalar> *deltaU;
	GenDistrVector<Scalar> *deltaF;

	// Extra Distributed Vectors for Nonlinear (not used for feti-dp)
	GenDistrVector<Scalar> *zz;
	GenDistrVector<Scalar> *rCompare;
	GenDistrVector<Scalar> *uu;

	// Added Vectors
	GenVector<Scalar> *alpha;
	GenVector<Scalar> *beta;
	GenVector<Scalar> *gamma;
	GenVector<Scalar> *working;
	GenVector<Scalar> *fc;
	GenVector<Scalar> *uc;
	GenVector<Scalar> *duc;
	GenVector<Scalar> *e;

	// Contact
	GenDistrVector<Scalar> *lambda_copy;
	GenDistrVector<Scalar> *p_copy;
	GenDistrVector<Scalar> *r_copy;
	GenDistrVector<Scalar> *Fp_copy;
	GenDistrVector<Scalar> *du_copy;
	GenVector<Scalar> *uc_copy;
	GenVector<Scalar> *duc_copy;
	GenDistrVector<Scalar> *gc;
	GenDistrVector<Scalar> *gf;

public:
	GenFetiWorkSpace(DistrInfo& interface, DistrInfo& local, int isNonlinear,
	                 int numrbms, int numcrns);
	GenFetiWorkSpace(DistrInfo& interface, DistrInfo& local, DistrInfo& wet, int ngrbms, int numC, bool contact);
	~GenFetiWorkSpace();

	GenDistrVector<Scalar>& ret_r()       { return *r; }
	GenDistrVector<Scalar>& ret_lambda()  { return *lambda; }
	GenDistrVector<Scalar>& ret_w()       { return *w; }
	GenDistrVector<Scalar>& ret_y()       { return *y; }
	GenDistrVector<Scalar>& ret_z()       { return *z; }
	GenDistrVector<Scalar>& ret_p()       { return *p; }
	GenDistrVector<Scalar>& ret_pr()      { return *pr; }
	GenDistrVector<Scalar>& ret_Fp()      { return *Fp; }
	//GenDistrVector<Scalar>& ret_Fr()      { return *Fr; }
	GenDistrVector<Scalar>& ret_du()      { return *du; }
	GenDistrVector<Scalar>& ret_uzero()   { return *uzero; }
	GenDistrVector<Scalar>& ret_wrk1()    { return *wrk1; }
	GenDistrVector<Scalar>& ret_wrk2()    { return *wrk2; }
	GenDistrVector<Scalar>& ret_deltaU()  { return *deltaU; }
	GenDistrVector<Scalar>& ret_deltaF()  { return *deltaF; }
	GenDistrVector<Scalar>& ret_fr()      { return *fr;  }
	GenDistrVector<Scalar>& ret_fr2()     { return *fr2; }
	GenDistrVector<Scalar>& ret_ur()      { return *ur;  }
	GenDistrVector<Scalar>& ret_fw()      { return *fw;  }
	GenDistrVector<Scalar>& ret_zz()      { return *zz; }
	GenDistrVector<Scalar>& ret_rCompare(){ return *rCompare; }
	GenDistrVector<Scalar>& ret_uu()      { return *uu; }
	GenDistrVector<Scalar>& ret_lambda_copy() { return *lambda_copy; }
	GenDistrVector<Scalar>& ret_p_copy()  { return *p_copy; }
	GenDistrVector<Scalar>& ret_r_copy()  { return *r_copy; }
	GenDistrVector<Scalar>& ret_gc()      { return *gc; }
	GenDistrVector<Scalar>& ret_gf()      { return *gf; }
	GenVector<Scalar>& ret_alpha()        { return *alpha; }
	GenVector<Scalar>& ret_beta()         { return *beta;  }
	GenVector<Scalar>& ret_gamma()        { return *gamma; }
	GenVector<Scalar>& ret_working()      { return *working; }
	GenVector<Scalar>& ret_fc()           { return *fc;  }
	GenVector<Scalar>& ret_uc()           { return *uc; }
	GenVector<Scalar>& ret_duc()          { return *duc; }
	GenVector<Scalar>& ret_e()            { return *e; }
	// extra functions for contact
	void save();
	void restore();
	void clean_up();
	void zeroPointers();
};

template <typename Scalar>
class FetiBaseClass : public GenParallelSolver<Scalar> {
protected:
	FetiBaseClass(std::vector<FetiSub<Scalar> *> subdomains, int numThreads, bool verboseFlag);
	void resetOrthoSet();
	/** \brief Create a guess for lambda0.
	 *
	 * @param r RHS for which a guess is desired.
	 * @param lambda0 The best guess lambda.
	 * @return Whether a prediction was made.
	 */
	bool predict(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0) const;
	void scatterHalfInterface(int iSub, Scalar *v1, GenDistrVector<Scalar> *v2) const;
	void gatherHalfInterface(int iSub, const GenDistrVector<Scalar> *v1,
	                         const GenDistrVector<Scalar> *v2,
	                         Scalar *v3, Scalar *v4) const;
	void rebuildInterface(int iSub, GenDistrVector<Scalar> &v) const;

	int nlPreCondition(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Pr) const;
	void orthogonalize(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) const;
	void orthoAddCG(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp) const;
	void setAndStoreInfo(int iter, double finalPrimal2, double finalDual2 ) const;
	// GMRES functions
	void initGMRES(GenDistrVector<Scalar> &p);
	double orthoAddGMRES(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
	void GMRESSolution(GenDistrVector<Scalar> &p);
	// GCR functions
	int predictGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0);
	void orthogonalizeGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fr,
	                      GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
	void orthoAddGCR(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar FpFp);
	void fSplit(int iSub, GenDistrVector<Scalar> &force) const;
	void distributeForce(GenDistrVector<Scalar> &force) const;
	void distributeForce(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw) const;
	void makeRbmPat();
	double preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag = true) const;
	int halfOffset(int iSub) const;
	void multKbb(int iSub, const GenDistrVector<Scalar>& v, GenDistrVector<Scalar> &interfvec,
	             GenDistrVector<Scalar>& deltaU, GenDistrVector<Scalar> &deltaF, bool &errorFlag) const;
	void sendDeltaF(int iSub, GenDistrVector<Scalar>& deltaF) const;
	void normDeltaF(int iSub, double * subDots, GenDistrVector<Scalar>* deltaF) const;
public:
	DistrInfo &interfInfo() { return interface; }
	DistrInfo &localInfo()  { return internalDI; }
	void sendScale(int isub);
	void interfaceDiff(int iSub, GenDistrVector<Scalar> &v) const;
	void interfSend(int iSub, GenDistrVector<Scalar> &dv1) const;

	void reSolve(GenDistrVector<Scalar> &) override;
	Timings& getTimers() override { return times; }
	double getSolutionTime() const override { return times.solve; }
protected:
	std::vector<FetiSub<Scalar> *> subdomains;
	int nsub = 0, glNumSub = 0;
	FetiInfo *fetiInfo = nullptr;
	int verboseFlag;
	FSCommunicator *fetiCom = nullptr;
	int myCPU, numCPUs;
	const Connectivity *subToSub = nullptr, *mpcToSub = nullptr, *mpcToSub_primal = nullptr;
	Connectivity *edgeToSub = nullptr, *subToEdge = nullptr;
	Connectivity *coarseConnect = nullptr;  // first level coarse prob. connectivity
	compStruct *renum = nullptr;
	compStruct renumber;
	SimpleNumberer *eqNums = nullptr;
	SimpleNumberer *PFcNums = nullptr;
	int gOffset;
	int mOffset;
	// TODO Find ownership (and of many other variables)
	const Connectivity *cpuToSub = nullptr;
	int glNumMpc = 0, glNumMpc_primal = 0;
	int *glSubToLoc = nullptr;
	double epsilon2;
	int maxiter;
	FSCommPattern<Scalar> *vPat = nullptr;
	FSCommPattern<Scalar> *rbmPat = nullptr;
	FSCommPattern<int> *sPat = nullptr;
	FSCommPattern<Scalar> *wiPat = nullptr;
	GenFetiOp<Scalar> **fetiOps = nullptr;
	GenFetiWorkSpace<Scalar> *wksp = nullptr;
	int halfSize = 0; // number of  rbms, half interface size
	DistrInfo internalDI, interface;

	mutable GenCGOrthoSet<Scalar> *oSetCG = nullptr; // Workspace
	mutable GenGMRESOrthoSet<Scalar> *oSetGMRES = nullptr;
	mutable GenGCROrthoSet<Scalar> *oSetGCR = nullptr;
	mutable GenFetiOpControler<Scalar> *opControl = nullptr;
	mutable int numSystems = 0; // Nonlinear additions
	mutable Timings times;
};

template<typename Scalar>
FetiBaseClass<Scalar>::FetiBaseClass(std::vector<FetiSub<Scalar> *> subdomains, int numThreads, bool verboseFlag) :
		nsub(subdomains.size()), internalDI(nsub),
		interface(nsub), times(numThreads,nsub), verboseFlag(verboseFlag)
{
	this->subdomains = std::move(subdomains);

}


#endif //FEM_FETIBASECLASS_H
