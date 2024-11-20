#ifndef _SOLVER_CNTL_
#define _SOLVER_CNTL_
#include <Feti.d/FetiInfo.h>
#include <map>

enum class SolverSelection {
	Direct = 0,
	Iterative = 1,
	Feti = 2,
	BlockDiag = 3,
	FetiLib = 4
};

inline
bool isFeti(SolverSelection selection) {
	return selection == SolverSelection::Feti || selection == SolverSelection::FetiLib;
}

struct SolverCntl {
public:
	SolverCntl() {
		mumps_icntl[3] = 0; // supress diagnostic output
		mumps_mineq = 0; // 053014 JAT
		mumps_stride = 1; // 040715 JAT
		goldfarb_tol = 1.0;
		goldfarb_check = false;
		precond = 0;
		tol = 1.0e-8;
		maxit = 1000;
		iterType = 0;
		iterSubtype = 3;
		maxvecsize = 0;
		printMatLab        = false;
		printMatLabFile    = "";
		verbose = 1;
		ilu_droptol = 1e-11;
	}
	SolverSelection type = SolverSelection::Direct;  //!< Selected solver type
	int subtype = 1;  //!< Solver subtype 1 is direct sparse ... 9 is mumps  10 is diag
	double trbm = 1e-16;         //!< algebraic rbm tolerance
	double trbm2 = 1e-16;        //!< algebraic rbm tolerance used for sparse/skyline when GRBM is activated
	int sparse_renum = 1;  //!< renumbering scheme for BLKSparseMatrix: 0 = esmond MMD, 1 = metis ND (default)
	int sparse_maxsup = 100;
	int sparse_defblk = 30;
	bool pivot = false;  //!< true if pivoting is to be used in spooles/mumps solvers
	bool unsymmetric = false;
	bool scaled = false; // true if scaling is to be used in skyline solver
	int spooles_scale = 0; // true if scaling is to be used in spooles solver
	/// \brief Used when pivoting is enabled, all entries in L and U have magnitude less or equal to tau.
	double spooles_tau = 100;
	int    spooles_seed = 532196; //!< see Solvers.d/Spooles.C for description
	int    spooles_maxsize = 64;  //!< see Solvers.d/Spooles.C for description
	int    spooles_maxdomainsize = 24;  //!< see Solvers.d/Spooles.C for description
	double spooles_maxzeros = 0.04;  //!< see Solvers.d/Spooles.C for description
	int    spooles_msglvl = 0;  //!< see Solvers.d/Spooles.C for description
	int    spooles_renum = 0;  //!< see Solvers.d/Spooles.C for description
	std::map<int, int> mumps_icntl;
	std::map<int, double> mumps_cntl;
	int mumps_mineq; // 053014 JAT
	int mumps_stride; // 040715 JAT
	double goldfarb_tol;
	bool goldfarb_check;
	int iterType; // 0 = CG, 1 = GMRES, 2 = GCR, 4 = BCG, 5 = CR
	int iterSubtype; // matrix storage
	int precond;  // preconditioner 0 = none, 1 = jacobi
	int maxit;    // maximum number of iterations
	double tol;   //!< \brief Tolerance for convergence
	int maxvecsize;  // for pcg # of krylov vectors to store default = 0
	FetiInfo fetiInfo;
	bool printMatLab;
	const char * printMatLabFile;
	int verbose;
	double ilu_droptol;
};

#endif
