#ifndef _FETI_INFO_H_
#define _FETI_INFO_H_

// maxit    = maximum number of FETI iterations
// tol      = global tolerance
// maxortho = maximum number of reorthogonalization directions

// Preconditioner types:
// precno     = 0       no preconditioner
// precno     = 1       lumped
// precno     = 2       dirichlet               (default)

// Solver type and parameters for coarse problem: FETI-1: GtG, FETI-DP: Kcc^* 
// coarse_cntl

// Solver type and parameters for subdomain local problem
// local_cntl

// Projector types:
// nonLocalQ = 0	basic projector		(default)
// nonLocalQ = 1	Q(preconditioner type)  projector

// Scaling types:
// scaling = 1		kscaling (stiffness)
// scaling = 2		tscaling (topology)     (default)

// FETI Version
// version = 0		FETI 1			(default)
// version = 1		FETI 2

// Whether to solve coarse problem or not
// only valid for dynamics problems
// noCoarse = 0		coarse problem used	(default)
// noCoarse = 1		coarse problem not used

// Tolerance for coarse problem in FETI-1 and aux coarse problem in FETI-DP
// grbm_tol = 1.0E-6 				(default)

// Number of iterations to print the error
// printNumber = -1     no printing
// printNumber =  1	every iteration
// printNumber =  5     every 5 iterations
// printNumber =  n     every n iterations

// contactPrintFlag = 2         print status change and also other info
// contactPrintFlag = 1   	print status change info for FETI-DPC (x+- etc)
// contactPrintFlag = 0		don't print status change info (default)

// Nonlinear FETI information
// nTang = rebuild Tangent Stiffness matrix every N iterations
//         
// nPrec = rebuild FETI preconditioner every N iterations 
//         NOTE: if the user specifies to rebuild FETI every
//               n iterations and the preconditioner every N
//               iterations, then only rebuild preconditioner
//               when FETI is being rebuilt.

// Nonlinear Krylov information
// nlPrecFlg  = 0       do not use Krylov preconditioning
// nlPrecFlg  = 1       use Krylov preconditioning
// nlPrecFlg  = 2       use Krylov preconditioning per load step

// We need to keep track of the number of load steps while using FETI
//
// numLoadSteps = Number of Load Steps
//              = 0 to start and is incremented at new Load step

// Which corner degrees of freedom to clamp.
// corners = allCorners6    clamp all corner dofs (default)
//         = allCorners3
//         = noEndCorners3
//         = noEndCorners6

// FSI corner types: 
// fsi_corner = 0;    no corners on fluid/structure wet interface
// fsi_corner = 1;    all fluid wet subdomain interface nodes are corners 
// fsi_corner = 2;    all fluid and structure wet subdomain interface nodes are corners 
// fsi_corner = 3;    all fluid and a few structure wet subdomain interface nodes are corners

// f_projector = 0	don't project the force (default)
//	       = 1	use averaging projector for force
// 	       = 2	use trimming projector for force
// e_projector = 0	don't project estar (default)
//	       = 1	use averaging projector for estar
// u_projector = 0	don't project u
//	       = 1	project u (default)

// MPC Preconditioner types: (for Rixen method)
// mpc_precno     = 0       no preconditioner
// mpc_precno     = 1       diagonal                  
// mpc_precno     = 2       global 
// mpc_precno     = 3       topo block diagonal 
// mpc_precno     = 4       sub block diagonal
// mpc_precno     = 5       mortar edge block diagonal
// mpc_precno     = 6       auto select (default)

// MPC type
// mpcflag        = 0	    ignore mpcs
// mpcflag	  = 1  	    "dual" Rixen method (mpc lagrange multipliers) (default)
// mpcflag	  = 2	    "primal" include mpcs in coarse problem
// mpcflag        = 3       mixed dual and primal

// MPC Scaling types:
// mpc_scaling = 1          kscaling (stiffness) 
// mpc_scaling = 2          tscaling (topology) (default)

// Solver type and parameters for CC^t matrices
// cct_cntl

//HB: overlap level for mortar block CC^t approximate solve 
// mpcBlkOverlap = 0 (default)

// For FETI-H,
// - Construction of the mass interface matrix : lumpedinterface
// - Number of directions for the coarse grid : numcgm

// For FETI-DPH
// outerloop = 0                    use CG solver
// outerloop = 1                    use GMRES solver
// outerloop = 2                    use GCR solver
// outerloop = 3                    use CGAL solver

// numdir    = 0 (default)          number of wave directions added in Q matrix
// orthotol  = 1.0E-02 (default)    relative tolerance value in orthogonalizing Q matrix
// orthotol2 = 0.0     (default)    absolute tolerance value in orthogonalizing Q matrix

class SolverCntl;
extern SolverCntl default_cntl;

struct FetiInfo {

    // Data members
    int    maxit = 1000; //!< \brief Maximum number of FETI iterations.
    double tol = 1.0E-6; //!< \brief Relative tolerance.
    double absolute_tol = 0.0; //!< \brief Absolute tolerance.
    double stagnation_tol = 1.0E-6; //!< \brief Relative stagnation tolerance.
    double absolute_stagnation_tol; //!< \brief Absolute stagnation tolerance.
    double grbm_tol = 1.0E-6;
    double crbm_tol = 1.0E-12;
    double cct_tol = 1.0E-16;
    int    uproj = 1;
    int    maxortho = maxit; // \brief Maximum number of reorthogonalization directions.
    int    noCoarse = 0;
    int    nonLocalQ = 0;
    int    nQ = 0;
    int    primalFlag = 0; // whether to output primal residual
    bool    printMatLab = false;

    int    printNumber = 10;

    // Nonlinear Data members
    int    nPrec = 1;
    int    nTang = 1;
    int    nlPrecFlg = 0;
    int    numLoadSteps = 0;

    enum Preconditioner { noPrec=0, lumped, dirichlet, identity } precno = dirichlet;
    enum PreconditionerType { nonshifted, shifted } prectype = nonshifted;
	enum MpcPreconditioner { noMpcPrec=0, diagCCt, globalCCt, blockDiagCCt,
		subBlockDiagCCt, superBlockDiagCCt, autoSelectCCt } mpc_precno = globalCCt;
    enum MpcBlock { subBlock, topoBlock, mortarBlock } mpc_block = topoBlock;
    int mpcflag = 1;
    enum Solvertype { skyline, sparse, blocksky, llt, ldlt, cholmod,
                      umfpack, superlu, spooles, mumps, diagonal, dbsgal, eisgal, goldfarb, splu, ssqr, spqr };
    SolverCntl *local_cntl = &default_cntl;
    SolverCntl *coarse_cntl = &default_cntl;
    SolverCntl *auxcoarse_cntl = &default_cntl;
    SolverCntl *cct_cntl = &default_cntl;
    SolverCntl *kii_cntl = &default_cntl;
    enum Scaling { noscaling=0, kscaling=1, tscaling=2 };
    Scaling scaling = Scaling::tscaling,
		    mpc_scaling = Scaling::tscaling,
		    fsi_scaling = Scaling::tscaling;
    enum Version { feti1, feti2, feti3, fetidp } version = Version::fetidp;
    /** \details if this is true then reassemble and apply scaling to f for every system, not just the first
     * affects the convergence criteria for nonlinear and dynamics since relative primal error will be
     * defined using norm of rescaled f which is typically much smaller than the norm of the unscaled f */
    bool rescalef = true;
	enum Feti2Version { fullCoarse, sparseCoarse } feti2version = sparseCoarse;
    enum Type { linear, nonlinear, eigen } type = linear;
    enum CornerType { allCorners6, allCorners3, noEndCorners6,
                      noEndCorners3, interface3, interface6, ThreeD, noCorners } corners = noEndCorners3;
    enum AugmentType { none, Gs, Edges, WeightedEdges } augment = Edges;
    enum AugmentImplementation { Constraint, Primal } augmentimpl = Constraint;
    int isEdgeAugmentationOn() const {
      return (((augment == Edges) || (augment == WeightedEdges)) && ((nGs > 0) || (numdir > 0))) ? 1 : 0; 
    }
    enum RbmType { translation, rotation, all,
                   averageTran, averageRot, averageAll, None,
                   pressure, temperature } rbmType = all;
    double nullSpaceFilterTol = 0.0;

    // FETI-H
    double tolcgm = 1e-3;
    int numcgm = 0; // number of coarse grid modes
    int spaceDimension = 3;
    int krylovtype = 8;
    int lumpedinterface = 0; // 0 - default (consistent) 1 - lumped
    int saveMemCoarse = 0;

    int nGs = 6;

    int maxiter()      { return maxit;       }
    int maxorth()      { return maxortho;    }
    int nPrecond()     { return nPrec;       } 
    double tolerance() { return tol;         }
    int numPrint()     { return printNumber; }

    enum OuterloopType { CG, GMRES, GCR, CGAL }  outerloop = CG;
    enum WaveType      { solid, shell, fluid, any } waveType = any;
    /** \details  note: only use uniform if domain is homogeneous
     * and entirely solid or shell or fluid */
    enum WaveMethod    { averageK, averageMat, uniform } waveMethod = averageMat;
    int numdir = 0;
    double orthotol = 1.0E-2;
    double orthotol2 = 0.0;
    bool dph_flag = false;
    int contactPrintFlag = 0;

    bool rebuildcct = true;
    bool geometric_gap = false;
    int mpcBlkOverlap = 0; //0=no interaction, 1=1st order interactions, 2=1st & 2nd order interactions, etc.
    double gamma = 1.0;
    bool bmpc = false;
    bool cmpc = false;
    double linesearch_tau = 2.0/3.0;
    int linesearch_maxit = 100;
    bool c_normalize = false;

    bool useMRHS = true;
    bool gmresResidual = false; //HB: to force computing the "primal residual" at each GMRES iteration;
    bool wetcorners = false;
    bool splitLocalFsi = true;
    bool pickAnyCorner = true;
    bool pick_unsafe_corners = true; // if this is true (default) unsafe nodes will be selected as corners
                              // should only be set to false if pivoting is enabled for local solver

    bool fsi_element = false;
    bool mpc_element = false; // true means add to element set & decompose
    int fsi_corner = 2;
    bool complex_hermitian = false;
    double dual_proj_tol = 0.0;
    double primal_proj_tol = 0.0;
    double dual_plan_tol = 0.0;
    double primal_plan_tol = 0.0;
    int dual_plan_maxit = 20;
    int primal_plan_maxit = 20;
};

#endif
