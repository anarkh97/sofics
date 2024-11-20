#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <iostream>
#include <cassert>
#include <set>
#include <map>

#include <Utils.d/resize_array.h>
#include <Utils.d/SolverInfo.h>
#include <Utils.d/CompositeInfo.h>
#include <Utils.d/MyComplex.h>
#include <Timers.d/Timing.h>
#include <Timers.d/MatrixTimers.h>
#include <Driver.d/HData.h>
#include <Parser.d/AuxDefs.h>
#include <Sfem.d/Sfem.h>
#include <Utils.d/OutputInfo.h>
#include <Utils.d/SensitivityInfo.h>
#include <Utils.d/ShapeSensitivityData.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Rom.d/GalerkinProjectionSolver.h>
#include <Rom.d/EiGalerkinProjectionSolver.h>

class MortarHandler;
class MFTTData;
class SS2DTData;
class ControlInterface;
class ControlInfo;
class DofSetArray;
class ConstrainedDSA;
class DOFMap;
template <class Scalar> class GenNBSparseMatrix;
typedef GenNBSparseMatrix<double> NBSparseMatrix;
template <class Scalar> class GenDBSparseMatrix;
typedef GenDBSparseMatrix<double> DBSparseMatrix;
typedef GenDBSparseMatrix<DComplex> DBComplexSparseMatrix;
template <typename Scalar, typename SolverClass> class GenEiSparseMatrix;
class Connectivity;
template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<double> CuCSparse;
template <class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
template <class Scalar> class GenDynamMat;
typedef GenDynamMat<double> DynamMat;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenDecDomain;
template <class Scalar> class GenDistrDomain;
template <class Scalar> class GenSandiaDomain;
class FlExchanger;
class Rbm;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
typedef GenSolver<DComplex> ComplexSolver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
typedef GenSparseMatrix<DComplex> ComplexSparseMatrix;
template <class Scalar, class AnyVector> class KrylovProjector;
template <class AnyVector> class Preconditioner;
class GeomState;
class DistrGeomState;
template <class Scalar> class DistrBlockVector;
class IntFullM;
class ControlLawInfo;
class StaticTimers;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class FSCommunicator;
struct ModeData;

class SurfaceEntity;

extern Sfem *sfem;

// ... Structure used to store problem Operators buildSkyOps
// ... i.e. Only a Solver is needed for a static problem
template<class Scalar>
struct AllOps
{
	GenSolver<Scalar> *sysSolver;  // system solver: to solve (coeM*M+coeC*C+coeK*K)x = b
	GenSparseMatrix<Scalar> *spm;
	GenSolver<Scalar> *prec;       // preconditioner
	GenSparseMatrix<Scalar> *spp;

	GenSparseMatrix<Scalar> *Msolver;  // for assembling mass solver: to solve Mx = b
	GenSparseMatrix<Scalar> *K;    // stiffness matrix
	GenSparseMatrix<Scalar> *M;    // mass matrix
	GenSparseMatrix<Scalar> *C;	 // damping matrix
	GenSparseMatrix<Scalar> *Kuc;	 // constrained to unconstrained stiffness
	GenSparseMatrix<Scalar> *Muc;	 // constrained to unconstrained mass matrix
	GenSparseMatrix<Scalar> *Cuc;	 // constrained to unconstrained damping matrix
	GenSparseMatrix<Scalar> *Kcc;  // constrained to constrained stiffness matrix
	GenSparseMatrix<Scalar> *Mcc;	 // constrained to constrained mass matrix
	GenSparseMatrix<Scalar> *Ccc;  // constrained to constrained damping matrix
	GenSparseMatrix<Scalar> **C_deriv;    // derivatives of damping matrix for higher order sommerfeld
	GenSparseMatrix<Scalar> **Cuc_deriv;    // derivatives of constrained to unconstrained damping matrix for higher order sommerfeld
	int n_Kderiv;
	GenSparseMatrix<Scalar> **K_deriv;    // derivatives of K for rubber or non-linear with frequency
	GenSparseMatrix<Scalar> **Kuc_deriv;    // derivatives of Kuc for rubber or non-linear with frequency
	int num_K_arubber;
	GenSparseMatrix<Scalar> **K_arubber_l;    // lambda part of K for rubber materials
	GenSparseMatrix<Scalar> **K_arubber_m;    // mu part of K for rubber materials
	GenSparseMatrix<Scalar> **Kuc_arubber_l;    // lambda part of Kuc for rubber materials
	GenSparseMatrix<Scalar> **Kuc_arubber_m;    // mu part of Kuc for rubber materials

	GenVector<Scalar> *rhs_inpc;
	// Constructor
	AllOps() { sysSolver = 0; spm = 0; prec = 0; spp = 0;
		Msolver = 0; K = 0; M = 0; C = 0; Kuc = 0; Muc = 0;
		Cuc = 0; Kcc = 0; Mcc = 0; Ccc = 0; C_deriv = 0;
		Cuc_deriv = 0; K_deriv = 0; Kuc_deriv = 0; K_arubber_l = 0;
		K_arubber_m = 0; Kuc_arubber_l = 0; Kuc_arubber_m = 0; rhs_inpc = 0; }

	void zero() {if(K) K->zeroAll();
		if(M) M->zeroAll();
		if(C) C->zeroAll();
		if(Kuc) Kuc->zeroAll();
		if(Muc) Muc->zeroAll();
		if(Cuc) Cuc->zeroAll();
		if(Kcc) Kcc->zeroAll();
		if(Mcc) Mcc->zeroAll();
		if(Ccc) Ccc->zeroAll();
// RT: 053113 : not finished
		if (C_deriv) if (C_deriv[0]) C_deriv[0]->zeroAll();
		if (Cuc_deriv) if (Cuc_deriv[0]) Cuc_deriv[0]->zeroAll();
		if (K_deriv) for(int i=0;i<n_Kderiv;i++)
				if (K_deriv[i]) K_deriv[i]->zeroAll();
		if (Kuc_deriv) for(int i=0;i<n_Kderiv;i++)
				if (Kuc_deriv[i]) Kuc_deriv[i]->zeroAll();
		if (K_arubber_l) for(int i=0;i<num_K_arubber;i++)
				if (K_arubber_l[i]) K_arubber_l[i]->zeroAll();
		if (K_arubber_m) for(int i=0;i<num_K_arubber;i++)
				if (K_arubber_m[i]) K_arubber_m[i]->zeroAll();
		if (Kuc_arubber_l) for(int i=0;i<num_K_arubber;i++)
				if (Kuc_arubber_l[i]) Kuc_arubber_l[i]->zeroAll();
		if (Kuc_arubber_m) for(int i=0;i<num_K_arubber;i++)
				if (Kuc_arubber_m[i]) Kuc_arubber_m[i]->zeroAll();
	}
};

template<class Scalar>
struct AllSensitivities
{
#ifdef USE_EIGEN3
	double weight;           // total weight of the structure
	Eigen::Matrix<double, Eigen::Dynamic, 1> *dwrAggregatedStressVM;
	Eigen::Matrix<double, Eigen::Dynamic, 1> *dwrStressVM;
	Eigen::Matrix<double, Eigen::Dynamic, 1> *dwrDisp;
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *residual; // residual of the forward problem
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *weightWRTthick;    // derivatives of weight wrt thickness
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *weightWRTshape;    // derivatives of weight wrt shape variables
	GenSparseMatrix<Scalar> *vonMisesWRTthickSparse;             // derivatives of von Mises stress wrt thickness
	GenSparseMatrix<Scalar> *vonMisesWRTdispSparse;       // derivatives of von Mises stress wrt displacement
	GenSparseMatrix<Scalar> *vonMisesWRTshapeSparse;      // derivatives of von Mises stress wrt shape varibales
	GenSparseMatrix<Scalar> **stiffnessWRTthickSparse;    // derivatives of stiffness wrt thickness
	GenSparseMatrix<Scalar> **stiffnessWRTshapeSparse;    // derivatives of stiffness wrt shape variables
	GenSparseMatrix<Scalar> **dKucdthickSparse;           // derivatives of constrained stiffness wrt thickness
	GenSparseMatrix<Scalar> **dKucdshapeSparse;           // derivatives of constrained stiffness wrt shape variables
	GenSparseMatrix<Scalar> **linearstaticWRTthickSparse; // derivative of linear static structural formulation wrt thickness
	GenSparseMatrix<Scalar> **linearstaticWRTshapeSparse; // derivative of linear static structural formulation wrt shape variables
	GenSparseMatrix<Scalar> **dispWRTthickSparse;         // derivative of displacement wrt thickness
	GenSparseMatrix<Scalar> **dispWRTshapeSparse;         // derivative of displacement wrt shape variables

	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vonMisesWRTthick;      // derivatives of von Mises stress wrt thickness
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vonMisesWRTdisp;       // derivatives of von Mises stress wrt displacement
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vonMisesWRTshape;      // derivatives of von Mises stress wrt shape varibales
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vonMisesWRTmach;       // derivatives of von Mises stress wrt Mach number
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vonMisesWRTalpha;      // derivatives of von Mises stress wrt angle of attack
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vonMisesWRTbeta;       // derivatives of von Mises stress wrt yaw angle
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **dKucdthick;           // derivatives of constrained stiffness wrt thickness
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **dKucdshape;           // derivatives of constrained stiffness wrt shape variables
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *stressWeight;          // weight used to average stress sensitivity
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **linearstaticWRTthick; // derivative of linear static structural formulation wrt thickness
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **linearstaticWRTshape; // derivative of linear static structural formulation wrt shape variables
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *aggregatedVonMisesWRTthick;         // derivative of KS function of von Mises stresses wrt thickness
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *aggregatedVonMisesWRTshape;         // derivative of KS function of von Mises stresses wrt shape variable
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *aggregatedVonMisesWRTdisp;          // derivative of KS function of von Mises stresses wrt displacement
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *lambdaFluidQuantity;                // dual sensitivity of fluid quantities (e.g., lift, drag...)
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *lambdaAggregatedStressVM;           // dual sensitivity of KS function of von Mises stresses
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> **lambdaStressVM;                    // dual sensitivity of von Mises stress at a specified node
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> **lambdaDisp;                        // dual sensitivity of displacement at a specified node
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **gdispWRTthick;        // derivative of global displacement wrt thickness
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **gdispWRTshape;        // derivative of global displacement wrt shape variables
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *gdispWRTmach;          // derivative of global displacement wrt Mach number
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *gdispWRTalpha;         // derivative of global displacement wrt angle of attack
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *gdispWRTbeta;          // derivative of global displacement wrt yaw angle

	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **dispWRTthick;         // derivative of displacement wrt thickness
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **dispWRTshape;         // derivative of displacement wrt shape variables
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *dispWRTmach;           // derivative of displacement wrt Mach number
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *dispWRTalpha;          // derivative of displacement wrt angle of attack
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *dispWRTbeta;           // derivative of displacement wrt yaw angle

	// Constructor
	AllSensitivities() { residual = 0; weight = 0;        weightWRTshape = 0;              weightWRTthick = 0;
		vonMisesWRTthickSparse = 0;      dKucdthickSparse = 0;            vonMisesWRTshapeSparse = 0;
		vonMisesWRTdispSparse = 0;       stiffnessWRTthickSparse = 0;     dKucdshapeSparse = 0;
		lambdaAggregatedStressVM = 0;    dwrAggregatedStressVM = 0;       linearstaticWRTthickSparse = 0;
		linearstaticWRTshapeSparse = 0;  dispWRTthickSparse = 0;          dispWRTshapeSparse = 0;
		stiffnessWRTshapeSparse = 0;     dispWRTmach = 0;                 dispWRTalpha = 0;
		dispWRTbeta = 0;                 lambdaDisp = 0;                  dwrDisp = 0;
		lambdaStressVM = 0;              dwrStressVM = 0;                 vonMisesWRTthick = 0;
		dKucdthick = 0;                  vonMisesWRTshape = 0;            vonMisesWRTalpha = 0;
		vonMisesWRTbeta = 0;             lambdaFluidQuantity = 0;         aggregatedVonMisesWRTthick = 0;
		aggregatedVonMisesWRTshape = 0;  aggregatedVonMisesWRTdisp = 0;   vonMisesWRTdisp = 0;
		stressWeight = 0;                dKucdshape = 0;                  linearstaticWRTthick = 0;
		linearstaticWRTshape = 0;        dispWRTthick = 0;                dispWRTshape = 0;
		gdispWRTthick = 0;               gdispWRTshape = 0;               gdispWRTmach = 0;
		gdispWRTalpha = 0;               gdispWRTbeta = 0; }

	void zero(int numShapeVars=0, int numThicknessGroups=0) {
		if(weightWRTthick) weightWRTthick->setZero();
		if(weightWRTshape) weightWRTshape->setZero();
		if(vonMisesWRTthick) {  vonMisesWRTthick->setZero();   vonMisesWRTthickSparse->zeroAll(); }
		if(vonMisesWRTdisp)  {  vonMisesWRTdisp->setZero();    vonMisesWRTdispSparse->zeroAll();  }
		if(vonMisesWRTshape) {  vonMisesWRTshape->setZero();   vonMisesWRTshapeSparse->zeroAll(); }
		if(vonMisesWRTalpha) vonMisesWRTalpha->setZero();
		if(vonMisesWRTalpha) vonMisesWRTalpha->setZero();
		if(aggregatedVonMisesWRTthick) aggregatedVonMisesWRTthick->setZero();
		if(aggregatedVonMisesWRTshape) aggregatedVonMisesWRTshape->setZero();
		if(aggregatedVonMisesWRTdisp) aggregatedVonMisesWRTdisp->setZero();
		if(lambdaFluidQuantity) lambdaFluidQuantity->setZero();
		if(lambdaStressVM) lambdaStressVM->setZero();
		if(lambdaDisp)     lambdaDisp->setZero();
		if(residual) residual->setZero();
		if(dwrStressVM) dwrStressVM->setZero();
		if(dwrDisp) dwrDisp->setZero();
		if(dwrAggregatedStressVM) dwrAggregatedStressVM->setZero();
		if(stressWeight)     {  stressWeight->setZero();  }
		if(dispWRTmach)      {  dispWRTmach->setZero();              }
		if(dispWRTalpha)     {  dispWRTalpha->setZero();             }
		if(dispWRTbeta)      {  dispWRTbeta->setZero();              }
		if(gdispWRTmach)     {  gdispWRTmach->setZero();              }
		if(gdispWRTalpha)    {  gdispWRTalpha->setZero();             }
		if(gdispWRTbeta)     {  gdispWRTbeta->setZero();              }
		if(stiffnessWRTthickSparse) for(int i=0; i<numThicknessGroups; ++i) {  stiffnessWRTthickSparse[i]->zeroAll();  }
		if(stiffnessWRTshapeSparse) for(int i=0; i<numShapeVars; ++i) {  stiffnessWRTshapeSparse[i]->zeroAll();   }
		if(dKucdthick) for(int i=0; i<numThicknessGroups; ++i) { dKucdthick[i]->setZero();               dKucdthickSparse[i]->zeroAll();   }
		if(dKucdshape) for(int i=0; i<numShapeVars; ++i) { dKucdshape[i]->setZero();                     dKucdshapeSparse[i]->zeroAll();  }
		if(linearstaticWRTthick) for(int i=0; i<numThicknessGroups; ++i) { linearstaticWRTthick[i]->setZero();  linearstaticWRTthickSparse[i]->zeroAll();  }
		if(linearstaticWRTshape) for(int i=0; i<numShapeVars; ++i) { linearstaticWRTshape[i]->setZero(); linearstaticWRTshapeSparse[i]->zeroAll();   }
		if(dispWRTthick) for(int i=0; i<numThicknessGroups; ++i) { dispWRTthick[i]->setZero();           dispWRTthickSparse[i]->zeroAll();  }
		if(dispWRTshape) for(int i=0; i<numShapeVars; ++i) { dispWRTshape[i]->setZero();                 dispWRTshapeSparse[i]->zeroAll();  }
		if(gdispWRTthick) for(int i=0; i<numThicknessGroups; ++i) { gdispWRTthick[i]->setZero();         }
		if(gdispWRTshape) for(int i=0; i<numShapeVars; ++i) { gdispWRTshape[i]->setZero();               }
	}
#endif
};

// Structure to store discrete masses
#include <Driver.d/DMassData.h> // TG : struct DMassData moved here

// Structure to store dof renumbering information
struct Renumber {
	int *order;
	int *renumb;
	Renumber() { order = renumb = 0; }
};

struct PrevFrc {
	int lastTIndex;
	double lastFluidTime;
	Vector lastFluidLoad;

	PrevFrc(int neq) : lastFluidLoad(neq, 0.0) { lastTIndex = -1; }
};

struct AdjacencyLists {
	std::vector<std::pair<int,std::vector<int> > > surfp;
	std::vector<std::pair<int,std::vector<int> > > cdnm;
	std::vector<std::pair<int,std::vector<int> > > cdnf;
	std::vector<std::pair<DMassData*,std::vector<int> > > dimass;
	std::set<int> crot;
};

struct DispNode {
	int nodeID;
	int dofs[6];
	int numdofs;
};

/** Class representing a structure and containing all auxiliary data-structures
 *
 */
class Domain : public HData {
protected:
	MatrixTimers *matrixTimers;// timers to time factoring and assembly
	Solver *solver;            // this domains solver
	int numnodes;		// number of nodes
	int numnodesFluid;		// number of nodes for fluid
	CoordSet &nodes;   	// All the nodes
	int numele;		// number of elements
	Elemset packedEset;	// The set of compressed elements
	int numdofs; 		// the total number of degrees of freedom
	int numSensitivity;  // the total number of sensitivity types
	std::vector<int> thicknessGroups;
	std::vector<int> stressNodes;
	std::vector<DispNode> dispNodes;

	// BC related data members
	int numDirichlet;		// number of dirichlet bc
	int numDirichletFluid;	// number of dirichlet bc in fluid
	int numDispDirichlet;      // number of displacement dirichlet bc
	BCond* dbc;		// set of those dirichlet bc
	BCond* dbcFluid;		// set of those dirichlet bc in fluid
	int numLMPC; 		// total number of Linear Multi-Point Constraints (both real & complex)
	ResizeArray<LMPCons *> lmpc;  // set of LMPCs
	int numFSI; 		// total number of Fluid-Structure interaction Constraints
	ResizeArray<LMPCons *> fsi ;  // set of FSIs
	int numCTC;                // total number of contact constraints
	int numNeuman;		// number of Neuman bc
	BCond* nbc;		// set of Neuman bc
	int numNeumanModal;
	BCond* nbcModal;
	int numIDis;		// number of Initial displacements
	BCond *iDis;		// set of those initial displacements
	int numIDisModal;
	BCond *iDisModal;
	int numIDis6;		// number of Initial displacements (6 column)
	BCond* iDis6;              // set of those intitial displacements
	int numIVel;		// number of initial velocities
	BCond *iVel;		// set of those initial velocities
	int numIVelModal;
	BCond *iVelModal;

	DofSetArray *dsa = nullptr;		// Dof set array
	DofSetArray *dsaFluid = nullptr;	// Dof set array for fluid
	ConstrainedDSA *c_dsa = nullptr;	// Constrained dof set array
	ConstrainedDSA *c_dsaFluid = nullptr;// Constrained dof set array for fluid
	ConstrainedDSA *MpcDSA = nullptr;
	ConstrainedDSA *g_dsa = nullptr;
	Connectivity *allDOFs = nullptr;     // all dof arrays for each node
	Connectivity *allDOFsFluid = nullptr;     // all dof arrays for each node
	int maxNumDOFs; 		// maximum number of dofs an element has
	int maxNumDOFsFluid; 	// maximum number of dofs a fluid element has
	int maxNumNodes;           // maximum number of nodes an element has
	int maxNumNodesFluid;      // maximum number of nodes a fluid element has
	Connectivity *elemToNode, *nodeToElem;
	std::unique_ptr<Connectivity> nodeToNode;
	Connectivity *fsiToNode, *nodeToFsi, *nodeToNodeDirect;
	Connectivity *elemToNodeFluid, *nodeToElemFluid, *nodeToNodeFluid;
	Connectivity *nodeToNode_sommer; // for higher order sommerfeld
	SolverInfo sinfo;		// solver information structure
	ControlLawInfo *claw = nullptr;      // contains user defined routine information
	DMassData *firstDiMass;	// Discrete Mass
	ShapeSensitivityData shapeSenData;
	int nDimass;               // number of DMASS
	double *gravityAcceleration;   // (gx,gy,gz)
	double gravitySloshing;    // g, added for sloshing problem
	compStruct renumb;         // renumbered nodes per structural component
	compStruct renumbFluid;    // renumbered nodes per fluid component
	compStruct renumb_nompc;
	std::map<std::pair<int,int>,double> loadfactor;
	std::map<std::pair<int,int>,int> loadfactor_mftt;
	std::map<std::pair<int,int>,int> loadfactor_hftt;
	std::map<int,bool> loadfactor_temp;
	std::map<int,bool> loadfactor_grav;
	std::map<int,MFTTData*> mftval; // Mechanics Force Time Table
	std::map<int,MFTTData*> hftval; // Heat Fluxes Time Table
	int numHFTT;                    // number of HFTTs
	ResizeArray<MFTTData *> ymtt;         // Young's modulus vs. temperature table
	int numYMTT;                          // number of YM Temp tables
	ResizeArray<MFTTData *> sdetaft;      // RT: Structural damping versus frequency table
	int numSDETAFT;                       // number of SDETAF tables
#ifdef USE_EIGEN3
	ResizeArray<GenMFTTData<Eigen::Vector4d> *> rubdaft;      // RT: Rubber damping versus frequency table
#endif
	int numRUBDAFT;                       // number of RUBDAF tables
     ResizeArray<MFTTData *> ctett;        // Coefficient of thermal expansion vs. temperature table
	int numCTETT;                         // number of CTE Temp tables
     ResizeArray<MFTTData *> ss1dt;        // Stress vs. strain 1-dimensional table
     int numSS1DT;                         // number of SS 1-dimensional tables
     ResizeArray<SS2DTData *> ss2dt;       // Stress vs. strain 2-dimensional table
     int numSS2DT;                         // number of SS 2-dimensional tables
     ResizeArray<MFTTData *> ysst;         // Yield stress vs. effective plastic strain table
	int numYSST;                          // number of YSS tables
     ResizeArray<MFTTData *> yssrt;        // Yield stress scale factor vs. effective plastic strain rate table
	int numYSSRT;                         // number of YSSRT tables
     ResizeArray<MFTTData *> ymst;         // Young's modulus vs. strain table
     int numYMST;                          // number of YMS tables
	FlExchanger *flExchanger;  // Fluid Exchanger
	FILE *outFile;

	// for computing stress results in function getStressStrain(...)
	Vector *stress;
	Vector *weight;
	Vector *elDisp;
	Vector *elTemp;
	Vector *elstress;
	Vector *elweight;
	// for computing stress results in function getPrincipalStress(...)
	FullM *p_stress;
	FullM *p_elstress;
	Vector *stressAllElems; // stores stresses of all the elements : used Sfem
	int sizeSfemStress;

	// for compute energies
	double Wext;
	double Waero;
	double Wdmp;
	double pWela;
	double pWkin;
	double pWdis;
	double modalWela;
	double modalWkin;
	Vector *previousExtForce;
	Vector *previousAeroForce;
	Vector *previousDisp;
	Vector *previousCq;

	Vector *heatflux;
	Vector *elheatflux;

	Vector *fluidDispSlosh;
	Vector *elFluidDispSlosh;
	Vector *elPotSlosh;
	Vector *elFluidDispSloshAll;

	double savedFreq;
	bool sowering;
	bool output_match_in_top;

	int numThicknessGroups;  // number of thickness groups
	int numShapeVars;        // number of shape variables
	int numSensitivityQuantityTypes; // number of sensitivity quantities
	int numStressNodes;      // number of requested nodes for von mises stress sensitivity
	int numDispNodes;        // number of requested nodes for displacement sensitivity
	int numTotalDispDofs;         // number of requested dofs at each requested node for displacement sensitivity

	double* aggregatedStress;
	double* aggregatedStressDenom;
	bool aggregatedFlag;
	int aggregatedFileNumber;

	void writeTopFileElementSets(ControlInfo *cinfo, int * nodeTable, int* nodeNumber, int topFlag);

	Elemset elems_copy; // Needed for SFEM
	Elemset elems_fullcopy; // copy of the full elementset, needed for SFEM

	std::vector<AdjacencyLists> elemAdj;
	std::vector<int> followedElemList;

	std::set<int> newDeletedElements; // list of elements that have been deleted during the current time step
	std::vector<std::pair<double,int> > outDeletedElements; // used for "elemdele" output
	Connectivity **nodeToFaceElem;

	FSCommunicator *com;
	bool *thgreleFlag;
	int *thpaIndex;

public:
	bool runSAwAnalysis; // if true, analysis will be run first then compute sensitivity
	// if false, no analysis will be run before computing sensitivity
	SensitivityInfo *senInfo;  // sensitivity information structure array
	// Implements nonlinear dynamics postprocessing for file # fileId
	void postProcessingImpl(int fileId, GeomState*, Vector&, Vector&, double, int, double *, double *,
	                        Corotator **, double *acceleration = 0, double *acx = 0, GeomState *refState = 0,
	                        Vector *reactions = 0, SparseMatrix *M = 0, SparseMatrix *C = 0);

	Domain(int iniSize = 16);
	Domain(Domain &, int nele, const int *ele, int nnodes, const int *nodes);
	Domain(Domain &, Elemset *elems, CoordSet *nodes);  // PJSA: for new sower
	virtual ~Domain();

	double getSavedFreq()  { return savedFreq; }
	void setSavedFreq(double freq)  { savedFreq = freq; }
	double *temprcvd;          // temperature received by structure from
	// heat solution
	int numContactPairs;  // used for timing file
	char * optinputfile;

	int* glWetNodeMap;

	int*  umap_add;            // mapping for coupling matrix assembly
	int*  umap;                // mapping for coupling matrix
	int** pmap;                // mapping for coupling matrix
	int   nuNonZero;
	int*  npNonZero;
	double ** C_condensed;

	// functions for controlling printing to screen
	void setVerbose() { outFile = stderr; }
	void setSilent()  { outFile = 0;      }
	void setOutputMatchInTop(bool b) {output_match_in_top = b;};
	void readInModes(int modal_id, ModeData &modeData);
	void readInShapeDerivatives(char* shapeDerFileName);
	void setSowering(bool b) { sowering = b;}
	bool getSowering() { return sowering;}
	void make_bc(int *bc, double *bcx);
	void make_bc(int *bc, DComplex *bcxC);
	void make_bc(gsl::span<int> bc, gsl::span<double> bcx) {
		make_bc(bc.data(), bcx.data());
	}
	void make_bc(gsl::span<int> bc, gsl::span<DComplex> bcx) {
		make_bc(bc.data(), bcx.data());
	}
	void makeAllDOFs();
	void makeAllDOFsFluid();
	void setNumShapeVars(int _numS) { numShapeVars = _numS; }
	void setThicknessGroup(int d) { thicknessGroups.push_back(d-1); numThicknessGroups++; }
	void setStressNodes(int d) { stressNodes.push_back(d-1); numStressNodes++; }
	void setDispNode(int n, int d1) { DispNode dn;  dn.nodeID = n;  dn.dofs[0] = d1; dn.numdofs = 1;
		dispNodes.push_back(dn); numDispNodes++; numTotalDispDofs += 1; }
	void setDispNode(int n, int d1, int d2) { DispNode dn;  dn.nodeID = n;  dn.dofs[0] = d1; dn.dofs[1] = d2;
		dn.numdofs = 2;  dispNodes.push_back(dn); numDispNodes++; numTotalDispDofs += 2;}
	void setDispNode(int n, int d1, int d2, int d3) { DispNode dn;  dn.nodeID = n;  dn.dofs[0] = d1; dn.dofs[1] = d2;
		dn.dofs[2] = d3; dn.numdofs = 3;  dispNodes.push_back(dn);
		numDispNodes++; numTotalDispDofs += 3; }
	void setDispNode(int n, int d1, int d2, int d3, int d4) { DispNode dn;  dn.nodeID = n;  dn.dofs[0] = d1; dn.dofs[1] = d2;
		dn.dofs[2] = d3; dn.dofs[3] = d4; dn.numdofs = 4;
		dispNodes.push_back(dn); numDispNodes++; numTotalDispDofs += 4; }
	void setDispNode(int n, int d1, int d2, int d3, int d4, int d5) { DispNode dn;  dn.nodeID = n;  dn.dofs[0] = d1;
		dn.dofs[1] = d2; dn.dofs[2] = d3; dn.dofs[3] = d4;
		dn.dofs[4] = d5; dn.numdofs = 5;  dispNodes.push_back(dn);
		numDispNodes++; numTotalDispDofs += 5; }
	void setDispNode(int n, int d1, int d2, int d3, int d4, int d5, int d6) { DispNode dn;  dn.nodeID = n;  dn.dofs[0] = d1;
		dn.dofs[1] = d2; dn.dofs[2] = d3; dn.dofs[3] = d4;
		dn.dofs[4] = d5; dn.dofs[5] = d6; dn.numdofs = 6;
		dispNodes.push_back(dn); numDispNodes++;
		numTotalDispDofs += 6; }
	std::vector<int> *getStressNodes() { return &stressNodes; }
	std::vector<DispNode> *getDispNodes() { return &dispNodes; }

	void setIncludeStressNodes(bool *);
	bool checkIsInStressNodes(int,int &);
	void createKelArray(FullSquareMatrix *& kel);
	void createKelArray(FullSquareMatrix *& kel,FullSquareMatrix *& mel);
	void createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray, FullSquareMatrix *&cArray);
	void getElementForces( GeomState &geomState, Corotator **allCorot,
	                       int fileNumber, int forceIndex, double time);
	void getStressStrain(GeomState &geomState, Corotator **allCorot,
	                     int fileNumber, int stressIndex, double time,
	                     GeomState *refState = NULL);
	void getPrincipalStress(GeomState &geomState, Corotator **allCorot,
	                        int fileNumber, int stressIndex, double time,
	                        GeomState *refState = NULL);
	void updateStates(GeomState *refState, GeomState& geomState, Corotator **allCorot, double time);
	void updateWeightedElemStatesOnly(const std::map<int, double> &weights, GeomState *refState,
	                                  GeomState &geomState, Corotator **corotators, double time);
	void updateWeightedElemStatesOnly(const std::set<int> &weightedElems, GeomState *refState,
	                                  GeomState &geomState, Corotator **corotators, double time);
	void getElemStiffAndForce(const GeomState &geomState, double time,
	                          const GeomState *refState, const Corotator &elemCorot,
	                          double *elemForce, FullSquareMatrix &elemStiff);
	void getElemStiffAndForce(const GeomState &geomState, double time,
	                          const Corotator &elemCorot,
	                          double *elemForce, FullSquareMatrix &elemStiff);
	void getStiffAndForce(GeomState &u, Vector &elementInternalForce,
	                      Corotator **allCorot, FullSquareMatrix *kel,
	                      Vector &residual, double lambda = 1.0, double time = 0.0,
	                      GeomState *refState = NULL, Vector *reactions = NULL,
	                      FullSquareMatrix *mel = NULL, FullSquareMatrix *cel = NULL);
	void makeElementAdjacencyLists();
	std::vector<AdjacencyLists>& getElementAdjacencyLists() { return elemAdj; }
	std::vector<int>& getFollowedElemList() { return followedElemList; }
	void getFollowerForce(GeomState &u, Vector &elementInternalForce,
	                      Corotator **allCorot, FullSquareMatrix *kel,
	                      Vector &residual, double lambda = 1.0, double time = 0.0,
	                      Vector *reactions = NULL, bool compute_tangents = false);
	void getElemFollowerForce(int iele, GeomState &geomState, double *_f, int bufSize,
	                          Corotator *corotator, FullSquareMatrix &kel,
	                          double loadFactor, double time, bool compute_tangents,
	                          BlastLoading::BlastData *conwep);
	void getFictitiousForce(GeomState &geomState, Vector &elementForce, FullSquareMatrix *kel, Vector &residual,
	                        double time, GeomState *refState, Vector *reactions,
	                        FullSquareMatrix *mel, bool compute_tangents, Corotator **, FullSquareMatrix *cel = NULL);
	void getElemFictitiousForce(int iele, GeomState &geomState, double *f, FullSquareMatrix &kel,
	                            double time, GeomState *refState, FullSquareMatrix &mel, bool compute_tangents,
	                            Corotator *elemCorot = NULL, FullSquareMatrix *celArray = NULL);
#ifdef USE_EIGEN3
	void getNodeFictitiousForce(int inode, GeomState &geomState, double time, GeomState *refState, Eigen::Matrix3d &M,
	                            Eigen::Vector3d &f, Eigen::Matrix3d &K, bool compute_tangents);
#endif
	void getWeightedFictitiousForceOnly(const std::map<int, double> &weights, GeomState &geomState, Vector &elementForce,
	                                    FullSquareMatrix *kel, Vector &residual,
	                                    double time, GeomState *refState, Vector *reactions,
	                                    FullSquareMatrix *mel, bool compute_tangents);
	void getUnassembledFictitiousForce(GeomState &geomState, Vector &elementForce,
	                                   FullSquareMatrix *kel, Vector &unassemResidual,
	                                   double time, GeomState *refState, Vector *reactions,
	                                   FullSquareMatrix *mel, bool compute_tangents);
	void getUDEIMFictitiousForceOnly(const std::map<int, std::vector<int> > &weights, GeomState &geomState, Vector &elementForce,
	                                 FullSquareMatrix *kel, Vector &residual,
	                                 double time, GeomState *refState, Vector *reactions,
	                                 FullSquareMatrix *mel, bool compute_tangents);
	void transformElemStiffAndForce(const GeomState &geomState, double *elementForce,
	                                FullSquareMatrix &kel, int iele, bool compute_tangents);
	void transformElemStiffAndForce_S2E(const GeomState &geomState, double *elementForce,
	                                    FullSquareMatrix &kel, int iele, bool compute_tangents);
	void transformNodalMoment(const GeomState &geomState, double G[],
	                          double H[][3], int nnum, bool compute_tangents);
	void transformElemStiff(const GeomState &geomState, FullSquareMatrix &kel, int iele);
	void getWeightedStiffAndForceOnly(const std::map<int, double> &weights,
	                                  GeomState &u, Vector &elementInternalForce,
	                                  Corotator **allCorot, FullSquareMatrix *kel,
	                                  Vector &residual, double lambda, double time,
	                                  GeomState *refState, FullSquareMatrix *mel = NULL,
	                                  FullSquareMatrix *kelCopy = NULL);
	void getUnassembledStiffAndForceOnly(const std::map<int, std::vector<int> > &weights,
	                                     GeomState &u, Vector &elementInternalForce,
	                                     Corotator **allCorot, FullSquareMatrix *kel,
	                                     Vector &residual, int dispSize, double lambda, double time,
	                                     GeomState *refState, FullSquareMatrix *mel = NULL,
	                                     FullSquareMatrix *kelCopy = NULL);
	void getElemInternalForce(const GeomState &geomState, double time,
	                          const GeomState *refState, const Corotator &elemCorot,
	                          double *elemForce, FullSquareMatrix &elemStiff);
	void getElemInternalForce(const GeomState &geomState, double time,
	                          const Corotator &elemCorot,
	                          double *elemForce, FullSquareMatrix &elemStiff);
	void getInternalForce(GeomState &u, Vector &elementInternalForce,
	                      Corotator **allCorot, FullSquareMatrix *kel,
	                      Vector &residual, double lambda = 1.0, double time = 0.0,
	                      GeomState *refState = NULL, Vector *reactions = NULL,
	                      FullSquareMatrix *mel = NULL, FullSquareMatrix *cel = NULL);
	void getWeightedInternalForceOnly(const std::map<int, double> &weights,
	                                  GeomState &u, Vector &elementInternalForce,
	                                  Corotator **allCorot, FullSquareMatrix *kel,
	                                  Vector &residual, double lambda, double time,
	                                  GeomState *refState, FullSquareMatrix *mel = NULL,
	                                  FullSquareMatrix *kelCopy = NULL);
	void getUDEIMInternalForceOnly(const std::map<int, std::vector<int> > &weights,
	                               GeomState &u, Vector &elementInternalForce,
	                               Corotator **allCorot, FullSquareMatrix *kel,
	                               Vector &residual, int dispSize, double lambda, double time,
	                               GeomState *refState, FullSquareMatrix *mel = NULL,
	                               FullSquareMatrix *kelCopy = NULL);
	void getUnassembledNonLinearInternalForce(GeomState &u, Vector &elementInternalForce,
	                                          Corotator **allCorot, FullSquareMatrix *kel,
	                                          Vector &unassemResidual, std::map<int, std::pair<int,int> > &uDOFaDOFmap,
	                                          double lambda = 1.0, double time = 0.0, int tIndex = 0,
	                                          GeomState *refState = NULL, Vector *reactions = NULL,
	                                          FullSquareMatrix *mel = NULL,FullSquareMatrix *kelCopy = NULL);

	void applyResidualCorrection(GeomState &geomState, Corotator **corotators, Vector &residual, double rcoef = 1.0);
	void initializeMultipliers(GeomState &geomState, Corotator **corotators);
	void initializeParameters(bool flag=true);
	void updateMultipliers(GeomState &geomState, Corotator **corotators);
	void updateParameters(bool flag=true);
	double getError(Corotator **corotators, GeomState &gs);
	void getElementDisp(int iele, GeomState& geomState, Vector& disp);
	void getElementVelo(int iele, GeomState& geomState, Vector& velo);
	void getElementAccel(int iele, GeomState& geomState, Vector& accel);
	double getKineticEnergy(double *velocity, FullSquareMatrix *mel);
	double getStrainEnergy(GeomState *geomState, Corotator **allCorot);
	double getDissipatedEnergy(GeomState *geomState, Corotator **allCorot);
	// compute dissipated energy for elements related to specific attribute
	double* getDissipatedEnergy(GeomState *geomState, Corotator **allCorot, int groupId);
	void setModalEnergies(double Wele, double Wkin, double Wdmp);
	void computeEnergies(GeomState *geomState, Vector &force, double time, Vector *aeroForce, double *vel,
	                     Corotator **allCorot, SparseMatrix *M, SparseMatrix *C, double &Wela, double &Wkin,
	                     double &Wdis, double &error); // Nonlinear statics and dynamics
	void computeEnergies(Vector &disp, Vector &force, double time, Vector *aeroForce, Vector *vel,
	                     SparseMatrix *K, SparseMatrix *M, SparseMatrix *C, double &Wela, double &Wkin,
	                     double &error); // Linear dynamics
	void computeExtAndDmpEnergies(Vector &disp, Vector &force, double time, Vector *aeroForce,
	                              Vector *vel, SparseMatrix *C, Vector *folForce = NULL);
	double getWext() { return Wext; }
	double getWaero() { return Waero; }
	double getWdmp() { return Wdmp; }
	void handleElementDeletion(int iele, GeomState &geomState, double time, Corotator &elemCorot, double *elemForce = 0);

	void getGeometricStiffness(GeomState &u, Vector &elementInternalForce,
	                           Corotator **allCorot, FullSquareMatrix *&kel);
	void computeGeometricPreStress(Corotator **&allCorot, GeomState *&geomState,
	                               FullSquareMatrix *&kelArray, StaticTimers *times,
	                               FullSquareMatrix *&geomKelArray, FullSquareMatrix *&melArray,
	                               bool melFlag = false);
	ControlInterface* getUserSuppliedFunction();
	void setsizeSfemStress(int fileNumber);
	int getsizeSfemStress() { return sizeSfemStress; }
	double * getSfemStress(int fileNumber, double* dummystress); // dummystress is an arbitrary vector that tells that the stress is double*
	DComplex * getSfemStress(int fileNumber, DComplex* dummystress) {std::cerr <<"DComplex * Domain::getSfemStress not implemented" << std::endl; return 0;};
	void updateSfemStress(double* str, int fileNumber);
	void updateSfemStress(DComplex* str, int fileNumber) {std::cerr <<"Domain::updateSfemStress(DComplex*,.) not implemented" << std::endl;};
	// Functions to suport the parsing:
	void buildSensitivityInfo();
	void addSensitivity(OutputInfo &oInfo);
	int  addDMass(int, int, double, int jdof = -1);
	DMassData* getFirstDMassData() { return firstDiMass; }
	int getDMassNumber() { return nDimass; }
	int  setDirichlet(int,BCond *);
	int  setDirichletFluid(int,BCond *);
	int  addLMPC(LMPCons *, bool checkflag=true);
	void printLMPC();
	void printLMPC2();
	void normalizeLMPC();
	void setPrimalLMPCs(int &numDual, int &numPrimal);
	void checkLMPCs(Connectivity *nodeToSub);
	void printFSI(FILE* file=stderr);
	Connectivity * makeLmpcToNode();
	Connectivity * makeLmpcToNode_primal();
	void makeFsiToNode();
	Connectivity *getFsiToNode() { return fsiToNode; }
	int getNumShapeVars() { return numShapeVars; }
	int getNumSensitivityQuantityTypes() { return numSensitivityQuantityTypes; }
	int getNumStressNodes() { return numStressNodes; }
	int getNumThicknessGroups() { return numThicknessGroups; }
	int getNumDispNodes() { return numDispNodes; }
	int getTotalNumDispDofs() { return numTotalDispDofs; }
	int getNumFSI() { return numFSI; }
	void setNumFSI(int n) { numFSI = n; }
	ResizeArray<LMPCons *> &getFSI() { return fsi; }
	virtual int  setNeuman(int,BCond *);
	int  setNeumanModal(int, BCond *);
	int  setIDis6(int, BCond *);
	int  setIDisModal(int, BCond *);
	int  setIDis(int, BCond *);
	int  setIVel(int, BCond *);
	int  setIVelModal(int, BCond *);
	int  setIAcc(int, BCond *);
	void setLoadFactor(int, int, double);
	void setLoadFactorMFTT(int, int, int);
	std::map<std::pair<int,int>,int>& getLoadFactorMFTT() { return loadfactor_mftt; };
	void setLoadFactorHFTT(int, int, int);
	void setLoadFactorTemp(int, bool);
	void setLoadFactorGrav(int, bool);
	void checkCases();
	double getLoadFactor(int) const;
	int  setMFTT(MFTTData *, int);
	MFTTData * getDefaultMFTT() const;
	MFTTData * getMFTT(int) const;
	int getNumMFTT() const;
	int  setHFTT(MFTTData *, int i);
	MFTTData * getDefaultHFTT() const;
	MFTTData * getHFTT(int) const;
	int getNumHFTT() const;
	int  addYMTT(MFTTData *);
	int  addSDETAFT(MFTTData *);
	void updateSDETAF(StructProp* p, double omega);
#ifdef USE_EIGEN3
	int  addRUBDAFT(GenMFTTData<Eigen::Vector4d> *);
#endif
	void updateRUBDAFT(StructProp* p, double omega);
	void printYMTT();
	int  addCTETT(MFTTData *);
	std::pair<int, ResizeArray<MFTTData*>* >* getCTETT() { return new std::pair<int, ResizeArray<MFTTData*>* >(numCTETT,&ctett); };
	std::pair<int, ResizeArray<MFTTData*>* >* getYMTT() { return new std::pair<int, ResizeArray<MFTTData*>* >(numYMTT,&ymtt); };
	std::pair<int, ResizeArray<MFTTData*>* >* getYSST() { return new std::pair<int, ResizeArray<MFTTData*>* >(numYSST,&ysst); };
	std::pair<int, ResizeArray<MFTTData*>* >* getYSSRT() { return new std::pair<int, ResizeArray<MFTTData*>* >(numYSSRT,&yssrt); };
     std::pair<int, ResizeArray<MFTTData*>* >* getSS1DT() { return new std::pair<int, ResizeArray<MFTTData*>* >(numSS1DT,&ss1dt); };
     std::pair<int, ResizeArray<SS2DTData*>* >* getSS2DT() { return new std::pair<int, ResizeArray<SS2DTData*>* >(numSS2DT,&ss2dt); };
     std::pair<int, ResizeArray<MFTTData*>* >* getYMST() { return new std::pair<int, ResizeArray<MFTTData*>* >(numYMST,&ymst); };
	void printCTETT();
	int  addYSST(MFTTData *);
	int  addYSSRT(MFTTData *);
     int  addSS1DT(MFTTData *);
     int  addSS2DT(SS2DTData *);
     int  addYMST(MFTTData *);
	void computeTDProps();

	ShapeSensitivityData getShapeSensitivityData() { return shapeSenData; }
	int getNumSensitivities() { return numSensitivity; }

	void setUpData(int topFlag);
	SolverInfo  &solInfo() { return sinfo; }
	const SolverInfo  &solInfo() const { return sinfo; }
	MatrixTimers &getTimers() { return *matrixTimers; }
	void setGravity(double ax, double ay, double az);

	void setGravitySloshing(double gg);

	virtual int glToPackElem(int e);
	void ProcessSurfaceBCs(int topFlag);
	void setNewProperties(int s);
	void assignRandMat();
	void retrieveElemset();

	void computeWeightWRTthicknessSensitivity(int, AllSensitivities<double> &allSens);
	void computeWeightWRTShapeVariableSensitivity(int, AllSensitivities<double> &allSens);
	void computeStiffnessWRTthicknessSensitivity(int, AllSensitivities<double> &allSens);
	void computeStiffnessWRTShapeVariableSensitivity(int, AllSensitivities<double> &allSens);
	void makePreSensitivities(AllSensitivities<double> &allSens, double *);
	void makePreSensitivities(AllSensitivities<DComplex> &allSens, DComplex *);

	void subtractGravityForceSensitivityWRTthickness(int, AllSensitivities<double> &allSens);
	void subtractGravityForceSensitivityWRTShapeVariable(int, AllSensitivities<double> &allSens);
	void computeDisplacementWRTShapeVariableDirectSensitivity(int, GenSolver<double> *,
	                                                          AllSensitivities<double> &,
	                                                          GenSparseMatrix<double> *spm=0,
	                                                          GenSparseMatrix<double> *K=0);
	void computeDisplacementWRTShapeVariableAdjointSensitivity(int,
	                                                           AllSensitivities<double> &,
	                                                           GenSparseMatrix<double> *spm=0);
	void computeDisplacementWRTthicknessDirectSensitivity(int, GenSolver<double> *,
	                                                      AllSensitivities<double> &,
	                                                      GenSparseMatrix<double> *spm=0,
	                                                      GenSparseMatrix<double> *K=0);
	void computeDisplacementWRTthicknessAdjointSensitivity(int,
	                                                       AllSensitivities<double> &,
	                                                       GenSparseMatrix<double> *spm=0);
	void computeStressVMDualSensitivity(int, GenSolver<double> *,
	                                    AllSensitivities<double> &,
	                                    GenSparseMatrix<double> *spm=0,
	                                    GenSparseMatrix<double> *K=0);
	void computeNormalizedVonMisesStress(Vector&, double*, int, Vector&, bool normalized=true);
	void computeNormalizedNLVonMisesStress(GeomState &, GeomState *, Corotator **, int, Vector &, Vector &, bool normalized=true);
	void scaleToTrueVonMisesStress(Vector&);
	void computeAggregatedStressDenom(Vector &stress);
	void computeAggregatedStressVMDualSensitivity(int, GenSolver<double> *,
	                                              AllSensitivities<double> &,
	                                              GenSparseMatrix<double> *spm=0,
	                                              GenSparseMatrix<double> *K=0);
	void computeDisplacementDualSensitivity(int, GenSolver<double> *,
	                                        AllSensitivities<double> &,
	                                        GenSparseMatrix<double> *spm=0,
	                                        GenSparseMatrix<double> *K=0);
	void computeNLStaticWRTthicknessSensitivity(int, AllSensitivities<double> &allSens,
	                                            GeomState *refState, GeomState *geomState, Corotator **allCorot);
	void computeLinearStaticWRTthicknessSensitivity(int, AllSensitivities<double> &allSens, GenVector<double> *sol,
	                                                GeomState *refState, GeomState *geomState, Corotator **allCorot);
	void computeLinearStaticWRTShapeVariableSensitivity(int, AllSensitivities<double> &allSens, GenVector<double> *sol);
	void computeStressVMWRTthicknessDirectSensitivity(int, AllSensitivities<double> &allSens, GenVector<double> *sol, double *bcx,
	                                                  GeomState *refState, GeomState *geomState, Corotator **allCorot,
	                                                  bool isDynam = false);
	void computeNLStressVMWRTthicknessDirectSensitivity(int, AllSensitivities<double> &allSens, GeomState *geomState,
	                                                    Corotator **allCorot, bool isDynam = false);
	void computeNLStressVMWRTthicknessAdjointSensitivity(int, AllSensitivities<double> &allSens, GeomState *geomState,
	                                                     Corotator **allCorot, bool *includeStressNodes, bool isDynam = false);
	void computeAggregatedStressVMWRTShapeVariableSensitivity(int, AllSensitivities<double> &allSens,
	                                                          GenVector<double> *sol, double *bcx,
	                                                          bool isDynam = false);
	void computeAggregatedStressVMWRTthicknessSensitivity(int, AllSensitivities<double> &allSens,
	                                                      GenVector<double> *sol, double *bcx,
	                                                      bool isDynam = false);
	void computeAggregatedNLStressVMWRTthicknessSensitivity(int, AllSensitivities<double> &allSens,
	                                                        GeomState *geomState, GeomState *refState, Corotator **allCorot,
	                                                        Vector &, Vector &, Vector &, bool isDynam = false);
	void computeStressVMWRTthicknessAdjointSensitivity(int, AllSensitivities<double> &allSens,
	                                                   GenVector<double> *sol, double *bcx,
	                                                   bool *includeStressNodes, bool isDynam = false);
	void computeStressVMWRTdisplacementSensitivity(int, AllSensitivities<double> &allSens,
	                                               GenVector<double> *sol, double *bcx);
	void computeNLStressVMWRTdisplacementSensitivity(int, AllSensitivities<double> &allSens,
	                                                 GeomState *geomState, Corotator **allCorot);
	void computeAggregatedStressVMWRTdisplacementSensitivity(int, AllSensitivities<double> &allSens,
	                                                         GenVector<double> *sol, double *bcx);
	void computeAggregatedNLStressVMWRTdisplacementSensitivity(int, AllSensitivities<double> &allSens,
	                                                           GeomState *geomState, GeomState *refState, Corotator **allCorot,
	                                                           Vector &, Vector &, Vector &);
	void computeStressVMWRTShapeVariableDirectSensitivity(int, AllSensitivities<double> &allSens,
	                                                      GenVector<double> *sol, double *bcx,
	                                                      bool isDynam = false);
	void computeStressVMWRTShapeVariableAdjointSensitivity(int, AllSensitivities<double> &allSens,
	                                                       GenVector<double> *sol, double *bcx,
	                                                       bool *includeStressNodes, bool isDynam = false);
	void computeStressVMWRTMachNumberSensitivity(AllSensitivities<double> &allSens);
	void computeStressVMWRTangleOfAttackSensitivity(AllSensitivities<double> &allSens);
	void computeStressVMWRTyawAngleSensitivity(AllSensitivities<double> &allSens);
	void makePostSensitivities(GenSolver<double> *, GenSparseMatrix<double> *, AllSensitivities<double> &allSens,
	                           GenVector<double> *sol, double *, GenSparseMatrix<double> *K=0, bool isDynam = false,
	                           GeomState *rs=NULL, GeomState *gs=NULL, Corotator **allCorot = NULL, bool isNonLin = false);
	void makePostSensitivities(GenSolver<DComplex> *, GenSparseMatrix<DComplex> *, AllSensitivities<DComplex> &allSens,
	                           GenVector<DComplex> *sol, DComplex *, GenSparseMatrix<DComplex> *K=0, bool isDynam = false,
	                           GeomState *rs=NULL, GeomState *gs=NULL, Corotator **allCorot = NULL, bool isNonLin = false);
	void makeNLPostSensitivities(GenSolver<double> *, AllSensitivities<double> &allSens,
	                             GeomState *rs, GeomState *gs, Corotator **allCorot, bool isDynam = false);
	void makeNLPostSensitivities(GenSolver<DComplex> *, AllSensitivities<DComplex> &allSens,
	                             GeomState *rs, GeomState *gs, Corotator **allCorot, bool isDynam = false);
	void makeThicknessGroupElementFlag();

/** ... General build functions to replace the specialized build
  * ... functions and allow us to reuse the code in each problem
  * ... type (i.e. use makeSparseOps in statics, dynamics, eigen, etc.)
  */
	template<class Scalar>
	void buildPreSensitivities(AllSensitivities<Scalar> &ops, Scalar *);

	template<class Scalar>
	void buildPostSensitivities(GenSolver<Scalar> *sysSolver, GenSparseMatrix<Scalar> *, GenSparseMatrix<Scalar> *,
	                            AllSensitivities<Scalar> &ops, GenVector<Scalar> *sol, Scalar *, bool isDynam = false,
	                            GeomState *refState = NULL, GeomState *geomState = NULL, Corotator **allCorot = NULL, bool isNonLin = false);

	template<class Scalar>
	void buildNLPostSensitivities(GenSolver<Scalar> *sysSolver, AllSensitivities<Scalar> &ops,
	                              GeomState *refState, GeomState *geomState, Corotator **allCorot, bool isDynam=false);

	template<class Scalar>
	void buildOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef, double Ccoef,
	              Rbm *rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *melArray = 0,
	              FullSquareMatrix *celArray = 0, bool factor = true);

	template<class Scalar>
	void makeStaticOpsAndSolver(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
	                            double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
	                            Rbm *rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *melArray = 0,
	                            FullSquareMatrix *celArray = 0);

	template<class Scalar>
	void makeDynamicOpsAndSolver(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
	                             double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
	                             Rbm *rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *mel = 0,
	                             FullSquareMatrix *celArray = 0);

	template<class Scalar>
	void rebuildOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef, double Ccoef,
	                Rbm* rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *mel = 0,
	                FullSquareMatrix *celArray = 0, bool factor = true);

	template<class Scalar>
	void makeSparseOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
	                   double Ccoef, GenSparseMatrix<Scalar> *mat = 0,
	                   FullSquareMatrix *kelArray = 0, FullSquareMatrix *melArray = 0,
	                   FullSquareMatrix *celArray = 0);

	template<class Scalar>
	GenDBSparseMatrix<Scalar> *constructDBSparseMatrix(DofSetArray *dof_set_array=0,
	                                                   Connectivity *cn=0);

#ifdef USE_EIGEN3
	template<typename Scalar>
	GenEiSparseMatrix<Scalar, Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *
	constructEiSparse(DofSetArray *dof_set_array=0, Connectivity *cn=0, bool flag=true);
#endif

	template<typename Scalar, typename SolverClass>
	GenEiSparseMatrix<Scalar,SolverClass> *constructEiSparseMatrix(DofSetArray *dof_set_array=0,
	                                                               Connectivity *cn=0, bool flag=true);

	template<class Scalar>
	GenCuCSparse<Scalar> *constructCuCSparse(DofSetArray *dof_set_array=0);

	template<class Scalar>
	GenCuCSparse<Scalar> *constructCCSparse(DofSetArray *dof_set_array=0);

	template<class Scalar>
	GenBLKSparseMatrix<Scalar> *constructBLKSparseMatrix(DofSetArray*, Rbm *rbm = 0);

	template<class Scalar>
	GenNBSparseMatrix<Scalar> *constructNBSparseMatrix();

	Rbm              *constructRbm(bool printFlag = true);
	Rbm              *constructHzem(bool printFlag = true);
	Rbm              *constructSlzem(bool printFlag = true);

	template<class Scalar>
	void addGravityForce(GenVector<Scalar>& force);

	template<class Scalar>
	void addGravityForceSensitivity(GenVector<Scalar>& forceSen);

	template<class Scalar>
	void addPressureForce(GenVector<Scalar>& force, int which = 2, double time = 0.0);

	template<class Scalar>
	void addAtddnbForce(GenVector<Scalar>& force, int which = 2, double time = 0.0);

	template<class Scalar>
	void addAtdrobForce(GenVector<Scalar>& force, int which = 2, double time = 0.0);

	template<class Scalar>
	void addThermalForce(GenVector<Scalar>& force);

	template<class Scalar>
	void addMpcRhs(GenVector<Scalar>& force, double t = 0);

	void buildPrescDisp(Vector &v, double lambda);
	void buildPrescDisp(Vector &v, double t, double dt);

	template<class Scalar>
	void buildRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *kuc = 0);

	template<class Scalar>
	void computeReactionForce(GenVector<Scalar> &fc, GenVector<Scalar> &Vu,
	                          GenSparseMatrix<Scalar> *kuc, GenSparseMatrix<Scalar> *kcc = 0);
	void computeReactionForce(Vector &fc, Vector &Du, Vector &Vu, Vector &Au,
	                          double *bcx, double *vcx, double *acx,
	                          SparseMatrix *_kuc, SparseMatrix *_kcc,
	                          SparseMatrix *_cuc, SparseMatrix *_ccc,
	                          SparseMatrix *_muc, SparseMatrix *_mcc);
	void computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
	                          FullSquareMatrix *kel, double lambda, GeomState *refState);
	void computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
	                          FullSquareMatrix *kel, double time, GeomState *refState,
	                          Vector &Vu, Vector &Au, double *vcx, double *acx,
	                          SparseMatrix *_cuc, SparseMatrix *_ccc,
	                          SparseMatrix *_muc, SparseMatrix *_mcc);
	bool reactionsReqd(double time, int step);

	template<class Scalar>
	void buildFreqSweepRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *muc,
	                            GenSparseMatrix<Scalar> **cuc_deriv,
	                            GenSparseMatrix<Scalar> **kuc_deriv,
	                            int iRHS, double omega);
	template<class Scalar>
	void buildDeltaK(double w0, double w, GenSparseMatrix<Scalar> *deltaK,
	                 GenSparseMatrix<Scalar> *deltaKuc);

	template<class Scalar>
	void buildRHSForce(GenVector<Scalar> &force,GenVector<Scalar> &tmp,
	                   GenSparseMatrix<Scalar> *kuc,
	                   GenSparseMatrix<Scalar> *muc,
	                   GenSparseMatrix<Scalar> **cuc_deriv,
	                   GenSparseMatrix<Scalar> **kuc_deriv,
	                   GenSparseMatrix<Scalar> **kuc_arubber_l,
	                   GenSparseMatrix<Scalar> **kuc_arubber_m,
	                   double omega, double delta_omega,
	                   GeomState *gs=0);

	void initNodalTemperatures();
	double * getNodalTemperatures();

	// Main program control functions.
	void arcLength();
	void createCorotators(Corotator **allCorot);
	void createContactCorotators(Corotator **allCorot, FullSquareMatrix *kArray, FullSquareMatrix *mArray);
	void preProcessing();
	FILE * openFile(char *fileName, const char *extension);
	void printStatistics(bool domain_decomp);

	// static & freq response post processing function
	template<class Scalar>
	void postProcessing(GenVector<Scalar> &sol, Scalar *bcx, GenVector<Scalar> &force,
	                    int ndflag = 0, int index = 0, double time = 0, double eigV = 0.0,
	                    GenSparseMatrix<Scalar> *kuc = NULL, GenSparseMatrix<Scalar> *kcc = NULL);

	// sensitivity post-processing function
	template<class Scalar>
	void sensitivityPostProcessing(AllSensitivities<Scalar> &allSens, GenVector<Scalar> *sol = 0, Scalar *bcx=0,
	                               GeomState *gs=0, GeomState *rs=0, Corotator **ac=0);

	template<class Scalar>
	void sensitivityPostProcessing(AllSensitivities<Scalar> &allSens, GenDistrVector<Scalar> *sol, Scalar *bcx=0,
	                               GeomState *gs=0, GeomState *rs=0, Corotator **ac=0) {}

	template<class Scalar>
	void sensitivityPostProcessing(AllSensitivities<Scalar> &allSens, DistrBlockVector<Scalar> *sol, Scalar *bcx=0,
	                               GeomState *gs=0, GeomState *rs=0, Corotator **ac=0) {}

	// Nonlinear post processing function
	void postProcessing(GeomState *geomState, Vector &force, Vector &aeroForce, double time = 0.0,
	                    int step = 0, double *velocity = 0, double *vcx = 0, Corotator **allCorot = 0,
	                    double *acceleration = 0, double *acx = 0, GeomState *refState = 0,
	                    Vector *reactions = 0, SparseMatrix *M = 0, SparseMatrix *C = 0);

	// Pita Nonlinear post processing function
	void pitaPostProcessing(int timeSliceRank, GeomState *geomState, Vector &force, Vector &aeroForce, double time = 0.0,
	                        int step = 0, double *velocity = 0, double *vcx = 0, Corotator **allCorot = 0,
	                        double *acceleration = 0, double *acx = 0, GeomState *refState = 0,
	                        Vector *reactions = 0, SparseMatrix *M = 0, SparseMatrix *C = 0);

	// Dynamic functions and thermal functions
	void getHeatFlux(Vector &tsol, double *bcx, int fileNumber, int hgIndex,
	                 double time=0);
	void getHeatFlux(ComplexVector &tsol, DComplex *bcx, int fileNumber, int hgIndex, double time=0)
	{ std::cerr << " *** WARNING: Domain::getHeatFlux(Complex) is not implemented \n"; }
	void getTrussHeatFlux(Vector &tsol, double *bcx, int fileNumber,
	                      int hgIndex, double time=0);
	void getTrussHeatFlux(ComplexVector &tsol, DComplex *bcx, int fileNumber, int hgIndex, double time=0)
	{ std::cerr << " *** WARNING: Domain::getTrussHeatFlux(Complex) is not implemented \n"; }
	template <class Scalar>
	void computeConstantForce(GenVector<Scalar>& constantForce, GenSparseMatrix<Scalar>* kuc = 0);
	template <class Scalar>
	void addConstantForceSensitivity(GenVector<Scalar>& constantForce, GenSparseMatrix<Scalar>* kuc = 0);
	template <class Scalar>
	void computeExtForce4(GenVector<Scalar>& force, const GenVector<Scalar>& constantForce, double t,
	                      GenSparseMatrix<Scalar> *kuc = 0, ControlInterface *userSupFunc = 0,
	                      GenSparseMatrix<Scalar> *cuc = 0, double tm = 0, GenSparseMatrix<Scalar> *muc = 0);
	template <class Scalar>
	void computeExtForce(GenVector<Scalar>& force, double t, GenSparseMatrix<Scalar> *kuc = 0, ControlInterface *userSupFunc = 0,
	                     GenSparseMatrix<Scalar> *cuc = 0, double tm = 0, GenSparseMatrix<Scalar> *muc = 0);
	void computeExtForce(Vector &f, double t, int tIndex,
	                     SparseMatrix *kuc, Vector &prev_f);
	void computeUnamplifiedExtForce(GenVector<double>& fcon, int loadsetid);

	int  probType() { return sinfo.probType; }

	template<typename DynamMatType>
	double computeStabilityTimeStep(DynamMatType&);
	double computeStabilityTimeStepROM(GenFullSquareMatrix<double>&);
	double computeStabilityTimeStep(FullSquareMatrix *kelArray, FullSquareMatrix *melArray, GeomState *geomState, int &eid);

	void initDispVeloc(Vector& d_n, Vector& v_n, Vector& a_n, Vector &v_p, const char* = "");
	void initDispVelocOnTimeSlice (Vector& d_n, Vector& v_n, int sliceRank); // PITA: Use user-provided initial seeds
	void initTempVector(Vector& d_n, Vector& v_n, Vector& v_p);
	void writeRestartFile(double time, int timeIndex,
	                      Vector &d_n, Vector &v_n, Vector &v_p, double Fref = 0.0, const char* = "");
	void writeRestartFile(double time, int timeIndex, Vector &v_n, Vector &a_n,
	                      GeomState *geomState, const char* = "");
	void readRestartFile(Vector &d_n, Vector &v_n,
	                     Vector &a_n, Vector &v_p, double *bcx,
	                     double *vcx, GeomState &geomState, const char* = "");
	void getOrAddDofForPrint(bool ad, Vector& d_n, double* bcx, int iNode,
	                         double *xdata, int *dofx, double *ydata=0, int *dofy=0, double *zdata=0, int *dofz=0);
	void addVariationOfShape_StructOpt(int iNode, CoordSet *nodescopy, double &x, double &y, double &z);
	void aeroSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState = 0);
	void aeroheatSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState = 0);
	void thermohSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState = 0);
	void buildAeroelasticForce(Vector &f, PrevFrc& prevFrc, int tIndex, double t, double gamma, double alphaf, GeomState* geomState = 0);
	void buildAeroheatFlux(Vector &f, Vector &prev_f, int tIndex, double t);
	void thermoeComm();
	void dynamOutput(int, double, double*, DynamMat&, Vector&, Vector &, Vector&, Vector&, Vector&, Vector &, double*, double* = 0);
	void pitaDynamOutput(int, double* bcx, DynamMat&, Vector&, Vector &, Vector&, Vector&, Vector&, Vector &,
	                     double* vcx, double* acx, int sliceRank, double time);

protected:
	void dynamOutputImpl(int, double* bcx, DynamMat&, Vector&, Vector &, Vector&, Vector&, Vector&, Vector &, double*, double*, double, int, int);

public:
	void tempdynamOutput(int, double*, DynamMat&, Vector&, Vector&, Vector&,
	                     Vector&);

     double computeStructureMass(bool printFlag = true, int groupId = 0);
	double computeFluidMass();
	double getStructureMass();
	int returnLocalDofNum(int, int);

	// Condition number estimated routines
	double computeConditionNumber(DynamMat&);

	//   Output Related functions
	template<class Scalar>
	int processOutput(OutputInfo::Type &type, GenVector<Scalar> &sol, Scalar *bcx, int iInfo,
	                  double time, double freq = 0, int printFlag = 0);

	template<class Scalar>
	int processDispTypeOutputs(OutputInfo &oinfo, Scalar (*sol)[11], int iInfo,
	                           int numNodes, double time, double freq = 0, int printFlag = 0);
	void getStressStrain(Vector &sol, double *bcx, int fileNumber,
	                     int strInd, double time = 0, int printFlag =0);
	void getStressStrain(ComplexVector &sol, DComplex *bcx, int fileNumber,
	                     int strInd, double time = 0, int printFlag =0);
	void getPrincipalStress(Vector &sol, double *bcx, int fileNumber,
	                        int strInd, double time = 0);
	void getPrincipalStress(ComplexVector &sol, DComplex *bcx, int fileNumber,
	                        int strInd, double time = 0) { std::cerr << "Domain::getPrincipalStress is not implemented for complex\n"; }
	void getElementForces(Vector &sol, double *bcx, int fileNumber,
	                      int forceIndex, double time = 0);
	void getElementForces(ComplexVector &sol, DComplex *bcx, int fileNumber,
	                      int forceIndex, double time = 0) { std::cerr << "Domain::getElementForces is not implemented for complex\n"; }
	void getElementAttr(int fileNumber, int typ, double time=0.0);
	void getKtimesU(Vector &dsp, double *bcx, Vector &ext_f, double eta,
	                FullSquareMatrix *kelArray=0);
	void getWeightedKtimesU(const std::map<int, double> &weights,
	                        Vector &dsp, double *bcx, Vector &ext_f, double eta,
	                        FullSquareMatrix *kelArray=0);
	void getUnassembledKtimesU(const std::map<int, std::vector<int> > &weights,
	                           Vector &dsp, double *bcx, Vector &ext_f, double eta,
	                           FullSquareMatrix *kelArray=0);
	void getElemKtimesU(int iele, int numEleDOFs, Vector &dsp, double *elForce,
	                    FullSquareMatrix *kelArray, double *karray);
	void getSloshDispAll(Vector &tsol, double *bcx, int fileNumber, double time);
	void getSloshDispAll(ComplexVector &tsol, complex<double> *bcx, int fileNumber, double time) { std::cerr << "getSloshDispAll(complex) not implemented\n"; }
	void getSloshDisp(Vector &tsol, double *bcx, int fileNumber, int hgIndex, double time);
	void getSloshDisp(ComplexVector &tsol, complex<double> *bcx, int fileNumber, int hgIndex, double time) { std::cerr << "getSloshDisp(complex) not implemented\n"; }
	double getKineticEnergy(Vector &sol, SparseMatrix *gMass);

	template<class Scalar>
	void scaleDisp(Scalar *u);
	template<class Scalar>
	void scaleInvDisp(Scalar *u);
	template<class Scalar>
	void scaleDisp(Scalar *u, double alpha);
	template<class Scalar>
	void forceContinuity(Scalar *u) {}
	template<class Scalar>
	void forceAssemble(Scalar *u) {}

	template<class Scalar>
	int mergeDistributedDisp(Scalar (*xyz)[11], Scalar *u, Scalar *bcx = 0, Scalar (*xyz_loc)[11] = NULL);
#ifdef USE_EIGEN3
	template<class Scalar>
	void mergeDistributedDispSensitivity(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *,
	                                     Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *);
	template<class Scalar>
	void mergeAdjointDistributedDispSensitivity(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *,
	                                            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *,
	                                            Scalar *, GenVector<Scalar>*);
#endif
	template<class Scalar>
	void forceDistributedContinuity(Scalar *u, Scalar (*xyz)[11]);//DofSet::max_known_nonL_dof

	Connectivity makeSommerToNode();
	Connectivity *prepDirectMPC();
	// renumbering functions
	Renumber  getRenumbering();
	Renumber* getRenumberingFluid();

	void makeNodeToNode_sommer();

	// Eigen solver
	void eigenOutput(Vector& eigenValues, VectorSet& eigenVectors, double* bcx = 0, int convEig = 0);
#ifdef USE_EIGEN3
	void eigenQROutput(Eigen::MatrixXd& Xmatrix, Eigen::MatrixXd& Qmatrix, Eigen::MatrixXd& Rmatrix);
#endif
	void setEigenValue(double _lbound, int _nshifts, int _maxArnItr = 0);
	void setEigenValues(double _lbound, double _ubound, int _neigps = 0, int _maxArnItr = 0);

	Solver *getSolver() { return solver; }

	template<class Scalar>
	void getSolverAndKuc(AllOps<Scalar> &allOps, FullSquareMatrix *kelArray, Rbm *rbm, bool factorize=true);

	void make_constrainedDSA();
	void make_constrainedDSA(const int *bc);

	void make_constrainedDSA(gsl::span<const int> bc) {
		make_constrainedDSA(bc.data());
	}
	void make_constrainedDSA(int fake);

	// returns the number of dof
	int numdof() const { return dsa ? dsa->size() : -1; }

	// returns the number of unconstrained dof
	int numUncon() const {
		return c_dsa ? c_dsa->size() : dsa ? dsa->size() : 0;
	}

	int coordVecSize() const {
		return dsa ? dsa->numNodes()*6 : 0;
	}

	// returns the number of unconstrained Fluid dof
	int numUnconFluid() const {
		return c_dsaFluid ? c_dsaFluid->size() : dsaFluid->size();
	}
	// returns a pointer to the dirichlet boundary condtions
	BCond* getDBC() { return dbc; }

	// returns the number of dirichlet bc
	int  nDirichlet() const { return numDirichlet; }
	int  nDirichletFluid() const { return numDirichletFluid; }
	int  nDispDirichlet() const { return numDispDirichlet; }
	void setNumDispDirichlet(int n) { numDispDirichlet = n; }

	// returns the number of initial displacements
	int numInitDisp() const { return numIDis;  }
	int numInitDispModal() const { return numIDisModal; }
	int numInitDisp6() const { return numIDis6; }

	// returns a pointer to the initial displacement boundary condtions
	BCond* getInitDisp()  { return iDis;  }
	BCond* getInitDispModal() { return iDisModal; }
	BCond* getInitDisp6() { return iDis6; }

	// returns a pointer to the initial velocities
	int    numInitVelocity() const { return numIVel; }
	BCond* getInitVelocity() { return iVel; }
	int    numInitVelocityModal() const { return numIVelModal; }
	BCond* getInitVelocityModal() { return iVelModal; }

	// returns the number of neumann bc
	int  nNeumann() const { return numNeuman; }
	int  nNeumannModal() const { return numNeumanModal; }

	// returns a pointer to the neumann boundary condtions
	BCond* getNBC() { return nbc; }
	BCond* getNBCModal() { return nbcModal; }

	// returns the number of nodes
	int  numNodes() const { return (nodeToNode) ? nodeToNode->csize() : numnodes; }
	int  numNode() const { return numnodes; }
	void setNumNodes(int n)  { numnodes = n; }  // includes virtual nodes
	int  numGlobalNodes() const { return numnodes; }

	// returns the number of elements
	int  numElements() const { return numele; }
	void setNumElements(int n) { numele = n; }
	int addElem(int ele, int type, int nnd, int *nd)
	{ packedEset.elemadd(ele, type, nnd, nd); return ele; }

	// returns the maximum possible number of elements for simulations in which the
	// number elements may change from one load/time step to another
	int maxNumElements() const { return numele - contactSurfElems.size() + maxContactSurfElems; }

	// returns the number of dofs
	void setNumDofs(int n) { numdofs = n; }
	int  numDofs() const { return numdofs; }

	// returns the packed element set (allows element numbering gaps)
	Elemset& getElementSet() { return packedEset; }

	// returns the value of the gravity Force flag
	int  gravityFlag();

	// returns the value of the pressure force flag
	int  pressureFlag();

	// returns the value of the contact force flag
	int  tdenforceFlag() { return int(nMortarCond > 0 && sinfo.newmarkBeta == 0.0 && sinfo.tdenforceFlag); } // TD enforcement (contact/tied surfaces with ACME) used for explicit dynamics

	int  thermalFlag();

	int  radiationFlag() { return sinfo.radiationFlag; }

	// returns the maximum number of dofs per element
	int  maxNumDOF() const { return maxNumDOFs; }

	// returns a pointer to the Connectivity allDOFs
	Connectivity *getAllDOFs() { return allDOFs; }
	void deleteAllDOFs() { delete allDOFs; allDOFs = 0; }

	CoordSet& getNodes() { return nodes; }
	const CoordSet& getNodes() const { return nodes; }

	ConstrainedDSA * getCDSA() { return c_dsa; }
	DofSetArray *    getDSA()  { return dsa;   }
	ConstrainedDSA * getCDSAFluid() { return c_dsaFluid; }
	DofSetArray *    getDSAFluid()  { return dsaFluid;   }

	// function that returns composite layer info
	LayInfo *getLayerInfo(int num);

	// function to get fluid exchanger
	FlExchanger * getFileExchanger() { return flExchanger; }
	void makeNodeTable(int topFlag);
	void makeTopFile(int topFlag);

	template<class Scalar> friend class GenDecDomain;
	template<class Scalar> friend class GenDistrDomain;
	friend class HData;

	// for output of tensor data

	double crossScale;

	double *** CPoint;
	double **  MidPoint;

	char *getBasename( char * fname);

	void getCompositeData(int iInfo,double time);

	void aeroPreProcess(Vector&, Vector&, Vector&, Vector &,double*,double*);
	void aeroSensitivityPreProcess(Vector&, Vector&, Vector&, Vector &,double*,double*);
	void sendDisplacements(Vector&, Vector&, Vector&, Vector&, double*, double *);
	void thermoePreProcess();

	void aeroHeatPreProcess(Vector&, Vector&, Vector&, double *bcx );
	void thermohPreProcess(Vector&, Vector&, Vector&, double *bcx );

	const Connectivity *getNodeToNode();
	Connectivity *getNodeToElem() { return nodeToElem; }
	Connectivity* getNodeToNode_sommer() { return nodeToNode_sommer;}
	ConstrainedDSA *makeCDSA(int nbc, BCond *bcs);
	Elemset* getEset() { return &packedEset; }
	void setNumnodes(int n) { numnodes = n; }

	void addNode(int nd, double *xyz) { nodes.nodeadd(nd,xyz); }

	// PJSA: mpc stuff
protected:
	bool haveNodes;
public:
	int getNumLMPC() { return numLMPC; }
	void setNumLMPC(int n) { numLMPC = n; }
	LMPCons* getLMPC(int i) { return lmpc[i]; }
	LMPCons* getFsi(int i) { return fsi[i]; }
	ResizeArray<LMPCons *>* getLMPC() { return &lmpc; }
	int addNodalCTC(int n1, int n2, double nx, double ny, double nz,
	                double normalGap = 0.0, int _mode = -1, int lagrangeMult = -1, double penalty = 0.0);
	int getNumCTC();
	void addDirichletLMPCs(int _numDirichlet, BCond *_dbc);
	void deleteAllLMPCs();
	void deleteSomeLMPCs(mpc::ConstraintSource s);
	void UpdateContactSurfaceElements(GeomState *);

	// HB: mortar stuff (EXPERIMENTAL)
protected:
	int nSurfEntity;                          // I know, this should be in GeoSource !!!
	ResizeArray<SurfaceEntity*> SurfEntities; //
	int nMortarLMPCs;                         // total number of Mortar LMPCs generated
	Connectivity* mortarToMPC;                //
	std::vector<int> contactSurfElems;
	int maxContactSurfElems;
	std::set<int> aeroEmbeddedSurfaceId;  //KW: Ids of wet surfaces
public:
	int AddSurfaceEntity(SurfaceEntity*);
	int AddSurfaceEntity(SurfaceEntity*, int isurf);
	void PrintSurfaceEntities();

	int AddAeroEmbedSurfaceId(int Id);
	std::set<int> & GetAeroEmbedSurfaceId() { return aeroEmbeddedSurfaceId; }

	int nMortarCond;
	int nContactSurfacePairs;
	ResizeArray<MortarHandler*> MortarConds;

	int AddMortarCond(MortarHandler*);
	void PrintMortarConds();
	void DeleteMortarConds();

	void SetMortarPairing();
	void SetUpSurfaces(CoordSet* cs = 0);
	void UpdateSurfaceTopology(int numSub = 0, SubDomain **sd = 0);
	void UpdateSurfaces(GeomState *, int config_type = 1);
	void UpdateSurfaces(DistrGeomState *geomState, int config_type, SubDomain **sd);
	void MakeNodalMass(SparseMatrix *M, SparseMatrix *Mcc);
	void MakeNodalMass(SubDOp *M, SubDomain **sd);

	void InitializeDynamicContactSearch(int numSub = 0, SubDomain **sd = 0);
	void PerformDynamicContactSearch(double dt_old, double dt);
	void AddContactForces(double dt_old, double dt, Vector &f);
	void AddContactForces(double dt_old, double dt, DistrVector &f);

	void InitializeStaticContactSearch(MortarHandler::Interaction_Type t, int numSub = 0, SubDomain **sd = 0);
	void ReInitializeStaticContactSearch(MortarHandler::Interaction_Type t, int numSub = 0, SubDomain **sd = 0);
	void UpdateSurfaces(MortarHandler::Interaction_Type t, GeomState *geomState);
	void UpdateSurfaces(MortarHandler::Interaction_Type t, DistrGeomState *geomState, SubDomain **sd);
	void PerformStaticContactSearch(MortarHandler::Interaction_Type t);
	void ExpComputeMortarLMPC(MortarHandler::Interaction_Type t, int nDofs = 0, int* Dofs = 0);

	void ComputeMortarLMPC(int nDofs = 0, int* Dofs = 0);
	void computeMatchingWetInterfaceLMPC();

	void CreateMortarToMPC();

	Connectivity* GetMortarToMPC();
	int GetnMortarConds() { return nMortarCond; }
	int GetnContactSurfacePairs() { return nContactSurfacePairs; }
	int GetnMortarLMPCs();
	MortarHandler* GetMortarCond(int i) { return MortarConds[i]; }
	ResizeArray<SurfaceEntity*>* viewSurfEntities() { return(&SurfEntities); }
	SurfaceEntity* GetSurfaceEntity(int i) { return SurfEntities[i]; }
	void setNumSurfs(int nSurfs) { nSurfEntity = nSurfs; }
	int getNumSurfs() { return(nSurfEntity); }

#ifdef HB_ACME_FFI_DEBUG
	void WriteFFITopFile(FILE* file);
#endif
protected:
	void initialize();
	void initializeNumbers();

	// FETI-DPH acoustics
	template<class Scalar>
	void assembleSommer(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops = 0);
	template<class Scalar>
	void computeSommerDerivatives(double HH, double KK, int curvatureFlag, int *dofs, FullSquareMatrix &ms,
	                              DComplex **bt2nMatrix, double kappa, double ss, AllOps<Scalar> *ops);
	template<class Scalar>
	void assembleATDROB(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops = 0, double Kcoef = 0.0);
	template<class Scalar>
	void updateMatrices(AllOps<Scalar> *ops, GenSparseMatrix<Scalar> *K, int *dofs, int *dofs_mdds,
	                    FullSquareMatrix *reEl, FullSquareMatrix *imEl,double Kcoef = 0.0);
	template<class Scalar>
	void updateDampingMatrices(AllOps<Scalar> *ops, int *dofs, FullSquareMatrix *reEl,
	                           FullSquareMatrix *imEl, double ss, int n);
	struct WetInterface {
		int fluidSurfaceID;
		int structureSurfaceID;
	};
	ResizeArray<WetInterface *> *wetInterfaces;
	int nWetInterface;
	bool firstOutput;

public:
	int isFluidElement(int i);
	int isStructureElement(int i);
	bool isComplex() {
		return ((numComplexDirichlet > 0) || (numComplexNeuman > 0)
		        || (numSommer > 0) || (numComplexLMPC > 0) ||
		        sinfo.hasDamping() || PMLFlag || packedEset.hasDamping() );
	}
	bool isHomogeneous();
	void addWetInterface(int fluidSurfaceID, int structureSurfaceID, double normal_tol = 0.1, double tangential_tol = 0.001);
	int* getAllWetInterfaceNodes(int &count);
	double getFrequencyOrWavenumber();
     void computeAverageProps(int &structure_element_count, int &fluid_element_count, double &global_average_E,
                              double &global_average_nu, double &global_average_rhof);
	void computeCoupledScaleFactors();
	void getInterestingDofs(DofSet &ret, int glNode);

	double** getCMatrix();
	void multCV(const Vector&, Vector&);
	void trMultCV(const Vector&, Vector&);

	ControlLawInfo* getClaw() { return  claw;}
     void setClaw(ControlLawInfo* _claw) { claw = _claw; }

	virtual double densProjCoeff(int dof) { return 1.0; }
	virtual void densProjectStiffness(GenFullSquareMatrix<double>& kel, int num) { /* do nothing */ }
	virtual void densProjectStiffnessC(GenFullSquareMatrix<DComplex>& kel, int num) { /* do nothing */ }
	void transformMatrix(GenFullSquareMatrix<double>& kel, int num);
	void transformMatrixInv(GenFullSquareMatrix<double>& kel, int num);
	void transformVector(Vector &vec, int iele);
	void transformNeumVector(Vector &vec, int iele);
	void transformVector(ComplexVector &vec, int iele);
	void transformElementSensitivityInv(GenFullM<double> *vec, int iele);
	void transformVectorInv(Vector &vec, int iele);
	void transformVectorInv(ComplexVector &vec, int iele);
	void transformVector(double *data, int inode, bool hasRot);
	void transformVector(complex<double> *data, int inode, bool hasRot);
	void transformElementSensitivityInv(double *data, int inode, int numNodes, bool hasRot);
	void transformVectorInv(double *data, int inode, bool hasRot);
	void transformVectorInv(complex<double> *data, int inode, bool hasRot);
	void transformStressStrain(FullM &mat, int iele, OutputInfo::FrameType oframe = OutputInfo::Local);
	void transformStressStrain(FullMC &mat, int iele, OutputInfo::FrameType oframe = OutputInfo::Local);
	void transformMatrix(double *data, int inode, bool sym = true);
	void transformMatrix(complex<double> *data, int inode, bool sym = true);
	void transformMatrixInv(double *data, int inode, bool sym = true);
	void transformMatrixInv(complex<double> *data, int inode, bool sym = true);

	int getMaxNumNodes() { return maxNumNodes; }

	void initSfem();

	ConstrainedDSA *makeMaps(DofSetArray *dsa, ConstrainedDSA *cdsa, DOFMap *baseMap, DOFMap *eqMap);

	/** Replaces all 6-DOF rigid elements that share a node by a single element */
	void collapseRigid6();

	// new controls
	void updateUsddInDbc(double* userDefineDisp, int* map = 0);
	void updateUsdfInNbc(double* userDefineForce, int* map = 0, double* weight = 0);
	void updateActuatorsInNbc(double* actuatorsForce, int* map = 0, double* weight = 0);

	int outFlag; // if this is set to 1 then the output file should use the compressed numbering (i.e. with gaps removed)
	// compatible with top files generated using -T command line argument
	// the default is 0, for compatibility with top files generated using -t command line argument
	int *nodeTable;
	int exactNumNodes;

	void assembleNodalInertiaTensors(FullSquareMatrix *mel);
	std::set<int> & getNewDeletedElements() { return newDeletedElements; }
	std::vector<std::pair<double,int> > & getDeletedElements() { return outDeletedElements; }
#ifdef USE_EIGEN3
protected:
	Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> Jn; // array of nodal inertia tensors for each node including contributions
	// from both elements and DIMASS. Used for nonlinear explicit ROM only
#endif
};

#endif
