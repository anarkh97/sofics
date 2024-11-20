#ifndef _MD_NL_DYNAMIC_H_
#define _MD_NL_DYNAMIC_H_

#include <Paral.d/MultiDomainBase.h>

#include <cstddef>
#include <map>
#include <vector>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif
#include <Rom.d/ModalGeomState.h>

class Domain;
template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class Scalar> class GenParallelSolver;
typedef GenParallelSolver<double> ParallelSolver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class Corotator;
class StaticTimers;
class ControlInterface;
class ControlLawInfo;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class DistrGeomState;
template <class Scalar> class GenMultiDomainPostProcessor;
typedef GenMultiDomainPostProcessor<double> MultiDomainPostProcessor;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
template <class Scalar> class GenMDDynamMat;
typedef GenMDDynamMat<double> MDDynamMat;
class DistrInfo;
class DistFlExchanger;
class LinesearchInfo;
template <class Scalar> class GenFetiDPSolver;
template <typename Scalar> struct AllSensitivities;
struct SensitivityInfo;

// Multiple Domain Nonlinear Dynamic problem descriptor

class MDNLDynamic : public MultiDomainBase
{
  protected:
    Domain *domain;
    DecDomain  *decDomain;
    ParallelSolver *solver;
    GenFetiDPSolver<double> *fetiSolver;
    int     totIter;             // counter of iterations

    MDDynamMat *allOps;
    SubDOp *M;
    SubDOp *C;
    SubDOp *Kuc, *Muc, *Cuc;
    SubDOp *Mcc, *Ccc;
    Corotator ***allCorot;       // element corotators per subdomain
    DistrVector *localTemp;

    FullSquareMatrix **kelArray; // element stiffness matrices per subdomain
    FullSquareMatrix **melArray; // element mass matrices per subdomain
    FullSquareMatrix **celArray; // element damping matrices per subdomain

    DistrVector *reactions;

    double t0;                   // initial time
    double totalTime;		 // total simulation time
    double dt0;			 // initial time step size
    double remainder;		 // remaining time after max steps
    int    maxStep;		 // maximum number of time steps

    double tolerance;		 // newton convergence criteria tolerance
    double firstRes;		 // first iteration residual norm
    double secondRes;            // second iteration residual norm
    double firstDv;		 // first iteration displacement increment norm
    double firstEng;             // first iteration energy increment norm
    double externalForceNorm;    // current external force norm

    int numSystems;              // number of linear systems solved

    ControlInterface *userSupFunc;
    ControlLawInfo *claw;
    int **clawDofs;

    StaticTimers *times;	 // timing information

    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;

    // aero data
    DistFlExchanger *distFlExchanger;
    DistrVector *prevFrc;
    int prevIndex;
    double prevTime;
    DistrVector *aeroForce;
    DistrVector *nodalTemps;

    double Kcoef_p;

    std::map<std::pair<int,int>, double> *mu, *muCopy; // lagrange multipliers for the contact surfaces
    std::vector<double> *lambda; // lagrange multipliers for all the other constraints

    bool updateCS;
    AllSensitivities<double> *allSens;

#ifdef USE_EIGEN3
    Eigen::MatrixXd VtMV;
#endif

 public:

    // Constructor
    MDNLDynamic(Domain *d);
    virtual ~MDNLDynamic();

    void clean();

    MultiDomainPostProcessor *getPostProcessor();
    void getInitialTime(int &initTimeIndex, double &initTime);
    void readRestartFile(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                         DistrVector &v_p, DistrGeomState &geomState); 

    int getInitState(DistrVector &d, DistrVector& v, DistrVector &a, DistrVector &v_p);
    void updatePrescribedDisplacement(DistrGeomState *geomState);

    DistrInfo& solVecInfo();
    DistrInfo& elemVecInfo();
    DistrInfo& sysVecInfo();

    double getTolerance();

    void computeTimeInfo();

    double getDelta() const     { return dt0/2; }
    double getDt() const        { return dt0;   }
    double getLastDt() const    { return remainder; }
    int    getMaxStep() const   { return maxStep;   }
    double getTotalTime() const { return totalTime; }
    int    getMaxit();
    double getDeltaLambda();

    void getConstForce(DistrVector &gravityForce);
    void getExternalForce(DistrVector &externalForce, DistrVector &constantForce,
                          int tIndex, double time, DistrGeomState *geomState, 
                          DistrVector &elementInternalForce, DistrVector &aeroF, double localDelta);

    void getIncDisplacement(DistrGeomState *geomState, DistrVector &du, DistrGeomState *refState,
                            bool zeroRot);

    double formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                           DistrVector& residual, DistrVector& rhs, DistrGeomState *geomState, double localDelta);

    void formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual,
                          DistrVector& rhs, DistrGeomState &, double mid, double localDelta);

    void formRHSinitializer(DistrVector &fext, DistrVector &velocity, DistrVector &elementInternalForce,
                            DistrGeomState &geomState, DistrVector &rhs, DistrGeomState *refState = NULL);

    virtual void preProcess();

    void processLastOutput();
    ParallelSolver *getSolver();

    DistrGeomState* createGeomState();
    DistrGeomState* copyGeomState(DistrGeomState* geomState);

    int getNumStages();
    int checkConvergence(int iter, double rN, DistrVector& residual, DistrVector& dv, double time);
    void updateContactSurfaces(DistrGeomState& geomState, DistrGeomState *refState);
    void updateStates(DistrGeomState *refState, DistrGeomState& geomState, double time);

    // getStiffAndForce forms element stiffness matrices and/or
    // returns the residual force = external - internal forces
    virtual double getStiffAndForce(DistrGeomState& geomState, DistrVector& residual, DistrVector& elementInternalForce,
                            double midtime=-1, DistrGeomState *refState = NULL, bool forceOnly = false);

    // reBuild assembles new dynamic stiffness matrix
    void reBuild(DistrGeomState& geomState, int iter, double localDelta, double t);

    void printTimers(double timeLoop);
    void dynamOutput(DistrGeomState* geomState, DistrVector& velocity, DistrVector &vp,
                     double time, int timestep, DistrVector& force, DistrVector &aeroF, DistrVector &acceleration,
                     DistrGeomState *refState);
    void dynamCommToFluid(DistrGeomState* geomState, DistrGeomState* bkGeomState,
                          DistrVector& velocity, DistrVector& bkVelocity,
                          DistrVector& vp, DistrVector& bkVp, int step, int parity,
                          int aeroAlg, double time);

    void getConstraintMultipliers(DistrGeomState &geomState);
    double getResidualNorm(DistrVector &vec, DistrGeomState &, double);

    int getAeroAlg();
    int getThermoeFlag();
    int getThermohFlag();
    int getAeroheatFlag();
    void getNewmarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);

    void initializeParameters(int step, DistrGeomState *geomState);
    void updateParameters(DistrGeomState *geomState);
    bool checkConstraintViolation(double &err, DistrGeomState *geomState);

    LinesearchInfo& linesearch();
    bool getResizeFlag();
    void resize(DistrGeomState *refState, DistrGeomState *geomState, DistrGeomState *stepState, DistrVector *stateIncr,
                DistrVector &v, DistrVector &a, DistrVector &vp, DistrVector &force);

    void sensitivityAnalysis(DistrGeomState*, DistrGeomState*) {}
    void sensitivityAnalysis(DistrModalGeomState*, DistrModalGeomState *){}
    void preProcessSA() {}
    void postProcessNLSA(DistrGeomState*, DistrGeomState*) {}
    AllSensitivities<double> *getAllSensitivities() { return allSens; }
    SensitivityInfo *getSensitivityInfo();
    int getNumSensitivities();

  protected:
    Domain *getDomain() { return domain; }
    DecDomain *getDecDomain() { return decDomain; }

  private:
    void makeSubDofs(int isub);
    void makeSubCorotators(int isub);
    void makeSubElementArrays(int isub);
    void subGetExternalForce(int isub, DistrVector& f, DistrVector& constantForce, double tf, double tm);
    void subGetStiffAndForce(int isub, DistrGeomState &geomState, DistrVector &res,
                             DistrVector &elemIntForce, double t, DistrGeomState *refState, bool forceOnly);
    void subUpdatePrescribedDisplacement(int isub, DistrGeomState& geomState);
    void addConstraintForces(int isub, DistrVector& rhs, DistrGeomState &geomState);
    void subGetConstraintMultipliers(int isub, DistrGeomState &geomState);
    void subUpdateGeomStateUSDD(int isub, DistrGeomState *geomState, double *userDefineDisp,
                                double *userDefineVel, double *userDefineAcc);
    void makeSubClawDofs(int isub);
    void subKucTransposeMultSubtractClaw(int iSub, DistrVector& residual, double *userDefineDisp);
    void subExtractControlDisp(int isub, DistrGeomState &geomState, double *ctrdsp);
    int aeroPreProcess(DistrVector &, DistrVector &, DistrVector &, DistrVector &, DistrGeomState *geomState = 0);
    void thermoePreProcess();
    void thermohPreProcess(DistrVector &);
    void aeroheatPreProcess(DistrVector &, DistrVector &, DistrVector &);
    void subDynamCommToFluid(int isub, DistrVector& v, DistrGeomState* distrGeomState,
                             DistrGeomState* bkDistrGeomState, int parity, int aeroAlg);
    void subDynamCommToFluidAeroheat(int isub, DistrVector& v, DistrGeomState* distrGeomState);
    void subUpdateStates(int isub, DistrGeomState *refState, DistrGeomState *geomState, double time);
    void subReadRestartFile(int i, DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                            DistrVector &v_p, DistrGeomState &geomState);
    void subWriteRestartFile(int i, double &t, int &index, DistrVector &vel_n, DistrVector &acc_n, DistrGeomState &geomState);
    void subGetReactionForce(int i, DistrGeomState &geomState, DistrGeomState &refState, DistrVector &vel_n,
                             DistrVector &acc_n, double &time);

    virtual bool factorWhenBuilding() const;
    void deleteSubCorotators(int isub);
    void deleteSubElementArrays(int isub);

    void subInitializeMultipliers(int isub, DistrGeomState& geomState);
    void subInitializeParameters(int isub);
    void subUpdateMultipliers(int isub, DistrGeomState& geomState);
    void subUpdateParameters(int isub);
};

#endif
