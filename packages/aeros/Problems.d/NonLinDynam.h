#ifndef _NON_LIN_DYNAM_H_
#define _NON_LIN_DYNAM_H_

#include <Problems.d/SingleDomainBase.h>
#include <Math.d/Vector.h>
#include <Utils.d/OutputInfo.h>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

class Domain;
class Rbm;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class StaticTimers;
class GeomState;
class Corotator;
class ControlInterface;
class SDDynamPostProcessor;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <typename T> class SysState;
class PrevFrc;
template <typename T> struct AllOps;
class ControlLawInfo;
template <typename T> class GenFSFullMatrix;
typedef GenFSFullMatrix<double> FSFullMatrix;
class LinesearchInfo;
struct SensitivityInfo;
template<class Scalar> struct AllSensitivities;

class NLDynamPostProcessor
{
public:
  virtual ~NLDynamPostProcessor() {}
  virtual void dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState, 
                                Vector& velocity, Vector& bkVelocity,
                                Vector& vp, Vector& bkVp, int step, int parity,
                                int aeroAlg, double time) = 0;
  virtual void dynamOutput(GeomState *, GenVector<double> &, GenVector<double> &, double, int, GenVector<double> &, GenVector<double> &,
                           GenVector<double> &, GeomState *) const = 0;
};

// Virtual methods to allow derived class PitaNonLinDynamic in Pita.d/PitaNonLinDynam.d
class NonLinDynamic : public NLDynamPostProcessor, public SingleDomainBase {

  protected:
    Domain *domain;
    double *bcx;	// displacement prescribed values
    double *vcx;        // velocity prescribed values
    double *acx;        // acceleration prescribed values
    Solver *solver;
    SparseMatrix *spm;
    Solver *prec;
    SparseMatrix *spp;
    FILE   *res;        // file to store residuals
    int     totIter;     // counter of iterations
    /*int    *dofTypeArray; */ 

    int *clawDofs;      // array containing cdsa dof numbers for usdd nodes

    SparseMatrix *K;    // Stiffness matrix
    SparseMatrix *M;    // Mass matrix
    SparseMatrix *C;    // Damping matrix
    SparseMatrix *Kuc;
    SparseMatrix *Muc, *Mcc, *Cuc, *Ccc;
    Corotator **allCorot;
    Vector localTemp;

    FullSquareMatrix *kelArray; // array of element stiffness matrices
    FullSquareMatrix *celArray; // array of element damping matrices
    FullSquareMatrix *melArray; // array of element mass matrices

    PrevFrc *prevFrc;   // previous Aeroelastic force at time step t(n-1)
    double t0;          // initial time
    double totalTime;   // total time
    double dt0;         // initial time step size
    int maxStep;        // maximum number of time steps

    double tolerance;   // convergence criteria tolerance
    double firstRes;    // first iteration residual norm
    double secondRes;    // second iteration residual norm
    double firstDv;     // first iteration displacement increment norm
    double firstEng;    // first iteration energy increment norm
    double firstForceNorm;
    double firstMomenNorm;
    double externalForceNorm; //current external force norm

    int numSystems;     // number of linear systems solved

    StaticTimers *times; // Timing information class

    ControlInterface *userSupFunc;
    ControlLawInfo *claw;
    AllSensitivities<double> *allSens;

    void extractControlData(Vector& d_n, Vector& v_n, Vector& a_n,
                            double* ctrdis, double* ctrvel, double* ctracc);
    void extractControlDisp(GeomState *, double *);

    void addCtrl(Vector& externalForce, double *controlForce);
    void addUserForce(Vector & externalForce, double *userForce);

    double resN;
    Vector *reactions;
    bool factor;
    bool updateCS;
#ifdef USE_EIGEN3
    Eigen::MatrixXd VtMV;
#endif

 public:
    // Constructor
    NonLinDynamic(Domain *d);
    virtual ~NonLinDynamic();

    SDDynamPostProcessor * getPostProcessor();
    virtual const NLDynamPostProcessor & defaultPostProcessor() const;
    void getInitialTime(int &initTimeIndex, double &initTime);
    void readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                         Vector &v_p, GeomState &geomState);
    void setBC(double *userDefineDisplacement, double *userDefineVel, double *userDefineAcc);

    virtual int getInitState(Vector& d, Vector& v, Vector& a, Vector &v_p);
    void updateUserSuppliedFunction(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p, double initialTime);
    void updatePrescribedDisplacement(GeomState *geomState);

    virtual int solVecInfo() const;
    int  sysVecInfo();
    int  elemVecInfo();

    double getTolerance();

    void   computeTimeInfo();

    double getDelta() const     { return dt0/2;     }
    double getDt() const        { return dt0;        }
    int    getMaxStep() const   { return maxStep;   }
    double getTotalTime() const { return totalTime; }
    int    getMaxit();
    double getDeltaLambda();

    bool getZeroRot() const;

    virtual void getConstForce(Vector& constantForce);

    void getExternalForce(Vector& externalForce, Vector& constantForce, int tIndex, double time,
                          GeomState* geomState, Vector& elementInternalForce, Vector& aeroF, double localDelta);

    void getIncDisplacement(GeomState *geomState, Vector &du, GeomState *refState, bool zeroRot);

    double formRHScorrector(Vector& inc_displac, Vector &velocity, Vector& acceleration,
                            Vector &residual, Vector &rhs, GeomState *geomState, double localDelta);

    void formRHSpredictor(Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs, GeomState &, double mid, double localDelta);

    void formRHSinitializer(Vector &fext, Vector &velocity, Vector &elementInternalForce, GeomState &geomState, Vector &rhs, GeomState *refState = NULL);

    virtual void preProcess(double Kcoef = 0, double Mcoef = 1, double Ccoef = 0);

    virtual void openResidualFile();

    void processLastOutput();
    Solver *getSolver();

    const SparseMatrix* getMassMatrix() { return M; }

    GeomState* createGeomState();
    GeomState* copyGeomState(GeomState* geomState);

    int getNumStages();
    int checkConvergence(int iter, double rN, Vector& residual, Vector& dv, double time);

    void updateContactSurfaces(GeomState& geomState, GeomState *refState);
    virtual void updateStates(GeomState *refState, GeomState& geomState, double time);

    // getStiffAndForce forms element stiffness matrices and/or
    // returns the residual force = external - internal forces
    double getStiffAndForce(GeomState& geomState, Vector& residual, Vector& elementInternalForce,
                            double midtime=-1, GeomState *refState = NULL, bool forceOnly = false);

    AllSensitivities<double> *getAllSensitivities() { return allSens; }
    SensitivityInfo *getSensitivityInfo();
    int getNumSensitivities();

  private:
    // Overridable implementation of getStiffAndForce
    virtual void getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                            Corotator **allCorot, FullSquareMatrix *kelArray,
                                            Vector &residual, double lambda, double time, GeomState *refState,
                                            FullSquareMatrix *melArray, bool forceOnly);

  public:
    // reBuild assembles new dynamic stiffness matrix
    void reBuild(GeomState& geomState, int iter, double localDelta, double t);

    void printTimers(double timeLoop);
    void dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState,
                          Vector& velocity, Vector& bkVelocity,
                          Vector& vp, Vector& bkVp, int step, int parity,
                          int aeroAlg, double time);
    void dynamOutput(GeomState* geomState, Vector& velocity, Vector &vp,
                     double time, int timestep, Vector& force, Vector &aeroF, Vector &acceleration,
                     GeomState *refState) const;
    void getConstraintMultipliers(GeomState &geomState) { /* deliberately empty */ }
    double getResidualNorm(const Vector &rhs, GeomState &geomState, double localDelta);

    int getAeroAlg();
    int getThermoeFlag();
    int getThermohFlag();
    int getAeroheatFlag();
    void getNewmarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);

    void initializeParameters(int step, GeomState *geomState);
    void updateParameters(GeomState *geomState);
    bool checkConstraintViolation(double &err, GeomState *geomState);

    LinesearchInfo& linesearch();
    bool getResizeFlag();
    void resize(GeomState *refState, GeomState *geomState, GeomState *stepState, Vector *stateIncr,
                Vector &v, Vector &a, Vector &vp, Vector &force) {} // XXX

    void preProcessSA();
    void postProcessSA(Vector &sol);

    void postProcessNLSA(GeomState *, GeomState *);
    void sensitivityAnalysis(GeomState *, GeomState *);
private:
    virtual bool factorWhenBuilding() const;
    void clean();
    void computeLambdaNLStressVM(int); 
    void computeLambdaAggregatedStress();
    void computeLambdaDisp(int , int);
    void setUpPODSolver(OutputInfo::Type);
};

inline const NLDynamPostProcessor &
NonLinDynamic::defaultPostProcessor() const
{
  return *this;
}

class DummyNLDynamPostProcessor {
public:
  virtual void dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState,
                          Vector& velocity, Vector& bkVelocity,
                          Vector& vp, Vector& bkVp, int step, int parity,
                          int aeroAlg){}
  virtual void dynamOutput(GeomState *, GenVector<double> &, GenVector<double> &, double, int, GenVector<double> &, GenVector<double> &,
                           GenVector<double> &) const {}
};

#endif
