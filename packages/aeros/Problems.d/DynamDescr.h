#ifndef _DYNAM_DESCR_H_
#define _DYNAM_DESCR_H_

#include <map>
#include <Driver.d/Domain.h>
#include <Problems.d/SingleDomainBase.h>

template <class Scalar> class GenDynamMat;
typedef GenDynamMat<double> DynamMat;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
class ControlInterface;
class StaticTimers;
class GeomState;
class Corotator;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenFSFullMatrix;
typedef GenFSFullMatrix<double> FSFullMatrix;
class PrevFrc;
class FlExchanger;
template <typename T> class SysState;
class Domain;
class ControlLawInfo;
template <typename Scalar> class GenSolver;
template <typename Scalar> struct AllSensitivities;
struct SensitivityInfo;
// Single Domain Dynamic Post Processor Class

class SDDynamPostProcessor {
  protected:
    Domain *domain;
    double *bcx;
    double *vcx, *acx;
    StaticTimers *times;
    GeomState *geomState;
    Corotator **allCorot;
    FullSquareMatrix *melArray;
    Vector *reactions;
    Vector *dummy;

  public:
    SDDynamPostProcessor(Domain *d, double *_bcx, double *_vcx, double *_acx,
                         StaticTimers *_times, GeomState *_geomState = 0,
                         Corotator **_allCorot = 0, FullSquareMatrix *_melArray = 0,
                         Vector *_reactions = 0);

    ~SDDynamPostProcessor();

    // File management
    void openOutputFiles();
    void openOutputFilesForPita(int timeSliceRank);
    void closeOutputFiles();
    void closeOutputFilesForPita(int timeSliceRank);

    // Perform output (linear dynamics only)
    void dynamOutput(int timeStepIndex, double time, DynamMat & dMat,
                     Vector & externalForce, Vector * aeroForce,
                     SysState<Vector> & systemState);
    void pitaDynamOutput(int timeStepIndex, DynamMat & dMat,
                         Vector & externalForce, Vector * aeroForce,
                         SysState<Vector> & systemState,
                         int sliceRank, double time);
    double getKineticEnergy(Vector & vel, SparseMatrix * gMass);

  private:
    void fillBcxVcx(double time);
};


// Single Domain Dynamic Problem Description

class SingleDomainDynamic : public SingleDomainBase
{
  protected:
    Domain *domain;
    int    *bc;
    double *bcx;	// displacement bc values
    double *userDefineDisp;
    double *vcx;	// velocity bc values
    double *acx;        // acceleration bc values
    StaticTimers *times;
    int    *nodeOrder;
    SparseMatrix *kuc, *muc, *cuc;
    ControlInterface *userSupFunc;
    ControlLawInfo *claw;
    AllSensitivities<double> *allSens;

    // members for nonlinear eigen problem
    FullSquareMatrix *kelArray, *melArray;
    Corotator **allCorot;
    GeomState *geomState, *refState;
    double t0; // initial time
    PrevFrc *prevFrc;
    PrevFrc *prevFrcBackup;

    double *alfa;
    int maxmode;      // number of eigenvectors for modal decomposition analysis
    double **eigmodes;  
    double **tPhiM;   // Stores Transpose(Phi)*M

    FlExchanger *flExchanger;
    Vector *reactions;
    bool firstSts;

 protected:
    // extract gets an acceleration, velocity and disp vector according to the claw sensors
    void extractControlData(SysState<Vector> &, double *, double *, double*);
    // addCtrl adds the force from the control
    void addCtrl(Vector &f, double *controlForce);
    void addCtrl(Vector &, double *, double *);
    void addUserForce(Vector&f, double *userDefineForce);
    void setBC(double *userDefineDisp, double *userDefineVel, double *userDefineAcc);
 public:
    SingleDomainDynamic(Domain *d); 
    ~SingleDomainDynamic();

    int* boundary() { return bc;}
    double* boundaryValue() { return bcx; }
    double* boundaryVeloc() { return vcx; }
    double* boundaryAccel() { return acx; }

    int solVecInfo();
    int masterSolVecInfo();
    int bcInfo();
    int dbcVecInfo();
    int getTimeIntegration();
    int getFilterFlag();
    
    Domain *getDomain() {return domain;}

    void getTimes(double &dt, double &t);
    void getNewMarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);
    void getQuasiStaticParameters(double &maxVel, double &delta);
    SensitivityInfo *getSensitivityInfo() { return domain->senInfo; }
    int getNumSensitivities() { return domain->getNumSensitivities(); }
    void getRayleighCoef(double &alpha);
    void printFullNorm(Vector &){};
    void getInitState(SysState<Vector> & currentState);
    void getInitialTime(int &tIndex, double &initialTime); 
    double getInitialForceNorm(); 
    void getConstForce(Vector &gravityForce);
    void addConstForceSensitivity(Vector &gravityForceSen);
    void getSteadyStateParam(int &steadyFlag, int &steadyMin, int &steadMax,
                             double &steadyTol); 
    void getSensitivityStateParam(double &sensitivityTol, double &ratioSensitivityTol);

    void getContactForce(Vector &d_n, Vector &dinc, Vector &ctc_f, double d_n_p, double dt, double dt_old);
    void computeExtForce2(SysState<Vector> &, Vector &, Vector &, int t_index,
                         double t, Vector * aero_f=0, double gamma=0.5, double alphaf=0.5, double *pt_dt=0);
    void getAeroelasticForceSensitivity(int t_index, double t, Vector * aero_f=0, double gamma=0.5, double alphaf=0.5);
    void preProcess();
    void preProcessSA();
    void postProcessSA(DynamMat *dMat, Vector &sol);
    void sensitivityPostProcessing(Vector *sol);
    void processLastOutput();
    DynamMat * buildOps(double coeM, double coeC, double coeK);
    GenSolver<double> *getSolver();
    AllSensitivities<double> *getAllSensitivities() { return allSens; }
    SDDynamPostProcessor *getPostProcessor();
    void printTimers(DynamMat *, double);
    double betaDamp() const;
    double alphaDamp() const; 
    void setDamping( double betaDamp, double alphaDamp );
	
    SparseMatrix * getpK(DynamMat * dynOps);
    SparseMatrix * getpM(DynamMat * dynOps);
    SparseMatrix * getpC(DynamMat * dynOps);

    // Central Difference only related subroutines
    void computeStabilityTimeStep(double& dt, DynamMat& dMat);

    void updateState(double dt_n_h, Vector& v_n_h, Vector& d_n);
    void pull_back(Vector& f);
    void push_forward(Vector& a);

    // Mode Decomposition parameters and subroutines
    int getModeDecompFlag();
    void modeDecompPreProcess(SparseMatrix*M);
    void modeDecomp(double t, int tIndex, Vector& d_n);

    void getInternalForce(Vector&, Vector&, double t, int tIndex);
    void getUnassembledNonLinearInternalForce(Vector &, Vector &, std::map<int, std::pair<int,int> > &,
                                              FullSquareMatrix *kelCopy, double t, int tIndex);

    // Aeroelastic problems related subroutines
    void computeTimeInfo();
    int aeroPreProcess(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
    int aeroSensitivityPreProcess(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
    int sendDisplacements(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
    void sendNumParam(int numParam, int actvar, double steadyTol);
    void getNumParam(int &numParam);
    void sendRelativeResidual(double relres);
    int cmdCom(int cmdFlag);
    int getAeroAlg();
    void aeroSend(double time, Vector& d, Vector& v, Vector& a, Vector& v_p);
    void a5TimeLoopCheck(int& parity, double& t, double dt);
    void a5StatusRevise(int parity, SysState<Vector>& curState, SysState<Vector>& bkState);

    // Thermoelastic problems related subroutines
    void thermoePreProcess(Vector& d_n, Vector& v_n, Vector& v_p);
    int getThermoeFlag();
    void thermohPreProcess(Vector& d_n, Vector& v_n, Vector& v_p);
    int getThermohFlag();

    // Aeroheat
    void aeroHeatPreProcess(Vector& d_n, Vector& v_n, Vector& v_p);
    int getAeroheatFlag();

    // Non-linear quasi-static
    void solveAndUpdate(Vector &force, Vector &dinc, Vector &d, double relaxFac, double time);
};

#endif
