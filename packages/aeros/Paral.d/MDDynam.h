#ifndef _MD_DYNAM_DESCR_H_
#define _MD_DYNAM_DESCR_H_

#include <Paral.d/MultiDomainBase.h>
#include <Driver.d/Domain.h>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<double> CuCSparse;
class StaticTimers;
class DistFlExchanger;
class ControlInterface;

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenDistrVectorSet;
typedef GenDistrVectorSet<double> DistrVectorSet;
template <typename Scalar> class GenFullSquareMatrix;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
class Domain;
class Rbm;
template <class Scalar> class MultiDomainRbm;

template <typename T> class SysState;
template <typename T> class GenParallelSolver;
template <typename Scalar> struct AllSensitivities;

class ControlLawInfo;
class Corotator;
class DistrInfo;
class DistrGeomState;
struct SensitivityInfo;

template<class Scalar>
class GenMDDynamMat {
 public:
   union {
     GenParallelSolver<Scalar> *dynMat;
     GenParallelSolver<Scalar> *sysSolver;
   };
   GenSubDOp<Scalar> *spMat;

   GenParallelSolver<Scalar> *prec;       // preconditioner
   GenSubDOp<Scalar> *spp;

   GenParallelSolver<Scalar> *Msolver;
   GenSubDOp<Scalar> *K;
   GenSubDOp<Scalar> *Kuc;
   GenSubDOp<Scalar> *Kcc;
   GenSubDOp<Scalar> *C;
   GenSubDOp<Scalar> *Cuc;
   GenSubDOp<Scalar> *Ccc;
   GenSubDOp<Scalar> *M;
   GenSubDOp<Scalar> *Muc;
   GenSubDOp<Scalar> *Mcc;
   GenSubDOp<Scalar> **C_deriv;
   GenSubDOp<Scalar> **Cuc_deriv;
   int num_K_deriv;
   GenSubDOp<Scalar> **K_deriv;
   GenSubDOp<Scalar> **Kuc_deriv;
   int num_K_arubber;
   GenSubDOp<Scalar> **K_arubber_l;
   GenSubDOp<Scalar> **Kuc_arubber_l;
   GenSubDOp<Scalar> **K_arubber_m;
   GenSubDOp<Scalar> **Kuc_arubber_m;
   Rbm* rigidBodyModes;

   GenMDDynamMat() { dynMat = 0; spMat = 0; prec = 0; spp = 0; Msolver = 0; K = 0; Kuc = 0; Kcc = 0; C = 0; Cuc = 0; Ccc = 0; M = 0; Muc = 0; Mcc = 0; 
                     C_deriv = 0; Cuc_deriv = 0; 
                     K_deriv = 0; Kuc_deriv = 0; num_K_deriv = 0;
                     K_arubber_l = 0; K_arubber_m = 0;
                     Kuc_arubber_l = 0; Kuc_arubber_m = 0;
                     num_K_arubber = 0;
                     rigidBodyModes = 0; };
};

typedef GenMDDynamMat<double> MDDynamMat;

class MultiDomDynPostProcessor 
{
  protected:
    DistFlExchanger *distFlExchanger;
 
    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;
    GenDecDomain<double> *decDomain;
    StaticTimers *times;
    DistrVector *nodalTemps;
    DistrGeomState *geomState;
    Corotator ***allCorot;
    GenFullSquareMatrix<double> **melArray;
    GenDistrVector<double> *reactions;

  public:
    MultiDomDynPostProcessor(DecDomain *d, StaticTimers* _times, DistrGeomState *_geomState = 0,
                             Corotator ***_allCorot = 0, GenFullSquareMatrix<double> **_melArray = 0,
                             GenDistrVector<double> *_reactions = 0) {
      decDomain = d;
      times = _times;
      geomState = _geomState;
      allCorot = _allCorot;
      melArray = _melArray;
      reactions = _reactions;
    }
    MultiDomDynPostProcessor(DecDomain *d, DistFlExchanger *_distFlExchanger, StaticTimers* _times,
                             DistrGeomState *_geomState = 0, Corotator ***_allCorot = 0,
                             GenFullSquareMatrix<double> **_melArray = 0, GenDistrVector<double> *_reactions = 0) {
      decDomain = d;
      distFlExchanger = _distFlExchanger;
      times = _times;
      geomState = _geomState;
      allCorot = _allCorot;
      melArray = _melArray;
      reactions = _reactions;
    }
    void setPostProcessor(DistFlExchanger *);
    void setUserDefs(double **, double **);
    void setNodalTemps(DistrVector*);
    void dynamOutput(int, double, MDDynamMat &, DistrVector &, DistrVector *aeroF, SysState<DistrVector> &, DistrVector *resF = 0);
    double getKineticEnergy(DistrVector & vel, SubDOp * gMass) { return 0.0; }

  private:
    void subUpdateGeomStateUSDD(int isub, double *userDefineDisp, DistrGeomState *geomState,
                                double *userDefineVel, double *userDefineAcc);
};

class MultiDomainDynam : public MultiDomainBase
{
protected:
    DecDomain *decDomain;
    Domain *domain;
    AllSensitivities<double> *allSens;

private:
    CuCSparse **cucs;
    StaticTimers *times;

    // control law data
    ControlInterface *userSupFunc;
    ControlLawInfo *claw;

protected:
    GenFullSquareMatrix<double> **kelArray;
    GenFullSquareMatrix<double> **melArray;
    Corotator ***allCorot;
    DistrGeomState *geomState, *refState;
    MDDynamMat *dynMat;

private:
    MultiDomDynPostProcessor *mddPostPro;

    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;

    // aero data
    DistFlExchanger *distFlExchanger;
    DistrVector *aeroForce;
    DistrVector *prevFrc;
    int prevIndex;
    double prevTime;
    DistrVector *prevFrcBackup;
    int prevIndexBackup;
    double prevTimeBackup;

    // thermoe/thermoh data
    DistrVector* nodalTemps;

    // reaction forces
    DistrVector *reactions;

  public:
    MultiDomainDynam(Domain *d);
    virtual ~MultiDomainDynam();
    MDDynamMat * buildOps(double, double, double);
    MultiDomDynPostProcessor *getPostProcessor();

    const DistrInfo &solVecInfo() const;
    const DistrInfo &masterSolVecInfo() const;
    DistrInfo &bcInfo();

    int getTimeIntegration();
    int getFilterFlag();
    int* boundary();
    double* boundaryValue();
    Domain* getDomain();
    AllSensitivities<double> *getAllSensitivities() { return allSens; }
    SensitivityInfo *getSensitivityInfo() { return domain->senInfo; }
    int getNumSensitivities() { return domain->getNumSensitivities(); }
    void getTimes(double &dt, double &t);
    void getNewMarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);
    void getQuasiStaticParameters(double &maxVel, double &delta);
    void getInitState(SysState<DistrVector> &);
    void printFullNorm(DistrVector &){};
    void getInitialTime(int &tIndex, double &initialTime);
    double getInitialForceNorm();
    void getSensitivityStateParam(double &sensitivityTol,double &ratioSensitivityTol);
    void getSteadyStateParam(int &steadyFlag, int &steadyMin, int &steadMax,
                             double &steadyTol); 
    void getConstForce(DistrVector &);
    void addConstForceSensitivity(DistrVector &);
    void getContactForce(DistrVector &d_n, DistrVector &dinc, DistrVector &ctc_f, double t_n_p, double dt, double dt_old);
    void computeExtForce2(SysState<DistrVector> &, DistrVector &, 
                          DistrVector &, int tIndex, double t,
                          DistrVector * aero_f=0,
                          double gamma=0.5, double alphaf=0.5);
    void getGravityForce(DistrVector &);
    void getUnamplifiedExtForce(DistrVector &, int);
    void getAeroelasticForceSensitivity(int t_index, double t, DistrVector * aero_f=0, double gamma=0.5, double alphaf=0.5);

    void getRHS(DistrVector &);
    void preProcess();
    void preProcessSA();
    void postProcessSA(MDDynamMat *, DistrVector &sol);
    void sensitivityPostProcessing(DistrVector *sol);
    void processLastOutput();
    void printTimers(MDDynamMat *, double);

    void getRayleighCoef(double& alpha);

    SubDOp* getpK(MDDynamMat* dynOps);
    SubDOp* getpM(MDDynamMat* dynOps);
    SubDOp* getpC(MDDynamMat* dynOps);

    // Central Difference only related subroutines
    virtual void computeStabilityTimeStep(double&, MDDynamMat&);

    void updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n);
    void pull_back(DistrVector& f);
    void push_forward(DistrVector& a);

    // Mode Decomposition parameters and subroutines
    int getModeDecompFlag();
    void modeDecompPreProcess(SparseMatrix* M);
    void modeDecomp(double t, int tIndex, DistrVector& d_n);

    void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex);
    void getFollowerForce(DistrVector &f, double t, int tIndex);

    // Aeroelastic problems related subroutines
    void computeTimeInfo();
    int aeroPreProcess(DistrVector &, DistrVector &, DistrVector &, DistrVector &); 
    int aeroSensitivityPreProcess(DistrVector &, DistrVector &, DistrVector &, DistrVector &); 
    int sendDisplacements(DistrVector &, DistrVector &, DistrVector &, DistrVector &); 
    void sendNumParam(int numParam, int actvar, double steadyTol) {}
    void getNumParam(int &numParam) {}
    void sendRelativeResidual(double relres) {}
    int cmdCom(int cmdFlag);
    int getAeroAlg();
    void aeroSend(double time, DistrVector& d, DistrVector& v, DistrVector& a, DistrVector& v_p);
    void a5TimeLoopCheck(int&, double&, double);
    void a5StatusRevise(int, SysState<DistrVector>&, SysState<DistrVector>&);

    // Thermoelastic problems related subroutines
    void thermoePreProcess(DistrVector&, DistrVector&, DistrVector&);
    int getThermoeFlag();
    void thermohPreProcess(DistrVector&, DistrVector&, DistrVector&);
    int getThermohFlag();

    // Aeroheat
    void aeroHeatPreProcess(DistrVector&, DistrVector&, DistrVector&);
    int getAeroheatFlag();

    // Non-linear quasi-static
    void solveAndUpdate(DistrVector &force, DistrVector &dinc, DistrVector &d, double relaxFac, double time);
   
  private:
    void subGetInternalForce(int isub, DistrVector &res, double &t, int &tIndex);
    void subGetKtimesU(int isub, DistrVector &d, DistrVector &f);
    void subGetFollowerForce(int isub, DistrVector &res, double &t, int &tIndex);
    void makeSubCorotators(int isub);
    void makeSubElementArrays(int isub);
    void initSubPrescribedDisplacement(int isub);
    void subUpdateGeomStateUSDD(int isub, double *userDefineDisp, DistrGeomState *geomState,
                                double *userDefineVel, double *userDefineAcc);
    void subUpdateUsrDefDispsAndVels(int isub, double *userDefineDisp, double *userDefineVel);
    void subExplicitUpdate(int isub, DistrVector &d, DistrGeomState *geomState);
    void subGetGravityForce(int isub, DistrVector &);
    void subGetUnamplifiedExtForce(int isub, DistrVector &, int);
    void subUpdateStates(int isub, DistrGeomState *refState, DistrGeomState *geomState, double time);
};

#endif
