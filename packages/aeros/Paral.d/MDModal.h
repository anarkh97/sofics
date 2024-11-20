#ifndef _MD_MODAL_H_
#define _MD_MODAL_H_

#include <Problems.d/ModalBase.h>
#include <Paral.d/MDModalBase.h>
#include <Driver.d/Domain.h>

#include <Hetero.d/DistFlExchange.h>

template<class VecType> class SysState;
template <typename Scalar> struct AllSensitivities;
struct SensitivityInfo;

class MultiDomainModal : public MDModalBase {
/* problem descriptor for linear dynamics problems using the eigenmdoes of
     non-zero frequency as the solution space
*/
private:

  ModalOps &modalOps;
  double Wdmp;
  Vector *previousCq;
  Vector *previousDisp;

public:

  MultiDomainModal(Domain *d);
  ~MultiDomainModal();

  void projectForce(DistrVector &fullV, Vector& modalV);
  void expand(const Vector &modalV, DistrVector& fullV);

  void preProcess();
  void preProcessSA();
  void postProcessSA(ModalOps *,Vector &sol);
  void processLastOutput();

  MultiDomainModal* getPostProcessor() { return this; }
  void getTimes(double &dt, double &tmax);
  void printFullNorm(Vector &){};
  void getInitState(SysState<Vector> &state);

  double getInitialForceNorm();

  int solVecInfo() { return numModes; }
  int masterSolVecInfo() { return numModes; }
  int bcInfo() {return 0;}

  void getConstForce(Vector &constF);
  void addConstForceSensitivity(Vector &gravityForceSen);
  void getSteadyStateParam(int &flag, int &min, int &max, double &tol);
  void getSensitivityStateParam(double &tol, double &ratioTol);
  int getTimeIntegration() { return domain->solInfo().timeIntegration; }
  void getNewMarkParameters(double &beta, double &gamma,
    double &alphaf, double  &alpham);
  Domain *getDomain() { return domain; };
  AllSensitivities<double> *getAllSensitivities();
  void getInitialTime(int &tIndex, double &time);
  void getRayleighCoef(double &alpha) { alpha =  domain->solInfo().alphaDamp; }

  ModalOps* buildOps(double kcoef, double ccoef, double mcoef);

  void computeStabilityTimeStep(double &dt, ModalOps &);
  void updateState(double, Vector &, Vector &){ /* leave blank */ }
  void push_forward(Vector &){ /* leave blank */ }
  void pull_back(Vector &){ /* leave blank */ }
  void getQuasiStaticParameters(double &maxVel, double &delta);
  SensitivityInfo *getSensitivityInfo() { return domain->senInfo; }
  int getNumSensitivities() { return domain->getNumSensitivities(); }
  int getFilterFlag() { return domain->solInfo().filterFlags; }
  void project(Vector &v) { /* leave blank */ }

  void getContactForce(Vector&, Vector&, Vector &extF, double, double, double) { extF.zero(); };
  void computeExtForce2(SysState<Vector>& state, Vector &extF,
                        Vector &constF, int tIndex, double time, Vector *aeroF = 0,
                        double gamma = 0.5, double alphaf = 0.5);
  void getAeroelasticForceSensitivity(int t_index, double t, Vector * aero_f=0, double gamma=0.5, double alphaf=0.5);

  void getInternalForce(Vector &d, Vector &f, double t, int tIndex);
  double getElasticEnergy(Vector &d);
  double getKineticEnergy(Vector &v);
  double getDampingEnergy(Vector &d, Vector &v, double time);

  void printTimers(ModalOps *, double) { /* leave blank */}

  int getModeDecompFlag(){ return domain->solInfo().modeDecompFlag; }
  void modeDecomp(double time, int tIndex, Vector& dsp){ /* leave blank */ }

  void dynamOutput(int tIndex, double time, ModalOps &mOps, Vector &extF,
                   Vector *aeroF, SysState<Vector> &state);

  DistFlExchanger *flExchanger;

  void computeTimeInfo() {};
  int aeroPreProcess(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
  int aeroSensitivityPreProcess(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
  int sendDisplacements(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p);
  void sendNumParam(int numParam, int actvar, double steadyTol) {}
  void getNumParam(int &numParam) {}
  void sendRelativeResidual(double relres) {}
  int cmdCom(int cmdFlag) { return -1; }

  void a5TimeLoopCheck(int& parity, double& t, double dt);
  void a5StatusRevise(int parity, SysState<Vector>& curState, SysState<Vector>& bkState);
  int getAeroAlg() { return domain->solInfo().aeroFlag; }
  void aeroSend(double time, Vector& d, Vector& v, Vector& a, Vector& v_p) {}

  void thermoePreProcess(Vector& d_n, Vector& v_n, Vector& v_p) {}
  int getThermoeFlag() { return -1; }
  void thermohPreProcess(Vector& d_n, Vector& v_n, Vector& v_p) {}
  int getThermohFlag() { return domain->solInfo().thermohFlag; }

  void aeroHeatPreProcess(Vector& d_n, Vector& v_n, Vector& v_p) {}
  int getAeroheatFlag() { return -1; }

  void solveAndUpdate(Vector &, Vector &, Vector &, double, double){ /* leave blank */ }

  void sensitivityPostProcessing(Vector *){/*leave blank */}
};

#endif
