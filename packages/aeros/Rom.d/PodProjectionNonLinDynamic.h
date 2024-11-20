#ifndef ROM_PODPROJECTIONNONLINDYNAMIC_H
#define ROM_PODPROJECTIONNONLINDYNAMIC_H

#include <Driver.d/GeoSource.h>
#include <Problems.d/DynamDescr.h>
#include <Problems.d/NonLinDynam.h>
#include <Driver.d/StateUpdater.h>
#include <Rom.d/ModalGeomState.h>
#include <Corotational.d/GeomState.h>
#include "VecBasis.h"
#include "PodProjectionSolver.h"

#include <memory>

namespace Rom {

//template <typename Scalar, template<typename> class GenVecType> class GenPodProjectionSolver;
template <int> class VecNodeDofConversion;
typedef VecNodeDofConversion<6> VecNodeDof6Conversion;
typedef GenModalGeomState<Vector> ModalGeomState; 

class SDDynamPodPostProcessor : public SDDynamPostProcessor
{
  public:
    SDDynamPodPostProcessor(Domain *d, double *bcx, double *vcx, double *acx,
                            StaticTimers *times, GeomState *geomState = 0,
                            Corotator **allCorot = 0);
    ~SDDynamPodPostProcessor();

    // Perform output
    void dynamOutput(int timeStepIndex, double time, DynamMat & dMat,
                     Vector & externalForce, Vector * aeroForce,
                     SysState<Vector> & systemState);
    void printPODSize(int);
    void makeSensorBasis(VecBasis *);

  private:
    OutputInfo *oinfo;
    int numOutInfo;
    int podSize;

    void buildSensorNodeVector();
    void printSensorValues(GenVector<double> &, OutputInfo *, double *, double *);
    double GetPrescribedSensorValue(int, int, double *);

    VecBasis *SensorBasis;
    GenVector<double> *DispSensorValues;
    GenVector<double> *AccSensorValues;
    GenVector<double> *VelSensorValues;

    bool DispSensor;
    bool AccSensor;
    bool VelSensor;

    std::vector<int> nodeVector;
};

class PodProjectionNonLinDynamic : public NonLinDynamic {
public:
  explicit PodProjectionNonLinDynamic(Domain *);
  virtual ~PodProjectionNonLinDynamic();

  virtual int solVecInfo() const;
  virtual void preProcess();

  void readRestartFile(Vector &, Vector &, Vector &, Vector &, ModalGeomState &);
  int getInitState(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p);
  void updatePrescribedDisplacement(ModalGeomState *);
  void getConstForce(Vector &);
  void getExternalForce(Vector &, Vector &, int, double, ModalGeomState *, Vector &, Vector &, double);
  void getIncDisplacement(ModalGeomState *, Vector &, ModalGeomState *, bool);
  double formRHScorrector(Vector &, Vector &, Vector &, Vector &, Vector &, ModalGeomState *, double);
  void formRHSpredictor(Vector &, Vector &, Vector &, Vector &, ModalGeomState &, double, double);
  void formRHSinitializer(Vector &, Vector &, Vector &, ModalGeomState &, Vector &, ModalGeomState * = NULL);
  ModalGeomState* createGeomState();
  ModalGeomState* copyGeomState(ModalGeomState *);
  virtual void updateStates(ModalGeomState *, ModalGeomState &, double);
  virtual double getStiffAndForce(ModalGeomState &, Vector &, Vector &, double = -1, ModalGeomState * = NULL, bool = false);

  virtual void reBuild(ModalGeomState &, int, double, double);
  void dynamCommToFluid(ModalGeomState *, ModalGeomState *, Vector &, Vector &, Vector &, Vector &, int, int, int, double);
  void dynamOutput(ModalGeomState *, Vector &, Vector &, double, int, Vector &, Vector &, Vector &, ModalGeomState *) const;
  void getConstraintMultipliers(ModalGeomState &);
  void initializeParameters(int step, ModalGeomState *);
  void updateParameters(ModalGeomState *);
  bool checkConstraintViolation(double &err, ModalGeomState *);

  // Hiding NonLinDynamic::getSolver
  GenPodProjectionSolver<double> *getSolver();
  const GenPodProjectionSolver<double> *getSolver() const;

  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

  // Hooks in NLDynamSolver
  virtual double getResidualNorm(const Vector &, ModalGeomState &, double);
  int checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time); // relies on function hiding

  void resize(ModalGeomState *refState, ModalGeomState *geomState, ModalGeomState *stepState, Vector *stateIncr,
              Vector &v, Vector &a, Vector &vp, Vector &force) {}

  // Local bases
  int selectLocalBasis(Vector &q);
  void initLocalBasis(Vector &q0);
  void setLocalBasis(ModalGeomState *refState, ModalGeomState *geomState, Vector &q_n, Vector &v, Vector &a);
  virtual void setLocalReducedMesh(int j) {}
  void readLocalBasesCent(const VecNodeDof6Conversion &vecNodeDof6Conversion);
  void readLocalBasesAuxi();
  void projectLocalBases(int i, int j, Vector &q);

  void postProcessNLSA(ModalGeomState *, ModalGeomState *) {}
  void sensitivityAnalysis(ModalGeomState *, ModalGeomState *) {}

protected:
  class Impl;
  GenPodProjectionSolver<double> *solver_;
  GeomState *geomState_Big, *refState_Big;
  SDDynamPodPostProcessor *podPostPro;
  Vector *d0_Big, *v0_Big;
  GenFullSquareMatrix<double> K_reduced;
  int localBasisId;
#ifdef USE_EIGEN3
  Eigen::Array<Eigen::MatrixXd,Eigen::Dynamic,Eigen::Dynamic> VtV;
  Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> d;
  Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,Eigen::Dynamic> w;
  Eigen::MatrixXd uc;
#endif

private:
  virtual bool factorWhenBuilding() const; // Overriden

  void saveMidTime(double); 
  void saveDelta(double);
  void saveStateSnapshot(const GeomState &);
  void saveVelocitySnapshot(const Vector &);
  void saveAccelerationSnapshot(const Vector &);
  void handleResidualSnapshot(const Vector &);
  void expandForce(Vector &fr, Vector &f) const;
  void reduceDisp(Vector &d, Vector &dr) const;

protected:
  std::unique_ptr<Impl> impl_;
  std::unique_ptr<Impl> sttImpl_;
  std::unique_ptr<Impl> velImpl_;
  std::unique_ptr<Impl> accImpl_;
  std::unique_ptr<Impl> resImpl_;
  std::unique_ptr<Impl> jacImpl_;
  
private:
  friend class Updater;
  friend class Impl;

  // Disallow copy and assignment
  PodProjectionNonLinDynamic(const PodProjectionNonLinDynamic &);
  PodProjectionNonLinDynamic &operator=(const PodProjectionNonLinDynamic &);
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class PodProjectionNonLinDynamic::Updater : public IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState> {
public:

  static double integrate(int iter, PodProjectionNonLinDynamic *pbd, ModalGeomState *refState, ModalGeomState *geomState,
                          GenVector<double> *du, GenVector<double> &residual,
                          GenVector<double> &elementInternalForce, GenVector<double> &gRes, GenVector<double> &vel_n,
                          GenVector<double> &accel, double midTime, bool forceOnly = false) {
    pbd->saveMidTime(midTime);

    return IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState>::integrate(
        iter, pbd, refState, geomState, du, residual, elementInternalForce, gRes, vel_n, accel, midTime, forceOnly);
  }
  
  static void midpointIntegrate(PodProjectionNonLinDynamic *pbd, GenVector<double> &velN,
                                double delta, ModalGeomState *refState, 
                                ModalGeomState *geomState, GenVector<double> *dummy1,
                                GenVector<double> &dummy2, GenVector<double> &dummy3,
                                GenVector<double> &dummy4, GenVector<double> &acceleration, bool zeroRot) {
    pbd->saveDelta(delta);

    Vector qN = refState->q; // Local bases
    IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);

/* XXX snapshot collection from Model II is currently not supported
    pbd->saveStateSnapshot(*geomState);
    pbd->saveVelocitySnapshot(velN);
    pbd->saveAccelerationSnapshot(acceleration);
*/
    // Local bases
    pbd->setLocalBasis(refState, geomState, qN, velN, acceleration);

  } 

 static double formRHScorrector(PodProjectionNonLinDynamic *pbd, GenVector<double> &inc_displac,
                                GenVector<double> &vel_n, GenVector<double> &accel,
                                GenVector<double> &residual, GenVector<double> &rhs,
                                ModalGeomState *geomState, double delta) {
    const double result = IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState, delta);

    pbd->handleResidualSnapshot(rhs);

    return result;
  }
};

} /* end namespace Rom */

#endif /* ROM_PODPROJECTIONNONLINDYNAMIC_H */
