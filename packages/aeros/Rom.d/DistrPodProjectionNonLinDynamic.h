#ifndef ROM_DISTRPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTRPODPROJECTIONNONLINDYNAMIC_H

#include <Paral.d/MDNLDynam.h>

#include <Corotational.d/DistrGeomState.h>
#include <Feti.d/DistrVector.h>
#include <Driver.d/StateUpdater.h>

#include "ModalGeomState.h"
#include "DistrVecBasis.h"
#include "EiGalerkinProjectionSolver.h"
#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"

#include <memory>
#include <vector>

class Domain;

namespace Rom {

// Specialization of the non-linear dynamics problem enabling the collection the snapshots
class DistrPodProjectionNonLinDynamic : public MDNLDynamic {
public:

  // Helper class to be used as template parameter in NLDynamSolver
  class Updater;

  explicit DistrPodProjectionNonLinDynamic(Domain *);
  virtual ~DistrPodProjectionNonLinDynamic();

  virtual DistrInfo &solVecInfo();

  // Required additional pre-processing
  virtual void preProcess();
 
  void readRestartFile(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                       DistrVector &v_p, DistrModalGeomState &geomState);
  int getInitState(DistrVector &d, DistrVector& v, DistrVector &a, DistrVector &v_p);
  void updatePrescribedDisplacement(DistrModalGeomState *geomState);
  void getConstForce(DistrVector &constantForce);
  void getExternalForce(DistrVector &rhs, DistrVector &constantForce,
                        int tIndex, double time, DistrModalGeomState *geomState, 
                        DistrVector &elementInternalForce, DistrVector &aeroF, double localDelta);

  void getIncDisplacement(DistrModalGeomState *geomState, DistrVector &du, DistrModalGeomState *refState, bool zeroRot);

  double formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                          DistrVector& residual, DistrVector& rhs, DistrModalGeomState *geomState, double localDelta);

  void formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual,
                        DistrVector& rhs, DistrModalGeomState &, double mid, double localDelta);

  void formRHSinitializer(DistrVector &fext, DistrVector &velocity, DistrVector &elementInternalForce,
                          DistrModalGeomState &geomState, DistrVector &rhs, DistrModalGeomState *refState = NULL);
  
  
  DistrModalGeomState* createGeomState();
  DistrModalGeomState* copyGeomState(DistrModalGeomState* geomState);

  virtual void updateStates(DistrModalGeomState *refState, DistrModalGeomState& geomState, double time);

  virtual double getStiffAndForce(DistrModalGeomState& geomState, DistrVector& residual, DistrVector& elementInternalForce,
                            double midtime=-1, DistrModalGeomState *refState = NULL, bool forceOnly = false);

  void reBuild(DistrModalGeomState& geomState, int iter, double localDelta, double t);

  void dynamCommToFluid(DistrModalGeomState* geomState, DistrModalGeomState* bkGeomState,
                        DistrVector& velocity, DistrVector& bkVelocity,
                        DistrVector& vp, DistrVector& bkVp, int step, int parity,
                        int aeroAlg, double time);

  void dynamOutput(DistrModalGeomState* geomState, DistrVector& velocity, DistrVector &vp,
                   double time, int timestep, DistrVector& force, DistrVector &aeroF, DistrVector &acceleration,
                   DistrModalGeomState *refState);

  void getConstraintMultipliers(DistrModalGeomState &geomState);

  void initializeParameters(int step, DistrModalGeomState *geomState);
  void updateParameters(DistrModalGeomState *geomState);
  bool checkConstraintViolation(double &err, DistrModalGeomState *geomState);

  // hide call to MDNLDynamic::getSolver()
  GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> > *getSolver();
  const GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> > *getSolver() const;

  // Hook in MNLDynamSolver
  double getResidualNorm(DistrVector &vec, DistrModalGeomState &, double);
//  int checkConvergence(int iter, double rN, DistrVector& residual, DistrVector& dv, double time);

  void resize(DistrModalGeomState *refState, DistrModalGeomState *geomState, DistrModalGeomState *stepState, DistrVector *stateIncr,
              DistrVector &v, DistrVector &a, DistrVector &vp, DistrVector &force){}

  //Local bases
  int selectLocalBasis(DistrVector &q);
  void initLocalBasis(DistrVector &q0);
  void setLocalBasis(DistrModalGeomState *refState, DistrModalGeomState *geomState, DistrVector &q_n, DistrVector &v, DistrVector &a);
  virtual void setLocalReducedMesh(int j) {}
  void readLocalBasesCent(DistrVecNodeDof6Conversion &converter, DistrNodeDof6Buffer &buffer, std::vector<int> &locBasisVec);
  void readLocalBasesAuxi();
  void projectLocalBases(int i, int j, DistrVector &q);

protected:
  GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> > *solver_;
  DistrGeomState *geomState_Big, *refState_Big;
  MultiDomDynPodPostProcessor *podPostPro;
  DistrVector *d0_Big, *v0_Big;
  GenFullSquareMatrix<double> K_reduced;
  int localBasisId;
  bool resetFromClean;
  DistrInfo reducedInfo;
  DistrVecBasis  centroids;
  DistrVecBasis projectionBasis_;
  DistrVecBasis dualProjectionBasis_; 
  std::vector<std::map<int,double> > dualMapVectors_; //first -> slave node global ID, second -> mu value
#ifdef USE_EIGEN3
  // data structures for fast basis switching
  Eigen::Array<Eigen::MatrixXd,Eigen::Dynamic,Eigen::Dynamic> VtV;
  Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> d;
  Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,Eigen::Dynamic> w;
#endif

private:
  friend class Updater;

  virtual bool factorWhenBuilding() const; // Override
  
  void expandForce(DistrVector &fr, DistrVector &f);
  void reduceDisp(DistrVector &d, DistrVector &dr); 
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class DistrPodProjectionNonLinDynamic::Updater : public IncrUpdater<DistrPodProjectionNonLinDynamic, GenDistrVector<double>, DistrModalGeomState> {
public:
  static double integrate(int iter, DistrPodProjectionNonLinDynamic *pbd, DistrModalGeomState *refState, DistrModalGeomState *geomState,
                          GenDistrVector<double> *du, GenDistrVector<double> &residual,
                          GenDistrVector<double> &elementInternalForce, GenDistrVector<double> &gRes, GenDistrVector<double> &vel_n,
                          GenDistrVector<double> &accel, double midTime, bool forceOnly=false) {

    return IncrUpdater<DistrPodProjectionNonLinDynamic, GenDistrVector<double>, DistrModalGeomState>::integrate(
        iter, pbd, refState, geomState, du, residual, elementInternalForce, gRes, vel_n, accel, midTime, forceOnly);
  }

  static void midpointIntegrate(DistrPodProjectionNonLinDynamic *pbd, GenDistrVector<double> &velN,
                                double delta, DistrModalGeomState *refState, 
                                DistrModalGeomState *geomState, GenDistrVector<double> *dummy1,
                                GenDistrVector<double> &dummy2, GenDistrVector<double> &dummy3,
                                GenDistrVector<double> &dummy4, GenDistrVector<double> &acceleration, bool zeroRot) {

    // local basis crap goes here
    DistrVector qN = refState->q;
    IncrUpdater<DistrPodProjectionNonLinDynamic, GenDistrVector<double>, DistrModalGeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);
   
    pbd->setLocalBasis(refState, geomState, qN, velN, acceleration);
 
  }

  

};

} /* end namespace Rom */

#endif /* ROM_DISTRPODPROJECTIONNONLINDYNAMIC_H */
