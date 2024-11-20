#ifndef ROM_SNAPSHOTNONLINDYNAMIC_H
#define ROM_SNAPSHOTNONLINDYNAMIC_H

#include <Problems.d/NonLinDynam.h>

#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>
#include <Driver.d/StateUpdater.h>

#include <memory>

class Domain;

namespace Rom {

// Specialization of the non-linear dynamics problem enabling the collection the snapshots
class SnapshotNonLinDynamic : public NonLinDynamic {
public:

  explicit SnapshotNonLinDynamic(Domain *);

  //Hooks in NLDynamSolver
  int checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time); //function hiding

  // Required additional pre- and post-processing
  virtual void preProcess();
  void postProcess(); 
  
  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

protected:
  // Interface to implementation
  class Impl {
  public:
    virtual void lastMidTimeIs(double) = 0;
    virtual void lastDeltaIs(double) = 0;
    virtual void stateSnapshotAdd(const GeomState &) = 0;
    virtual void internalStateSnapshotAdd(const GeomState &) = 0;
    virtual void velocSnapshotAdd(const Vector &veloc) = 0;
    virtual void accelSnapshotAdd(const Vector &accel) = 0;
    virtual void handleResidualSnapshot(const Vector &res) = 0;
    virtual void handleJacobianSnapshot() = 0;
    virtual void postProcess() = 0;

    virtual ~Impl() {}

  protected:
    Impl() {}

  private:
    // Disallow copy and assigment
    Impl(const Impl &);
    Impl &operator=(const Impl &);
  };

private:
  // Snapshot collection 
  void saveMidTime(double t);
  void saveDelta(double delta);
  void saveStateSnapshot(const GeomState &state);
  void saveInternalStateSnapshot(const GeomState &state);
  void saveVelocitySnapshot(const GeomState &state, const Vector &veloc);
  void saveAccelerationSnapshot(const GeomState &state, const Vector &accel);
  void handleResidualSnapshot(const Vector &snap); 


//  BasisType outputBasisType_;
  std::unique_ptr<Impl> stateImpl_;
  std::unique_ptr<Impl> internalStateImpl_;
  std::unique_ptr<Impl> velocImpl_;
  std::unique_ptr<Impl> accelImpl_;
  std::unique_ptr<Impl> resImpl_;
  std::unique_ptr<Impl> jacImpl_;

  friend class Updater;
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class SnapshotNonLinDynamic::Updater : public IncrUpdater<SnapshotNonLinDynamic, GenVector<double>, GeomState> {
public:
  static double integrate(int iter, SnapshotNonLinDynamic *pbd, GeomState *refState, GeomState *geomState,
                          GenVector<double> *du, GenVector<double> &residual,
                          GenVector<double> &elementInternalForce, GenVector<double> &gRes, GenVector<double> &vel_n,
                          GenVector<double> &accel, double midTime, bool forceOnly = false) {
    pbd->saveMidTime(midTime);

    return IncrUpdater<SnapshotNonLinDynamic, GenVector<double>, GeomState>::integrate(
        iter, pbd, refState, geomState, du, residual, elementInternalForce, gRes, vel_n, accel, midTime, forceOnly);
  }

  static void midpointIntegrate(SnapshotNonLinDynamic *pbd, GenVector<double> &velN,
                                double delta, GeomState *refState, 
                                GeomState *geomState, GenVector<double> *dummy1,
                                GenVector<double> &dummy2, GenVector<double> &dummy3,
                                GenVector<double> &dummy4, GenVector<double> &acceleration, bool zeroRot) {
    pbd->saveDelta(delta); 

    IncrUpdater<SnapshotNonLinDynamic, GenVector<double>, GeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);
    
    pbd->saveStateSnapshot(*geomState);
    pbd->saveInternalStateSnapshot(*geomState);
    pbd->saveVelocitySnapshot(*geomState, velN);
    pbd->saveAccelerationSnapshot(*geomState, acceleration);
  }

  static double formRHScorrector(SnapshotNonLinDynamic *pbd, GenVector<double> &inc_displac,
                                 GenVector<double> &vel_n, GenVector<double> &accel,
                                 GenVector<double> &residual, GenVector<double> &rhs,
                                 GeomState *geomState, double delta) {
    const double result = IncrUpdater<SnapshotNonLinDynamic, GenVector<double>, GeomState>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState, delta);

    pbd->handleResidualSnapshot(rhs);

    return result;
  }

};

} /* end namespace Rom */

#endif /* ROM_SNAPSHOTNONLINDYNAMIC_H */
