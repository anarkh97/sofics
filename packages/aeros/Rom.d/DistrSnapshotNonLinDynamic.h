#ifndef ROM_DISTRSNAPSHOTNONLINDYNAMIC_H
#define ROM_DISTRSNAPSHOTNONLINDYNAMIC_H

#include <Paral.d/MDNLDynam.h>

#include <Corotational.d/DistrGeomState.h>
#include <Feti.d/DistrVector.h>
#include <Driver.d/StateUpdater.h>

#include <memory>

class Domain;

namespace Rom {

// Specialization of the non-linear dynamics problem enabling the collection the snapshots
class DistrSnapshotNonLinDynamic : public MDNLDynamic {
public:
  explicit DistrSnapshotNonLinDynamic(Domain *);

  // Required additional pre- and post-processing
  virtual void preProcess();
  void postProcess() { impl_->postProcess(); }
  
  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

  void updateParameters(DistrGeomState *geomState) { impl_->muvarSnapshotAdd(*geomState); MDNLDynamic::updateParameters(geomState); }

protected:
  // Interface to implementation
  class Impl {
  public:
    virtual void lastMidTimeIs(double) = 0;
    virtual void lastDeltaIs(double) = 0;
    virtual void stateSnapshotAdd(const DistrGeomState &) = 0;
    virtual void velocSnapshotAdd(const DistrVector &) = 0;
    virtual void accelSnapshotAdd(const DistrVector &) = 0;
    virtual void dsvarSnapshotAdd(const DistrGeomState &) = 0;
    virtual void muvarSnapshotAdd(const DistrGeomState &) = 0;
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
  void saveMidTime(double t) { impl_->lastMidTimeIs(t); }
  void saveDelta(double delta) { impl_->lastDeltaIs(delta); }
  void saveStateSnapshot(const DistrGeomState &state) { impl_->stateSnapshotAdd(state); impl_->dsvarSnapshotAdd(state); }
  void saveVelocSnapshot(DistrGeomState &state, const DistrVector &veloc);
  void saveAccelSnapshot(DistrGeomState &state, const DistrVector &accel);
 
  std::unique_ptr<Impl> impl_;

  friend class Updater;
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class DistrSnapshotNonLinDynamic::Updater : public IncrUpdater<DistrSnapshotNonLinDynamic, GenDistrVector<double>, DistrGeomState> {
public:
  static double integrate(int iter, DistrSnapshotNonLinDynamic *pbd, DistrGeomState *refState, DistrGeomState *geomState,
                          GenDistrVector<double> *du, GenDistrVector<double> &residual,
                          GenDistrVector<double> &elementInternalForce, GenDistrVector<double> &gRes, GenDistrVector<double> &vel_n,
                          GenDistrVector<double> &accel, double midTime, bool forceOnly=false) {
    pbd->saveMidTime(midTime);

    return IncrUpdater<DistrSnapshotNonLinDynamic, GenDistrVector<double>, DistrGeomState>::integrate(
        iter, pbd, refState, geomState, du, residual, elementInternalForce, gRes, vel_n, accel, midTime, forceOnly);
  }

  static void midpointIntegrate(DistrSnapshotNonLinDynamic *pbd, GenDistrVector<double> &velN,
                                double delta, DistrGeomState *refState, 
                                DistrGeomState *geomState, GenDistrVector<double> *dummy1,
                                GenDistrVector<double> &dummy2, GenDistrVector<double> &dummy3,
                                GenDistrVector<double> &dummy4, GenDistrVector<double> &acceleration, bool zeroRot) {
    pbd->saveDelta(delta);

    IncrUpdater<DistrSnapshotNonLinDynamic, GenDistrVector<double>, DistrGeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);
    
    pbd->saveStateSnapshot(*geomState);
    pbd->saveVelocSnapshot(*geomState, velN);
    pbd->saveAccelSnapshot(*geomState, acceleration);
  }
};

} /* end namespace Rom */

#endif /* ROM_DISTRSNAPSHOTNONLINDYNAMIC_H */
