#ifndef ROM_CHECKNONLINDYNAMIC_H
#define ROM_CHECKNONLINDYNAMIC_H

#include <Problems.d/NonLinDynam.h>

#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>
#include <Driver.d/StateUpdater.h>

#include <memory>

class Domain;

namespace Rom {

class CheckNonLinDynamic : public NonLinDynamic {
public:
  explicit CheckNonLinDynamic(Domain *);
  ~CheckNonLinDynamic();

  // Required additional preprocessing
  virtual void preProcess();
  
  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

private:
  void handleStateSnap(const GeomState &state);
 
  class Impl;
  std::unique_ptr<Impl> impl_;

  friend class Updater;
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class CheckNonLinDynamic::Updater : public IncrUpdater<CheckNonLinDynamic, GenVector<double>, GeomState> {
public:
  static void midpointIntegrate(CheckNonLinDynamic *pbd, GenVector<double> &velN,
                                double delta, GeomState *refState, 
                                GeomState *geomState, GenVector<double> *dummy1,
                                GenVector<double> &dummy2, GenVector<double> &dummy3,
                                GenVector<double> &dummy4, GenVector<double> &acceleration, bool zeroRot) {
    IncrUpdater<CheckNonLinDynamic, GenVector<double>, GeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);
   
    pbd->handleStateSnap(*geomState);
  }
};

} /* end namespace Rom */

#endif /* ROM_CHECKNONLINDYNAMIC_H */
