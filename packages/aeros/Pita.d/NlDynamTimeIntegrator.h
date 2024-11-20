#ifndef PITA_NLDYNAMTIMEINTEGRATOR_H
#define PITA_NLDYNAMTIMEINTEGRATOR_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"
#include "DynamTimeIntegrator.h"
#include "NlDynamOps.h"

class GeomState;

namespace Pita {

class PitaNonLinDynamic;

class NlDynamTimeIntegrator : public DynamTimeIntegrator {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDynamTimeIntegrator);

  virtual void timeStepSizeIs(Seconds dt);
  virtual void initialConditionIs(const DynamState & initialState, Seconds initialTime = Seconds(0.0));
  virtual void currentTimeInc(Seconds increment);
  virtual void timeStepCountInc(TimeStepCount steps = TimeStepCount(1));

  static NlDynamTimeIntegrator::Ptr New(PitaNonLinDynamic * pbDesc) {
    return new NlDynamTimeIntegrator(pbDesc); 
  }

  NlDynamOps::Ptr nlDynamOpsNew() const { return NlDynamOps::New(probDesc_); }
  
  const PitaNonLinDynamic * probDesc() const { return probDesc_; }
  
  // Accessors for output
  GeomState * geomState() const { return geomState_; }
  const GenVector<double> & externalForce() const { return externalForce_; }
  const GenVector<double> & acceleration() const { return acceleration_; }

protected:
  explicit NlDynamTimeIntegrator(PitaNonLinDynamic * pbDesc);
  ~NlDynamTimeIntegrator();

private:
  void integrate(int stepCount);
  void updateInitialCondition(const DynamState & initialState, Seconds initialTime);
  void updateDt(double timeStep);
  void updateAcceleration();
  
  GenVector<double> & displacement() { return internalState_.displacement(); }
  const GenVector<double> & displacement() const { return internalState_.displacement(); }

  GenVector<double> & velocity() { return internalState_.velocity(); }
  const GenVector<double> & velocity() const { return internalState_.velocity(); }

  PitaNonLinDynamic * probDesc_;
  
  GeomState * geomState_;
  GeomState * stepState_;
  GeomState * refState_; 

  DynamState internalState_;
  GenVector<double> acceleration_;

  GenVector<double> externalForce_;
  GenVector<double> gravityForce_;
  GenVector<double> aeroForce_;
  GenVector<double> elementInternalForce_;
  
  GenVector<double> incDisplac_;
  GenVector<double> previousVelocity_; 
  GenVector<double> residual_;
  GenVector<double> rhs_;
  GenVector<double> stateIncr_;
  
  double beta_, gamma_, alphaf_, alpham_;
  double localDt_;
  double localDelta_;
  int maxNumIter_;
  
  double currTime_;
  double midTime_; // Time for equilibrium enforcement
  int currStep_;
};
  
} // end namespace Pita

#endif /* PITA_NLDYNAMTIMEINTEGRATOR_H */
