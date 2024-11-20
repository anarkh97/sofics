#ifndef PITA_OLD_NL_DYNAM_TIME_INTEGRATOR_H
#define PITA_OLD_NL_DYNAM_TIME_INTEGRATOR_H

class NonLinDynamic;

#include <Problems.d/NonLinDynam.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>

class NLDynamTimeIntegrator
{
public:
  typedef GenVector<double> VecType;
  
  explicit NLDynamTimeIntegrator(NonLinDynamic & pbDesc);
  ~NLDynamTimeIntegrator();  

  // Initial conditions 
  void setCurrentDisplacement(const VecType &);
  void getCurrentDisplacement(VecType &) const; 
  void setCurrentVelocity(const VecType &);
  void getCurrentVelocity(VecType &) const; 

  // Time integration parameters
  double currentTime() const;
  void currentTimeIs(double); 
  double timeStep() const;
  void timeStepIs(double); 
  int currentTimeStepNumber() const;
  void currentTimeStepNumberIs(int);

  // PostProcessor
  const NLDynamPostProcessor & postProcessor() const;
  void postProcessor(NLDynamPostProcessor &); 

  // Integrator
  void integrate(int numSteps = 1);

private:
  NonLinDynamic & probDesc;
  const NLDynamPostProcessor * postProcessor_;
  GeomState * geomState;
  GeomState * stepState;
  GeomState * refState; 
  VecType velocity;
  VecType inc_displac;
  VecType gravityForce;
  VecType elementInternalForce;
  VecType residual;
  VecType rhs;
  VecType aeroForce;
  VecType external_force;
  VecType stateIncr;
  VecType dummyVp; 
  VecType acceleration;
  double currTime;
  double midTime;
  int currStep;
  double localDt;
  double localDelta;
  int numStages;
  int maxNumIter;
  double dlambda;
  double beta, gamma, alphaf, alpham;
};

inline void
NLDynamTimeIntegrator::getCurrentDisplacement(NLDynamTimeIntegrator::VecType & disp) const
{
  geomState->get_inc_displacement(disp, *refState, false);
}
                                                                                                                                                                                                     
inline void
NLDynamTimeIntegrator::setCurrentDisplacement(const NLDynamTimeIntegrator::VecType & disp)
{
  *geomState = *refState;
  geomState->update(disp);
  *stepState = *geomState;
}

inline void
NLDynamTimeIntegrator::setCurrentVelocity(const NLDynamTimeIntegrator::VecType & vel)
{
  velocity = vel;
}

inline void
NLDynamTimeIntegrator::getCurrentVelocity(NLDynamTimeIntegrator::VecType & vel) const
{
  vel = velocity;
}

inline double
NLDynamTimeIntegrator::currentTime() const
{
  return currTime;
}

inline double
NLDynamTimeIntegrator::timeStep() const
{
  return localDt;
}

inline void
NLDynamTimeIntegrator::timeStepIs(double newStep)
{
  localDt = newStep;
  localDelta = 0.5 * newStep;
  midTime = currTime + localDelta;
}

inline int
NLDynamTimeIntegrator::currentTimeStepNumber() const
{
  return currStep;
}

inline void
NLDynamTimeIntegrator::currentTimeStepNumberIs(int newTimeStep)
{
  if (currStep == newTimeStep) {
    return;
  }
  currStep = newTimeStep;
}

inline const NLDynamPostProcessor &
NLDynamTimeIntegrator::postProcessor() const
{
  return *postProcessor_;
}

inline void
NLDynamTimeIntegrator::postProcessor(NLDynamPostProcessor & pp)
{
  postProcessor_ = &pp;
}

#endif
