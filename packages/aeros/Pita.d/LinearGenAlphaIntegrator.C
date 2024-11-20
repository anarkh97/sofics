#include "LinearGenAlphaIntegrator.h"
#include <Problems.d/DynamDescr.h>
#include <Driver.d/DynamProbType.h>

namespace Pita {

// ---------------------------------------
// LinearGenAlphaIntegrator implementation
// ---------------------------------------

// Constructors

LinearGenAlphaIntegrator::LinearGenAlphaIntegrator(LinearDynamOps::Manager * dOpsMgr,
                                                   const GeneralizedAlphaParameter & param,
                                                   ExternalForceStatus efs) :
  AffineDynamTimeIntegrator(dOpsMgr->probDesc()->solVecInfo(), efs),
  dynamOpsManager_(dOpsMgr),
  probDesc_(dOpsMgr->probDesc()),
  dynamOps_(NULL),
  rhoInfinity_(param.rhoInfinity()),
  beta_(param.beta()),
  gamma_(param.gamma()),
  alphaf_(param.alphaf()),
  alpham_(param.alpham()),
  constForce_(probDesc_->solVecInfo(), 0.0),
  externalForce_(probDesc_->solVecInfo(), 0.0),
  aeroForce_(probDesc_->solVecInfo(), 0.0),
  currentDisplacement_(probDesc_->solVecInfo()),
  currentVelocity_(probDesc_->solVecInfo()),
  currentAcceleration_(probDesc_->solVecInfo()),
  nextDisplacement_(probDesc_->solVecInfo()),
  nextVelocity_(probDesc_->solVecInfo()),
  nextAcceleration_(probDesc_->solVecInfo()),
  previousVelocity_(initialState().vectorSize(), 0.0),
  rhs_(probDesc_->solVecInfo(), 0.0),
  temp_(probDesc_->solVecInfo(), 0.0),
  temp2_(probDesc_->solVecInfo(), 0.0)
{
  getDynamOps(param);
  probDesc()->getConstForce(constForce_);
  currentDisplacement_ = initialState().displacement();
  currentVelocity_ = initialState().velocity();
  currentAcceleration_ = VectorType(initialState().vectorSize(), 0.0);
  setTimeStepSize(param.timeStepSize());
}

// Public overriden functions

void
LinearGenAlphaIntegrator::timeStepSizeIs(Seconds dt) {
  if (dt != timeStepSize()) {
    getDynamOps(GeneralizedAlphaParameter(dt, rhoInfinity_));
    setTimeStepSize(dt);
  }
}

void
LinearGenAlphaIntegrator::initialConditionIs(const DynamState & initState, Seconds initTime) {
  currentDisplacement_ = initState.displacement();
  currentVelocity_ = initState.velocity();
  currentAcceleration_ = VectorType(initState.vectorSize(), 0.0);
  previousVelocity_ = VectorType(initState.vectorSize(), 0.0);

  setInitialTime(initTime);
  setInitialState(initState);
  setTimeStepCount(TimeStepCount(0));
  setCurrentTime(initTime);
  setCurrentState(initState);
  performNotification(&NotifieeConst::onInitialCondition);
}

void
LinearGenAlphaIntegrator::currentTimeInc(Seconds timeIncrement) {
  TimeStepCount stepCount(static_cast<unsigned int>(std::ceil(timeIncrement.value() / timeStepSize().value())));
  integrate(stepCount);
}

void
LinearGenAlphaIntegrator::timeStepCountInc(TimeStepCount steps) {
  integrate(steps);
}

// Public specific functions

void
LinearGenAlphaIntegrator::rhoInfinityIs(double r) {
  GeneralizedAlphaParameter param(timeStepSize(), r);
  getDynamOps(param);
  rhoInfinity_ = r;
}


// External force computation

inline
void
LinearGenAlphaIntegrator::computeExternalForceImpl(Seconds forceEvalTime, SysState<VectorType> & currentState) {
  probDesc_->computeExtForce2(currentState, externalForce_, constForce_, timeStepCount().value(), forceEvalTime.value(), &aeroForce_, gamma_, alphaf_);
}

// Private functions

// Update the dynamic operators 
void
LinearGenAlphaIntegrator::getDynamOps(const GeneralizedAlphaParameter & param) {
  // Retrieve the dynamic operators 
  dynamOps_ = dynamOpsManager_->dynOps(param);
  if (!dynamOps_) {
    // Corresponding dynamic operators do not exist, create them
    dynamOps_ = dynamOpsManager_->dynOpsNew(param);
  }
}

// Perform time integration
void
LinearGenAlphaIntegrator::integrate(TimeStepCount s) {

  unsigned int stepCount = s.value();
  SysState<VectorType> currentSysState(currentDisplacement_, currentVelocity_, currentAcceleration_, previousVelocity_);
  double dt = timeStepSize().value();

  for (unsigned int step = 0; step < stepCount; ++step) {
    Seconds endStepTime = currentTime() + timeStepSize(); 

    // ... Construct force vector, includes time-independent constant force

    // ... Compute external force at time t+dt*(1-alphaf)
    Seconds externalForceTime = currentTime() + (alphaf_ * timeStepSize());
    computeExternalForce(externalForceTime, currentSysState);

    // ... Construct R.H.S. vector
    // ... d_n_h = ((1-alpham)/(1-alphaf))*d_n 
    //           + dt*(1-alpham)*v_n 
    //           + dt*dt*((1-alpham)/2 - beta)*a_n
    // (if alphaf=alpham=0.5,beta=0.25: d_n_h = d_n + dt*0.5*v_n + zero*a_n)

    // First: d_n_h = ((1-alpham)/(1-alphaf))*d_n + dt*(1-alpham)*v_n 
    //d_n_h.linC( ((1.0-alpham)/(1.0-alphaf)), d_n, dt*(1.0-alpham), v_n );
   
    temp_.linC(((1.0 - alpham_) / (1.0 - alphaf_)), currentDisplacement_, dt * (1.0 - alpham_), currentVelocity_); 

    // Second: d_n_h = d_n_h + dt*dt*((1-alpham)/2 - beta)*a_n
    double coef = 0.5 * (1.0 - alpham_) - beta_;
    if (coef != 0.0) {
      temp_.linAdd(dt * dt * coef, currentAcceleration_);
    }

    // Third: Multiply by Mass Matrix M
    dynamOps_->massMatrix()->mult(temp_, rhs_);

    // Accumulate in rhs vector: rhs = Md_n_h + beta*dt*dt*ext_f
    rhs_.linAdd(beta_ * dt * dt, externalForce_);

    if (dynamOps_->dampingMatrix()) {
      // ... d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma(1-alphaf))*v_n
      // ...       - dt*dt*dt*0.5(1-alphaf)*(2*beta - gamma)*a_n
      // (if alphaf=alpham=gamma=0.5,beta=0.25: d_n_h = dt*0.5*d_n - zero*v_n - zero*a_n)

      // ... d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma(1-alphaf))*v_n
      temp_.linC(dt * gamma_, currentDisplacement_, dt * dt * (gamma_ * (1.0 - alphaf_) - beta_), currentVelocity_);

      // ... d_n_h = d_n_h - dt*dt*dt*0.5(1-alphaf)*(2*beta - gamma)*a_n
      temp_.linAdd(dt * dt * dt * 0.5 * (alphaf_ - 1.0) * (2.0 * beta_ - gamma_), currentAcceleration_);

      // Multiply by Damping Matrix C
      dynamOps_->dampingMatrix()->mult(temp_, temp2_);

      // Accumulate in rhs vector 
      rhs_ += temp2_;
    }

    dynamOps_->dynamicMassSolver()->reSolve(rhs_); // Now rhs contains (d_n+1-alphaf)

    // call projector for RBMs in case of rbmfilter level 2
    if (probDesc()->getFilterFlag() == 2)
      probDesc()->project(rhs_);
    
    // one time step forward
    //d_n_p.linC(rhs,(-1.0*alphaf),d_n);
    //d_n_p *= (1.0/(1-alphaf));
    //(temp - alphaf * state0.disp) / (1.0 - alphaf);
    nextDisplacement_.linC(1.0 / (1.0 - alphaf_), rhs_, alphaf_ / (alphaf_ - 1.0), currentDisplacement_);

    //state1.accel = (1.0 / (beta * dt^2)) * (state1.disp - state0.disp) + (-1.0 / (dt * beta)) * state0.velo + (1.0 - 0.5 / beta) * state0.accel;
    temp_.linC(nextDisplacement_, -1.0, currentDisplacement_);
    nextAcceleration_.linC(- 1.0 / (dt * beta_), currentVelocity_, (1.0 - 0.5 / beta_), currentAcceleration_);
    nextAcceleration_.linAdd((1.0 / (beta_ * dt * dt)), temp_);

    //state1.velo = ((dt * (1.0 - gamma)) * state0.accel + (dt * gamma) * state1.accel) + state0.velo;
    nextVelocity_.linC(dt * (1.0 - gamma_), currentAcceleration_, dt * gamma_, nextAcceleration_);
    nextVelocity_ += currentVelocity_;

    // Now swap v_n_p -> v_n and d_n_p -> d_n
    previousVelocity_ = currentVelocity_;
    currentDisplacement_.swap(nextDisplacement_);
    currentVelocity_.swap(nextVelocity_);
    currentAcceleration_.swap(nextAcceleration_);

    // FORCE PRINTING AT LAST ITERATION
    //if (t+1.01*dt > tmax)  probDesc->processLastOutput();

    //postProcessor->dynamOutput( tIndex, dynOps, ext_f, aeroForce, curState);

    //Update state
    setCurrentState(DynamState(currentDisplacement_, currentVelocity_));
    setCurrentTime(endStepTime);
    setTimeStepCount(timeStepCount() + TimeStepCount(1));

    performNotification(&NotifieeConst::onCurrentCondition); 
  }
}


// Specific implementations
LinearGenAlphaIntegratorImpl::LinearGenAlphaIntegratorImpl(
    LinearDynamOps::Manager * dOpsMgr,
    const GeneralizedAlphaParameter & param) :
  LinearGenAlphaIntegrator(dOpsMgr, param, NONHOMOGENEOUS)
{}

void
LinearGenAlphaIntegratorImpl::externalForceStatusIs(LinearGenAlphaIntegrator::ExternalForceStatus efs) {
  if (efs != NONHOMOGENEOUS) {
    throw Fwk::RangeException();
  }
}

void
LinearGenAlphaIntegratorImpl::computeExternalForce(
    Seconds forceEvalTime,
    SysState<VectorType> & currentState) {
  computeExternalForceImpl(forceEvalTime, currentState);
}

HomogeneousGenAlphaIntegrator::HomogeneousGenAlphaIntegrator(
    LinearDynamOps::Manager * dOpsMgr,
    const GeneralizedAlphaParameter & param) :
  LinearGenAlphaIntegrator(dOpsMgr, param, HOMOGENEOUS)
{}

void
HomogeneousGenAlphaIntegrator::externalForceStatusIs(LinearGenAlphaIntegrator::ExternalForceStatus efs) {
  if (efs != HOMOGENEOUS) {
    throw Fwk::RangeException();
  }
}

void
HomogeneousGenAlphaIntegrator::computeExternalForce(
    Seconds forceEvalTime,
    SysState<VectorType> & currentState) { 
  // Do nothing
}

AffineGenAlphaIntegrator::AffineGenAlphaIntegrator(
    LinearDynamOps::Manager * dOpsMgr,
    const GeneralizedAlphaParameter & param) :
  LinearGenAlphaIntegrator(dOpsMgr, param, NONHOMOGENEOUS)
{}

void
AffineGenAlphaIntegrator::externalForceStatusIs(AffineGenAlphaIntegrator::ExternalForceStatus efs) {
  zeroExternalForce();
  setExternalForceStatus(efs);
}

void
AffineGenAlphaIntegrator::computeExternalForce(
    Seconds forceEvalTime,
    SysState<VectorType> & currentState) { 
  if (externalForceStatus() == NONHOMOGENEOUS) {
    computeExternalForceImpl(forceEvalTime, currentState);
  }
}
} // end namespace Pita
