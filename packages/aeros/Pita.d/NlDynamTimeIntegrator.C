#include "NlDynamTimeIntegrator.h"
#include "PitaNonLinDynam.h"
#include <Corotational.d/GeomState.h>
#include <Solvers.d/Solver.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

extern int verboseFlag;
extern int totalNewtonIter;

using std::fprintf;

namespace Pita {

NlDynamTimeIntegrator::NlDynamTimeIntegrator(PitaNonLinDynamic * pbDesc) :
  DynamTimeIntegrator(pbDesc->solVecInfo()),
  probDesc_(pbDesc),
  geomState_(pbDesc->createGeomState()),
  stepState_(pbDesc->createGeomState()),
  refState_(pbDesc->createGeomState()),
  internalState_(pbDesc->solVecInfo()),
  acceleration_(pbDesc->solVecInfo()),
  externalForce_(pbDesc->solVecInfo()),
  gravityForce_(pbDesc->solVecInfo()),
  aeroForce_(pbDesc->solVecInfo()),
  elementInternalForce_(pbDesc->elemVecInfo()),
  incDisplac_(pbDesc->solVecInfo()),
  previousVelocity_(pbDesc->solVecInfo()),
  residual_(pbDesc->solVecInfo()),
  rhs_(pbDesc->solVecInfo()),
  stateIncr_(pbDesc->solVecInfo()),
  beta_(), gamma_(), alphaf_(), alpham_(),
  localDt_(pbDesc->getDt()), localDelta_(pbDesc->getDelta()),
  maxNumIter_(pbDesc->getMaxit()),
  currTime_(), midTime_(), currStep_()
{
  probDesc_->getNewmarkParameters(beta_, gamma_, alphaf_, alpham_);
  pbDesc->getConstForce(gravityForce_);
  setTimeStepSize(Seconds(localDt_));
  probDesc_->getInitState(displacement(), velocity(), acceleration_, previousVelocity_);
  double initTime;
  probDesc_->getInitialTime(currStep_, initTime);
  setTimeStepCount(TimeStepCount(currStep_));
  updateInitialCondition(internalState_, Seconds(initTime));
}

NlDynamTimeIntegrator::~NlDynamTimeIntegrator() {
  delete refState_;
  delete stepState_;
  delete geomState_;
}

void
NlDynamTimeIntegrator::timeStepSizeIs(Seconds dt) {
  if (dt != timeStepSize()) { 
    updateDt(dt.value());
    setTimeStepSize(dt);
  }
}

void
NlDynamTimeIntegrator::initialConditionIs(const DynamState & initialState, Seconds initialTime) {
  updateInitialCondition(initialState, initialTime);
  setTimeStepCount(TimeStepCount(0));
  performNotification(&NotifieeConst::onInitialCondition);
}

void
NlDynamTimeIntegrator::currentTimeInc(Seconds increment) {
  int stepCount = static_cast<int>(std::floor(increment.value() / localDt_));
  double remainder = increment.value() - stepCount * localDt_;
  integrate(stepCount);
  updateDt(remainder);
  integrate(1);
  updateDt(timeStepSize().value());
}

void
NlDynamTimeIntegrator::timeStepCountInc(TimeStepCount steps) {
  integrate(steps.value());
}


// Private methods

void
NlDynamTimeIntegrator::updateInitialCondition(const DynamState & initialState, Seconds initialTime) {
  currTime_ = initialTime.value();
  rhs_.zero();
  aeroForce_.zero();
  *geomState_ = *refState_;
  geomState_->update(initialState.displacement());
  *stepState_ = *geomState_;
  velocity() = initialState.velocity();
  setInitialTime(initialTime);
  setInitialState(initialState);
  setCurrentTime(initialTime);
  setCurrentState(initialState);
  updateAcceleration();
}

void
NlDynamTimeIntegrator::updateDt(double timeStep) {
  localDt_ = timeStep;
  localDelta_ = 0.5 * timeStep;
}

void
NlDynamTimeIntegrator::updateAcceleration() {
  if (probDesc()->getInitialAcceleration()) {
    if(verboseFlag) {
      fprintf(stderr, " ... Computing consistent initial acceleration ...\n");
    }
    // a_0 = M^{-1} * (fext(t_0) - fint(u_0) - C * v_0)
    probDesc_->getExternalForce(externalForce_, gravityForce_, -1, currTime_, geomState_, elementInternalForce_, aeroForce_, localDelta_);
    probDesc_->formRHSinitializer(externalForce_, velocity(), elementInternalForce_, *geomState_, acceleration_);
    probDesc_->getSolver()->reSolve(acceleration_);
  } else {
    // a_0 = 0
    acceleration_.zero();
  }
}

void
NlDynamTimeIntegrator::integrate(int steps) {
  // Time-marching
  for (int lastStep = currStep_ + steps; currStep_ < lastStep; ++currStep_) {
    midTime_ = currTime_ + (1 - alphaf_) * localDt_; // Equilibrium enforced at t_{n + (1 - alpha_f)}
    currTime_ += localDt_;
    if (verboseFlag) {
      fprintf(stderr, "\n**** Begin Time Step %d - Time %e ****\n", currStep_ + 1, currTime_);
    }
    probDesc_->getExternalForce(externalForce_, gravityForce_, currStep_, midTime_, geomState_, elementInternalForce_, aeroForce_, localDelta_);
    stateIncr_.zero();

    GenVector<double> & velocity_n = velocity(); // Local reference is safe as long as internalState_ is not shared again

    // Newton iterations
    int convergenceStatus;
    double initialResNorm = 0.0, lastResNorm;
    for (int iter = 0; iter < maxNumIter_; ++iter, ++totalNewtonIter) {
      residual_ = externalForce_;
      probDesc_->getStiffAndForce(*geomState_, residual_, elementInternalForce_, midTime_); // Update elementary tangential stiffness & residual force (force imbalance) at equilibrium time
      probDesc_->reBuild(*geomState_, iter, localDelta_, midTime_); // Assemble and factor Mdyn
      geomState_->get_inc_displacement(incDisplac_, *stepState_, probDesc()->getZeroRot()); // Compute incremental displacement
      lastResNorm = probDesc_->formRHScorrector(incDisplac_, velocity_n, acceleration_, residual_, rhs_, geomState_, localDelta_); // Assemble residual (in rhs) 
      if (iter == 0) {
        initialResNorm = lastResNorm;
      }
      residual_ = rhs_;
      probDesc_->getSolver()->reSolve(rhs_); // Displacement increment: rhs <- Mdyn^{-1} * rhs
      convergenceStatus = probDesc_->checkConvergence(iter, lastResNorm, residual_, rhs_, currTime_);
      stateIncr_ = rhs_;
      if (convergenceStatus == 1) {
        break;
      } else if (convergenceStatus == -1) {
        fprintf(stderr, " ... WARNING: Solution diverging\n");
      }
      geomState_->update(stateIncr_);
    }
    if (convergenceStatus == 0) {
      fprintf(stderr,
              " *** WARNING: At time %f Newton solver did not reach convergence after %d iterations (residual = %9.3e, final = %9.3e, target = %9.3e)\n",
              currTime_, maxNumIter_, initialResNorm, lastResNorm, probDesc_->getTolerance());
    }

    // Update internal state
    previousVelocity_ = velocity_n;
    geomState_->midpoint_step_update(velocity_n, acceleration_, localDelta_, *stepState_, beta_, gamma_, alphaf_, alpham_, probDesc()->getZeroRot());
    acceleration_.linC(-(1-gamma_)/gamma_, acceleration_, -1/(localDt_*gamma_), previousVelocity_, 1/(localDt_*gamma_), velocity_n);

    // Update attributes 
    setCurrentTime(Seconds(currTime_));
    geomState_->get_inc_displacement(displacement(), *refState_, false);
    
    // Publish results and perform notification
    // Allows side-effects, such as post-processing
    // Invalidates reference velocity_n
    setCurrentTime(Seconds(currTime_));
    setCurrentState(internalState_);
    setTimeStepCount(TimeStepCount(currStep_ + 1));
    performNotification(&NotifieeConst::onCurrentCondition);
  }
}

} // end namespace Pita
