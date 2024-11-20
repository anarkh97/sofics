#include "NLDynamTimeIntegrator.h"

#include <Problems.d/NonLinDynam.h>
#include <Solvers.d/Solver.h>

#include <cstdio>
#include <cstdlib>

using std::fprintf;
using std::exit;

extern int verboseFlag;
extern int totalNewtonIter;

// Constructor & destructor
// ------------------------
NLDynamTimeIntegrator::NLDynamTimeIntegrator(NonLinDynamic & pbDesc) :
  probDesc(pbDesc),
  postProcessor_(&(pbDesc.defaultPostProcessor())),
  geomState(pbDesc.createGeomState()),
  stepState(pbDesc.createGeomState()),
  refState(pbDesc.createGeomState()),
  velocity(pbDesc.solVecInfo()),
  acceleration(pbDesc.solVecInfo()),
  inc_displac(pbDesc.solVecInfo()),
  gravityForce(pbDesc.solVecInfo()),
  elementInternalForce(pbDesc.elemVecInfo()),
  residual(pbDesc.solVecInfo()),
  rhs(pbDesc.solVecInfo()),
  aeroForce(pbDesc.solVecInfo()),
  external_force(pbDesc.solVecInfo()),
  stateIncr(pbDesc.solVecInfo()),
  dummyVp(pbDesc.solVecInfo()),
  localDt(pbDesc.getDt()),
  localDelta(pbDesc.getDelta()),
  numStages(pbDesc.getNumStages()),
  maxNumIter(pbDesc.getMaxit()),
  dlambda(1.0 / pbDesc.getNumStages())
{
  probDesc.getConstForce(gravityForce);
  VecType initialDisplacement(pbDesc.solVecInfo());
  /*int aeroAlg =*/ probDesc.getInitState(initialDisplacement, velocity, acceleration, dummyVp);
  setCurrentDisplacement(initialDisplacement);
  double initialTime;
  probDesc.getInitialTime(currStep, initialTime);
  currentTimeIs(initialTime);  
  probDesc.getNewmarkParameters(beta, gamma, alphaf, alpham);
}

NLDynamTimeIntegrator::~NLDynamTimeIntegrator()
{
  delete refState;
  delete stepState;
  delete geomState;
}

// Public methods
// --------------
void
NLDynamTimeIntegrator::currentTimeIs(double newTime)
{
  currTime = newTime;
  midTime = newTime + localDt;
  rhs.zero();
  aeroForce.zero();
}

void
NLDynamTimeIntegrator::integrate(int numSteps)
{
  if (currStep == 0) {
    postProcessor().dynamOutput(geomState, velocity, dummyVp, currTime, -1, external_force, aeroForce, acceleration, refState);
  }
  
  // Begin time-marching
  for (int lastStep = currStep + numSteps; currStep < lastStep; ++currStep)
  {
    currTime += localDt;
    midTime = currTime - localDelta;
    if (verboseFlag)
    {
       fprintf(stderr,"\n**** Begin Time Step %d - Time %e ****\n", currStep + 1, currTime);
    }
    probDesc.getExternalForce(external_force, gravityForce, currStep, midTime, geomState, elementInternalForce, aeroForce, localDelta);
    stateIncr.zero();

    // Newton iterations
    int converged;
    double currentRes, resN;
    for (int iter = 0; iter < maxNumIter; ++iter, ++totalNewtonIter)
    {
      residual = external_force;
      probDesc.getStiffAndForce(*geomState, residual, elementInternalForce, midTime); // Update elementary matrices & residual force
      probDesc.reBuild(*geomState, iter, localDelta, midTime); // Assemble [Kt] and factor ([M] + delta^2 * [Kt])
      geomState->get_inc_displacement(inc_displac, *stepState, probDesc.getZeroRot()); // Compute incremental displacement
      resN = probDesc.formRHScorrector(inc_displac, velocity, acceleration, residual, rhs, geomState, localDelta); // rhs = delta^2 * residual - [M] (inc_displac - delta * velocity) 
      if (iter == 0)
      {
        currentRes = resN;
      }
      residual = rhs;
      probDesc.getSolver()->reSolve(rhs); // Solve ([M] + delta^2 [Kt]) du = rhs
      converged = probDesc.checkConvergence(iter, resN, residual, rhs, currTime);
      stateIncr = rhs;
      if (converged == 1)
      {
        break;
      }
      else if (converged == -1)
      {
        fprintf(stderr," ... WARNING, Solution diverged\n");
      }
      geomState->update(stateIncr);
    }

    if (converged == 0)
    {
      fprintf(stderr," *** WARNING: Newton solver did not reach convergence after %d iterations (res = %e, target = %e)\n", maxNumIter, currentRes, probDesc.getTolerance());
    }

    geomState->midpoint_step_update(velocity, acceleration, localDelta, *stepState, beta, gamma, alphaf, alpham, probDesc.getZeroRot());

    postProcessor().dynamOutput(geomState, velocity, dummyVp, currTime, currStep, external_force, aeroForce, acceleration, refState);
  }
}

