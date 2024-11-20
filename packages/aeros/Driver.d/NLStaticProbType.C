#include <cstdio>
#include <limits>
#include <Timers.d/GetTime.h>

extern int verboseFlag;
extern int totalNewtonIter;
extern SolverInfo &solInfo;

template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
           class StateUpdate >
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >::~NLStaticSolver()
{
  if(geomState) delete geomState;
}

template < class OpSolver, 
           class VecType,  
           class PostProcessor, 
           class ProblemDescriptor, 
           class GeomType,
           class StateUpdate >
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >::solve()
{
 // Set up nonlinear problem 
 probDesc->preProcess(false);

 // Allocate appropriate vectors to store external force and residual
 VecType force(probDesc->solVecInfo());
 VecType residual(probDesc->solVecInfo());
 VecType totalRes(probDesc->solVecInfo());

 VecType elementInternalForce(probDesc->elemVecInfo());

 residual.zero();
 totalRes.zero();
 elementInternalForce.zero();

 // Get right hand side (external force)
 probDesc->getRHS(force);

 // Initialize geometric state of problem
 geomState = probDesc->createGeomState();
 stateIncr = StateUpdate::initInc(geomState, &residual);
 
 bool failSafe = solInfo.getNLInfo().failsafe;
 refState = (solInfo.soltyp == 2 && !failSafe) ? 0 : StateUpdate::initRef(geomState);

 double lambda = 0.0;

 // Incremental step parameter for load / prescribed displacement control 
 double deltaLambda0 = probDesc->getDeltaLambda0(); 

 // Newton iteration loop
 double maxLambda = probDesc->getMaxLambda();

 // Output structure initial configuration
 if(deltaLambda0 != maxLambda)
   probDesc->staticOutput(geomState, lambda, force, totalRes, refState);

 int numIter = 0;

 int numSteps = int(0.5+(maxLambda / deltaLambda0));

 int step = 1;
 double deltaLambda;
 bool failed = false;
 int numConverged = 0;
 int p = step, q = 1;

 for( ; step <= numSteps || failed; ) {

   deltaLambda = deltaLambda0/q;

   filePrint(stderr, " --------------------------------------\n");
   filePrint(stderr, " ... Newton : Start Step #%d --- Lambda = %e\n", step, lambda+deltaLambda);

   if(failSafe && failed) {
     StateUpdate::copyState(refState, geomState);
   } else
   if(solInfo.soltyp != 2) StateUpdate::copyState(geomState, refState);

   probDesc->updatePrescribedDisplacement(geomState, lambda+deltaLambda);

   probDesc->initializeParameters(step,geomState); // for augmented Lagrangian (a) Lagrange multipliers are initialized at step 1 only
                                                   // unless SolverInfo::reinit_lm is true, (b) penalty parameter is always reset at each step
   int converged;
   bool feasible;
   double err, resN;
   for(int i = 0; i < solInfo.num_penalty_its; ++i) {

     // call newton iteration with load step lambda
     try {
       converged = newton(force, residual, totalRes,
                          elementInternalForce, probDesc,
                          refState, geomState, stateIncr, numIter, resN, lambda+deltaLambda, step);
     }
     catch(std::runtime_error& e) {
       if(!failSafe) std::cerr << "\nexception: " << e.what() << std::endl;
       converged = 0;
       break;
     }

     if(converged == 1) {
       filePrint(stderr," ... Newton : Step #%d, Iter #%d converged after %d iterations\n",
                 step, i+1, numIter+1);
     }
     else if(converged == -1) {
       if(!failSafe) {
         filePrint(stderr," ... Newton : Step #%d, Iter #%d diverged after %d iterations\n",
                   step, i+1, numIter+1);
         filePrint(stderr," ... Newton : analysis interrupted by divergence\n");
       }
       break;
     } 
     else if(converged == 0 && !failSafe && solInfo.getNLInfo().stepUpdateK != std::numeric_limits<int>::max()) {
       filePrint(stderr," *** WARNING: Newton solve did not converge after %d iterations (res = %e, target = %e)\n",
                 numIter, resN, probDesc->getTolerance());
     }

     if((failed = (failSafe && converged != 1 && resN > solInfo.getNLInfo().failsafe_tol && q < (std::numeric_limits<int>::max() >> 1)))) {
       // if a Newton solve fails to converge, terminate constraint enforcement iterations
       break;
     }

     // check constraint violation error
     feasible = probDesc->checkConstraintViolation(err, geomState);

     // update lagrange multipliers and/or penalty parameters
     if((!feasible && i+1 < solInfo.num_penalty_its) || (solInfo.lm_update_flag == 1)) {
       probDesc->updateParameters(geomState);
     }

     filePrint(stderr, " ... Newton : End Step #%d, Iter #%d --- Max Steps = %d, Max Iters = %d\n",
               step, i+1, numSteps, solInfo.num_penalty_its);
     if(err > std::numeric_limits<double>::epsilon()) filePrint(stderr," ... Constraint violation: %8.2e ...\n", err);
     filePrint(stderr," --------------------------------------\n");
   
     if(feasible) break;

   } // end of constraint enforcement iteration loop

   if(failed) {
     // if a Newton solve fails to converge, decrease the load step and try again
     filePrint(stderr," ... Failsafe: bisecting load step  ...\n");
     p *= 2;
     q *= 2;
     numConverged = 0;
     continue;
   }

   if(solInfo.soltyp != 2) probDesc->updateStates(refState, *geomState, lambda+deltaLambda);

   p++;
   numConverged++;
   if(numConverged >= 2 && (q-p)%2 == 0 && q >= 2) { p /= 2; q /= 2; numConverged = 0; } // increase the load step

   // increment load parameter
   lambda += deltaLambda;

   // Output current load step results
   if(p%q == 0) { // finished a whole step
     step = p/q;
     lambda = (step-1)*deltaLambda0;
     double time = (deltaLambda0 == maxLambda) ? 0.0 : lambda;
     probDesc->staticOutput(geomState, time, force, totalRes, refState);
   }

   // Exit loop in the case of divergence
   if(converged == -1) break;
 }

 if(refState) delete refState;
 if(stateIncr) delete stateIncr;

 probDesc->printTimers();
}


template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
	   class StateUpdate>
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >
::arclength()
{
 // WARNING THE CREATION OF REFSTATE MUST BE DONE IF THIS IS TO BE USED
 // WITH TOTALUPDATE MODE
 // Set up  nonlinear Problem
 probDesc->preProcess();

 // Allocate Vectors
 VecType force(probDesc->solVecInfo());
 VecType residual(probDesc->solVecInfo());
 VecType totRes(probDesc->solVecInfo());
 VecType arcLenResid(probDesc->solVecInfo());
 VecType pVec(probDesc->solVecInfo());

 VecType elementInternalForce(probDesc->elemVecInfo());

 // HB
 residual.zero();
 totRes.zero();
 elementInternalForce.zero();

 // Get right hand side
 probDesc->getRHS(force);

 // Compute initial force norm
 double forceNorm = force.norm();

 // Initialize geometric state of Problem
 GeomType *u0 = probDesc->createGeomState();
 GeomType *u  = probDesc->createGeomState();
 stateIncr = StateUpdate::initInc(u, &residual);

 // ... Output initial configuration
 probDesc->staticOutput(u, 0.0, force, totRes, u0);

 // ... DEFINE deltaLambda0
 double deltaLambda0 = probDesc->getDeltaLambda0();

 // ... DEFINE minimum delta Lambda and maximum delta Lambda
 // lmin <= deltaLamba <= lmax
 double lfactor = probDesc->getScaleFactor();
 double lmin = deltaLambda0 / lfactor;
 double lmax = deltaLambda0 * lfactor;

 // ... DEFINE initial lambda0
 double lambda0 = 0.0;

 // ... COMPUTE FIRST LOAD CONTROL PARAMETER
 double lambda = lambda0 + deltaLambda0;

 int numIter = 0;
 int step    = 1;
 double resN;
 newton(force, residual, totRes, elementInternalForce, probDesc, u0,
        u, stateIncr, numIter, resN, lambda, step);

 // ... Declare Vector dU
 VecType dU(probDesc->solVecInfo());

 // ... Define initial dU = u - u0
 u->diff(*u0, dU);

 // ... Define deltaS   (arc length distance)
 double deltaS = dU.norm();
 filePrint(stderr," -> deltaS = %e\n",deltaS);
 double deltaSmax = 10*deltaS;
 double deltaSmin = 0.1*deltaS;

 // ... COMPUTE W = SMOOTHING PARAMETER
 //double w = (deltaS*deltaS) / (deltaLambda0*deltaLambda0);
 // HB: it seems to work but for me this is wrong !!!
 //deltaS = deltaS*deltaS;
 double w = (deltaS*deltaS) / (deltaLambda0*deltaLambda0);
 w = 0;

 // ... DEFINE NU = load control parameter multiplier
 double nu = 1.0;
 int extMax = solInfo.getNLInfo().extMax;
 int extMin = solInfo.getNLInfo().extMin;

 // ... DEFINE MAXIMUM NUMBER OF TRAJECTORY ITERATIONS
 int maxNumTrajectory = 2000;

// ---------------------------------------------------

 double deltaLambda;
 int iter, numExtIter;

  // ... COMPUTE TRAJECTORY (EQUILIBRIUM) PATH
 for(iter=0; iter<maxNumTrajectory; ++iter) {

   probDesc->staticOutput(u, lambda, force, totRes, u0);

   // Compute: dU = nu*(u - u0);
   u->diff(*u0, dU);
   filePrint(stderr," -> ||dU|| = %e\n",dU.norm());
   //dU *= nu;

   *u0 = *u;

   // ... COMPUTE NEXT DELTALAMBDA
   //deltaLambda = nu*(lambda - lambda0);
   deltaLambda = (lambda - lambda0);

   // ... CHECK MAGNITUDE OF DELTALAMBDA
   if(deltaLambda > lmax ) deltaLambda = lmax;
   if(deltaLambda < lmin ) deltaLambda = lmin;

   // ... SET lambda0 = lambda
   lambda0 = lambda;

   // ... COMPUTE DELTAS
   double deltaS0 = deltaS;
   deltaS *= nu;
   if(deltaS > deltaSmax) deltaS = deltaSmax;
   if(deltaS < deltaSmin) deltaS = deltaSmin;

   filePrint(stderr,"**** Equilibrium #%d deltaLambda = %10.6f lambda = %10.6f\n"
                   , step, deltaLambda, lambda);
   step++;

   predictorStep(*u, *u0, dU, lambda, deltaLambda, deltaS, deltaS0, w,
                 force, residual,totRes, elementInternalForce,  pVec, step);

   //u->diff(*u0, dU);
   filePrint(stderr," -> ||dU|| = %e\n", dU.norm());

   // ... CALL EXTENDED NEWTON
   extendedNewton(*u, *u0, dU, lambda, deltaLambda, deltaS, w, numExtIter,
                  force, residual, totRes, arcLenResid, forceNorm,
                  elementInternalForce, pVec, step);

   probDesc->updateStates(u0, *u, lambda);

   if( std::abs(lambda) >= probDesc->getMaxLambda() ) break;

   //...DETERMINE CONTROL PARAMETER BASED ON # OF ITERATIONS IN EXTENDED NEWTON
   nu = 1.0;
   //HB: disable this for debug
   //nu = sqrt(4./numExtIter);
   //nu = pow(4./numExtIter,0.75);
   if(numExtIter < extMin) nu = 2.0;
   if(numExtIter > extMax) nu = 0.5;
 }
 //HB
 filePrint(stderr,"**** Equilibrium #%d deltaLambda = %10.6f lambda = %10.6f\n"
                 , step, deltaLambda, lambda);


 // ... DEFINE INITIAL GUESS
 u->interp((1.0 - lambda0)/(lambda - lambda0), *u, *u0);
 lambda = 1.0;

 // ... CALL NEWTON FOR FINAL SOLUTION
 *u0 = *u;
 step++;
 newton(force, residual, totRes, elementInternalForce, probDesc, u0, u, stateIncr, numIter, resN, lambda, step);

 probDesc->updateStates(u0, *u, lambda);

 // CALL POST PROCESSING OF DISPLACEMENTS
 probDesc->staticOutput(u, lambda, force, totRes, u0);
}

template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
	   class StateUpdate>
int
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >
::newton(VecType& force, VecType& residual, VecType &totalRes,
         VecType& elementInternalForce, ProblemDescriptor *probDesc,
         typename StateUpdate::RefState *refState,
         GeomType* geomState, typename StateUpdate::StateIncr *stateIncr,
         int &numIter, double &residualNorm, double lambda, 
         int step)
{
  // Accumulate time spent in solving and geomstate update for one step
  double timeSolve   = 0.0;
  double timeStiff   = 0.0;
  double timeRebuild = 0.0;

  int maxit = probDesc->getMaxit(), rebuildFlag;
  bool useTolInc = (solInfo.getNLInfo().tolInc != std::numeric_limits<double>::infinity()
                 || solInfo.getNLInfo().absTolInc != std::numeric_limits<double>::infinity());
  double normDv;
  
  // Zero the state increment
  StateUpdate::zeroInc(stateIncr);
  
  // Main Newton Iteration Loop
  VecType residualCopy(probDesc->solVecInfo());
  int iter, converged;
  for(iter = 0; iter < maxit; ++iter) {

    // residual = lambda*force;
    residual.linC(force, lambda);
 
    // Update geomState then compute current tangent stiffness and residual force (including follower force contributions)
    timeStiff -= getTime();
    StateUpdate::integrate(iter, probDesc, refState, geomState, stateIncr,
                           residual, elementInternalForce, totalRes, lambda);
    timeStiff += getTime();
#ifdef PRINT_TIMERS
    filePrint(stderr,"  Rebuild Element Stiffness & Internal Force time = %13.4f s\n", 
              timeStiff/1000.0);
#endif    

    // Rebuild tangent stiffness matrix when necessary
    if(solInfo.mpcDirect) {
      timeRebuild -= getTime();
      rebuildFlag = probDesc->reBuild(iter, step, *geomState);
      timeRebuild += getTime();
#ifdef PRINT_TIMERS
      filePrint(stderr,"  Rebuild Tangent Stiffness Matrix time = %13.4f s\n",
                timeRebuild/1000.0);
#endif
    }

    // note the residual norm needs to be computed after reBuild for "constraints direct"
    residualNorm = probDesc->getResidualNorm(residual, *geomState);

    // If the convergence criteria does not involve the solution increment, then 
    // check for convergence now (to avoid potentially unnecessary solve)
    if(useTolInc || !(converged = probDesc->checkConvergence(iter, normDv, residualNorm)) ) {

      // Copy residual if necessary before it gets overwritten
      if(probDesc->linesearch().type != 0) residualCopy = residual;

      // Rebuild tangent stiffness matrix when necessary
      if(!solInfo.mpcDirect) {
        timeRebuild -= getTime();
        rebuildFlag = probDesc->reBuild(iter, step, *geomState);
        timeRebuild += getTime();
#ifdef PRINT_TIMERS
        filePrint(stderr,"  Rebuild Tangent Stiffness Matrix time = %13.4f s\n",
                  timeRebuild/1000.0);
#endif
      }

      if(rebuildFlag) {
        filePrint(stderr," ... Newton : Iter #%d --- Rebuild Tangent Stiffness (res = %e)\n", iter+1, residualNorm);
      }
      else filePrint(stderr," ... Newton : Iter #%d (res = %e)\n", iter+1, residualNorm);

      totalNewtonIter++;

      // Solve current system Kt*u = residual, overwrite residual with u
      timeSolve -= getTime();
      probDesc->getSolver()->reSolve(residual);
      timeSolve += getTime();
#ifdef PRINT_TIMERS
      filePrint(stderr,"  Solve Incremental Displacement %13.4f s\n",
                timeSolve/1000.0);
#endif

      if(probDesc->linesearch().type != 0) {
        // Optional adjustment of the step length
        StateUpdate::linesearch(probDesc, refState, geomState, stateIncr,
                                residualCopy, elementInternalForce, totalRes, lambda, force, residual);
      }
      StateUpdate::updateIncr(stateIncr, residual);

      // Update state here if the maximum number of iterations is reached
      if(iter == maxit-1) StateUpdate::updateState(probDesc, geomState, *stateIncr);

      // Compute incremental displacement norm
      normDv = residual.norm();
    }

    // If the converged criteria does involve the solution increment, then
    // check for convergence now
    if(useTolInc) {
      converged = probDesc->checkConvergence(iter, normDv, residualNorm);
    }
    else if(converged) {
      filePrint(stderr," ... Newton : Iter #%d (res = %e)\n", iter+1, residualNorm);
    }

#ifdef DEBUG_NEWTON
    probDesc->staticOutput(geomState, double(iter), force, totalRes, refState);
#endif

    // If converged or diverged, break out of loop
    if(converged) break;
  }

  // return with the number of iterations newton took to converge/diverge
  numIter = iter;

  return converged;
}

// HB
template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
           class StateUpdate>
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >
::predictorStep(GeomType &u, GeomType &u0, VecType &dU, double &lambda, double &deltaLambda,
                double &deltaS, double &deltaS0, double w,
                VecType& force, VecType& residual, VecType &totRes, VecType& elementInternalForce, 
                VecType& duds, int step)
{
  filePrint(stderr," ### GET IN PREDICTOR STEP\n");
  //HB: be aware that for multi-domain problem most of the dot product here are not "consistent"
  //    in the sense that they won't give the same answer as for the same non multi-domains problem
  //    (see the way a dot product is computed for a GenDistrVector). But in the case were the 
  //    distributed vector has the same value on its share nodes, the computed dot product can be
  //    interpreted as a scaled/weighted dot product: this is why it should be ok to use it here
  //    But be carefull, if we use this with a residual-like vector (i.e. a "disassembled-like vector") ...
  probDesc->getRHS(force); //PJSA nonlinear external force computed in getStiffAndForce probDesc->getRHS(force, &u);
  residual.zero();
  double residualNorm = probDesc->getStiffAndForce(u, residual, elementInternalForce, totRes, lambda, &u0);
  int rebuildFlag = probDesc->reBuild(0, step, u);
  if(rebuildFlag) { filePrint(stderr," ... ||res|| = %e\n", residualNorm); }

  duds.zero(); 
  probDesc->getSolver()->solve(force, duds);

  double m = 2*(dU*duds);
  double n = 2*w*deltaLambda*force.sqNorm();
  //double dlds = 1./(m+n);
  double dlds = 2.*deltaS0/(m+n);
  filePrint(stderr," -> m = %e, n = %e, dlds = %e\n",m,n,dlds);
  duds *= (dlds*deltaS);
  filePrint(stderr," -> ||duds||    = %e\n",duds.norm());

  u.update(duds);
  deltaLambda = dlds*deltaS;
  lambda += deltaLambda;
  filePrint(stderr," -> deltaLambda = %e\n",deltaLambda);
  u.diff(u0,dU);
  filePrint(stderr," -> ||dU||      = %e\n",dU.norm());
  //dU = duds;
}

template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
	   class StateUpdate>
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >
::extendedNewton(GeomType &u, GeomType &u0, VecType &dU, double &lambda, double deltaLambda, 
                 double &deltaS, double w, int &numExtIter,
                 VecType& force, VecType& residual,  VecType &totRes,
                 VecType& arcLenResid, 
                 double forceNorm, VecType& elementInternalForce, VecType& pVec, int step)
{

  //HB
  filePrint(stderr," #####################################\n");
  filePrint(stderr," ### GET IN EXTENDED-NEWTON SOLVER ###\n");
  filePrint(stderr," #####################################\n");
  // Compute new lambda
  double lambda0 = lambda-deltaLambda;
  //lambda += deltaLambda;

  // Compute first norm
  double firstNorm = (lambda*lambda)*forceNorm;
  //HB
  filePrint(stderr,"--- firstNorm = %e\n",firstNorm);


  if(firstNorm == 0.0) return;

  //predictorStep(u, lambda, deltaLambda, deltaS, w,
  //              force, residual, pVec, step);


  //double Dlbda = lambda-lambda0;
  double Dlbda = deltaLambda;
  double r0= (dU*dU) + w*(Dlbda*Dlbda) - (deltaS*deltaS);
  //r0 = 1.;

  int maxExtIter = probDesc->getMaxit();
  for(numExtIter = 0; numExtIter < maxExtIter; ++numExtIter) {

    filePrint(stderr," ------------------------------------\n");
    filePrint(stderr," ### Extended-Newton iteration %d ###\n",numExtIter);
    filePrint(stderr," ------------------------------------\n");
    // HB 
    //probDesc->getRHS(force, &u);

    // Compute residual = lambda*force
    residual.linC(force, lambda);

    // Compute stiffness and residual force
    double residualNorm = probDesc->getStiffAndForce(u, residual, elementInternalForce, totRes, lambda, &u0);

    // rebuild tangent stiffness matrix when necessary
    // step # should be the second argument in probDesc->reBuild() 
    //probDesc->reBuild(numExtIter, 1, u);
    int rebuildFlag = probDesc->reBuild(numExtIter, step, u); 

    if(rebuildFlag) { filePrint(stderr," ... ||res|| = %e\n", residualNorm); }

    //HB
    Dlbda = lambda-lambda0;
    double r = (dU*dU) + w*(Dlbda*Dlbda) - (deltaS*deltaS);
    filePrint(stderr," r = %e, |r/r0| = %e,|dl| = %e, ||dU|| = %e\n",r, fabs(r/r0), fabs(Dlbda), sqrt(dU*dU));

    arcLenResid = residual;

    filePrint(stderr," ### First solve ###\n");
    probDesc->getSolver()->reSolve(residual);

    filePrint(stderr," ### Second solve ###\n");
    probDesc->getSolver()->solve(force, pVec);

    //HB
    //double mu = -1.0*( dU * residual ) / (w*deltaLambda + dU * pVec);
    double mu = -0.5*( r + 2.*(dU*residual)) / (w*Dlbda + dU*pVec);

    // Compute: pVec = mu*pVec + residual
    pVec.linC(mu, pVec, 1.0, residual);
    //HB
    dU.linC(dU,1.0,pVec);
    //u.diff(*u0,dU);

    // Update geometrical state
    u.update(pVec);

    // Update lambda
    lambda += mu;
    // HB
    filePrint(stderr,"              lambda = %f, mu = %15.6e\n",lambda,mu);

    // Compute arcLenResid = arcLenResid + mu*force
    arcLenResid.linC(arcLenResid, mu, force);

    // double residualNorm = arcLenResid.norm();
    double normDv       = pVec.norm();

    // KHP: same problem as before, norm of force contains 
    //      lagrange multiplier forces and thus norm does not go to zero.
    //int converged = probDesc->checkConvergence(numExtIter, normDv, 
    //                residual.norm());
    // HB
    int converged = probDesc->checkConvergence(numExtIter, normDv, residualNorm);

    totalNewtonIter++;
    if (converged) break;
  }

  filePrint(stderr,"... Number of External Iterations = %d\n", numExtIter);

}
