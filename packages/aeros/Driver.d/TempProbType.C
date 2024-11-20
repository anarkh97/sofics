#include <Driver.d/TempProbType.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Math.d/TVectorSet.h>
#include <Utils.d/DistHelper.h>

#include <limits>

extern int verboseFlag;

template<
     class DynOps,             // Data Structure for K, M and dynMat
     class VecType,
     class PostProcessor,
     class ProblemDescriptor>
void
TempSolver<
     DynOps,
     VecType,
     PostProcessor,
     ProblemDescriptor>
::solve()
{

   // Take care of all the renunmbering, connectivities, dofsets
   // and boundary condition. Suppose to be common for
   // all the different dynamic scheme.

   probDesc->preProcess();
 
   postProcessor = probDesc->getPostProcessor();

   VecType d_n(probDesc->solVecInfo());
   VecType v_n(probDesc->solVecInfo());
   VecType v_p(probDesc->solVecInfo());

   // Set up the initial conditions

   TempState<VecType> curState(d_n, v_n, v_p);
   probDesc->tempInitState(curState);
  
   if(probDesc->getAeroheatFlag() >= 0) probDesc->aeroHeatPreProcess(d_n, v_n, v_p);

   if(probDesc->getThermohFlag() >= 0) probDesc->thermohPreProcess(d_n, v_n, v_p);

   // Get time loop information

   double tmax, dtemp;
   probDesc->getTempTimes(dtemp, tmax);

   // Get Time Integration Scheme

   int algType = probDesc->getTimeIntegration();

   switch ( algType )
   {
     // Newmark
     case 0:
       NewmarkTempLoop(curState, dtemp, tmax);
       break;

     // Quasi-Static
     case 1:
       filePrint(stderr, " ... Quasistatics                   ...\n");
       probDesc->getQuasiStaticParameters(maxVel, qsbeta);
       probDesc->getSteadyStateParam(steadyFlag,steadyMin,steadyMax,steadyTol);
       quasistaticLoop(curState, dtemp, tmax);
   }

} 

// -------------------
//
//  NEWMARK TIME LOOP
//
// -------------------
template<
     class DynOps,
     class VecType,
     class PostProcessor,
     class ProblemDescriptor>
void
TempSolver<
     DynOps,
     VecType,
     PostProcessor,
     ProblemDescriptor>
::NewmarkTempLoop(TempState<VecType>& curState, double dt, double tmax)
{ 
  double gamma;
  probDesc->getTempNeum(gamma);
  if(gamma == 0)
    filePrint(stderr, " ... Explicit Forward Euler Method  ...\n"
                      " ... i.e. α = 0                     ...\n");
  else if(gamma == 0.5)
    filePrint(stderr, " ... Implicit Midpoint Rule         ...\n"
                      " ... i.e. α = ½                     ...\n");
  else if(gamma == 1.0)
    filePrint(stderr, " ... Implicit Backward Euler Method ...\n"
                      " ... i.e. α = 1                     ...\n");
  else
    filePrint(stderr, " ... Imp. Generalized Midpoint Rule ...\n"
                      " ... with α = %5.3f                 ...\n", gamma);

  VecType &d_n = curState.getDisp();
  VecType &v_n = curState.getVeloc();
  VecType &v_n_p = curState.getPrevVeloc();

  // Some vectors for Newmark time loop
  VecType rhs(probDesc->solVecInfo());
  VecType tmp(probDesc->solVecInfo());
  VecType ext_f(probDesc->solVecInfo());
  VecType prev_f(probDesc->solVecInfo());

  prev_f.zero() ;

  // ... Built LHS Matrix (M + 0*D +gamma*dt*K)
  DynOps dynOps = probDesc->buildOps(1.0, 0., gamma*dt);

  // Project in initial conditions in case of filter is activated
  if(probDesc->getHzemFlag()) {
    probDesc->tempProject(d_n);
    probDesc->tempProject(v_n);
  }

  // Get initial Time
  int n = 0;
  double t = 0.0;
  probDesc->getInitialTime(n, t);

  // compute the external force at the initial time
  probDesc->computeExtForce(ext_f, t, -1, prev_f);

  // Compute the initial temperature gradient: v^0 = M^{-1}*(f^0 - K*d^0) if requested and possible and not unnecessary
  if(domain->solInfo().iacc_switch && (dynOps.Msolver || gamma == 0.0) && !geoSource->getCheckFileInfo()->lastRestartFile) {
    if(verboseFlag) filePrint(stderr," ... Computing initial first time derivative of temperature ...\n");
    dynOps.K->mult(d_n, tmp);
    v_n.linC(ext_f, -1.0, tmp);
    if(gamma == 0.0) dynOps.dynMat->reSolve(v_n);
    else dynOps.Msolver->reSolve(v_n);
  }

  // Output the initial state: d^0, v^0, f^0
  postProcessor->tempdynamOutput(n, dynOps, ext_f, curState);

  // ... BEGIN MAIN TIME-LOOP
  double s0 = -getTime(), s1 = -51, s2 = 0;
  char ch[4] = { '|', '/', '-', '\\' };

  bool coupled = (probDesc->getAeroheatFlag() >= 0 || probDesc->getThermohFlag() >= 0);
  int printNumber = (domain->solInfo().printNumber > 0) ? domain->solInfo().printNumber : std::numeric_limits<int>::max();
  if(!coupled && printNumber < std::numeric_limits<int>::max()) filePrint(stderr, " ⌈\x1B[33m Time Integration Loop In Progress: \x1B[0m⌉\n");

  for( ; t < tmax-0.01*dt; t += dt, s2 = s0+getTime()) {

    if(!coupled && (s2-s1 > printNumber)) {
      s1 = s2;
      filePrint(stderr, "\r ⌊\x1B[33m %c t = %9.3e Δt = %8.2e %3d%% \x1B[0m⌋",
                ch[int(s1/250.)%4], t+dt, dt, int((t+dt)/(tmax-0.01*dt)*100));
    }

    // Mode decomposition of displacement
    if(probDesc->getModeDecompFlag()) probDesc->modeDecomp(t, n, d_n);

    // Compute external flux: f^{n+1/2}
    probDesc->computeExtForce(ext_f, t+dt/2.0, n, prev_f);

    if(gamma != 0.0) { // IMPLICIT
      // Solve for temperature: d^{n+1/2} = (M + gamma*dt*K)^{-1}(gamma*dt*f^{n+1/2} + M*(d^n+dt/2*(1-2*gamma)*v^n))
      if(gamma != 0.5) {
        tmp.linC(d_n, dt/2.0*(1.0-2.0*gamma), v_n);
        dynOps.M->mult(tmp, rhs);
      }
      else dynOps.M->mult(d_n, rhs);
      rhs.linAdd(gamma*dt, ext_f);
      dynOps.dynMat->reSolve(rhs);

      if(probDesc->getHzemFlag() == 2) probDesc->tempProject(rhs);

      // Extrapolate temperature solution to t^{n+1} : d^{n+1} = 2*d^{n+1/2} - d^n
      d_n.linC(2.0, rhs, -1.0, d_n);

      // Compute the first time derivative of temperature at t^{n+1}: v^{n+1} = 2/(gamma*dt)*(d^{n+1/2} - d^n) - (1-gamma)/(gamma)*v^n
      v_n_p.linC(2.0/(gamma*dt), d_n, -2.0/(gamma*dt), rhs);
      if(gamma != 1.0) v_n_p.linAdd(-(1.0-gamma)/gamma, v_n);
    }
    else { // EXPLICIT if gamma = 0.0 and M is diagonal
      // Solve for first time derivative of temperature: v^{n+1/2} = (M + gamma*dt*K)^{-1}(f^{n+1/2} - K*(d^n + dt/2*(1-2*gamma)*v^n))
      //                                                           = M^{-1}(f^{n+1/2} - K*(d^n+dt/2*v^n)) if gamma == 0
      tmp.linC(-1.0, d_n, -dt/2.0, v_n);
      if(domain->solInfo().isNonLin()) probDesc->getInternalForce(tmp, rhs);
      else dynOps.K->mult(tmp, rhs);
      rhs += ext_f;
      dynOps.dynMat->reSolve(rhs);

      if(probDesc->getHzemFlag() == 2) probDesc->tempProject(rhs);

      // Extrapolate the first time derivative of temperature to t^{n+1}: v^{n+1} = 2*v^{n+1/2} - v^n
      v_n_p.linC(2.0, rhs, -1.0, v_n);

      // Compute the temperature at t^{n+1}: d^{n+1} = d^n + dt*(2*gamma*v^{n+1/2} + (1-2*gamma)*v^n)
      //                                             = d^n + dt*v^n if gamma == 0
      d_n.linAdd(dt, v_n);
    }

    // Now swap v_n_p <-> v_n
    v_n.swap(v_n_p);

    // Increment time index
    n++;

    // Output d^n, v^n, f^{n-1/2}
    postProcessor->tempdynamOutput(n, dynOps, ext_f, curState);

  }
  if(!coupled && printNumber < std::numeric_limits<int>::max())
    filePrint(stderr, "\r ⌊\x1B[33m   t = %9.3e Δt = %8.2e 100%% \x1B[0m⌋\n", t, dt);

#ifdef PRINT_TIMERS
  if(verboseFlag) filePrint(stderr, " ... Total Loop Time = %.2e s   ...\n", s2/1000.0);
#endif
}

// ----------------------------
// ... Quasi-Static Time Loop
// ----------------------------

template< class DynOps,        
          class VecType,
          class PostProcessor, 
          class ProblemDescriptor>
void
TempSolver< DynOps, VecType, PostProcessor, ProblemDescriptor >
::quasistaticLoop(TempState<VecType>& curState, double dtemp, double tmax)
{
   VecType &d_n = curState.getDisp();

   // Some vectors

   VecType deltad_n_p(probDesc->solVecInfo());
   VecType deltad_n(probDesc->solVecInfo());
   VecType rhs(probDesc->solVecInfo());
   VecType ext_f(probDesc->solVecInfo());
   VecType prev_f(probDesc->solVecInfo());
   int i;

   int tIndex   = 0;
   int iSteady  = 0;

   // double startTime=0;
  
   double forceRef;
   prev_f.zero();
   deltad_n.zero();

   double relaxFac = maxVel;
  
   // Build K
   DynOps dynOps = probDesc->buildOps(0, 0, 1.);

   VecSet<VecType> RBM(dynOps.dynMat->numRBM(), probDesc->solVecInfo());
   //if(dynOps.dynMat->numRBM() > 0)
   //   fprintf(stderr," ... Number of ZEM        =   %d     ...\n", dynOps.dynMat->numRBM());

   // Case  K singular
   Vector alpha(dynOps.dynMat->numRBM());
   alpha = 0.0;
   if(dynOps.dynMat->numRBM() > 0) {
     dynOps.dynMat->getRBMs( &RBM[0]);
     for(i = 0; i < dynOps.dynMat->numRBM(); ++i) {
       RBM[i] *= 1.0/RBM[i].norm();
       alpha[i] = d_n*RBM[i];
     }
   }

   // Output state of model
   postProcessor->tempdynamOutput(tIndex, dynOps, ext_f, curState);

   //-----------------------------------------------------------
   // ... BEGIN MAIN TIME-LOOP
   //-----------------------------------------------------------

   for(tIndex = 1; tIndex <= steadyMax; tIndex++) {

     // ... Compute external force at time t
     probDesc->computeExtForce(ext_f, (double) tIndex, tIndex, prev_f);

    // ... build force reference norm
    if(tIndex == 1) forceRef = ext_f.norm();

    // ... build internal force
    probDesc->getInternalForce(d_n,rhs);

    // ... check for convergence
    double relres = norm(rhs-ext_f) / forceRef;
    if(relres <= steadyTol) iSteady = 1;

    fprintf(stderr," ... Pseudo-Step = %d  Rel. Res. = %10.4e\n", tIndex, relres);

    // command communication with fluid or structure for coupled aerothermal or thermoelastic
    if(probDesc->getAeroheatFlag() >= 0 /*|| probDesc->getThermohFlag() >= 0*/) {
      if(tIndex == steadyMax  && !iSteady) {
        probDesc->cmdComHeat(1);          
        break;
      }
      else {
        iSteady = probDesc->cmdComHeat(iSteady);
      }
    }

    // ... stop quasi-transient simulation if converged
    if(iSteady) {
      fprintf(stderr," --------------------------------------\n");
      fprintf(stderr," ... Quasistatic Analysis Converged After %d Steps \n",tIndex);
      fprintf(stderr," --------------------------------------\n");
      break;
    }

    // ... save load vector
    rhs=ext_f;
    if(dynOps.dynMat->numRBM() > 0) {
      for(i = 0; i < dynOps.dynMat->numRBM(); ++i) {
        fprintf(stderr, "Alpha going from %e to ", alpha[i]);
        alpha[i] += qsbeta * (ext_f*RBM[i]);
        fprintf(stderr, "%e\n",alpha[i]);
      }
    }

    for(i = 0; i < dynOps.dynMat->numRBM(); ++i)
      rhs.linAdd(-(rhs*RBM[i]), RBM[i]);

    // ... solve System for current load, rhs now contains solution
    dynOps.dynMat->reSolve(rhs);

    for(i = 0; i < dynOps.dynMat->numRBM(); ++i)
      rhs.linAdd(-(rhs*RBM[i]), RBM[i]);

    // ... add R*alphaR if K is singular
     
    // ... compute displacement increment;
    deltad_n_p.linC(rhs, -1.0, deltad_n);

    // ... apply relaxation factor
    deltad_n_p *= relaxFac;

    // ... update solution
    deltad_n += deltad_n_p;

    d_n = deltad_n;
    if(dynOps.dynMat->numRBM() > 0) {
      for(i = 0; i < dynOps.dynMat->numRBM(); ++i)
        d_n.linAdd(alpha[i], RBM[i]);
    }

    // ... output current solution and send displacements to fluid
    postProcessor->tempdynamOutput( tIndex, dynOps, ext_f, curState );

  }

   // ... check for convergence
   if(!iSteady) {
     fprintf(stderr," --------------------------------------\n");
     fprintf(stderr," ... Quasistatic Analysis Did Not Converge After %d Steps \n",tIndex);
     fprintf(stderr," --------------------------------------\n");
   }

}

//------------------------------------------------------------------------------
