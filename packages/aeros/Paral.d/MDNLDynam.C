#include <algorithm>
#include <cstdio>
#include <Element.d/MpcElement.d/MpcElement.h>
#include <Threads.d/Paral.h>
#include <Driver.d/Domain.h>
#include <Paral.d/MDNLDynam.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Driver.d/SubDomain.h>
#include <Threads.d/PHelper.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Math.d/Vector.h>
#include <Utils.d/DistHelper.h>
#include <Paral.d/MDStatic.h>
#include <Paral.d/MDDynam.h>
#include <Paral.d/MDOp.h>
#include <Rom.d/EiGalerkinProjectionSolver.h>
#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif
#include <Hetero.d/DistFlExchange.h>
#include <Utils.d/ModeData.h>
#include <Control.d/ControlInterface.h>
#include <Corotational.d/DistrGeomState.h>
#include <Solvers.d/ParallelSolver.h>
#include <Feti.d/Feti.h>
#include <Paral.d/MDDynam.h>
#include <Solvers.d/MultiDomainRbm.h>
#include <Driver.d/SysState.h>

extern int verboseFlag;
extern ModeData modeData;

// ***************************************************************
// *                                                             *
// *  Purpose: Multiple Domain implementation of nonlinear       *
// *           dynamics				                 *	
// *                                                             *
// *  Coded by: Kendall H. Pierson and Philip J. S. Avery        *
// ***************************************************************

void
MDNLDynamic::getConstForce(DistrVector &v)
{
  times->formRhs -= getTime();

  MultiDomainOp mdop(&MultiDomainOp::getConstForce, decDomain->getAllSubDomains(), &v, Kuc);
  threadManager->execParal(decDomain->getNumSub(), &mdop);

  times->formRhs += getTime();
}

int
MDNLDynamic::getInitState(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n, DistrVector &v_p)
{
  MultiDomainOp mdop(&MultiDomainOp::getInitState, decDomain->getAllSubDomains(),
                     &d_n, &v_n, &a_n, &v_p);
  threadManager->execParal(decDomain->getNumSub(), &mdop);

  if(sinfo.filterFlags) {
    project(d_n);
    project(v_n);
    project(a_n);
    project(v_p);
  }

  if(claw && userSupFunc && claw->numSensor) {
    double *ctrdisp = new double[claw->numSensor];
    double *ctrvel  = new double[claw->numSensor];
    double *ctracc  = new double[claw->numSensor];
#ifdef DISTRIBUTED
    for(int i=0; i<claw->numSensor; ++i) ctrdisp[i] = ctrvel[i] = ctracc[i] = std::numeric_limits<double>::min();
#endif
    decDomain->extractControlData(d_n, v_n, a_n, ctrdisp, ctrvel, ctracc);
#ifdef USE_MPI
    structCom->globalMax(claw->numSensor, ctrdisp);
    structCom->globalMax(claw->numSensor, ctrvel);
    structCom->globalMax(claw->numSensor, ctracc);
#endif
    userSupFunc->init(ctrdisp, ctrvel, ctracc);
    delete [] ctrdisp; delete [] ctrvel; delete [] ctracc;
  }

  int aeroAlg = domain->solInfo().aeroFlag;
  if(aeroAlg >= 0)
    aeroPreProcess(d_n, v_n, a_n, v_p);

  if(domain->solInfo().thermoeFlag >= 0)
    thermoePreProcess();

  if(domain->solInfo().thermohFlag >= 0)
    thermohPreProcess(d_n);

  if(domain->solInfo().aeroheatFlag >= 0)
    aeroheatPreProcess(d_n, v_n, v_p);

  return aeroAlg;
}

void
MDNLDynamic::updatePrescribedDisplacement(DistrGeomState *geomState)
{
  times->timePresc -= getTime();
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdatePrescribedDisplacement, *geomState);
  times->timePresc += getTime();
}

void
MDNLDynamic::subUpdatePrescribedDisplacement(int isub, DistrGeomState& geomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  sd->updatePrescribedDisp(geomState[isub]);
}

void
MDNLDynamic::getIncDisplacement(DistrGeomState *geomState, DistrVector &du, DistrGeomState *refState,
                                bool zeroRot)
{
  if(domain->GetnContactSurfacePairs()) {
    du.resize(solVecInfo());
    refState->resize(decDomain);
  }
  geomState->get_inc_displacement(du, *refState, zeroRot); 
}

void
MDNLDynamic::formRHSinitializer(DistrVector &fext, DistrVector &velocity, DistrVector &elementInternalForce, 
                                DistrGeomState &geomState, DistrVector &rhs, DistrGeomState *refState)
{
  // rhs = (fext - fint - Cv)
  rhs = fext;
  elementInternalForce.zero();
  getStiffAndForce(geomState, rhs, elementInternalForce, domain->solInfo().initialTime, refState);
  if(domain->solInfo().order == 2 && C) {
    C->mult(velocity, *localTemp);
    rhs.linC(rhs, -1.0, *localTemp);
  }
}

double
MDNLDynamic::formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                              DistrVector& residual, DistrVector& rhs, DistrGeomState *geomState, double localDelta)
{
  times->correctorTime -= getTime();
  if(domain->GetnContactSurfacePairs()) {
    velocity.conservativeResize(solVecInfo());
    acceleration.conservativeResize(solVecInfo());
    rhs.resize(solVecInfo());
  }

  if(domain->solInfo().order == 1) {
    M->mult(inc_displacement, rhs);
    rhs.linC(localDelta, residual, -1.0, rhs);
  }
  else {
    double beta, gamma, alphaf, alpham, dt = 2*localDelta;
    getNewmarkParameters(beta, gamma, alphaf, alpham);
    if(domain->solInfo().quasistatic) rhs = 0.;
    else {
      localTemp->linC(-(1-alpham)/(1-alphaf), inc_displacement, dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
      M->mult(*localTemp, rhs);
    }
    if(C) {
      localTemp->linC(-dt*gamma, inc_displacement, -dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(*localTemp, rhs);
    }
    rhs.linAdd(dt*dt*beta, residual);
  }

  times->correctorTime += getTime();
  return rhs.norm();
}

void
MDNLDynamic::formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual,
                              DistrVector& rhs, DistrGeomState &geomState,
                              double midtime, double localDelta)
{
  times->predictorTime -= getTime();

  if(claw && userSupFunc) {
    if(claw->numUserDisp > 0) {

      // allocate memory for user defined motion
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineDispLast = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      double *userDefineAcc = new double[claw->numUserDisp];
      for(int i=0; i<claw->numUserDisp; ++i) {
        userDefineVel[i] = 0;
        userDefineAcc[i] = 0;
      }

      // get user defined motion
      userSupFunc->usd_disp(midtime, userDefineDisp, userDefineVel, userDefineAcc);
      userSupFunc->usd_disp(midtime-localDelta, userDefineDispLast, userDefineVel, userDefineAcc);

      // update state
      execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, &geomState, userDefineDisp,
                  userDefineVel, userDefineAcc);

      // get delta disps
      for(int j = 0; j < claw->numUserDisp; j++)
        userDefineDisp[j] -= userDefineDispLast[j];

      // update force residual with KUC
      if(Kuc) 
        execParal(decDomain->getNumSub(), this, &MDNLDynamic::subKucTransposeMultSubtractClaw, residual, userDefineDisp);

      delete [] userDefineDisp; delete [] userDefineDispLast; delete [] userDefineVel; delete [] userDefineAcc;
    }
  }

  if(domain->solInfo().order == 1)
    rhs.linC(residual, localDelta);
  else {
    double beta, gamma, alphaf, alpham, dt = 2*localDelta;
    getNewmarkParameters(beta, gamma, alphaf, alpham);
    // rhs = dt*dt*beta*residual + (dt*(1-alpham)*M - dt*dt*(beta-(1-alphaf)*gamma)*C)*velocity
    //       + (dt*dt*((1-alpham)/2-beta)*M - dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2*C)*acceleration
    localTemp->linC(dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
    M->mult(*localTemp, rhs);
    if(C) {
      localTemp->linC(-dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(*localTemp, rhs);
    }
    rhs.linAdd(dt*dt*beta, residual);
  }

  times->predictorTime += getTime();
}

void
MDNLDynamic::subKucTransposeMultSubtractClaw(int isub, DistrVector& residual, double *userDefineDisp)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  if(subClaw) {
    if(subClaw->numUserDisp) {
      double *subUserDefineDisp = new double[subClaw->numUserDisp];
      for(int i=0; i<subClaw->numUserDisp; ++i) {
        subUserDefineDisp[i] = userDefineDisp[sd->getUserDispDataMap()[i]];
      }
      (*Kuc)[isub]->transposeMultSubtractClaw(subUserDefineDisp, residual.subData(isub), subClaw->numUserDisp, clawDofs[isub]);
      delete [] subUserDefineDisp;
    }
  }
}

void
MDNLDynamic::computeTimeInfo()
{
  // Time integration information

  // Get total time and time step size and store them 
  totalTime = domain->solInfo().tmax;
  dt0        = domain->solInfo().getTimeStep();

  // Compute maximum number of steps
  maxStep = (int) ( (totalTime+0.49*dt0)/dt0 );

  // Compute time remainder
  double remainder = totalTime - maxStep*dt0;
  if(std::abs(remainder)>0.01*dt0){
    domain->solInfo().tmax = totalTime = maxStep*dt0;
    filePrint(stderr, " Warning: Total time is being changed to : %e\n", totalTime);
  }
}

MDNLDynamic::MDNLDynamic(Domain *d)
 : MultiDomainBase(d->solInfo())
{
  domain = d;
#ifdef DISTRIBUTED
  decDomain = new GenDistrDomain<double>(domain);
#else
  decDomain = new GenDecDomain<double>(domain);
#endif
  numSystems = 0;
  secondRes = 0.0;
  claw = 0; 
  clawDofs = 0;
  userSupFunc = 0;
  aeroForce = 0;
  kelArray = 0;
  melArray = 0;
  celArray = 0;
  times = 0;
  mu = 0; muCopy = 0; lambda = 0;
  solver = 0;
  allOps =0;
  M = 0;
  C = 0;
  Kuc = 0;
  Muc = 0;
  Cuc = 0;
  Mcc = 0;
  Ccc = 0;
  allCorot = 0;
  localTemp = 0;
  reactions = 0;
  usrDefDisps = 0;
  usrDefVels = 0;
  prevFrc = 0;
  distFlExchanger = 0;
  updateCS = false;
}

MDNLDynamic::~MDNLDynamic()
{
  clean();
  if(mu) delete [] mu;
  if(lambda) delete [] lambda;
  if(times) delete times;

  if(aeroForce) delete aeroForce;
  if(prevFrc) delete prevFrc;
  if(distFlExchanger) delete distFlExchanger;

  if(decDomain) delete decDomain;
}

void
MDNLDynamic::clean()
{
  if(solver) { delete solver; solver = 0; }
  if(allOps) { delete allOps; allOps = 0; }
  if(M) { delete M; M = 0; }
  if(C) { delete C; C = 0; }
  if(Kuc) { delete Kuc; Kuc = 0; }
  if(Muc) { delete Muc; Muc = 0; }
  if(Cuc) { delete Cuc; Cuc = 0; }
  if(Mcc) { delete Mcc; Mcc = 0; }
  if(Ccc) { delete Ccc; Ccc = 0; }
  if(localTemp) delete localTemp;
  if(allCorot) {
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::deleteSubCorotators);
    delete [] allCorot;
    allCorot = 0;
  }
  if(kelArray || melArray || celArray) {
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::deleteSubElementArrays);
    if(kelArray) { delete [] kelArray; kelArray = 0; }
    if(melArray) { delete [] melArray; melArray = 0; }
    if(celArray) { delete [] celArray; celArray = 0; }
  }
  if(usrDefDisps) {
    for(int i=0; i<decDomain->getNumSub(); ++i) delete [] usrDefDisps[i];
    delete [] usrDefDisps;
    usrDefDisps = 0;
  }
  if(usrDefVels) {
    for(int i=0; i<decDomain->getNumSub(); ++i) delete [] usrDefVels[i];
    delete [] usrDefVels;
    usrDefVels = 0;
  }
  if(reactions) { delete reactions; reactions = 0; }
  if(clawDofs) {
    for(int i=0; i<decDomain->getNumSub(); ++i) if(clawDofs[i]) delete [] clawDofs[i];
    delete [] clawDofs;
  }
}

void
MDNLDynamic::initializeParameters(int step, DistrGeomState *geomState)
{
  if(step == 1 || domain->solInfo().reinit_lm) {
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::subInitializeMultipliers, *geomState);
  }
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subInitializeParameters);
  domain->initializeParameters();
}

void
MDNLDynamic::subInitializeMultipliers(int isub, DistrGeomState& geomState)
{
  decDomain->getSubDomain(isub)->initializeMultipliers(*(geomState[isub]), allCorot[isub]);
}

void
MDNLDynamic::subInitializeParameters(int isub)
{
  decDomain->getSubDomain(isub)->initializeParameters(false);
}

void
MDNLDynamic::updateParameters(DistrGeomState *geomState)
{
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateMultipliers, *geomState);
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateParameters);
  domain->updateParameters();
}

void
MDNLDynamic::subUpdateMultipliers(int isub, DistrGeomState& geomState)
{
  decDomain->getSubDomain(isub)->updateMultipliers(*(geomState[isub]), allCorot[isub]);
}

void
MDNLDynamic::subUpdateParameters(int isub)
{
  decDomain->getSubDomain(isub)->updateParameters(false);
}

bool
MDNLDynamic::checkConstraintViolation(double &err, DistrGeomState *gs)
{
  err = 0;
  for(int isub=0; isub<decDomain->getNumSub(); ++isub)
    err = std::max(err, decDomain->getSubDomain(isub)->getError(allCorot[isub], *(*gs)[isub]));
#ifdef DISTRIBUTED
  if(structCom) err = structCom->globalMax(err);
#endif
  return (err <= domain->solInfo().penalty_tol);
}

int
MDNLDynamic::checkConvergence(int iteration, double normRes, DistrVector &residual, DistrVector& dv, double time)
{
  times->timeCheck -= getTime();
  if(domain->GetnContactSurfacePairs()) dv.conservativeResize(solVecInfo());

  // Note when useTolInc is false, this function is called before normDv is calculated
  bool useTolInc = (domain->solInfo().getNLInfo().tolInc != std::numeric_limits<double>::infinity() ||
                    domain->solInfo().getNLInfo().absTolInc != std::numeric_limits<double>::infinity());

  double normDv  = dv.norm();
  double normEnergy = residual*dv;

  if(iteration == 0) {
    firstRes = normRes;
    if(useTolInc) {
      firstDv  = normDv;
      firstEng = normEnergy;
    }
    else {
      normDv = 0; firstDv = 1;
      normEnergy = 0; firstEng = 1;
    }
  }
  if(iteration == 1) {
    secondRes = normRes;
    if(!useTolInc) {
      firstDv  = normDv;
      firstEng = normEnergy;
    }
  }

  double relRes = normRes/firstRes;
  double relDv  = normDv /firstDv;

  int converged = 0;

  // Check for convergence
  if(iteration > 0 && ((normRes <= tolerance*firstRes && normDv <= domain->solInfo().getNLInfo().tolInc*firstDv) 
     || (normRes < domain->solInfo().getNLInfo().absTolRes && normDv < domain->solInfo().getNLInfo().absTolInc)))
    converged = 1;

  // Check for divergence
  else if(iteration > 1 && (normRes >= 1.0e10 * firstRes && normRes > secondRes)) 
    converged = -1;

  if(verboseFlag) {
    filePrint(stderr, " Iteration # %d\n", iteration);
    if(useTolInc || iteration >= 1) {
      double relEng = normEnergy/firstEng;
      filePrint(stderr, " r      = %e dv      = %e energy      = %e\n"
                        " rel. r = %e rel. dv = %e rel. energy = %e\n",
                        normRes, normDv, normEnergy,
                        relRes, relDv, relEng);
    }
    else {
      filePrint(stderr, " r      = %e\n"
                        " rel. r = %e\n",
                        normRes, relRes);
    }
  }

  totIter++;

  // Store residual norm and dv norm for output
  times->norms[numSystems].normDv      = normDv;
  times->norms[numSystems].relativeDv  = relDv;
  times->norms[numSystems].normRes     = normRes;
  times->norms[numSystems].relativeRes = relRes;
  times->numSystems = numSystems;

  numSystems += 1;
  times->timeCheck += getTime();
  return converged;
}

void
MDNLDynamic::updateContactSurfaces(DistrGeomState& geomState, DistrGeomState *refState)
{
  if(fetiSolver || (!domain->solInfo().readInDualROB.empty() && domain->solInfo().activatePodRom)) {
    domain->UpdateSurfaces(MortarHandler::CTC, &geomState, decDomain->getAllSubDomains());
    domain->PerformStaticContactSearch(MortarHandler::CTC);
    domain->deleteSomeLMPCs(mpc::ContactSurfaces);
    domain->ExpComputeMortarLMPC(MortarHandler::CTC);
    domain->CreateMortarToMPC();
    decDomain->reProcessMPCs();
    if(fetiSolver)
      fetiSolver->reconstructMPCs(decDomain->mpcToSub_dual.get(),
        decDomain->mpcToMpc.get(), decDomain->mpcToCpu.get());
  }
  else {
    clean();
    //std::cout << "line number " << __LINE__ << " of " << __FILE__ << std::endl;
    domain->ReInitializeStaticContactSearch(MortarHandler::CTC, decDomain->getNumSub(), decDomain->getAllSubDomains());
    domain->UpdateSurfaces(MortarHandler::CTC, &geomState, decDomain->getAllSubDomains());
    if(domain->solInfo().trivial_detection) {
    } else {
      domain->PerformStaticContactSearch(MortarHandler::CTC);
      domain->ExpComputeMortarLMPC(MortarHandler::CTC);
    }
    paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::renumberElementsGlobal);
    geoSource->UpdateContactSurfaceElements(&geomState, *mu);
    domain->deleteSomeLMPCs(mpc::ContactSurfaces);
    domain->getElementSet().removeAll();
    geoSource->getElems(domain->getElementSet());
    domain->setNumElements(domain->getElementSet().last());
    decDomain->clean();
    //MDNLDynamic::preProcess();
    this->preProcess();
    if(muCopy) for(int i=0; i<decDomain->getNumSub(); ++i) muCopy[i] = mu[i];
    geomState.resize(decDomain, mu);
    if(refState) {
      refState->resize(decDomain);
    }

  }
}

double
MDNLDynamic::getStiffAndForce(DistrGeomState& geomState, DistrVector& residual,
                              DistrVector& elementInternalForce, double t, DistrGeomState *refState, bool forceOnly) 
{
  times->buildStiffAndForce -= getTime();

  // update the geomState according to the USDD prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp > 0) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel  = new double[claw->numUserDisp];
      double *userDefineAcc  = new double[claw->numUserDisp];
      for(int i=0; i<claw->numUserDisp; ++i) {
        userDefineVel[i] = 0;
        userDefineAcc[i] = 0;
      }
      userSupFunc->usd_disp(t, userDefineDisp, userDefineVel, userDefineAcc);
      execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, &geomState, userDefineDisp,
                  userDefineVel, userDefineAcc);
      delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
    }
  }

  // update the contact surfaces
  if(t != domain->solInfo().initialTime) {
    if(domain->GetnContactSurfacePairs()) {
      if(!domain->solInfo().piecewise_contact || updateCS) {
        updateContactSurfaces(geomState, refState);
        updateCS = false;
      }
      elementInternalForce.resize(elemVecInfo());
      residual.conservativeResize(solVecInfo());
      localTemp->resize(solVecInfo());
    }
    // set the gap for the linear constraints
    if(fetiSolver) decDomain->setConstraintGap(&geomState, refState, fetiSolver, t);
  }

  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subGetStiffAndForce, geomState,
              residual, elementInternalForce, t, refState, forceOnly);

  times->buildStiffAndForce += getTime();

  return sqrt(solver->getFNormSq(residual));
}

void
MDNLDynamic::subGetStiffAndForce(int isub, DistrGeomState &geomState,
                                 DistrVector &res, DistrVector &elemIntForce, double t,
                                 DistrGeomState *refState, bool forceOnly)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector residual(res.subData(isub), res.subLen(isub));
  // eIF = element internal force
  StackVector eIF(elemIntForce.subData(isub), elemIntForce.subLen(isub));
  GeomState *subRefState = (refState) ? (*refState)[isub] : 0;
  if(forceOnly) {
    sd->getInternalForce(*geomState[isub], eIF, allCorot[isub], kelArray[isub], residual, 1.0, t, subRefState,
                         (Vector *) NULL, melArray[isub], (celArray) ? celArray[isub] : NULL);
  }
  else {
    sd->getStiffAndForce(*geomState[isub], eIF, allCorot[isub], kelArray[isub], residual, 1.0, t, subRefState,
                         (Vector *) NULL, melArray[isub], (celArray) ? celArray[isub] : NULL);
  }
}

void
MDNLDynamic::subUpdateGeomStateUSDD(int isub, DistrGeomState *geomState, double *userDefineDisp,
                                    double *userDefineVel, double *userDefineAcc)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  if(subClaw) {
    if(subClaw->numUserDisp) {
      double *subUserDefineDisp = new double[subClaw->numUserDisp];
      double *subUserDefineVel = new double[subClaw->numUserDisp];
      double *subUserDefineAcc = new double[subClaw->numUserDisp];
      for(int i=0; i<subClaw->numUserDisp; ++i) {
        int globalIndex = sd->getUserDispDataMap()[i];
        subUserDefineDisp[i] = userDefineDisp[globalIndex];
        subUserDefineVel[i] = userDefineVel[globalIndex];
        subUserDefineAcc[i] = userDefineAcc[globalIndex];
      }
      (*geomState)[isub]->updatePrescribedDisplacement(subUserDefineDisp, subClaw, sd->getNodes(),
                                                       subUserDefineVel, subUserDefineAcc);
      delete [] subUserDefineDisp;
      delete [] subUserDefineVel;
      delete [] subUserDefineAcc;
    }
  }
}

DistrGeomState*
MDNLDynamic::createGeomState()
{
  times->timeGeom -= getTime();
  DistrGeomState *geomState = new DistrGeomState(decDomain);
  times->timeGeom += getTime();
  return geomState;
}

DistrGeomState*
MDNLDynamic::copyGeomState(DistrGeomState* geomState)
{
 times->timeGeom -= getTime();
 DistrGeomState *geomStateCopy = new DistrGeomState(*geomState);
 times->timeGeom += getTime();
 return geomStateCopy;
}

void
MDNLDynamic::reBuild(DistrGeomState& geomState, int iteration, double localDelta, double t)
{
 int step = (t+localDelta)/(2*localDelta);
 times->rebuild -= getTime();

 if((iteration % domain->solInfo().getNLInfo().updateK == 0 && (step-1) % domain->solInfo().getNLInfo().stepUpdateK == 0)
    || (t == domain->solInfo().initialTime)) {
   times->norms[numSystems].rebuildTang = 1;

   double Kcoef, Ccoef, Mcoef;
   if(t == domain->solInfo().initialTime) {
     Kcoef = 0;
     Ccoef = 0;
     Mcoef = 1;
   }
   else {
     if(verboseFlag) filePrint(stderr, " ... Rebuilding Tangent Stiffness for Step %d Iteration %d ...\n", step, iteration);
     double beta, gamma, alphaf, alpham, dt = 2*localDelta;
     getNewmarkParameters(beta, gamma, alphaf, alpham);
     Kcoef = (domain->solInfo().order == 1) ? localDelta : dt*dt*beta;
     Ccoef = (domain->solInfo().order == 1) ? 0.0 : dt*gamma;
     if(domain->solInfo().quasistatic) Mcoef = 0;
     else Mcoef = (domain->solInfo().order == 1) ? 1 : (1-alpham)/(1-alphaf);
   }
   GenMDDynamMat<double> ops;
   ops.sysSolver = solver;
   ops.Kuc = Kuc;

   decDomain->rebuildOps(ops, Mcoef, Ccoef, Kcoef, kelArray, melArray, celArray);

   Kcoef_p = Kcoef;

 } else
   times->norms[numSystems].rebuildTang = 0;

 times->rebuild += getTime();
}

void
MDNLDynamic::preProcess()
{
  // Structure used to store timers
  if(!times) times = new StaticTimers;

  times->memoryPreProcess -= threadManager->memoryUsed();

  totIter = 0;
  
  // Set the nonlinear tolerance
  tolerance = domain->solInfo().getNLInfo().tolRes;

  // Constructs all renumbering, connectivities and dofsets
  times->preProcess -= getTime();
  decDomain->preProcess();
  times->preProcess += getTime();

  // Make each subdomain's dofs
  times->makeDOFs -= getTime();
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubDofs);
  times->makeDOFs += getTime();

  // Make each subdomain's corotators
  times->corotatorTime -= getTime();
  allCorot = new Corotator**[decDomain->getNumSub()]; 
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubCorotators);
  times->corotatorTime += getTime();

  // Allocate each subdomain's array of stiffness matrices
  kelArray = new FullSquareMatrix*[decDomain->getNumSub()];

  // Allocate each subdomain's array of mass matrices
  melArray = new FullSquareMatrix*[decDomain->getNumSub()];

  // Allocate each subdomain's array of damping matrices
  if(domain->solInfo().alphaDamp != 0 || domain->solInfo().betaDamp != 0 || domain->getElementSet().hasDamping()) 
    celArray = new FullSquareMatrix*[decDomain->getNumSub()];

  // Compute the geometric rigid body modes if requested
  MultiDomainRbm<double> *rigidBodyModes = 0;
  if(!mu && (domain->solInfo().rbmflg || domain->solInfo().filterFlags)) {
    rigidBodyModes = decDomain->constructRbm();
  }

  // Allocate vector to store reaction forces
  if(!reactions) reactions = new DistrVector(*decDomain->pbcVectorInfo());

  // Now make those arrays
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubElementArrays);

  times->memoryPreProcess += threadManager->memoryUsed();

  // Construct FETI Solver
  times->getFetiSolverTime -= getTime();
  allOps = new MDDynamMat;
  double Kcoef = 0.0;
  double Mcoef = 1.0;
  double Ccoef = 0.0;
  decDomain->buildOps(*allOps, Mcoef, Ccoef, Kcoef, (Rbm **)NULL, kelArray, true, melArray, celArray, factorWhenBuilding());
  solver = (ParallelSolver *) allOps->dynMat;
  fetiSolver = dynamic_cast<GenFetiDPSolver<double> *>(solver);
  M = allOps->M;
  C = allOps->C;
  Kuc = allOps->Kuc;
  Muc = allOps->Muc;
  Cuc = allOps->Cuc;
  Mcc = allOps->Mcc;
  Ccc = allOps->Ccc;
  if(allOps->K) delete allOps->K;
  times->getFetiSolverTime += getTime();

  times->memoryPreProcess -= threadManager->memoryUsed();

 if(sinfo.filterFlags) projector_prep(rigidBodyModes, M);
 if(rigidBodyModes) delete rigidBodyModes;

  // Look if there is a user supplied routine for control
  claw = geoSource->getControlLaw();

  // create list of usdd node dofs mapped to cdsa dof numbers
  if(claw && claw->numUserDisp)  {
    clawDofs = new int * [decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubClawDofs);
  }
  // Check to see if there is a user supplied function
  // for displacements, forces or control law
  userSupFunc = domain->getUserSuppliedFunction();

  localTemp = new DistrVector(decDomain->solVecInfo());

  if(mu == 0) { // first time only
    domain->InitializeStaticContactSearch(MortarHandler::CTC, decDomain->getNumSub(), decDomain->getAllSubDomains());
    mu = new std::map<std::pair<int,int>,double>[decDomain->getNumSub()];
    if(domain->solInfo().getNLInfo().failsafe) muCopy = new std::map<std::pair<int,int>,double>[decDomain->getNumSub()];
    lambda = new std::vector<double>[decDomain->getNumSub()];
    updateCS = true;
  }

  times->memoryPreProcess += threadManager->memoryUsed();
}

void
MDNLDynamic::makeSubDofs(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  sd->makeAllDOFs();
}

void
MDNLDynamic::makeSubCorotators(int isub)
{
  SubDomain *sd  = decDomain->getSubDomain(isub);
  int numele     = sd->numElements();
  allCorot[isub] = new Corotator*[numele];
  sd->createCorotators(allCorot[isub]);
}

void
MDNLDynamic::deleteSubCorotators(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  if(allCorot[isub]) {
    for (int iElem = 0; iElem < sd->numElements(); ++iElem) {
      if(allCorot[isub][iElem] && (allCorot[isub][iElem] != dynamic_cast<Corotator*>(sd->getElementSet()[iElem])))
        delete allCorot[isub][iElem];
    }
    delete [] allCorot[isub];
  }
}

void
MDNLDynamic::makeSubElementArrays(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  if(celArray) sd->createKelArray(kelArray[isub], melArray[isub], celArray[isub]); 
  else sd->createKelArray(kelArray[isub], melArray[isub]);
}

void
MDNLDynamic::deleteSubElementArrays(int isub)
{
  if(kelArray && kelArray[isub]) delete [] kelArray[isub];
  if(melArray && melArray[isub]) delete [] melArray[isub];
  if(celArray && celArray[isub]) delete [] celArray[isub];
}

void
MDNLDynamic::makeSubClawDofs(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  int nClaw = subClaw->numUserDisp;
  if(nClaw > 0) {
    clawDofs[isub] = new int[nClaw];
    for(int j = 0; j < nClaw; ++j) {
      int dd = sd->getDSA()->locate(subClaw->userDisp[j].nnum, (1 << subClaw->userDisp[j].dofnum));
      clawDofs[isub][j] = sd->getCDSA()->invRCN(dd);
    }
  }
  else clawDofs[isub] = 0;
}

void
MDNLDynamic::getExternalForce(DistrVector& f, DistrVector& constantForce,
                              int tIndex, double t, DistrGeomState* geomState,
                              DistrVector& elementInternalForce, DistrVector& aero_f, double localDelta)
{
  times->formRhs -= getTime();

  // update nodal temperature for thermoe
  if(domain->solInfo().thermoeFlag >= 0 && tIndex >= 0) {
    distFlExchanger->getStrucTemp(nodalTemps->data());
    if(verboseFlag) filePrint(stderr, " ... [E] Received temperatures      ...\n");
    geomState->setNodalTemperatures(*nodalTemps);
  }

  double beta, gamma, alphaf, alpham;
  getNewmarkParameters(beta, gamma, alphaf, alpham);
  double dt = 2*localDelta;
  double t0 = domain->solInfo().initialTime;
  double tm = (t == t0) ? t0 : t + dt*(alphaf-alpham);
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subGetExternalForce,
              f, constantForce, t, tm);

  // add the USDF and ACTUATOR forces
  if(claw && userSupFunc) {
    if(claw->numUserForce > 0) {
      double *userDefineForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(t, userDefineForce);
      decDomain->addUserForce(f, userDefineForce);
      delete [] userDefineForce;
    }

    if(claw->numActuator > 0) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel  = new double[claw->numSensor];
      double *ctracc  = new double[claw->numSensor];
      double *ctrfrc  = new double[claw->numActuator];

      for(int i=0; i<claw->numSensor; ++i) ctrvel[i] = ctracc[i] = 0.0; // not supported
#ifdef DISTRIBUTED
      for(int i=0; i<claw->numSensor; ++i) ctrdisp[i] = std::numeric_limits<double>::min();
#endif
      for(int i=0; i<decDomain->getNumSub(); ++i) subExtractControlDisp(i, *geomState, ctrdisp);
#ifdef DISTRIBUTED
      structCom->globalMax(claw->numSensor, ctrdisp);
#endif

      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t);
      decDomain->addCtrl(f, ctrfrc);

      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }

  // AERO
  SolverInfo& sinfo = domain->solInfo();
  if(sinfo.aeroFlag >= 0 && tIndex >= 0) {

    domain->getTimers().receiveFluidTime -= getTime();
    aeroForce->zero();
    int iscollocated;
    double tFluid = distFlExchanger->getFluidLoad(*aeroForce, tIndex, t,
                                                  alphaf, iscollocated);
    if(verboseFlag) filePrint(stderr, " ... [E] Received fluid load        ...\n");

    if(sinfo.aeroFlag == 20) {
      if(prevIndex >= 0)
        aero_f.linC(0.5,*aeroForce,0.5,*prevFrc);
      else
        aero_f = *aeroForce;
    }
    else {
      if (iscollocated == 0) {
        if(prevIndex >= 0) {
          *aeroForce *= (1/gamma);
          aeroForce->linAdd(((gamma-1.0)/gamma),*prevFrc);
        }
      }

      double alpha = 1.0-alphaf;
      if(prevIndex < 0) alpha = 1.0;

      aero_f.linC(alpha, *aeroForce, (1.0-alpha), *prevFrc);
    }
    f += aero_f;

    *prevFrc = *aeroForce;
    prevTime = tFluid;
    prevIndex = tIndex;
    domain->getTimers().receiveFluidTime += getTime();
  }

  // AEROH
  if(sinfo.aeroheatFlag >= 0 && tIndex >= 0) {

    aeroForce->zero();
    double tFluid = distFlExchanger->getFluidFlux(*aeroForce, tIndex, t);
    if(verboseFlag) filePrint(stderr, " ... [T] Received fluid fluxes      ...\n");

    /*  Compute fluid flux at n+1/2, since we use midpoint rule in thermal */

    int useProjector = domain->solInfo().filterFlags;

    if(tIndex == 0)
      f += *aeroForce;
    else {
      if(useProjector) f = *aeroForce;
      else
        f.linAdd(0.5, *aeroForce, 0.5, *prevFrc);
    }

    *prevFrc = *aeroForce;
  }

  // apply rbmfilter projection
  if(sinfo.filterFlags) trProject(f);

  times->formRhs += getTime();
}

void
MDNLDynamic::subGetExternalForce(int isub, DistrVector& f, DistrVector& constantForce, double tf, double tm)
{
  StackVector localf(f.subData(isub), f.subLen(isub));
  StackVector localg(constantForce.subData(isub), constantForce.subLen(isub));

  SubDomain *sd = decDomain->getSubDomain(isub);

  SparseMatrix *localCuc = (Cuc) ? (*Cuc)[isub] : 0;
  sd->computeExtForce4(localf, localg, tf, (*Kuc)[isub], userSupFunc, localCuc, tm, (*Muc)[isub]);
}

void
MDNLDynamic::subExtractControlDisp(int isub, DistrGeomState &geomState, double *ctrdsp)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  if(subClaw && subClaw->numSensor > 0) {
    CoordSet &nodes = sd->getNodes();
    for(int i=0; i<subClaw->numSensor; ++i) {
      switch (subClaw->sensor[i].dofnum) {
        case 0:
          ctrdsp[sd->getSensorDataMap()[i]] = (*geomState[isub])[subClaw->sensor[i].nnum].x
                        - nodes[subClaw->sensor[i].nnum]->x;
          break;
        case 1:
          ctrdsp[sd->getSensorDataMap()[i]] = (*geomState[isub])[subClaw->sensor[i].nnum].y
                        - nodes[subClaw->sensor[i].nnum]->y;
          break;
        case 2:
          ctrdsp[sd->getSensorDataMap()[i]] = (*geomState[isub])[subClaw->sensor[i].nnum].z
                        - nodes[subClaw->sensor[i].nnum]->z;
          break;
        default:
          fprintf(stderr, "ERROR: Sensor dof %d not available in MDNLDynam::subExtractControlDisp\n",subClaw->sensor[i].dofnum+1);
      }
    }
  }
}

ParallelSolver *
MDNLDynamic::getSolver()
{
  return solver;
}

MultiDomainPostProcessor *
MDNLDynamic::getPostProcessor()
{
  return new MultiDomainPostProcessor(decDomain, solver);
}

void
MDNLDynamic::printTimers(double timeLoop)
{
  int i;
  long (*memory)=(long *) dbg_alloca(sizeof(long)*decDomain->getNumSub());
  for (i = 0; i < decDomain->getNumSub(); ++i)
    memory[i] = 0;

  MultiDomainPostProcessor *mdpp = getPostProcessor();

  execParal(decDomain->getNumSub(), mdpp,
           &MultiDomainPostProcessor::getMemoryPrec, memory);
  long totMemPrec = 0;
  for(i = 0; i < decDomain->getNumSub(); ++i)
    totMemPrec += memory[i];

  for(i = 0; i < decDomain->getNumSub(); ++i)
    memory[i] = 0;

  execParal(decDomain->getNumSub(), mdpp,
           &MultiDomainPostProcessor::getMemoryK, memory);

  delete mdpp;

  long totMemK = 0;
  for(i=0; i<decDomain->getNumSub(); ++i)
    totMemK += memory[i];

  Timings &fetiTimers = solver->getTimers();
  fetiTimers.preconditioner.addOverAll(totMemPrec, 0.0);
  fetiTimers.kMem.addOverAll(totMemK, 0.0);

#ifdef DISTRIBUTED

  double mem1 = (double) totMemPrec;
  if(structCom) mem1 = structCom->globalSum(mem1);
  totMemPrec = (long) mem1;

  mem1 = (double) totMemK;
  if(structCom) mem1 = structCom->globalSum(mem1);
  totMemK = (long) mem1;

#endif

  times->memoryK = totMemK;
  times->memoryPrecond = totMemPrec;

  times->timeTimers -= getTime();

  times->numSubdomain = decDomain->getNumSub();

  times->printTimers(domain, solver->getTimers(),
                     solver->getSolutionTime());

  times->timeTimers += getTime();
}

void
MDNLDynamic::dynamOutput(DistrGeomState *geomState, DistrVector &vel_n, DistrVector &vel_p, 
                         double time, int index, DistrVector &ext_force, DistrVector &aeroF, DistrVector &acc_n,
                         DistrGeomState *refState)
{
  if(claw && claw->numUserDisp) {
    double *userDefineDisp = new double[claw->numUserDisp];
    double *userDefineVel  = new double[claw->numUserDisp];
    double *userDefineAcc  = new double[claw->numUserDisp];
    for(int i=0; i<claw->numUserDisp; ++i) {
      userDefineVel[i] = 0;
      userDefineAcc[i] = 0;
    }
    userSupFunc->usd_disp(time,userDefineDisp,userDefineVel,userDefineAcc);
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, geomState, userDefineDisp,
                userDefineVel, userDefineAcc);
    paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC, userDefineDisp,
               userDefineVel, userDefineAcc, true);
    delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
  }

  if(domain->solInfo().nRestart > 0) {
#ifdef RESTART_NO_EXECPARAL
    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      StackVector vel_ni(vel_n.subData(i), vel_n.subLen(i));
      StackVector acc_ni(acc_n.subData(i), acc_n.subLen(i));
      int extlen = std::log10((double) sd->subNum()+1) + 1;
      char *ext = new char[extlen+2];
      sprintf(ext,"_%d",sd->subNum()+1);
      sd->writeRestartFile(time, index+1, vel_ni, acc_ni, (*geomState)[i], ext);
      delete [] ext;
    }
#else
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::subWriteRestartFile, time, index, vel_n, acc_n, *geomState);
#endif
  }

  if(domain->reactionsReqd(time, index+1)) {
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::subGetReactionForce, *geomState, *refState, vel_n, acc_n, time);
  }

  DistrVector d_n(decDomain->solVecInfo()); d_n.zero();
  SysState<DistrVector> distState(d_n, vel_n, acc_n, vel_p); 
  decDomain->postProcessing(geomState, ext_force, allCorot, time, &distState, &aeroF, refState, reactions, allOps);
}

void
MDNLDynamic::subWriteRestartFile(int i, double &time, int &index, DistrVector &vel_n, DistrVector &acc_n, DistrGeomState &geomState)
{
  SubDomain *sd = decDomain->getSubDomain(i);
  StackVector vel_ni(vel_n.subData(i), vel_n.subLen(i));
  StackVector acc_ni(acc_n.subData(i), acc_n.subLen(i));
  int extlen = (int)std::log10((double) sd->subNum()+1) + 1;
  char *ext = new char[extlen+2];
  sprintf(ext,"_%d",sd->subNum()+1);
  sd->writeRestartFile(time, index+1, vel_ni, acc_ni, geomState[i], ext);
  delete [] ext;
}

void
MDNLDynamic::subGetReactionForce(int i, DistrGeomState &geomState, DistrGeomState &refState, DistrVector &vel_n, 
                                 DistrVector &acc_n, double &time)
{
  SubDomain *sd = decDomain->getSubDomain(i);
  StackVector ri(reactions->subData(i), reactions->subLen(i));
  StackVector velocity(vel_n.subData(i), vel_n.subLen(i));
  StackVector acceleration(acc_n.subData(i), acc_n.subLen(i));
  SparseMatrix *Cuci = (Cuc) ? (*Cuc)[i] : 0;
  SparseMatrix *Ccci = (Ccc) ? (*Ccc)[i] : 0;
  sd->computeReactionForce(ri, geomState[i], allCorot[i], kelArray[i], time, refState[i], velocity,
                           acceleration, sd->getVcx(), sd->getAcx(), Cuci, Ccci, (*Muc)[i], (*Mcc)[i]);
 
  // TODO: the lagrange multipliers should probably be extrapolated to t^{n+1}
  //       check if the equality constraint forces are incremental
  sd->addCConstraintForces(geomState.mu[i], geomState.lambda[i], ri, 1/Kcoef_p);
}

void
MDNLDynamic::processLastOutput()  
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  domain->solInfo().lastIt = true;
  for(int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

void
MDNLDynamic::getInitialTime(int &initTimeIndex, double &initTime)
{
  initTimeIndex = domain->solInfo().initialTimeIndex;
  initTime      = domain->solInfo().initialTime;
}

void
MDNLDynamic::getNewmarkParameters(double &beta, double &gamma,
                                  double &alphaf, double &alpham)
{
  beta  = domain->solInfo().newmarkBeta;
  gamma = domain->solInfo().newmarkGamma;
  alphaf = domain->solInfo().newmarkAlphaF;
  alpham = domain->solInfo().newmarkAlphaM;
}

int
MDNLDynamic::getMaxit()
{
  return domain->solInfo().getNLInfo().maxiter;
}

double
MDNLDynamic::getDeltaLambda()
{
  return domain->solInfo().getNLInfo().dlambda;
}

int
MDNLDynamic::getNumStages()
{
  return int(0.2 + domain->solInfo().getNLInfo().maxLambda / domain->solInfo().getNLInfo().dlambda);
}

DistrInfo&
MDNLDynamic::solVecInfo()
{
  return decDomain->solVecInfo();
}

DistrInfo&
MDNLDynamic::elemVecInfo()
{ 
  return *decDomain->elementVectorInfo();
} 
  
DistrInfo&
MDNLDynamic::sysVecInfo()
{ 
  return decDomain->sysVecInfo();
}

double
MDNLDynamic::getTolerance()
{
 return std::max(tolerance*firstRes, domain->solInfo().getNLInfo().absTolRes);
}

double
MDNLDynamic::getResidualNorm(DistrVector &r, DistrGeomState &geomState, double)
{
 //returns: sqrt( (r+c^T*lambda)**2 + pos_part(gap)**2 )
 DistrVector w(r);
 execParal(decDomain->getNumSub(), this, &MDNLDynamic::addConstraintForces, w, geomState); // w = r + C^T*lambda
                  // note C = grad(gap) has already been updated in getStiffAndForce.
 return sqrt(solver->getFNormSq(w));
}

bool
MDNLDynamic::factorWhenBuilding() const {
  return false;
}

void
MDNLDynamic::addConstraintForces(int isub, DistrVector& vec, DistrGeomState& geomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector localvec(vec.subData(isub), vec.subLen(isub));
  sd->addConstraintForces(geomState.mu[isub], geomState.lambda[isub], localvec);  // C^T*lambda added to vec
}

void
MDNLDynamic::getConstraintMultipliers(DistrGeomState& geomState)
{
  // this function extracts the constraint multipliers from the subdomains,
  // where they are stored at the send of the feti solver's solve function
  // and stores them in geomState 
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subGetConstraintMultipliers, geomState);
}

void
MDNLDynamic::subGetConstraintMultipliers(int isub, DistrGeomState& geomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  geomState.mu[isub].clear();
  geomState.lambda[isub].clear();
  sd->getConstraintMultipliers(geomState.mu[isub], geomState.lambda[isub]);
}

void 
MDNLDynamic::dynamCommToFluid(DistrGeomState* geomState, DistrGeomState* bkGeomState,
                              DistrVector& velocity, DistrVector& bkVelocity,
                              DistrVector& vp, DistrVector& bkVp, int step, int parity,
                              int aeroAlg, double time) 
{  
  if(domain->solInfo().aeroFlag >= 0 && !domain->solInfo().lastIt) {

    domain->getTimers().sendFluidTime -= getTime();
    // update the geomState according to the USDD prescribed displacements
    if(claw && userSupFunc) {
      if(claw->numUserDisp > 0) {
        // Note: the approprate value of "time" passed into this function should be t^{n+½} for A6 and t^{n+1}
        // otherwise, where t^n denotes the time at the end of the current structure timestep. Note that the
        // predictor in FlExchanger::sendDisplacements is not applied to prescribed displacements; we directly
        // compute here the desired values of the prescribed displacements/velocities rather than predicting them.
        double *userDefineDisp = new double[claw->numUserDisp];
        double *userDefineVel  = new double[claw->numUserDisp];
        double *userDefineAcc  = new double[claw->numUserDisp];
        for(int i=0; i<claw->numUserDisp; ++i) {
          userDefineVel[i] = 0;
          userDefineAcc[i] = 0;
        }
        userSupFunc->usd_disp(time, userDefineDisp, userDefineVel, userDefineAcc);
        execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, geomState, userDefineDisp,
                    userDefineVel, userDefineAcc);
       delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
      } 
    }

    DistrVector d_n(decDomain->solVecInfo()); d_n.zero();
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::subDynamCommToFluid, d_n, geomState, bkGeomState, parity, aeroAlg);
    if(!parity && aeroAlg == 5) {
      velocity.linC(0.5, velocity, 0.5, bkVelocity);
      vp.linC(0.5, vp, 0.5, bkVp);
    }
    DistrVector acceleration(decDomain->solVecInfo()); acceleration.zero();

    SysState<DistrVector> state(d_n, velocity, acceleration, vp);

    distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
    if(verboseFlag) filePrint(stderr, " ... [E] Sent displacements         ...\n");
    domain->getTimers().sendFluidTime += getTime();
  }

  if(domain->solInfo().aeroheatFlag >= 0) {
    DistrVector d_n(decDomain->solVecInfo()); d_n.zero();
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::subDynamCommToFluidAeroheat, d_n, geomState);

    SysState<DistrVector> tempState(d_n, velocity, vp);

    distFlExchanger->sendTemperature(tempState);
    if(verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
  }

  if(domain->solInfo().thermohFlag >= 0) {
    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      for(int j = 0; j < sd->numNodes(); ++j) {
        nodalTemps->subData(i)[j] = (*(*geomState)[i])[j].x;
      }
    }

    distFlExchanger->sendStrucTemp(*nodalTemps);
    if(verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
  }
}

void
MDNLDynamic::subDynamCommToFluid(int isub, DistrVector& v, DistrGeomState* distrGeomState, 
                                 DistrGeomState* bkDistrGeomState, int parity, int aeroAlg)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector d_n(v.subData(isub), v.subLen(isub));
  GeomState* geomState = (*distrGeomState)[isub];
  double* bcx = usrDefDisps[isub];
  double* vcx = usrDefVels[isub];

  // Make d_n_aero from geomState
  ConstrainedDSA *c_dsa = sd->getCDSA();
  DofSetArray *dsa = sd->getDSA();
  CoordSet &nodes = sd->getNodes();
  int numNodes = sd->numNodes();

  for(int i = 0; i < numNodes; ++i) {

    int xloc  = c_dsa->locate(i, DofSet::Xdisp );
    int xloc1 =   dsa->locate(i, DofSet::Xdisp );

    if(xloc >= 0)
      d_n[xloc]  = ( (*geomState)[i].x - nodes[i]->x);
    else if (xloc1 >= 0) {
      bcx[xloc1] = ( (*geomState)[i].x - nodes[i]->x);
      vcx[xloc1] = (*geomState)[i].v[0];
    }

    int yloc  = c_dsa->locate(i, DofSet::Ydisp );
    int yloc1 =   dsa->locate(i, DofSet::Ydisp );

    if(yloc >= 0)
      d_n[yloc]  = ( (*geomState)[i].y - nodes[i]->y);
    else if (yloc1 >= 0) {
      bcx[yloc1] = ( (*geomState)[i].y - nodes[i]->y);
      vcx[yloc1] = (*geomState)[i].v[1];
    }

    int zloc  = c_dsa->locate(i, DofSet::Zdisp);
    int zloc1 =   dsa->locate(i, DofSet::Zdisp);

    if(zloc >= 0)
      d_n[zloc]  = ( (*geomState)[i].z - nodes[i]->z);
    else if (zloc1 >= 0) {
      bcx[zloc1] = ( (*geomState)[i].z - nodes[i]->z);
      vcx[zloc1] = (*geomState)[i].v[2];
    }
  }

  if(!parity && aeroAlg == 5) {
    Vector d_n2(v.subLen(isub), 0.0);
    GeomState* bkGeomState = (*bkDistrGeomState)[isub];
    for(int i = 0; i < numNodes; ++i) {

      int xloc  = c_dsa->locate(i, DofSet::Xdisp );
      if(xloc >= 0)
        d_n2[xloc]  = ( (*bkGeomState)[i].x - nodes[i]->x);

      int yloc  = c_dsa->locate(i, DofSet::Ydisp );
      if(yloc >= 0)
        d_n2[yloc]  = ( (*bkGeomState)[i].y - nodes[i]->y);

      int zloc  = c_dsa->locate(i, DofSet::Zdisp);
      if(zloc >= 0)
        d_n2[zloc]  = ( (*bkGeomState)[i].z - nodes[i]->z);
    }
    d_n.linC(0.5, d_n, 0.5, d_n2);
  }
}

void
MDNLDynamic::subDynamCommToFluidAeroheat(int isub, DistrVector& v, DistrGeomState* distrGeomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector d_n(v.subData(isub), v.subLen(isub));
  GeomState* geomState = (*distrGeomState)[isub];
  double* bcx = usrDefDisps[isub];

  // Make d_n from geomState
  ConstrainedDSA *c_dsa = sd->getCDSA();
  DofSetArray *dsa = sd->getDSA();
  CoordSet &nodes = sd->getNodes();
  int numNodes = sd->numNodes();

  for(int i = 0; i < numNodes; ++i) {

    int xloc  = c_dsa->locate(i, DofSet::Temp );
    int xloc1 =   dsa->locate(i, DofSet::Temp );

    if(xloc >= 0)
      d_n[xloc]  = (*geomState)[i].x;
    else if (xloc1 >= 0)
      bcx[xloc1] = (*geomState)[i].x;
  }
}

void 
MDNLDynamic::readRestartFile(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                             DistrVector &v_p, DistrGeomState &geomState)
{
  if(geoSource->getCheckFileInfo()->lastRestartFile) {
    filePrint(stderr, " ... Restarting From a Previous Run ...\n");
#ifdef RESTART_NO_EXECPARAL
    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      StackVector d_ni(d_n.subData(i), d_n.subLen(i));
      StackVector v_ni(v_n.subData(i), v_n.subLen(i));
      StackVector a_ni(a_n.subData(i), a_n.subLen(i));
      StackVector v_pi(v_p.subData(i), v_p.subLen(i));
      int extlen = std::log10((double) sd->subNum()+1) + 1;
      char *ext = new char[extlen+2];
      sprintf(ext,"_%d",sd->subNum()+1);
      sd->readRestartFile(d_ni, v_ni, a_ni, v_pi, sd->getBcx(), sd->getVcx(), *(geomState[i]), ext);
      delete [] ext;
    }
#else
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::subReadRestartFile, d_n, v_n, a_n, v_p, geomState);
#endif
    domain->solInfo().initialTimeIndex = decDomain->getSubDomain(0)->solInfo().initialTimeIndex;
    domain->solInfo().initialTime = decDomain->getSubDomain(0)->solInfo().initialTime;

    int aeroFlag = domain->solInfo().aeroFlag;
    if(aeroFlag >= 0 && geoSource->getCheckFileInfo()->hotRestart()) {
      double time = domain->solInfo().initialTime;
      double dt = domain->solInfo().getTimeStep();
      double t_aero = (aeroFlag == 6) ? time+dt/2 : time+dt;
      dynamCommToFluid(&geomState, &geomState, v_n, v_n, v_p, v_p, domain->solInfo().initialTimeIndex, 1,
                       domain->solInfo().aeroFlag, t_aero);
    }

    // update geomState and bcx/vcx/acx for time dependent prescribed displacements and their time derivatives
    // at the initial time
    if(claw && userSupFunc) {
      if(claw->numUserDisp > 0) {
        double *userDefineDisp = new double[claw->numUserDisp];
        double *userDefineVel  = new double[claw->numUserDisp];
        double *userDefineAcc  = new double[claw->numUserDisp];
        for(int i=0; i<claw->numUserDisp; ++i) {
          userDefineVel[i] = 0;
          userDefineAcc[i] = 0;
        }
        userSupFunc->usd_disp(domain->solInfo().initialTime, userDefineDisp, userDefineVel, userDefineAcc);
        execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, &geomState, userDefineDisp,
                    userDefineVel, userDefineAcc);
       delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
      }
    }

    updateStates(&geomState, geomState, domain->solInfo().initialTime);
  }
}

void
MDNLDynamic::subReadRestartFile(int i, DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                                DistrVector &v_p, DistrGeomState &geomState)
{
  SubDomain *sd = decDomain->getSubDomain(i);
  StackVector d_ni(d_n.subData(i), d_n.subLen(i));
  StackVector v_ni(v_n.subData(i), v_n.subLen(i));
  StackVector a_ni(a_n.subData(i), a_n.subLen(i));
  StackVector v_pi(v_p.subData(i), v_p.subLen(i));
  int extlen = (int)std::log10((double) sd->subNum()+1) + 1;
  char *ext = new char[extlen+2];
  sprintf(ext,"_%d",sd->subNum()+1);
  sd->readRestartFile(d_ni, v_ni, a_ni, v_pi, sd->getBcx(), sd->getVcx(), *(geomState[i]), ext);
  delete [] ext;
}

int
MDNLDynamic::aeroPreProcess(DistrVector &disp, DistrVector &vel,
                            DistrVector &accel, DistrVector &lastVel, DistrGeomState *distrGeomState)
{
  // get solver info
  SolverInfo& sinfo = domain->solInfo();

  // Initialize previous force data
  // Reexamine for the case of restart
  prevFrc = new DistrVector(solVecInfo());
  prevFrc->zero();
  prevIndex = -1;
  prevTime = 0;

  // Initialize the aeroforce vector
  aeroForce = new DistrVector(solVecInfo());

  if(sinfo.aeroFlag < 0)
    return 0;

  auto cpuToSub = geoSource->getCpuToSub();

  // get cpu id
#ifdef USE_MPI
  int myId = structCom->myID();
#else
  int myId = 0;
#endif

  int numLocSub = cpuToSub->num(myId);

  SubDomain **subdomain = decDomain->getAllSubDomains();

  // allocate for pointer arrays
  CoordSet **cs = new CoordSet *[numLocSub];
  Elemset **elemSet = new Elemset *[numLocSub];
  DofSetArray **cdsa = new DofSetArray *[numLocSub];
  DofSetArray **dsa = new DofSetArray *[numLocSub];
  usrDefDisps = new double *[numLocSub];
  usrDefVels = new double *[numLocSub];

  int iSub;
  for(iSub = 0; iSub < numLocSub; iSub++) {

    // assemble coordsets in this mpi
    cs[iSub] = &subdomain[iSub]->getNodes();

    // assemble element sets in this mpi
    elemSet[iSub] = &subdomain[iSub]->getElementSet();

    // assemble constrained and unconstrained dofset arrays in this mpi
    cdsa[iSub] = subdomain[iSub]->getCDSA();
    dsa[iSub] = subdomain[iSub]->getDSA();

    // allocate and initialize for the user defined disps and vels
    int numDofs = dsa[iSub]->size();
    usrDefDisps[iSub] = new double[numDofs];
    usrDefVels[iSub] = new double[numDofs];

    for(int iDof = 0; iDof < numDofs; iDof++) {
      usrDefDisps[iSub][iDof] = 0.0;
      usrDefVels[iSub][iDof] = 0.0;
    }
    if(distrGeomState) {
      GeomState* geomState = (*distrGeomState)[iSub];
      CoordSet &nodes = subdomain[iSub]->getNodes();
      int numNodes = subdomain[iSub]->numNodes();

      for(int i = 0; i < numNodes; ++i) {

        int xloc  = cdsa[iSub]->locate(i, DofSet::Xdisp );
        int xloc1 =  dsa[iSub]->locate(i, DofSet::Xdisp );

        if (xloc < 0 && xloc1 >= 0) {
          usrDefDisps[iSub][xloc1] = ( (*geomState)[i].x - nodes[i]->x);
          usrDefVels[iSub][xloc1] = (*geomState)[i].v[0];
        }

        int yloc  = cdsa[iSub]->locate(i, DofSet::Ydisp );
        int yloc1 =  dsa[iSub]->locate(i, DofSet::Ydisp );

        if (yloc < 0 && yloc1 >= 0) {
          usrDefDisps[iSub][yloc1] = ( (*geomState)[i].y - nodes[i]->y);
          usrDefVels[iSub][yloc1] = (*geomState)[i].v[1];
        }

        int zloc  = cdsa[iSub]->locate(i, DofSet::Zdisp);
        int zloc1 =  dsa[iSub]->locate(i, DofSet::Zdisp);

        if (zloc < 0 && zloc1 >= 0) {
          usrDefDisps[iSub][zloc1] = ( (*geomState)[i].z - nodes[i]->z);
          usrDefVels[iSub][zloc1] = (*geomState)[i].v[2];
        }
      }
    }
  }

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int flag = 0;

  // Check if aero forces are requested for output
  int iInfo;
  for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
    if(oinfo[iInfo].type == OutputInfo::AeroForce) {
      flag = 1;
      break;
    }
  }

  // create distributed fluid exchanger
  OutputInfo *oinfo_aero = (flag) ? oinfo+iInfo : NULL;
  std::set<int> &aeroEmbeddedSurfaceId = domain->GetAeroEmbedSurfaceId();
  if(aeroEmbeddedSurfaceId.size() != 0) {
    int iSurf = -1;
    for(int i = 0; i < domain->getNumSurfs(); i++)
      if(aeroEmbeddedSurfaceId.find((*domain->viewSurfEntities())[i]->ID()) != aeroEmbeddedSurfaceId.end()) {
        iSurf = i;
        break; // only allows one surface.
      }
    if(iSurf<0) {
      fprintf(stderr,"ERROR: Embedded wet surface not found! Aborting...\n");
      exit(-1);
    }
    distFlExchanger = new DistFlExchanger(cs, elemSet, (*domain->viewSurfEntities())[iSurf],
                                          &domain->getNodes(), domain->getNodeToElem(),
                                          decDomain->getElemToSub().get(), subdomain,
                                          cdsa, dsa, oinfo_aero, false);
  }
  else {
    distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa, oinfo_aero);
  }

  // negotiate with the fluid code
  distFlExchanger->negotiate();

  int restartinc = (sinfo.nRestart >= 0) ? (sinfo.nRestart) : 0;

  DistrVector dispAero(disp);

  if(sinfo.gepsFlg == 1) {
    // If we are in the first time step, and we initialized with
    // IDISP6, do not send IDISP6
    if(domain->numInitDisp() == 0 && sinfo.zeroInitialDisp != 1) {
      filePrint(stderr," ... DO NOT SEND IDISP6             ...\n");
    } else {
      filePrint(stderr," ... SENDING IDISP6                 ...\n");
      for(iSub = 0; iSub < numLocSub; iSub++) {
        BCond* iDis6 = subdomain[iSub]->getInitDisp6();
        for(int i = 0; i < subdomain[iSub]->numInitDisp6(); ++i) {
          int dof = cdsa[iSub]->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
          if(dof >= 0)
            dispAero[dof] += iDis6[i].val;
        }
      }
    }
  }

  SysState<DistrVector> state(dispAero, vel, accel, lastVel);

  if(sinfo.aeroFlag == 8) {
    distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.mppFactor,
                               restartinc, sinfo.isCollocated, sinfo.alphas);
    distFlExchanger->sendModeFreq(modeData.frequencies, modeData.numModes);
    if(verboseFlag) filePrint(stderr, " ... [E] Sent parameters and mode frequencies ...\n");
    distFlExchanger->sendModeShapes(modeData.numModes, modeData.numNodes,
                                    modeData.modes, state, sinfo.mppFactor);
    if(verboseFlag) filePrint(stderr, " ... [E] Sent mode shapes           ...\n");
  }
  else {
    distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.tmax, restartinc,
                               sinfo.isCollocated, sinfo.alphas);
    if(verboseFlag) filePrint(stderr, " ... [E] Sent parameters            ...\n");

    // initialize the Parity
    if(sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
      distFlExchanger->initRcvParity(1);
      distFlExchanger->initSndParity(1);
    } else {
      distFlExchanger->initRcvParity(-1);
      distFlExchanger->initSndParity(-1);
    }

    // send initial displacements
    if(!(geoSource->getCheckFileInfo()->lastRestartFile)) { // XXX check
      distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
      if(verboseFlag) filePrint(stderr, " ... [E] Sent initial displacements ...\n");
    }

    if(sinfo.aeroFlag == 1) { // Ping pong only
      filePrint(stderr, "Ping Pong Only requested. Structure code exiting\n");
    }
  }

  return sinfo.aeroFlag;
}

void
MDNLDynamic::thermoePreProcess()
{
  if(domain->solInfo().thermoeFlag >=0) {

    auto cpuToSub = geoSource->getCpuToSub();
    int myId = structCom->myID();
    int numLocSub = cpuToSub->num(myId);
    SubDomain **subdomain = decDomain->getAllSubDomains();

    // if sinfo.aeroFlag >= 0, flExchanger has already been initialize before,
    // thus, only when sinfo.aeroFlag < 0 is necessary.
    if(domain->solInfo().aeroFlag < 0) {

      // allocate for pointer arrays
      CoordSet **cs = new CoordSet *[numLocSub];
      Elemset **elemSet = new Elemset *[numLocSub];
      DofSetArray **cdsa = new DofSetArray *[numLocSub];
      DofSetArray **dsa = new DofSetArray *[numLocSub];
      usrDefDisps = new double *[numLocSub];
      usrDefVels = new double *[numLocSub];

      int iSub;
      for(iSub = 0; iSub < numLocSub; iSub++)  {

        // assemble coordsets in this mpi
        cs[iSub] = &subdomain[iSub]->getNodes();

        // assemble element sets in this mpi
        elemSet[iSub] = &subdomain[iSub]->getElementSet();

        // assemble constrained and unconstrained dofset arrays in this mpi
        cdsa[iSub] = subdomain[iSub]->getCDSA();
        dsa[iSub] = subdomain[iSub]->getDSA();

        // allocate and initialize for the user defined disps and vels
        int numDofs = dsa[iSub]->size();
        usrDefDisps[iSub] = new double[numDofs];
        usrDefVels[iSub] = new double[numDofs];

        for(int iDof = 0; iDof < numDofs; iDof++)  {
          usrDefDisps[iSub][iDof] = 0.0;
          usrDefVels[iSub][iDof] = 0.0;
        }
      }

      distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);
    }

    nodalTemps = new DistrVector(decDomain->ndVecInfo());
    for(int iSub = 0; iSub < numLocSub; iSub++) subdomain[iSub]->temprcvd = nodalTemps->subData(iSub);
    int buffLen = nodalTemps->size();

    distFlExchanger->thermoread(buffLen);

    distFlExchanger->getStrucTemp(nodalTemps->data()) ;
    if(verboseFlag) filePrint(stderr, " ... [E] Received initial temperatures ...\n");
  }
}

void
MDNLDynamic::thermohPreProcess(DistrVector& d)
{
  if(domain->solInfo().thermohFlag >=0) {

    auto cpuToSub = geoSource->getCpuToSub();
    int myId = structCom->myID();
    int numLocSub = cpuToSub->num(myId);
    SubDomain **subdomain = decDomain->getAllSubDomains();

    // if sinfo.aeroheatFlag >= 0, flExchanger has already been initialize before,
    // thus, only when sinfo.aeroheatFlag < 0 is necessary.
    if(domain->solInfo().aeroheatFlag < 0) {

      // allocate for pointer arrays
      CoordSet **cs = new CoordSet *[numLocSub];
      Elemset **elemSet = new Elemset *[numLocSub];
      DofSetArray **cdsa = new DofSetArray *[numLocSub];
      DofSetArray **dsa = new DofSetArray *[numLocSub];
      usrDefDisps = new double *[numLocSub];
      usrDefVels = new double *[numLocSub];

      int iSub;
      for(iSub = 0; iSub < numLocSub; iSub++)  {

        // assemble coordsets in this mpi
        cs[iSub] = &subdomain[iSub]->getNodes();

        // assemble element sets in this mpi
        elemSet[iSub] = &subdomain[iSub]->getElementSet();

        // assemble constrained and unconstrained dofset arrays in this mpi
        cdsa[iSub] = subdomain[iSub]->getCDSA();
        dsa[iSub] = subdomain[iSub]->getDSA();

        // allocate and initialize for the user defined disps and vels
        int numDofs = dsa[iSub]->size();
        usrDefDisps[iSub] = new double[numDofs];
        usrDefVels[iSub] = new double[numDofs];

        for(int iDof = 0; iDof < numDofs; iDof++)  {
          usrDefDisps[iSub][iDof] = 0.0;
          usrDefVels[iSub][iDof] = 0.0;
        }
      }

      distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);
    }

    nodalTemps = new DistrVector(decDomain->ndVecInfo());
    int buffLen = nodalTemps->size();

    distFlExchanger->thermoread(buffLen);

    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      for(int j = 0; j < sd->numNodes(); ++j) {
        int tloc  = sd->getCDSA()->locate(j, DofSet::Temp);
        int tloc1 = sd->getDSA()->locate(j, DofSet::Temp);
        double temp  = (tloc >= 0) ? d.subData(i)[tloc] : sd->getBcx()[tloc1];
        if(tloc1 < 0) temp = 0.0;
        nodalTemps->subData(i)[j] = temp;
      }
    }

    distFlExchanger->sendStrucTemp(*nodalTemps);
    if(verboseFlag) filePrint(stderr, " ... [T] Sent initial temperatures  ...\n");
  }

}

void
MDNLDynamic::aeroheatPreProcess(DistrVector &disp, DistrVector &vel, DistrVector &lastVel)
{
  // get solver info
  SolverInfo& sinfo = domain->solInfo();

  // Initialize previous force data
  // Reexamine for the case of restart
  prevFrc = new DistrVector(solVecInfo());
  prevFrc->zero();
  prevIndex = -1;
  prevTime = 0;

  // Initialize the aeroforce vector
  aeroForce = new DistrVector(solVecInfo());

  if(sinfo.aeroheatFlag < 0)
    return;

  auto cpuToSub = geoSource->getCpuToSub();

  // get cpu id
#ifdef USE_MPI
  int myId = structCom->myID();
#else
  int myId = 0;
#endif

  int numLocSub = cpuToSub->num(myId);

  SubDomain **subdomain = decDomain->getAllSubDomains();

  // allocate for pointer arrays
  CoordSet **cs = new CoordSet *[numLocSub];
  Elemset **elemSet = new Elemset *[numLocSub];
  DofSetArray **cdsa = new DofSetArray *[numLocSub];
  DofSetArray **dsa = new DofSetArray *[numLocSub];
  usrDefDisps = new double *[numLocSub];
  usrDefVels = new double *[numLocSub];

  int iSub;
  for(iSub = 0; iSub < numLocSub; iSub++)  {

    // assemble coordsets in this mpi
    cs[iSub] = &subdomain[iSub]->getNodes();

    // assemble element sets in this mpi
    elemSet[iSub] = &subdomain[iSub]->getElementSet();

    // assemble constrained and unconstrained dofset arrays in this mpi
    cdsa[iSub] = subdomain[iSub]->getCDSA();
    dsa[iSub] = subdomain[iSub]->getDSA();

    // allocate and initialize for the user defined disps
    int numDofs = dsa[iSub]->size();
    usrDefDisps[iSub] = new double[numDofs];

    for(int iDof = 0; iDof < numDofs; iDof++)  {
      usrDefDisps[iSub][iDof] = 0.0;
    }
  }

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int flag = 0;

  // Check if aero fluxes are requested for output
  int iInfo;
  for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
    if(oinfo[iInfo].type == OutputInfo::AeroForce) {
      flag = 1;
      break;
    }
  }

  // create distributed fluid exchanger
  if(flag)
    distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa, oinfo+iInfo);
  else
    distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);

  // negotiate with the fluid code
  distFlExchanger->negotiate();

  int restartinc = (sinfo.nRestart >= 0) ? (sinfo.nRestart) : 0;

  SysState<DistrVector> state(disp, vel, lastVel);

  distFlExchanger->sendTempParam(sinfo.aeroheatFlag, sinfo.getTimeStep(), sinfo.tmax, restartinc,
                                 sinfo.alphat);
  if(verboseFlag) filePrint(stderr, " ... [T] Sent parameters            ...\n");

  // send initial displacements
  distFlExchanger->sendTemperature(state);
  if(verboseFlag) filePrint(stderr, " ... [T] Sent initial temperatures  ...\n");
}

int
MDNLDynamic::getAeroAlg()
{
  return domain->solInfo().aeroFlag;
}

int
MDNLDynamic::getThermoeFlag()
{
  return domain->solInfo().thermoeFlag;
}

int
MDNLDynamic::getThermohFlag()
{
  return domain->solInfo().thermohFlag;
}

int
MDNLDynamic::getAeroheatFlag()
{
  return domain->solInfo().aeroheatFlag;
}

void
MDNLDynamic::updateStates(DistrGeomState *refState, DistrGeomState& geomState, double time)
{
  if(domain->solInfo().piecewise_contact) updateCS = true;
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateStates, refState, &geomState, time);
}

void
MDNLDynamic::subUpdateStates(int isub, DistrGeomState *refState, DistrGeomState *geomState, double time)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  GeomState *subRefState = (refState) ? (*refState)[isub] : 0;
  sd->updateStates(subRefState, *(*geomState)[isub], allCorot[isub], time);
}

LinesearchInfo&
MDNLDynamic::linesearch()
{
  return domain->solInfo().getNLInfo().linesearch;
}

bool
MDNLDynamic::getResizeFlag()
{
  return (domain->GetnContactSurfacePairs() > 0);
}

void
MDNLDynamic::resize(DistrGeomState *refState, DistrGeomState *geomState, DistrGeomState *stepState, DistrVector *stateIncr,
                    DistrVector &v, DistrVector &a, DistrVector &vp, DistrVector &force)
{
  if(domain->GetnContactSurfacePairs() > 0 && domain->solInfo().piecewise_contact) {
    if(refState) {
      refState->resize(decDomain);
    }
    if(geomState) {
      geomState->resize(decDomain, muCopy);
    }
    if(stepState) {
      geomState->resize(decDomain);
    }
    if(stateIncr) stateIncr->conservativeResize(solVecInfo());
    v.conservativeResize(solVecInfo());
    a.conservativeResize(solVecInfo());
    vp.conservativeResize(solVecInfo());
    force.conservativeResize(solVecInfo());
  }
}

SensitivityInfo*
MDNLDynamic::getSensitivityInfo() { return domain->senInfo; }

int
MDNLDynamic::getNumSensitivities() { return domain->getNumSensitivities(); }
