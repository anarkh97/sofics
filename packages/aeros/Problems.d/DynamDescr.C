#include <cstdlib>
#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>

#include <Problems.d/DynamDescr.h>

#include <Math.d/FullMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DiagMatrix.h>
#include <Solvers.d/SolverFactory.h>

#include <Utils.d/dofset.h>
#include <Solvers.d/Solver.h>
#include <Element.d/State.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>

#include <Control.d/ControlInterface.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/ModeData.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>

#include <Hetero.d/FlExchange.h>
#include <Driver.d/SysState.h>

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

#include <algorithm>

typedef FSFullMatrix FullMatrix;

extern ModeData modeDataMode;
extern int verboseFlag;

// SDDynamPostProcessor implementation

SDDynamPostProcessor::SDDynamPostProcessor(Domain *d, double *_bcx, double *_vcx, double *_acx,
                                           StaticTimers *_times, GeomState *_geomState,
                                           Corotator **_allCorot, FullSquareMatrix *_melArray,
                                           Vector *_reactions)
{ domain = d; bcx = _bcx; vcx = _vcx; acx = _acx; times = _times; geomState = _geomState;
  allCorot = _allCorot; melArray = _melArray; reactions = _reactions; dummy = 0; }

SDDynamPostProcessor::~SDDynamPostProcessor() {
  geoSource->closeOutputFiles();
  if(dummy) delete dummy;
}

void
SDDynamPostProcessor::openOutputFiles() {
  geoSource->openOutputFiles();
}

void
SDDynamPostProcessor::openOutputFilesForPita(int timeSliceRank) {
  geoSource->openOutputFilesForPita(timeSliceRank);
}

void
SDDynamPostProcessor::closeOutputFiles() {
  geoSource->closeOutputFiles();
}

void
SDDynamPostProcessor::closeOutputFilesForPita(int timeSliceRank) {
  geoSource->closeOutputFilesForPita(timeSliceRank);
}

double 
SDDynamPostProcessor::getKineticEnergy(Vector & vel, SparseMatrix * gMass) {
  return domain->getKineticEnergy(vel,gMass); 
}

void
SDDynamPostProcessor::dynamOutput(int tIndex, double time, DynamMat& dMat, Vector& ext_f, Vector *aeroForce, SysState<Vector> &state)
{
  startTimerMemory(times->output, times->memoryOutput);

  if(!aeroForce) { // check to avoid dereferencing a null pointer
    if(!dummy) dummy = new Vector(domain->numUncon(), 0.0);
    aeroForce = dummy;
  }

  this->fillBcxVcx(time);

  if(domain->solInfo().isNonLin() && domain->solInfo().nRestart > 0) {
    domain->writeRestartFile(time, tIndex, state.getVeloc(), state.getAccel(), geomState);
  } 

  domain->dynamOutput(tIndex, time, bcx, dMat, ext_f, *aeroForce, state.getDisp(), state.getVeloc(),
                      state.getAccel(), state.getPrevVeloc(), vcx, acx);

  // need to output the stresses for nonlinear using corotator functions for some elements (bt shell is an exception)
  // also rotation, angular velocity, angular acceleration, energy and reaction force outputs use geomState
  if(domain->solInfo().isNonLin()) {
    geomState->setVelocityAndAcceleration(state.getVeloc(), state.getAccel());
    int numOutInfo = geoSource->getNumOutInfo();
    OutputInfo *oinfo = geoSource->getOutputInfo();
    for(int iInfo = 0; iInfo < numOutInfo; ++iInfo) {
      if(oinfo[iInfo].isStressOrStrain() || oinfo[iInfo].isRotation()) {
        domain->postProcessingImpl(iInfo, geomState, ext_f, *aeroForce, time, tIndex, state.getVeloc().data(), vcx,
                                   allCorot, state.getAccel().data(), acx, geomState, reactions, dMat.M, dMat.C);
      }
    }
  }

  stopTimerMemory(times->output, times->memoryOutput);
}

void
SDDynamPostProcessor::pitaDynamOutput(int tIndex, DynamMat& dMat, Vector& ext_f, Vector *aeroForce, 
                                      SysState<Vector> &state, int sliceRank, double time)
{
  startTimerMemory(times->output, times->memoryOutput);

  this->fillBcxVcx(time);

  domain->pitaDynamOutput(tIndex, bcx, dMat, ext_f, *aeroForce, state.getDisp(), state.getVeloc(),
                          state.getAccel(), state.getPrevVeloc(), vcx, acx,
                          sliceRank, time);

  stopTimerMemory(times->output, times->memoryOutput);
}

// Update bcx for time dependent prescribed displacements and velocities
void
SDDynamPostProcessor::fillBcxVcx(double time)
{
  ControlLawInfo *claw = geoSource->getControlLaw();
  ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
  if(claw && claw->numUserDisp) {
    double *userDefineDisp = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
    double *userDefineVel  = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
    double *userDefineAcc  = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
    for(int i=0; i<claw->numUserDisp; ++i) {
      userDefineVel[i] = 0;
      userDefineAcc[i] = 0;
    }
    userSupFunc->usd_disp(time, userDefineDisp, userDefineVel, userDefineAcc);
    DofSetArray *dsa = domain->getDSA();
    for(int i = 0; i < claw->numUserDisp; ++i) {
      int dof = dsa->locate(claw->userDisp[i].nnum,1 << claw->userDisp[i].dofnum);
      if(dof >= 0) {
        bcx[dof] = userDefineDisp[i];
        vcx[dof] = userDefineVel[i];
        acx[dof] = userDefineAcc[i];
      }
    }
  }
} 

typedef FSFullMatrix FullMatrix;

SingleDomainDynamic::SingleDomainDynamic(Domain *d)
 : SingleDomainBase(d->solInfo())
{ 
  domain = d;
  if(domain->solInfo().sensitivity) allSens = new AllSensitivities<double>; 
  else allSens = 0;
  kelArray = 0; 
  melArray = 0;
  allCorot = 0;
  geomState = 0; 
  refState = 0;
  userDefineDisp = 0;
  bcx = 0;
  vcx = 0;
  acx = 0;
  claw = 0;
  userSupFunc = 0;
  times = 0;
  prevFrc = 0;
  prevFrcBackup = 0;

  flExchanger = domain->getFileExchanger();
  reactions = 0;
  firstSts = true;
}

SingleDomainDynamic::~SingleDomainDynamic()
{
  if(bcx) delete [] bcx;
  if(vcx) delete [] vcx;
  if(acx) delete [] acx;
  if(allSens) delete allSens;
  if(geomState) delete geomState;
  if(refState) delete refState;
  if(times) delete times;
  if(prevFrc) delete prevFrc;
  if(prevFrcBackup) delete prevFrcBackup;
  if(reactions) delete reactions;
  if(kelArray) { delete [] kelArray; kelArray = 0; }
  if(melArray) { delete [] melArray; melArray = 0; }
  if(allCorot) {

    for (int iElem = 0; iElem < domain->numElements(); ++iElem) {
      if(allCorot[iElem] && (allCorot[iElem] != dynamic_cast<Corotator*>(domain->getElementSet()[iElem])))
        delete allCorot[iElem];
    }

    delete [] allCorot;
    allCorot = 0;
  }
}

int
SingleDomainDynamic::getFilterFlag()
{
 int filterFlag = std::max(domain->solInfo().hzemFilterFlag, domain->solInfo().filterFlags);
 // note: only level 1 rbm filter is used for nonlinear
 return (domain->solInfo().isNonLin()) ? std::min(filterFlag,1) : filterFlag;
}

int
SingleDomainDynamic::solVecInfo()
{
 return domain->numUncon();
}

int
SingleDomainDynamic::masterSolVecInfo()
{
 return domain->numUncon();
}

int
SingleDomainDynamic::dbcVecInfo()
{
 return domain->numdof();
}

int
SingleDomainDynamic::bcInfo()
{
 return domain->nDirichlet();
}

int
SingleDomainDynamic::getTimeIntegration()
{
 return domain->solInfo().timeIntegration;
}

void
SingleDomainDynamic::getTimes(double &dt, double &tmax)
{
 dt   = domain->solInfo().getTimeStep(); // time step size
 tmax = domain->solInfo().tmax; // total time
 if(userSupFunc)
   userSupFunc->setDt(dt);
}

void
SingleDomainDynamic::getInitialTime(int &initTimeIndex, double &initTime)
{
 initTimeIndex = domain->solInfo().initialTimeIndex;
 initTime      = domain->solInfo().initialTime;
 t0 = initTime;
}

double
SingleDomainDynamic::getInitialForceNorm()
{
  return domain->solInfo().initExtForceNorm;
}

void
SingleDomainDynamic::getNewMarkParameters(double &beta, double &gamma,
                                          double &alphaf, double &alpham)
{
 beta  = domain->solInfo().newmarkBeta;
 gamma = domain->solInfo().newmarkGamma;
 alphaf = domain->solInfo().newmarkAlphaF;
 alpham = domain->solInfo().newmarkAlphaM;
}

void
SingleDomainDynamic::getQuasiStaticParameters(double &maxVel, double &delta) 
{
 maxVel  = domain->solInfo().qsMaxvel;
 delta  = domain->solInfo().delta;
}

void
SingleDomainDynamic::getRayleighCoef(double &alpha)
{
 alpha = domain->solInfo().alphaDamp;
}

void
SingleDomainDynamic::getSteadyStateParam(int &steadyFlag, int &steadyMin, 
                                         int &steadyMax, double &steadyTol)
{
 steadyFlag  = domain->solInfo().steadyFlag;
 steadyMin   = domain->solInfo().steadyMin;
 steadyMax   = domain->solInfo().steadyMax;
 steadyTol   = domain->solInfo().steadyTol;
}

void
SingleDomainDynamic::getSensitivityStateParam(double &sensitivityTol, double &ratioSensitivityTol)
{
 sensitivityTol = domain->solInfo().sensitivityTol;
 ratioSensitivityTol = domain->solInfo().ratioSensitivityTol;
}

void
SingleDomainDynamic::computeStabilityTimeStep(double& dt, DynamMat& dMat)
{
 // ... Compute Stability Time Step
 double sts;
 int eid;
 if(domain->solInfo().isNonLin())
   sts = domain->computeStabilityTimeStep(kelArray, melArray, geomState, eid);
 else
   sts = domain->computeStabilityTimeStep(dMat);

 if(sts == std::numeric_limits<double>::infinity()) {
   filePrint(stderr," **************************************\n");
   filePrint(stderr," Stability max. timestep could not be  \n");
   filePrint(stderr," determined for this model.            \n");
   if(domain->solInfo().isNonLin() && eid > -1) {
     filePrint(stderr," Element with inf. time step = %7d\n",eid+1);
   }
   filePrint(stderr," Specified time step is selected\n");
   filePrint(stderr," **************************************\n");
   domain->solInfo().stable = 0;
 }
 else {
   filePrint(stderr," **************************************\n");
   if (domain->solInfo().modifiedWaveEquation) {
     sts = 1.73205*sts;
     filePrint(stderr," CONDITIONALLY STABLE MODIFIED WAVE EQUATION\n");
   }
   else
     filePrint(stderr," CONDITIONALLY STABLE NEWMARK ALGORITHM\n");

   filePrint(stderr," --------------------------------------\n");
   filePrint(stderr," Specified time step      = %10.4e\n",dt);
   filePrint(stderr," Stability max. time step = %10.4e\n",sts);
   if(domain->solInfo().isNonLin()) {
     filePrint(stderr," Element with min. time step = %7d\n",eid+1);
   }
   filePrint(stderr," **************************************\n");
   if( (domain->solInfo().stable == 1 && sts < dt) || domain->solInfo().stable == 2 ) {
     dt = sts;
     filePrint(stderr," Stability max. time step is selected\n");
   } else
     filePrint(stderr," Specified time step is selected\n");
   filePrint(stderr," **************************************\n");
 }
 if(getAeroAlg() < 0 && !firstSts && domain->solInfo().printNumber > 0) filePrint(stderr, " ⌈\x1B[33m Time Integration Loop In Progress: \x1B[0m⌉\n");
 firstSts = false;

 domain->solInfo().setTimeStep(dt);
}

void
SingleDomainDynamic::getInitState(SysState<Vector> &inState)
{
  // initialize state with IDISP/IDISP6/IVEL/IACC or RESTART
  domain->initDispVeloc(inState.getDisp(),  inState.getVeloc(),
                        inState.getAccel(), inState.getPrevVeloc()); // IVEL, IDISP, IDISP6, restart for linear 

  // apply rbmfilter projection
  if(sinfo.filterFlags || sinfo.hzemFilterFlag) {
    project(inState.getDisp());
    project(inState.getVeloc());
    project(inState.getAccel());
    project(inState.getPrevVeloc());
  }

  if(geoSource->getCheckFileInfo()->lastRestartFile) { 
    filePrint(stderr, " ... Restarting From a Previous Run (SingleDomainDynamic) ...\n");
    if(domain->solInfo().isNonLin()) { // restart for nonlinear
      domain->readRestartFile(inState.getDisp(), inState.getVeloc(), inState.getAccel(),
                              inState.getPrevVeloc(), bcx, vcx, *geomState);
      domain->updateStates(geomState, *geomState, allCorot, domain->solInfo().initialTime);
    }
  }
  else if(geomState) {
    geomState->update(inState.getDisp());
    geomState->setVelocityAndAcceleration(inState.getVeloc(), inState.getAccel());
  }

  // if we have a user supplied function, give it the initial state at the sensors
  // .. first update bcx, vcx in case any of the sensors have prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      double *userDefineAcc = new double[claw->numUserDisp];
      for(int i=0; i<claw->numUserDisp; ++i) {
        userDefineVel[i] = 0;
        userDefineAcc[i] = 0;
      }
      userSupFunc->usd_disp(domain->solInfo().initialTime, userDefineDisp, userDefineVel, userDefineAcc);
      setBC(userDefineDisp, userDefineVel, userDefineAcc);
      delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
    }
    if(claw->numSensor) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel = new double[claw->numSensor];
      double *ctracc = new double[claw->numSensor];
      extractControlData(inState, ctrdisp, ctrvel, ctracc);
      userSupFunc->init(ctrdisp, ctrvel, ctracc, this);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc;
    }
  }

}

void
SingleDomainDynamic::extractControlData(SysState<Vector> &state, double *ctrdsp, 
                                        double *ctrvel, double *ctracc)
{
  // get SENSOR state
  Vector &dsp = state.getDisp();
  Vector &vlc = state.getVeloc();
  Vector &acc = state.getAccel();

  DofSetArray *cdsa = domain->getCDSA();
  DofSetArray *dsa = domain->getDSA();

  for(int i = 0; i < claw->numSensor; ++i) {
    int dof = cdsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
    if(dof >= 0) { // free
      ctrdsp[i] = dsp[dof];
      ctrvel[i] = vlc[dof];
      ctracc[i] = acc[dof];
    }
    else { // either constrained or non-existant
      int dof2 = dsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
      if(dof2 >= 0) { // constrained
        ctrdsp[i] = bcx[dof2];
        ctrvel[i] = vcx[dof2];
        ctracc[i] = acx[dof2];
      }
    }
  }
}

/*
void
SingleDomainDynamic::addCtrl(Vector &f, double *ctrfrc)
{
  // add ACTUATOR forces
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numActuator; ++i) {
    int dof = cdsa->locate(claw->actuator[i].nnum, 1 << claw->actuator[i].dofnum);
    if(dof >= 0) f[dof] += ctrfrc[i];
  }
}

void
SingleDomainDynamic::addUserForce(Vector&f, double *userDefineForce)
{
  // add USDF forces
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numUserForce; ++i) {
    int dof = cdsa->locate(claw->userForce[i].nnum, 1<<claw->userForce[i].dofnum);
    if(dof >= 0) f[dof] += userDefineForce[i];
  }
}
*/

void
SingleDomainDynamic::setBC(double *userDefineDisp, double* userDefineVel, double *userDefineAcc)
{
  // update time-dependent prescribed displacements and velocities
  DofSetArray *dsa = domain->getDSA();
  for(int i = 0; i < claw->numUserDisp; ++i) {
    int dof = dsa->locate(claw->userDisp[i].nnum,1 << claw->userDisp[i].dofnum);
    if(dof >= 0) {
      bcx[dof] = userDefineDisp[i];
      vcx[dof] = userDefineVel[i];
      acx[dof] = userDefineAcc[i];
    }
  }
}

void
SingleDomainDynamic::getConstForce(Vector &cnst_f)
{
  domain->computeConstantForce(cnst_f, kuc);
}

void
SingleDomainDynamic::addConstForceSensitivity(Vector &cnst_fSen)
{
  domain->addConstantForceSensitivity(cnst_fSen, kuc);
}

void
SingleDomainDynamic::getContactForce(Vector &d_n, Vector &dinc, Vector &ctc_f, double t_n_p, double dt, double dt_old)
{
  ctc_f.zero();
  if(t_n_p < domain->solInfo().tdenforceInitia || t_n_p >= domain->solInfo().tdenforceFinal) return;
  if(domain->tdenforceFlag()) {
    times->tdenforceTime -= getTime();

    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaceTopology(); // remove deleted elements
    domain->UpdateSurfaces(geomState, 1); // update to current configuration
    times->updateSurfsTime += getTime();

    // copy and update the current state (geomState) to the predicted state
    GeomState *predictedState = new GeomState(*geomState);
    if(domain->solInfo().isNonLin()) {
      predictedState->update(dinc,1);
    }
    else {
      Vector d_n_p(domain->numUncon());
      d_n_p = d_n + dinc;
      predictedState->explicitUpdate(domain->getNodes(), d_n_p);
    }

    // update the prescribed displacements to their correct value at the time of the predictor
    if(claw && userSupFunc && claw->numUserDisp) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      double *userDefineAcc = new double[claw->numUserDisp];
      for(int i=0; i<claw->numUserDisp; ++i) {
        userDefineVel[i] = 0;
        userDefineAcc[i] = 0;
      }
      userSupFunc->usd_disp(t_n_p, userDefineDisp, userDefineVel, userDefineAcc);
      predictedState->updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes(), userDefineVel, userDefineAcc);
      delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
    }

    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaces(predictedState, 2); // update to predicted configuration
    times->updateSurfsTime += getTime();

    times->contactSearchTime -= getTime();
    domain->PerformDynamicContactSearch(dt_old, dt);
    times->contactSearchTime += getTime();

    times->contactForcesTime -= getTime();
    domain->AddContactForces(dt_old, dt, ctc_f);
    times->contactForcesTime += getTime();

    delete predictedState;
    times->tdenforceTime += getTime();
  }
}

void
SingleDomainDynamic::updateState(double dt_n_h, Vector &v_n_h, Vector &d_n)
{
  if(domain->solInfo().isNonLin()) {
    *refState = *geomState; // (AN) update refState values
    Vector dinc(solVecInfo()); dinc = dt_n_h*v_n_h;
    geomState->update(dinc, 1);
    geomState->setVelocity(v_n_h);
    geomState->get_tot_displacement(d_n, false);
  }
}

void
SingleDomainDynamic::computeExtForce2(SysState<Vector> &state, Vector &ext_f, 
                                      Vector &cnst_f, int tIndex, double t,
                                      Vector *aero_f, double gamma, double alphaf, double *pt_dt)
{
  times->formRhs -= getTime();
  SolverInfo& sinfo = domain->solInfo();

  ext_f.zero();
  double *userDefineDisp = 0;
  double *userDefineVel = 0;
  double *userDefineAcc = 0;

  if(claw && userSupFunc) {
    if(claw->numUserDisp) { // USDD
      userDefineDisp = new double[claw->numUserDisp];
      userDefineVel = new double[claw->numUserDisp];
      userDefineAcc = new double[claw->numUserDisp];
      for(int i=0; i<claw->numUserDisp; ++i) {
        userDefineVel[i] = 0;
        userDefineAcc[i] = 0;
      }
      userSupFunc->usd_disp(t, userDefineDisp, userDefineVel, userDefineAcc);
      setBC(userDefineDisp, userDefineVel, userDefineAcc); // update bcx, vcx, acx
      domain->updateUsddInDbc(userDefineDisp);
    }
    if(claw->numUserForce) { // USDF
      double *userDefineForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(t, userDefineForce);
      domain->updateUsdfInNbc(userDefineForce);
      delete [] userDefineForce;
    }
    if(claw->numActuator > 0) { // SENSORS+ACTUATORS
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel = new double[claw->numSensor];
      double *ctracc = new double[claw->numSensor];
      double *ctrfrc = new double[claw->numActuator];

      extractControlData(state, ctrdisp, ctrvel, ctracc);
      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t, &state, &ext_f);
      domain->updateActuatorsInNbc(ctrfrc);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }

  // finish update of geomState. note that for nonlinear problems the unconstrained positiion and rotation
  // nodal variables have already been updated in updateDisplacement
  if(sinfo.isNonLin() || domain->tdenforceFlag()) {
    if(!sinfo.isNonLin()) {
      geomState->explicitUpdate(domain->getNodes(), state.getDisp());
    }
    if(userDefineDisp) geomState->updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes(),
                                                               userDefineVel, userDefineAcc);
  }

  if(sinfo.isNonLin() && domain->GetnContactSurfacePairs() && !domain->tdenforceFlag()) {
    domain->UpdateSurfaces(MortarHandler::CTC, geomState);
    domain->PerformStaticContactSearch(MortarHandler::CTC);
    domain->deleteSomeLMPCs(mpc::ContactSurfaces);
    domain->ExpComputeMortarLMPC(MortarHandler::CTC);
    domain->UpdateContactSurfaceElements(geomState);
    domain->makeAllDOFs();
    domain->createContactCorotators(allCorot, kelArray, melArray);
  }

  // THERMOE update nodal temperatures
  if(sinfo.thermoeFlag >= 0 && tIndex >= 0) {
    domain->thermoeComm();
    if(geomState) geomState->setNodalTemperatures(domain->getNodalTemperatures());
  }

  // add f(t) to cnst_f
  // for linear problems also add contribution of non-homogeneous dirichlet (DISP/TEMP/USDD etc)
  double dt = sinfo.getTimeStep();
  double alpham = sinfo.newmarkAlphaM;
  double t0 = sinfo.initialTime;
  double tm = (t == t0) ? t0 : t + dt*(alphaf-alpham);
  domain->computeExtForce4(ext_f, cnst_f, t, kuc, userSupFunc, cuc, tm, muc);
  if(userDefineDisp) delete [] userDefineDisp;
  if(userDefineVel) delete [] userDefineVel;
  if(userDefineAcc) delete [] userDefineAcc;

  // add aeroelastic forces from fluid dynamics code
  if(sinfo.aeroFlag >= 0 && tIndex >= 0 &&
     !(geoSource->getCheckFileInfo()->lastRestartFile && sinfo.aeroFlag == 20 && !sinfo.dyna3d_compat && tIndex == sinfo.initialTimeIndex)) {
    if(tIndex % sinfo.subcycle == 0) {
      domain->buildAeroelasticForce(*aero_f, *prevFrc, tIndex, t, gamma, alphaf);
    }
    ext_f += *aero_f;
  }

  // add aerothermal fluxes from fluid dynamics code
  if(sinfo.aeroheatFlag >= 0 && tIndex >= 0) 
    domain->buildAeroheatFlux(ext_f, prevFrc->lastFluidLoad, tIndex, t);

  // apply rbmfilter projection
  if(sinfo.filterFlags || sinfo.hzemFilterFlag)
    trProject(ext_f); 

  if(tIndex == 1)
    sinfo.initExtForceNorm = ext_f.norm();

  times->formRhs += getTime();

}

void
SingleDomainDynamic::getAeroelasticForceSensitivity(int tIndex, double t,
                                                    Vector *aero_f, double gamma, double alphaf)
{
  // get aeroelastic force sensitivity from fluid code
  domain->buildAeroelasticForce(*aero_f, *prevFrc, tIndex, t, gamma, alphaf);
}

void
SingleDomainDynamic::processLastOutput()
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  domain->solInfo().lastIt = true;
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

void
SingleDomainDynamic::preProcessSA()
{
  domain->buildPreSensitivities<double>(*allSens, bcx);
}

void
SingleDomainDynamic::postProcessSA(DynamMat *dMat, Vector &sol)
{
  domain->buildPostSensitivities<double>(dMat->dynMat, dMat->K, dMat->K, *allSens, &sol, bcx, true);
}

void
SingleDomainDynamic::sensitivityPostProcessing(Vector *sol)
{
  domain->sensitivityPostProcessing(*allSens, sol, bcx);
}

void
SingleDomainDynamic::preProcess()
{
  // Allocate space for Timers
  times = new StaticTimers;
  startTimerMemory(times->preProcess, times->memoryPreProcess);

  // Makes renumbering, connectivities and dofsets
  domain->preProcessing();

  // initialize bcx and vcx with time-independent prescribed displacements and velocities
  // note: bcx and vcx are updated with time-dependent prescribed displacements and velocities (USDD) later
  times->makeBCs -= getTime();
  int numdof = domain->numdof();
  int *bc = new int[numdof];
  bcx = new double[numdof];
  domain->make_bc(bc, bcx);
  delete [] bc;
  vcx = new double[numdof]; 
  acx = new double[numdof];
  for(int i=0; i<numdof; ++i) { acx[i] = vcx[i] = 0.0; }
  times->makeBCs += getTime();

  times->makeDOFs -= getTime();
  domain->make_constrainedDSA();
  domain->makeAllDOFs();
  times->makeDOFs += getTime();

  // Check for user supplied routines (control, force or displacement)
  claw = geoSource->getControlLaw();
  userSupFunc = domain->getUserSuppliedFunction();

  if((domain->numInitDisp6() > 0 && domain->solInfo().gepsFlg == 1) || domain->solInfo().isNonLin()) {
    FullSquareMatrix *geomKelArray=0;
    // this function builds corotators, geomstate and kelArray 
    // for linear+geps only it updates geomState with IDISP6 and computes the element stiffness matrices using this updated geomState
    bool melFlag = (domain->solInfo().isNonLin() && domain->solInfo().newmarkBeta == 0);
       // currently we only need to store the element mass matrices
       // for nonlinear explicit dynamics to compute the critical time step
    domain->computeGeometricPreStress(allCorot, geomState, kelArray, times, geomKelArray, melArray, melFlag);
  }
  else if(domain->tdenforceFlag())
    geomState = new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes(), &domain->getElementSet(),
                              domain->getNodalTemperatures());

  if(domain->solInfo().isNonLin() || domain->tdenforceFlag()) {
    // for nonlinear explicit we only need to initialize geomState with the constant constrained displacements (DISP).
    // the geomState is always updated before use with the current unconstrained displacments plus any time-dependent constrained displacements (USDD)
    if(domain->nDirichlet() > 0) { 
      geomState->updatePrescribedDisplacement(domain->getDBC(), domain->nDirichlet(), domain->getNodes()); 
    }
  }

  // call GeomStateCopy constructor to initialize refState
  // assign refState for the first time step
  refState = new GeomState(*geomState);

  if(domain->tdenforceFlag())
    domain->InitializeDynamicContactSearch();
  else if(domain->solInfo().isNonLin() && domain->GetnContactSurfacePairs())
    domain->InitializeStaticContactSearch(MortarHandler::CTC);

  prevFrc = new PrevFrc(domain->numUncon());
  prevFrcBackup = new PrevFrc(domain->numUncon());
  reactions = new Vector(domain->nDirichlet());

  stopTimerMemory(times->preProcess, times->memoryPreProcess);
}

// In general a symmetric matrix can be partitioned as:
//
// A = | A11   A12 |
//     | A12^t A22 |
//

DynamMat *
SingleDomainDynamic::buildOps(double coeM, double coeC, double coeK)
{
 AllOps<double> allOps;
 DynamMat *dMat = new DynamMat(true);

 domain->getTimers().constructTime -= getTime();
 allOps.K   = domain->constructDBSparseMatrix<double>();
 if(geoSource->getMRatio() != 0) {
#ifdef USE_EIGEN3
   if(domain->solInfo().svdPodRom) allOps.M = domain->constructEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >();
   else
#endif
   allOps.M = domain->constructDBSparseMatrix<double>();
 }
 else {
   allOps.M = new DiagMatrix(domain->getCDSA());
 }
 allOps.Muc = domain->constructCuCSparse<double>();
 allOps.Kuc = domain->constructCuCSparse<double>();
 allOps.Mcc = domain->constructCCSparse<double>();
 allOps.Kcc = domain->constructCCSparse<double>();

 // Rayleigh Damping coefficients
 double alpha = domain->solInfo().alphaDamp;
 double beta  = domain->solInfo().betaDamp;

 // Damping Matrix: C = alpha*M + beta*K + D
 if(alpha != 0.0 || beta != 0.0 || domain->getElementSet().hasDamping() || domain->solInfo().ATDARBFlag != -2.0) {
   allOps.C   = domain->constructDBSparseMatrix<double>();
   allOps.Cuc = domain->constructCuCSparse<double>();
   allOps.Ccc = domain->constructCCSparse<double>();
 }

 // to compute a^0 = M^{-1}(f_ext^0-f_int^0-Cu^0)
 if(getTimeIntegration() != 1 && (domain->solInfo().newmarkBeta != 0.0 && domain->solInfo().iacc_switch)
    && domain->solInfo().solvercntl->type == SolverSelection::Direct) { // not required for explicit
   SparseMatrix *spp; Solver *prec; // XXX
   SolverCntl *m_cntl = (domain->solInfo().solvercntl->type == SolverSelection::Direct) ? domain->solInfo().solvercntl : &default_cntl;
   dMat->Msolver = GenSolverFactory<double>::getFactory()->createSolver(domain->getNodeToNode(), domain->getDSA(), domain->getCDSA(), 
                                                                        *m_cntl, allOps.Msolver, (Rbm*) NULL, spp, prec);
 }
 domain->getTimers().constructTime += getTime();

 Rbm *rigidBodyModes = 0;

 int useRbmFilter = domain->solInfo().filterFlags && !domain->solInfo().isNonLin();
 int useGrbm = domain->solInfo().rbmflg; 
 int useHzemFilter = domain->solInfo().hzemFilterFlag && !domain->solInfo().isNonLin();
 int useHzem = domain->solInfo().hzemFlag;

 if(useGrbm || useRbmFilter) 
   rigidBodyModes = domain->constructRbm();
 else if(useHzem || useHzemFilter)
   rigidBodyModes = domain->constructHzem();

 if((getTimeIntegration() == 1) && (useGrbm || useHzem) && !domain->solInfo().isNonLin()) // only use rigidBodyModes for linear quasistatics
   domain->buildOps(allOps, coeK, coeM, coeC, rigidBodyModes, kelArray, melArray);
 else
   domain->buildOps(allOps, coeK, coeM, coeC, 0, kelArray, melArray);

 if(useRbmFilter)
   filePrint(stderr," ... RBM Filter Level %d Requested   ...\n", useRbmFilter);
 else if(useHzemFilter)
   filePrint(stderr," ... HZEM Filter Requested          ...\n");

 if(useRbmFilter || useHzemFilter)
   projector_prep(rigidBodyModes, allOps.M);
 
 // Modal decomposition preprocessing
 int decompFlag = domain->solInfo().modeDecompFlag;
 if(decompFlag) {
   filePrint(stderr," ... Modal decomposition requested  ...\n");
   modeDecompPreProcess(allOps.M);
 }

 dMat->K         = allOps.K;
 dMat->M         = allOps.M;
 dMat->C         = allOps.C;
 cuc = dMat->Cuc = allOps.Cuc;
 dMat->Ccc       = allOps.Ccc;
 muc = dMat->Muc = allOps.Muc;
 dMat->Mcc       = allOps.Mcc;
 kuc = dMat->Kuc = allOps.Kuc;
 dMat->Kcc       = allOps.Kcc;
 dMat->dynMat    = allOps.sysSolver;
 if(dMat->Msolver) {
   if(verboseFlag) filePrint(stderr," ... Factoring mass matrix for iacc ...\n");
   domain->getTimers().factor -= getTime();
   dMat->Msolver->factor();
   domain->getTimers().factor += getTime();
 }

 if(domain->tdenforceFlag()) domain->MakeNodalMass(allOps.M, allOps.Mcc); 

 return dMat;
}

void
SingleDomainDynamic::getInternalForce(Vector& d, Vector& f, double t, int tIndex)
{
  if(domain->solInfo().isNonLin()) {
    Vector residual(domain->numUncon(),0.0);
    Vector fele(domain->maxNumDOF());
    if(reactions) reactions->zero();
    // NOTE #1: for explicit nonlinear dynamics, geomState and refState are the same object -- AN no longer true
    // NOTE #2: by convention, the internal variables associated with a nonlinear constitutive relation are not updated
    //          when getStiffAndForce is called, so we have to call updateStates.
    if(domain->solInfo().newmarkBeta == 0 && domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
      domain->getStiffAndForce(*geomState, fele, allCorot, kelArray, residual, 1.0, t, refState,
                               reactions, melArray);
/* PJSA 10/12/2014 this is done in getStiffAndForce now because it needs to be done before handleElementDeletion.
      domain->updateStates(geomState, *geomState, allCorot);
*/
    }
    else {
      // refState store data at t_n; and is sent to calculate intenal forces for foam like materials

      //domain->getInternalForce(*geomState, fele, allCorot, kelArray, residual, 1.0, t, geomState,
      //                         reactions, melArray);
      domain->getInternalForce(*geomState, fele, allCorot, kelArray, residual, 1.0, t, refState,
                               reactions, melArray);
    }
    f.linC(-1.0,residual); // f = -residual
  }
  else {
    f.zero();
    domain->getKtimesU(d, bcx, f, 1.0, kelArray);  // note: although passed as an argument, the bcx contribution is not computed in this function
  }
}

void
SingleDomainDynamic::pull_back(Vector& f)
{
  if(domain->solInfo().isNonLin() && !domain->solInfo().galerkinPodRom && !domain->solInfo().getNLInfo().linearelastic) {
    // Transform both moments and forces to convected frame: f = [R^T  I ]*f
    //                                                           [ I  R^T]
    geomState->pull_back(f); 
  }
}

void
SingleDomainDynamic::push_forward(Vector &a)
{
  if(domain->solInfo().isNonLin() && !domain->solInfo().galerkinPodRom && !domain->solInfo().getNLInfo().linearelastic) {
    // Transform 2nd time-derivative of displacement to spatial frame: a = [R I]*a
    //                                                                     [I I]
    // Note: the angular accelerations are deliberately not transformed.
    geomState->push_forward(a);
  }
}

void
SingleDomainDynamic::getUnassembledNonLinearInternalForce(Vector& d, Vector& uf, std::map<int, std::pair<int,int> > &uDOFaDOFmap, FullSquareMatrix *kelCopy, double t, int tIndex)
{// this function is only used to construct unassembled force snapshots in the UDEIM ROM driver
  Vector unassemResidual(uf.size(),0.0);
  Vector fele(domain->maxNumDOF());

  domain->getUnassembledNonLinearInternalForce(*geomState, fele, allCorot, kelArray, unassemResidual, uDOFaDOFmap, 1.0, t, tIndex, geomState,
                          (Vector*) NULL, melArray,kelCopy);

  uf.linC(-1.0,unassemResidual); // uf = -residual

}

void
SingleDomainDynamic::computeTimeInfo()
{
  // Time integration information

  // Get total time and time step size and store them
  double totalTime = domain->solInfo().tmax;
  double dt        = domain->solInfo().getTimeStep();

  // Compute maximum number of steps
  int maxStep = (int) ( (totalTime+0.49*dt)/dt );

  // Compute time remainder
  double remainder = totalTime - maxStep*dt;
  if(std::abs(remainder)>0.01*dt){
    domain->solInfo().tmax = maxStep*dt;
    filePrint(stderr, " Warning: Total time is being changed to : %e\n", domain->solInfo().tmax);
  }
}

int 
SingleDomainDynamic::aeroPreProcess(Vector& d_n, Vector& v_n, 
                                    Vector& a_n, Vector& v_p)
{
  domain->aeroPreProcess(d_n, v_n, a_n, v_p, bcx, vcx);
  return domain->solInfo().aeroFlag;
}

int 
SingleDomainDynamic::aeroSensitivityPreProcess(Vector& d_n, Vector& v_n, 
                                               Vector& a_n, Vector& v_p)
{
  domain->aeroSensitivityPreProcess(d_n, v_n, a_n, v_p, bcx, vcx);
  return domain->solInfo().aeroFlag;
}

int
SingleDomainDynamic::sendDisplacements(Vector& d_n, Vector& v_n,
                                       Vector& a_n, Vector& v_p)
{
  domain->sendDisplacements(d_n, v_n, a_n, v_p, bcx, vcx);
  return domain->solInfo().aeroFlag;
}

void
SingleDomainDynamic::a5TimeLoopCheck(int& parity, double& t, double dt)
{
  if(domain->solInfo().aeroFlag == 5) {
     if(!parity) t -= dt;
     parity = ( parity ? 0 : 1 );
  }
}

void
SingleDomainDynamic::a5StatusRevise(int parity, SysState<Vector>& curState, 
                                    SysState<Vector>& bkState)
{
  if(domain->solInfo().aeroFlag == 5) {
     if(parity) { // restore

       *prevFrc = *prevFrcBackup;
       curState.getDisp()      = bkState.getDisp();
       curState.getVeloc()     = bkState.getVeloc();
       curState.getAccel()     = bkState.getAccel();
       curState.getPrevVeloc() = bkState.getPrevVeloc();

     } else {      // backup

       *prevFrcBackup = *prevFrc;
       bkState.getDisp()      = curState.getDisp();
       bkState.getVeloc()     = curState.getVeloc();
       bkState.getAccel()     = curState.getAccel();
       bkState.getPrevVeloc() = curState.getPrevVeloc();

     }
  }
}

void
SingleDomainDynamic::thermoePreProcess(Vector&, Vector&, Vector&)
{
  domain->thermoePreProcess();
  if(geomState) geomState->setNodalTemperatures(domain->getNodalTemperatures());
}

void
SingleDomainDynamic::aeroHeatPreProcess(Vector& d_n, Vector& v_n, Vector& v_p)
{
  domain->aeroHeatPreProcess(d_n, v_n, v_p, bcx);
}

void
SingleDomainDynamic::thermohPreProcess(Vector& d_n, Vector& v_n, Vector& v_p)
{
  domain->thermohPreProcess(d_n, v_n, v_p, bcx);
}

// Single Domain Post Processor Functions

SDDynamPostProcessor *
SingleDomainDynamic::getPostProcessor()
{
  return new SDDynamPostProcessor(domain, bcx, vcx, acx, times, geomState, allCorot, melArray, reactions);
}

void
SingleDomainDynamic::printTimers(DynamMat *dynamMat, double timeLoop)
{
  long memoryUsed = (dynamMat->dynMat)->size();
  double solveTime  = dynamMat->dynMat->getSolutionTime();

  times->printStaticTimers(solveTime, memoryUsed, domain, timeLoop);

  if(domain->solInfo().massFlag) {
    double mass = domain->computeStructureMass();
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," ... Structure mass = %e  ...\n",mass);
    filePrint(stderr," --------------------------------------\n");
  }
}

double
SingleDomainDynamic::betaDamp() const {
  return domain->solInfo().betaDamp;
}

double
SingleDomainDynamic::alphaDamp() const {
  return domain->solInfo().alphaDamp;
}

void 
SingleDomainDynamic::setDamping( double betaDamp, double alphaDamp )
{ 
  domain->solInfo().setDamping( betaDamp, alphaDamp ); 
}

int
SingleDomainDynamic::getModeDecompFlag()
{
 return domain->solInfo().modeDecompFlag;
}   

void
SingleDomainDynamic::modeDecompPreProcess(SparseMatrix *M)
{
  int eigsize;
 
  if(sinfo.mode_modal_id < 0) {

    // Read Eigenvectors from file EIGENMODES
    // ======================================

    BinFileHandler modefile("EIGENMODES" ,"r");
    modefile.read(&maxmode, 1);
 
    modefile.read(&eigsize, 1);
 
    eigmodes = new double*[maxmode];
    for (int i = 0; i<maxmode; ++i)
      eigmodes[i] = new double[eigsize];

    // Check if the problem sizes are identical

    int numdof = solVecInfo();

    if (eigsize != numdof)
      fprintf(stderr, " *** ERROR: EigenVector and problem sizes differ ...\n");

    for (int i = 0; i<maxmode; ++i)
      modefile.read(eigmodes[i], eigsize);
  }
  else {

    // Use vectors read from file specified with READMODE command
    // ==========================================================

    maxmode = modeDataMode.numModes;
    DofSetArray *cdsa = domain->getCDSA();
    eigsize = cdsa->size();
    eigmodes = new double*[maxmode];
    for (int i = 0; i<maxmode; ++i) {
      eigmodes[i] = new double[eigsize];
      std::fill(eigmodes[i], eigmodes[i]+eigsize, 0.);
      BCond y; y.nnum = i; y.val = 1.;
      modeDataMode.addMultY(1, &y, eigmodes[i], cdsa);
    }
  }

  tPhiM = new double*[maxmode];
    for (int i = 0; i<maxmode; ++i)
      tPhiM[i] = new double[eigsize];
 
  if(sinfo.mode_modal_id < 0 || sinfo.readInModes[sinfo.mode_modal_id].type != ModalParams::Inorm) {

    // Compute Transpose(Phi_i)*M once and for all
    // ===========================================

    if(verboseFlag) filePrint(stderr, " ... Preparing Transpose(Phi_i)*M for %d modes ...\n", maxmode);

    for (int i = 0; i<maxmode; ++i)
      M->mult(eigmodes[i], tPhiM[i]);  // taking advantage of symmetry of M and computing
                                       // M*Phi_i instead of transpose(Phi_i)*M
/*
    // Verify that Phi_i*M*Phi_i = 1 and Phi_i*M*Phi_j = 0
    // NOTE: the following assumes M is diagonal
    // ===================================================

    double PhiMPhi = 0;

    for (int i = 0; i < maxmode; ++i) {
      for (int j = 0; j < maxmode; ++j) {
        for (int k = 0; k < eigsize; ++k) {
          PhiMPhi = PhiMPhi + eigmodes[i][k]*M->diag(k)*eigmodes[j][k];
        }
        fprintf(stderr, "Phi_%d*M*Phi_%d = %19.11e\n",i,j, PhiMPhi) ;
        PhiMPhi = 0;
      }
      fprintf(stderr,"\n");
    }
*/
  }
  else {
    for (int i = 0; i<maxmode; ++i)
      std::copy(eigmodes[i], eigmodes[i]+eigsize, tPhiM[i]);
  }
}

void
SingleDomainDynamic::modeDecomp(double t, int tIndex, Vector& d_n)
{
// Compute Alpha and error only if their output file is specified,
// otherwise, it wouldn't make sense

   // Compute alfa_i=PhiDiag_i*d_n
 
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int i, j, k;
  alfa = 0;

  for (i=0; i < numOutInfo; i++) {
    if(oinfo[i].interval != 0 && tIndex % oinfo[i].interval == 0) {

      int w = oinfo[i].width;
      int p = oinfo[i].precision;

      switch(oinfo[i].type) {
        case OutputInfo::ModeAlpha: {

          if(!alfa) {
            alfa = new double[maxmode];

            for (k=0; k<maxmode; ++k)
              alfa[k] = 0;

            for (j = 0; j < maxmode; ++j)
              for (k = 0; k < d_n.size(); ++k)
                 alfa[j] += tPhiM[j][k]*d_n[k];
          }

          // Write alfa
          fprintf(oinfo[i].filptr, "%e  ", t);
          for(j=0; j<maxmode; ++j)
            fprintf(oinfo[i].filptr, "% *.*E ", w, p, alfa[j]);
          fprintf(oinfo[i].filptr, "\n");
     
          fflush(oinfo[i].filptr);
        } break;

        case OutputInfo::ModeError: {

          if(!alfa) {
            alfa = new double[maxmode];
  
            for (k=0; k<maxmode; ++k)
              alfa[k] = 0;
  
            for (j = 0; j < maxmode; ++j)
              for (k = 0; k < d_n.size(); ++k)
                alfa[j] += tPhiM[j][k]*d_n[k];
          }

          double sumerror = 0;
          double normerror = 0;
          double sumdisp = 0;
          double normdisp = 0;

          int ersize = d_n.size()-geoSource->internalNumNodes(); // don't include Lagrange multipliers

          double *sumalfa = new double[ersize];
          double *error   = new double[ersize];

          for (k=0; k < ersize; ++k) sumalfa[k] = 0;

          for (k=0; k < ersize; ++k)
            for (j=0; j < maxmode; ++j)
              sumalfa[k] += alfa[j]*eigmodes[j][k];


          for (j=0; j < ersize; ++j) {
            error[j] = d_n[j]-sumalfa[j];
            sumerror += error[j]*error[j];
            sumdisp += d_n[j]*d_n[j];
          }

          normdisp = sqrt(sumdisp);
          if (normdisp == 0.0)  normerror = 0.0;
          else normerror = sqrt(sumerror)/normdisp;

          // Write error
          fprintf(oinfo[i].filptr, "%e % *.*E\n", t, w, p, normerror);
          fflush(oinfo[i].filptr);
          delete [] sumalfa;
          delete [] error;
        } break;
   
        default: 
          break;
      }

    }
  }
  if(alfa) delete [] alfa;

}

SparseMatrix* 
SingleDomainDynamic::getpK(DynamMat* dynOps)
{
  return dynOps->K;
}

SparseMatrix* 
SingleDomainDynamic::getpM(DynamMat* dynOps)
{
  return dynOps->M;
}

SparseMatrix* getpC(DynamMat* dynOps)
{
  return dynOps->C;
}

void
SingleDomainDynamic::aeroSend(double time, Vector& d, Vector& v, Vector& a, Vector& v_p)
{
  if(claw && userSupFunc) {
    if(claw->numUserDisp) { // USDD
      // Note: the approprate value of "time" passed into this function should be t^{n+½} for A6 and C0, and
      // t^{n+1} otherwise, where t^n denotes the time at the end of the current structure timestep. Note that
      // the predictor in FlExchanger::sendDisplacements is not applied to prescribed displacements; we directly
      // compute here the desired values of the prescribed displacements/velocities rather than predicting them.
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      double *userDefineAcc = new double[claw->numUserDisp];
      for(int i=0; i<claw->numUserDisp; ++i) {
        userDefineVel[i] = 0;
        userDefineAcc[i] = 0;
      }
      userSupFunc->usd_disp(time, userDefineDisp, userDefineVel, userDefineAcc);
      setBC(userDefineDisp, userDefineVel, userDefineAcc); // update bcx, vcx, acx
      delete [] userDefineDisp;
      delete [] userDefineVel;
      delete [] userDefineAcc;
    }
  }
  startTimerMemory(times->output, times->memoryOutput);
  domain->aeroSend(d, v, a, v_p, bcx, vcx);
  stopTimerMemory(times->output, times->memoryOutput);
}

void 
SingleDomainDynamic::sendNumParam(int numParam, int actvar, double steadyTol)
{
  flExchanger->sendNumParam(numParam,actvar,steadyTol);
}

void 
SingleDomainDynamic::getNumParam(int &numParam)
{
  flExchanger->getNumParam(numParam);
}

void 
SingleDomainDynamic::sendRelativeResidual(double relres)
{
  flExchanger->sendRelativeResidual(relres);
}

int 
SingleDomainDynamic::cmdCom(int cmdFlag) 
{
  return flExchanger->cmdCom(cmdFlag);
}

int 
SingleDomainDynamic::getAeroAlg()
{ 
  return domain->solInfo().aeroFlag;
}

int
SingleDomainDynamic::getThermoeFlag()
{
  return domain->solInfo().thermoeFlag;
} 

int
SingleDomainDynamic::getAeroheatFlag()
{
  return domain->solInfo().aeroheatFlag;
}

int
SingleDomainDynamic::getThermohFlag()
{
  return domain->solInfo().thermohFlag;
}

#include <Problems.d/NonLinQStatic.h>
#include <Driver.d/NLStaticProbType.h>
void
SingleDomainDynamic::solveAndUpdate(Vector &force, Vector &dinc, Vector &d, double relaxFac, double time)
{
  int numElemStates = geomState->getTotalNumElemStates();
  if(!refState) {
    // For the first coupling cycle refState is the initial state as defined by either IDISP or restart, if specified.
    refState = new GeomState(*geomState);
  }
  else if(numElemStates == 0) {
    // In this case dlambda is only used for the first cycle.
    domain->solInfo().getNLInfo().dlambda = domain->solInfo().getNLInfo().maxLambda = 1.0;
  }

  NonLinQStatic nlstatic(domain, force, refState);
  NLStaticSolver<Solver,Vector,SingleDomainPostProcessor<double,Vector,Solver>,NonLinStatic,GeomState> nlsolver(&nlstatic);
  nlsolver.solve();

  nlsolver.getGeomState()->get_inc_displacement(dinc, *geomState, false);
  if(numElemStates == 0) {
    // In this case refState now stores the solution of the previous non-linear solve, and will be used as the initial
    // guess for the next coupling cycle's non-linear solve.
    *refState = *nlsolver.getGeomState();
  }

  dinc *= relaxFac;
  geomState->update(dinc);
  if(numElemStates != 0) {
    domain->updateStates(refState, *geomState, allCorot, time);
  }

  geomState->get_tot_displacement(d, false);
}

