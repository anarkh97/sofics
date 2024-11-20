#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <cstdlib>

#include <Driver.d/Domain.h>
#include <Element.d/MpcElement.d/MpcElement.h>
#include <Problems.d/NonLinDynam.h>
#include <Problems.d/DynamDescr.h>
#include <Solvers.d/Solver.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Control.d/ControlInterface.h>
#include <Corotational.d/TemperatureState.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/Corotator.h>
#include <Timers.d/GetTime.h>
#include <Math.d/FullMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>
#include <Element.d/State.h>
#include <Driver.d/SysState.h>
#include <Rom.d/PodProjectionSolver.h>
#include <Hetero.d/FlExchange.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <exception>

typedef FSFullMatrix FullMatrix;

extern int verboseFlag;

NonLinDynamic::NonLinDynamic(Domain *d) :
  SingleDomainBase(d->solInfo()),
  domain(d),
  bcx(0),
  vcx(0),
  acx(0),
  solver(NULL),
  spm(NULL),
  prec(NULL),
  spp(NULL),
  res(NULL),
  clawDofs(NULL),
  K(NULL),
  M(NULL),
  C(NULL),
  Kuc(NULL),
  allCorot(NULL),
  localTemp(),
  kelArray(NULL),
  celArray(NULL),
  melArray(NULL),
  prevFrc(NULL),
  secondRes(0.0),
  numSystems(0),
  times(NULL),
  userSupFunc(NULL),
  claw(NULL),
  Cuc(NULL),
  Ccc(NULL),
  Muc(NULL),
  Mcc(NULL),
  reactions(NULL),
  factor(false)
{
  if(domain->GetnContactSurfacePairs()) {
     domain->InitializeStaticContactSearch(MortarHandler::CTC);
     updateCS = true;
  }
  if(domain->solInfo().sensitivity) allSens = new AllSensitivities<double>;
  else allSens = 0;
}

NonLinDynamic::~NonLinDynamic()
{
  clean();
  if (res) {
    fclose(res); 
  }
  if(clawDofs) delete [] clawDofs;
  delete times;
  if(reactions) delete reactions;
  if(allSens) delete allSens;
}

void
NonLinDynamic::clean()
{
  if(prevFrc)  { delete prevFrc; prevFrc = 0; }
  if(bcx)      { delete [] bcx; bcx = 0; }
  if(vcx)      { delete [] vcx; vcx = 0; }
  if(acx)      { delete [] acx; acx = 0; }
  if(solver)   { delete solver; solver = 0; }
  if(prec)     { delete prec; prec = 0; }
  if(kelArray) { delete [] kelArray; kelArray = 0; }
  if(celArray) { delete [] celArray; celArray = 0; }
  if(melArray) { delete [] melArray; melArray = 0; }
  if(K)        { delete K; K = 0; }
  if(M)        { delete M; M = 0; }
  if(C)        { delete C; C = 0; }
  if(Kuc)      { delete Kuc; Kuc = 0; }
  if(Muc)      { delete Muc; Muc = 0; }
  if(Mcc)      { delete Mcc; Mcc = 0; }
  if(Cuc)      { delete Cuc; Cuc = 0; }
  if(Ccc)      { delete Ccc; Ccc = 0; }
  if(allCorot) {

    for (int iElem = 0; iElem < domain->numElements(); ++iElem) {
      if(allCorot[iElem] && (allCorot[iElem] != dynamic_cast<Corotator*>(domain->getElementSet()[iElem])))
        delete allCorot[iElem];
    }

    delete [] allCorot;
    allCorot = 0;
  }
}

void
NonLinDynamic::readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                               Vector &v_p, GeomState &geomState)
{
  if(geoSource->getCheckFileInfo()->lastRestartFile) {
    filePrint(stderr, " ... Restarting From a Previous Run (NonLinDynanmic) ...\n");

    domain->readRestartFile(d_n, v_n, a_n, v_p, bcx, vcx, geomState);

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

        geomState.updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes(),
                                               userDefineVel, userDefineAcc);

        setBC(userDefineDisp, userDefineVel, userDefineAcc);
        delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
      }
    }

    updateStates(&geomState, geomState, domain->solInfo().initialTime);
  }
}

int
NonLinDynamic::getInitState(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p)
{
  // initialize state with IDISP/IDISP6/IVEL or RESTART
  domain->initDispVeloc(d_n, v_n, a_n, v_p);

  // apply rbmfilter projection
  if(sinfo.filterFlags || sinfo.hzemFilterFlag) {
    project(d_n);
    project(v_n);
    project(a_n);
    project(v_p);
  }

  updateUserSuppliedFunction(d_n, v_n, a_n, v_p, domain->solInfo().initialTime);

  int aeroAlg = domain->solInfo().aeroFlag; // by default, non-aeroelastic computation
  // call aeroPreProcess (note: for restarts the initial state is now sent to the fluid in readRestartFile, not aeroPreProcess)
  if(aeroAlg >= 0)
    domain->aeroPreProcess(d_n, v_n, a_n, v_p, bcx, vcx);

  if(domain->solInfo().aeroheatFlag >= 0)
    domain->aeroHeatPreProcess(d_n, v_n, v_p, bcx);

  if(domain->solInfo().thermoeFlag >= 0)
    domain->thermoePreProcess();

  if(domain->solInfo().thermohFlag >= 0)
    domain->thermohPreProcess(d_n, v_n, v_p, bcx);

  return aeroAlg;
}

void
NonLinDynamic::updateUserSuppliedFunction(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p, double initialTime)
{
  // if we have a user supplied function, give it the initial state at the sensors
  // .. first update bcx, vcx, acx in case any of the sensors have prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      double *userDefineAcc = new double[claw->numUserDisp];
      userSupFunc->usd_disp(initialTime, userDefineDisp, userDefineVel, userDefineAcc);
      setBC(userDefineDisp, userDefineVel, userDefineAcc);
      delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
    }
    if(claw->numSensor) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel = new double[claw->numSensor];
      double *ctracc = new double[claw->numSensor];
      extractControlData(d_n, v_n, a_n, ctrdisp, ctrvel, ctracc); 
      userSupFunc->init(ctrdisp, ctrvel, ctracc);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc;
    }
  }
}

void
NonLinDynamic::extractControlData(Vector& d_n, Vector& v_n, Vector& a_n,
                                  double *ctrdsp, double *ctrvel, double *ctracc)
{
  // get SENSOR states
  DofSetArray *cdsa = domain->getCDSA();
  DofSetArray *dsa = domain->getDSA();

  for(int i = 0; i < claw->numSensor; ++i) {
    int dof = cdsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
    if(dof >= 0) { // free
      ctrdsp[i] = d_n[dof];
      ctrvel[i] = v_n[dof];
      ctracc[i] = a_n[dof];
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

void
NonLinDynamic::extractControlDisp(GeomState *geomState, double *ctrdsp)  
{
  CoordSet &nodes = domain->getNodes();
  for(int i = 0; i < claw->numSensor; ++i) {
    switch(claw->sensor[i].dofnum) {
      case 0:
        ctrdsp[i] = (*geomState)[claw->sensor[i].nnum].x - nodes[claw->sensor[i].nnum]->x;
        break;
      case 1:
        ctrdsp[i] = (*geomState)[claw->sensor[i].nnum].y - nodes[claw->sensor[i].nnum]->y;
        break;
      case 2:
        ctrdsp[i] = (*geomState)[claw->sensor[i].nnum].z - nodes[claw->sensor[i].nnum]->z;
        break;
      default:
        fprintf(stderr, "ERROR: Sensor dof %d not available in NonLinDynamic::extractControlDisp\n",claw->sensor[i].dofnum+1);
    }
  }
}

/*
void
NonLinDynamic::addCtrl(Vector& f, double *ctrfrc)
{
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numActuator; ++i) {
    int dof = cdsa->locate(claw->actuator[i].nnum,1 << claw->actuator[i].dofnum);
    if(dof >= 0) f[dof] += ctrfrc[i];
  }
}

void
NonLinDynamic::addUserForce(Vector& f, double *userDefinedForce)
{
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numUserForce; ++i) {
    int dof = cdsa->locate(claw->userForce[i].nnum,1<<claw->userForce[i].dofnum);
    if(dof >= 0) f[dof] += userDefinedForce[i];
  }
}
*/

void
NonLinDynamic::getConstForce(Vector& constantForce)
{
  domain->computeConstantForce(constantForce, Kuc);
}

void
NonLinDynamic::computeTimeInfo()
{
  // Time integration information

  // Get total time and time step size and store them 
  totalTime = domain->solInfo().tmax;
  dt0       = domain->solInfo().getTimeStep();

  // Compute maximum number of steps
  maxStep = (int) ( (totalTime+0.49*dt0)/dt0 );

  // Compute time remainder
  double remainder = totalTime - maxStep*dt0;
  if(std::abs(remainder)>0.01*dt0){
    domain->solInfo().tmax = totalTime = maxStep*dt0; 
    filePrint(stderr, " Warning: Total time is being changed to : %e\n", totalTime);
  }

  // set half time step size in user defined functions 
  if(userSupFunc)
    userSupFunc->setDt(dt0/2);
}

void
NonLinDynamic::updateContactSurfaces(GeomState& geomState, GeomState *refState)
{
  clean();
  domain->UpdateSurfaces(MortarHandler::CTC, &geomState);

  if(!domain->solInfo().trivial_detection) 
    domain->PerformStaticContactSearch(MortarHandler::CTC);

  domain->deleteSomeLMPCs(mpc::ContactSurfaces);

  if(domain->solInfo().trivial_detection) {
    double threshold = -1.0;
    Elemset &eset = geoSource->getPackedEsetConstraintElementIeq();
     for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) {
       Element *ele = eset[i];
       static_cast<MpcElement*>(ele)->update(refState, geomState, domain->getNodes(), 0.0);
       int n = ele->getNumMPCs();
       LMPCons **l = ele->getMPCs();
       for(int j = 0; j < n; ++j) {
//         if(l[j]->rhs < threshold) {
           l[j]->type = mpc::ContactSurfaces;
           domain->addLMPC(l[j]);
 //        }
       }
       delete [] l;
     }
  } else {
    domain->ExpComputeMortarLMPC(MortarHandler::CTC);
  }

  domain->UpdateContactSurfaceElements(&geomState);
  factor = false;
  //NonLinDynamic::preProcess();
  this->preProcess();
  geomState.resizeLocAndFlag(*domain->getCDSA());
  if(refState) refState->resizeLocAndFlag(*domain->getCDSA());
}

void
NonLinDynamic::updateStates(GeomState *refState, GeomState& geomState, double time)
{
  if(domain->solInfo().piecewise_contact) updateCS = true;  
  domain->updateStates(refState, geomState, allCorot, time);
}

double
NonLinDynamic::getStiffAndForce(GeomState& geomState, Vector& residual,
                                Vector& elementInternalForce, double t, GeomState *refState, bool forceOnly)
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

      geomState.updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes(),
                                             userDefineVel, userDefineAcc);

      setBC(userDefineDisp, userDefineVel, userDefineAcc);
      delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
    }
  }

  if(t != domain->solInfo().initialTime) {
    if(domain->GetnContactSurfacePairs()) {
      if(!domain->solInfo().piecewise_contact || updateCS) {
        updateContactSurfaces(geomState, refState);
        updateCS = false;
      }
      residual.conservativeResize(domain->getCDSA()->size());
      elementInternalForce.resize(domain->maxNumDOF());
      localTemp.resize(domain->getCDSA()->size());
    }
  }
  getStiffAndForceFromDomain(geomState, elementInternalForce, allCorot, kelArray, residual, 1.0, t, refState, melArray, forceOnly);

  times->buildStiffAndForce += getTime();
 
  // return residual force norm
  return residual.norm();
}

void
NonLinDynamic::getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                          Corotator **allCorot, FullSquareMatrix *kelArray,
                                          Vector &residual, double lambda, double time, GeomState *refState,
                                          FullSquareMatrix *melArray, bool forceOnly) {
  if(forceOnly) {
    domain->getInternalForce(geomState, elementInternalForce, allCorot, kelArray, residual, lambda, time, refState,
                             (Vector*) NULL, (domain->solInfo().quasistatic ? (FullSquareMatrix*) NULL : melArray), celArray);
  }
  else {
    domain->getStiffAndForce(geomState, elementInternalForce, allCorot, kelArray, residual, lambda, time, refState,
                             (Vector*) NULL, melArray, celArray);
  }
}

int
NonLinDynamic::getNumStages()
{
  return int(0.2+ domain->solInfo().getNLInfo().maxLambda / domain->solInfo().getNLInfo().dlambda);
}

int
NonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector& dv, 
                                double time)
{
  if(domain->GetnContactSurfacePairs()) dv.conservativeResize(domain->getCDSA()->size());

#ifdef PRINT_FORCENORMS
  ConstrainedDSA *cdsa = domain->getCDSA();
  double momenNorm = 0.0;
  double forceNorm = 0.0;
  int i;
  for(i=0; i<domain->numNodes(); ++i) {
    if(domain->solInfo().order == 1) {
      int dof = cdsa->locate(i, DofSet::Temp);
      if(dof >= 0)
        forceNorm += (residual[dof]*residual[dof]);
    }
    else {
      int dof = cdsa->locate(i, DofSet::Xdisp);
      if(dof >= 0)
        forceNorm += (residual[dof]*residual[dof]);
      dof = cdsa->locate(i, DofSet::Ydisp);
      if(dof >= 0)
        forceNorm += (residual[dof]*residual[dof]);
      dof = cdsa->locate(i, DofSet::Zdisp);
      if(dof >= 0)
        forceNorm += (residual[dof]*residual[dof]);
      dof = cdsa->locate(i, DofSet::Xrot);
      if(dof >= 0)
        momenNorm += (residual[dof]*residual[dof]);
      dof = cdsa->locate(i, DofSet::Yrot);
      if(dof >= 0)
        momenNorm += (residual[dof]*residual[dof]); 
      dof = cdsa->locate(i, DofSet::Zrot);
      if(dof >= 0)
        momenNorm += (residual[dof]*residual[dof]); 
    }
  }
  if(iteration == 0) {
    firstForceNorm = (forceNorm == 0.0) ? 1.0 : sqrt(forceNorm);
    firstMomenNorm = (momenNorm == 0.0) ? 1.0 : sqrt(momenNorm);
  }
  fprintf(stderr,"===> time %f force Norm %e Moment Norm %e\n",
          time,
          sqrt(forceNorm)/firstForceNorm,
          sqrt(momenNorm)/firstMomenNorm);
#endif

  // Note when useTolInc is false, this function is called before normDv is calculated
  bool useTolInc = (domain->solInfo().getNLInfo().tolInc != std::numeric_limits<double>::infinity() ||
                    domain->solInfo().getNLInfo().absTolInc != std::numeric_limits<double>::infinity());

  double normDv     = dv.norm();
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
  filePrint(res, "%d %19.12e %e %e %e %e\n", totIter, time, normRes, relRes, normDv, relDv);
  fflush(res);

  // Store residual norm and dv norm for output
  times->norms[numSystems].normDv      = normDv;
  times->norms[numSystems].relativeDv  = relDv;
  times->norms[numSystems].normRes     = normRes;
  times->norms[numSystems].relativeRes = relRes;
  times->numSystems = numSystems;

  numSystems += 1;

  return converged;
}

GeomState*
NonLinDynamic::createGeomState()
{
  if(domain->solInfo().soltyp == 2)
    return new TemperatureState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes());
  else
    return new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes(), &domain->getElementSet(),
                         domain->getNodalTemperatures());
}

GeomState*
NonLinDynamic::copyGeomState(GeomState* geomState)
{
  if(domain->solInfo().soltyp == 2)
    return new TemperatureState(* (TemperatureState *) geomState);
  else
    return new GeomState(*geomState);
}

// Rebuild dynamic mass matrix
void
NonLinDynamic::reBuild(GeomState& geomState, int iteration, double localDelta, double t)
{
 int step = (t+localDelta)/(2*localDelta);

 // note: localDelta = deltat/2
 times->rebuild -= getTime();

 // Rebuild every updateK iterations
 if((iteration % domain->solInfo().getNLInfo().updateK == 0 && (step-1) % domain->solInfo().getNLInfo().stepUpdateK == 0)
    || (t == domain->solInfo().initialTime))  {

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
     if(domain->solInfo().quasistatic) {
       Mcoef = 0;
       Kcoef = 1;
       Ccoef = 0;
     }
     else {
       Mcoef = (domain->solInfo().order == 1) ? 1 : (1-alpham)/(1-alphaf);
       Kcoef = (domain->solInfo().order == 1) ? dt*gamma : dt*dt*beta;
       Ccoef = (domain->solInfo().order == 1) ? 0 : dt*gamma;
     }
   }

   if(domain->solInfo().mpcDirect != 0) {
     if(solver) delete solver;
     if(prec) delete prec;
     if(Kuc) delete Kuc;
     if(M) delete M; if(Muc) delete Muc; if(Mcc) delete Mcc;
     if(C) delete C; if(Cuc) delete Cuc; if(Ccc) delete Ccc;
     factor = true;
     preProcess(Kcoef, Mcoef, Ccoef);
   }
   else {
     if(!(domain->solInfo().DEIMPodRom || domain->solInfo().UDEIMPodRom))
       spm->zeroAll();
     AllOps<double> ops;
     if(Kuc) { Kuc->zeroAll(); ops.Kuc = Kuc; }
     if(spp) { spp->zeroAll(); ops.spp = spp; }
#ifdef USE_EIGEN3
     if(domain->solInfo().printMatLab) {
       ops.K = domain->constructEiSparse<double>(); ops.K->zeroAll();
       ops.M = domain->constructEiSparse<double>(); ops.M->zeroAll();
     }
#endif
     if(domain->solInfo().galerkinPodRom && (domain->solInfo().useMassNormalizedBasis || domain->solInfo().modalDIMASS)) {
       domain->makeSparseOps<double>(ops, Kcoef, 0.0, Ccoef, spm, kelArray, melArray, celArray);
       if(domain->solInfo().useMassNormalizedBasis)
         dynamic_cast<Rom::PodProjectionSolver*>(solver)->addReducedMass(Mcoef);
       else {
#ifdef USE_EIGEN3
         dynamic_cast<Rom::PodProjectionSolver*>(solver)->addToReducedMatrix(VtMV, Mcoef);
#endif
       }
     }
     else {
       domain->makeSparseOps<double>(ops, Kcoef, Mcoef, Ccoef, spm, kelArray, melArray, celArray);
     }
#ifdef USE_EIGEN3
     if(domain->solInfo().printMatLab) {
       ops.M->printSparse(std::string(domain->solInfo().printMatLabFile)+".mass");
       ops.K->printSparse(std::string(domain->solInfo().printMatLabFile)+".stiffness");
       if(domain->solInfo().printMatLabExit) throw std::exception();
     }
#endif
     if(!verboseFlag) solver->setPrintNullity(false);
     domain->getTimers().factor -= getTime();
     solver->factor();
     if(prec) prec->factor();
     domain->getTimers().factor += getTime();
   }

 }
 times->rebuild += getTime();
}

int
NonLinDynamic::solVecInfo() const
{
  return domain->numUncon();
}

int
NonLinDynamic::sysVecInfo()
{
  return domain->numdof();
}

int
NonLinDynamic::elemVecInfo()
{
 return domain->maxNumDOF();
}

double
NonLinDynamic::getTolerance()
{
 return std::max(tolerance*firstRes, domain->solInfo().getNLInfo().absTolRes);
}

int
NonLinDynamic::getMaxit()
{
 return domain->solInfo().getNLInfo().maxiter;
}

double
NonLinDynamic::getDeltaLambda()
{
 return domain->solInfo().getNLInfo().dlambda;
}

bool
NonLinDynamic::getZeroRot() const {
  return domain->solInfo().zeroRot;
}

void
NonLinDynamic::getExternalForce(Vector& rhs, Vector& constantForce, int tIndex, double t, 
                                GeomState* geomState, Vector&, 
                                Vector &aeroForce, double localDelta)
{
  // ... BUILD THE EXTERNAL FORCE at t_{n+1-alphaf}
  times->formRhs -= getTime();

  // update USDF and ACTUATORS
  if(claw && userSupFunc) {
    if(claw->numUserForce > 0) {
      double *userDefinedForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(t, userDefinedForce);
      domain->updateUsdfInNbc(userDefinedForce);
      delete [] userDefinedForce;
    }

    if(claw->numActuator > 0) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel  = new double[claw->numSensor];
      double *ctracc  = new double[claw->numSensor];
      double *ctrfrc  = new double[claw->numActuator];

      for(int j = 0; j < claw->numSensor; j++) ctrvel[j] = ctracc[j] = 0.0; // TODO f(v,a) currently not supported

      // KHP: we need the state of the control sensors to pass to
      //      the user supplied control function
      extractControlDisp(geomState, ctrdisp);

      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t);
      domain->updateActuatorsInNbc(ctrfrc);

      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }

  // update THERMOE 
  if(domain->solInfo().thermoeFlag >= 0 && tIndex >= 0) {
    domain->thermoeComm();
    geomState->setNodalTemperatures(domain->getNodalTemperatures());
  }

  // add f(t) to constantForce (not including follower forces)
  double beta, gamma, alphaf, alpham;
  getNewmarkParameters(beta, gamma, alphaf, alpham);
  double dt = 2*localDelta;
  double t0 = domain->solInfo().initialTime;
  double tm = (t == t0) ? t0 : t + dt*(alphaf-alpham);
  domain->computeExtForce4(rhs, constantForce, t, Kuc, userSupFunc, Cuc, tm, Muc);

  // add aeroelastic forces from fluid dynamics code
  if(domain->solInfo().aeroFlag >= 0 && tIndex >= 0) {
    domain->buildAeroelasticForce(aeroForce, *prevFrc, tIndex, t, gamma, alphaf);
    rhs += aeroForce;
  }

  // add aerothermal fluxes from fluid dynamics code
  if(domain->solInfo().aeroheatFlag >= 0 && tIndex >= 0)
    domain->buildAeroheatFlux(rhs, prevFrc->lastFluidLoad, tIndex, t);

  // apply rbmfilter projection
  if(sinfo.filterFlags || sinfo.hzemFilterFlag)
    trProject(rhs);

  times->formRhs += getTime();
}

void
NonLinDynamic::getIncDisplacement(GeomState *geomState, Vector &du, GeomState *refState,
                                  bool zeroRot)
{
  geomState->get_inc_displacement(du, *refState, zeroRot);
}

void
NonLinDynamic::formRHSinitializer(Vector &fext, Vector &velocity, Vector &elementInternalForce, GeomState &geomState, Vector &rhs, GeomState *refState)
{
  // rhs = (fext - fint - Cv)
  rhs = fext;
  getStiffAndForceFromDomain(geomState, elementInternalForce, allCorot, kelArray, rhs, 1.0, domain->solInfo().initialTime, refState,
                             melArray, false);
  if(domain->solInfo().order == 2 && C) {
    C->mult(velocity, localTemp);
    rhs.linC(rhs, -1.0, localTemp);
  }
}

void
NonLinDynamic::formRHSpredictor(Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs, GeomState &geomState, 
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
      geomState.updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes(),
                                             userDefineDisp, userDefineAcc);
      setBC(userDefineDisp, userDefineVel, userDefineAcc);

      // get delta disps
      for(int j = 0; j < claw->numUserDisp; j++)
        userDefineDisp[j] -= userDefineDispLast[j];
 
      // update force residual with KUC
      if(Kuc)
        Kuc->transposeMultSubtractClaw(userDefineDisp, residual.data(), claw->numUserDisp, clawDofs);

      delete [] userDefineDisp; delete [] userDefineDispLast; delete [] userDefineVel; delete [] userDefineAcc;
    }
  }

  if(domain->solInfo().order == 1) 
    rhs.linC(localDelta, residual);
  else {
    double beta, gamma, alphaf, alpham, dt = 2*localDelta;
    getNewmarkParameters(beta, gamma, alphaf, alpham);
    // rhs = dt*dt*beta*residual + (dt*(1-alpham)*M - dt*dt*(beta-(1-alphaf)*gamma)*C)*velocity
    //       + (dt*dt*((1-alpham)/2-beta)*M - dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2*C)*acceleration
    localTemp.linC(dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
    M->mult(localTemp, rhs);
    if(C) {
      localTemp.linC(-dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(localTemp.data(), rhs.data());
    }
    rhs.linAdd(dt*dt*beta, residual);
  }

  times->predictorTime += getTime();
}

double
NonLinDynamic::formRHScorrector(Vector &inc_displacement, Vector &velocity, Vector &acceleration,
                                Vector &residual, Vector &rhs, GeomState *geomState, double localDelta)
{
  times->correctorTime -= getTime();
  if(domain->GetnContactSurfacePairs()) {
    velocity.conservativeResize(solVecInfo());
    acceleration.conservativeResize(solVecInfo());
    rhs.resize(solVecInfo());
    inc_displacement.conservativeResize(solVecInfo());
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
      localTemp.linC(-(1-alpham)/(1-alphaf), inc_displacement, dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
      M->mult(localTemp, rhs);
    }
    if(C) {
      localTemp.linC(-dt*gamma, inc_displacement, -dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(localTemp.data(), rhs.data());
    }
    rhs.linAdd(dt*dt*beta, residual);
  }
  times->correctorTime += getTime();
  return rhs.norm();
}

void
NonLinDynamic::processLastOutput()  {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  domain->solInfo().lastIt = true;
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

void
NonLinDynamic::preProcessSA()
{
  domain->buildPreSensitivities<double>(*allSens, bcx);
}

void
NonLinDynamic::preProcess(double Kcoef, double Mcoef, double Ccoef)
{
 // Allocate space for the Static Timers
 if(!times) times = new StaticTimers;

 this->openResidualFile();

 // Set the nonlinear tolerance
 tolerance = domain->solInfo().getNLInfo().tolRes;

 // Makes renumbering, connectivities and dofsets
 startTimerMemory(times->preProcess, times->memoryPreProcess);
 domain->preProcessing();

 int numdof = domain->numdof();

 int *bc = (int *) dbg_alloca(sizeof(int)*numdof);
 if(!bcx) bcx      = new double[numdof];

 // vcx stores the prescribed velocities, which are associated
 // with prescribed displacements. If a user defined displacement
 // is used, then vcx will contain the user defined velocities also.
 if(!vcx) vcx      = new double[numdof];
 if(!acx) acx      = new double[numdof];

 int i;
 for(i=0; i<numdof; ++i)
   acx[i] = vcx[i] = 0.0;

 BCond* iVel = domain->getInitVelocity();

 // Make the boundary conditions info
 domain->make_bc( bc, bcx );
 if(!reactions) reactions = new Vector(domain->nDirichlet());

 // Now, call make_constrainedDSA(bc) to 
 // built c_dsa that will incorporate all 
 // the boundary conditions info
 domain->make_constrainedDSA(bc);

 domain->makeAllDOFs();

 AllOps<double> allOps;

 allOps.M = domain->constructDBSparseMatrix<double>();
 allOps.Muc = domain->constructCuCSparse<double>();
 allOps.Mcc = domain->constructCCSparse<double>();

 if(domain->solInfo().alphaDamp != 0.0 || domain->solInfo().betaDamp != 0.0 || domain->getElementSet().hasDamping()) {
   allOps.C = domain->constructDBSparseMatrix<double>();
   allOps.Cuc = domain->constructCuCSparse<double>();
   allOps.Ccc = domain->constructCCSparse<double>();
 }

 if(domain->solInfo().getNLInfo().linearelastic) {
   allOps.K = domain->constructDBSparseMatrix<double>();
 }
 allOps.Kuc = domain->constructCuCSparse<double>();

 Rbm *rigidBodyModes = 0;

 int useGrbm = domain->solInfo().rbmflg;
 int useHzem = domain->solInfo().hzemFlag;

 if (useGrbm || sinfo.filterFlags)
   rigidBodyModes = domain->constructRbm();
 else if(useHzem || sinfo.hzemFilterFlag)
   rigidBodyModes = domain->constructHzem();

 // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
 if(!kelArray) {
   if(allOps.C) domain->createKelArray(kelArray, melArray, celArray);
   else domain->createKelArray(kelArray, melArray);
 }

 domain->buildOps<double>(allOps, Kcoef, Mcoef, Ccoef, (Rbm *) NULL, kelArray,
                          melArray, celArray, factorWhenBuilding()); // don't use Rbm's to factor in dynamics

 K      = allOps.K;
 Kuc    = allOps.Kuc;
 M      = allOps.M;
 C      = allOps.C;
 solver = allOps.sysSolver;
 spm    = allOps.spm;
 prec   = allOps.prec;
 spp    = allOps.spp;
 Muc    = allOps.Muc;
 Mcc    = allOps.Mcc;
 Cuc    = allOps.Cuc;
 Ccc    = allOps.Ccc;

 if(sinfo.filterFlags || sinfo.hzemFilterFlag)
   projector_prep(rigidBodyModes, M);

 if(rigidBodyModes) delete rigidBodyModes;

 if(!allCorot) {
   // ... ALLOCATE MEMORY FOR THE ARRAY OF COROTATORS
   allCorot = new Corotator *[domain->numElements()];

   // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
   domain->createCorotators(allCorot);
 }

 // Look if there is a user supplied routine for control
 claw = geoSource->getControlLaw();

 // create list of usdd node dofs mapped to cdsa dof numbers
 if(claw) {
   int nClaw = claw->numUserDisp;
   clawDofs = new int[nClaw];
   for (int j = 0; j < nClaw; ++j) {
     int dd = domain->getDSA()->locate(claw->userDisp[j].nnum, (1 << claw->userDisp[j].dofnum));
     clawDofs[j] = domain->getCDSA()->invRCN(dd);
   }
 }

 // Check to see if there is a user supplied function
 // for displacements, forces or control law
 userSupFunc = domain->getUserSuppliedFunction();
 if(!prevFrc) prevFrc = new PrevFrc(domain->numUncon());

 localTemp.initialize(domain->numUncon());

 stopTimerMemory(times->preProcess, times->memoryPreProcess);
}

void NonLinDynamic::openResidualFile()
{
  if(!res) {
#ifdef ANDROID
    res = fopen("/sdcard/residuals", "w");
#else
    res = fopen("residuals", "w");
#endif
    totIter = 0;
    fprintf(res,"Iteration Time           Residual\trel. res\tdv\t rel. dv\n");
  }
}

Solver *
NonLinDynamic::getSolver()
{
  return solver;
}

SDDynamPostProcessor* 
NonLinDynamic::getPostProcessor()
{
  return new SDDynamPostProcessor(domain, bcx, vcx, acx, times);
}

void
NonLinDynamic::printTimers(double timeLoop)
{
  long memoryUsed = solver->size();
  double solveTime = solver->getSolutionTime();

  times->printStaticTimers( solveTime, memoryUsed, domain, timeLoop );

  if(domain->solInfo().massFlag) {
    double mass = domain->computeStructureMass();
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," ... Structure mass = %e  ...\n",mass);
    filePrint(stderr," --------------------------------------\n");
  }
}

int
NonLinDynamic::getAeroAlg()
{
  return domain->solInfo().aeroFlag;
}

int
NonLinDynamic::getThermoeFlag()
{
  return domain->solInfo().thermoeFlag;
}

int
NonLinDynamic::getThermohFlag()
{
  return domain->solInfo().thermohFlag;
}

int
NonLinDynamic::getAeroheatFlag()
{
  return domain->solInfo().aeroheatFlag;
}

void 
NonLinDynamic::dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState,
                                Vector& velocity, Vector& bkVelocity,
                                Vector& vp, Vector& bkVp, int step, int parity, 
                                int aeroAlg, double time)
{
  times->output -= getTime();

  if(domain->solInfo().aeroFlag >= 0 && !domain->solInfo().lastIt) {
    domain->getTimers().sendFluidTime -= getTime();

    // update geomState bcx/vcx/acx for time dependent prescribed displacements and their time derivatives
    ControlLawInfo *claw = geoSource->getControlLaw();
    ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
    if(claw && claw->numUserDisp) {
      // Note: the approprate value of "time" passed into this function should be t^{n+½} for A6 and t^{n+1}
      // otherwise, where t^n denotes the time at the end of the current structure timestep. Note that the
      // predictor in FlExchanger::sendDisplacements is not applied to prescribed displacements; we directly
      // compute here the desired values of the prescribed displacements/velocities rather than predicting them.
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel  = new double[claw->numUserDisp];
      double *userDefineAcc  = new double[claw->numUserDisp];
      for(int i = 0; i < claw->numUserDisp; ++i) {
        userDefineVel[i] = 0;
        userDefineAcc[i] = 0;
      }
      userSupFunc->usd_disp(time, userDefineDisp, userDefineVel, userDefineAcc);
      setBC(userDefineDisp, userDefineVel, userDefineAcc);
      geomState->updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes(), userDefineVel, userDefineAcc);
      delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
    }

    // Make d_n_aero from geomState
    Vector d_n(domain->numUncon(), 0.0);
    geomState->get_tot_displacement(d_n);

    ConstrainedDSA *c_dsa = domain->getCDSA();
    DofSetArray *dsa = domain->getDSA();

    if(!parity && aeroAlg==5) { 
      Vector d_n2(domain->numUncon(), 0.0);
      bkGeomState->get_tot_displacement(d_n2);
      d_n.linC(0.5,d_n,0.5,d_n2);
      velocity.linC(0.5,velocity,0.5,bkVelocity);
      vp.linC(0.5,vp,0.5,bkVp);
    }
    Vector a_n( domain->numUncon(), 0.0 );

    State state( c_dsa, dsa, bcx, vcx, d_n, velocity, a_n, vp );

    domain->getFileExchanger()->sendDisplacements(state, -1, geomState);
    if(verboseFlag) filePrint(stderr, " ... [E] Sent displacements         ...\n");

    domain->getTimers().sendFluidTime += getTime();
  }

  if(domain->solInfo().aeroheatFlag >= 0) {
    domain->getTimers().sendFluidTime -= getTime();

    // Make d_n_aero from geomState
    ConstrainedDSA *c_dsa = domain->getCDSA();
    DofSetArray *dsa = domain->getDSA();

    // Note that d_n and a_n are vectors being allocated and de-allocated at
    // each time-step being executed.

    Vector d_n( domain->numUncon(), 0.0 );

    CoordSet &nodes = domain->getNodes();

    for(int i=0; i<domain->numNodes(); ++i) {

      int tloc  = c_dsa->locate(i, DofSet::Temp );
      int tloc1 =   dsa->locate(i, DofSet::Temp );

      if(tloc >= 0)
        d_n[tloc]  = (*geomState)[i].x;
      else if (tloc1 >= 0)
        bcx[tloc1] = (*geomState)[i].x;
    }

    State tempState(c_dsa, dsa, bcx, d_n, velocity, vp);
    domain->getFileExchanger()->sendTemperature(tempState);
    if(verboseFlag) filePrint(stderr," ... [T] Sent temperatures ...\n");

    domain->getTimers().sendFluidTime += getTime();
  }

  if(domain->solInfo().thermohFlag >= 0) {
    /* we have to send the vector of temperatures in NODAL order, not
       in DOF order (in which is d_n)! */

    Vector tempsent(domain->numNodes());

    for(int iNode=0; iNode<domain->numNodes(); ++iNode)
      tempsent[iNode] = (*geomState)[iNode].x;

    domain->getFileExchanger()->sendStrucTemp(tempsent);
    if(verboseFlag) filePrint(stderr," ... [T] Sent temperatures ...\n");
  }

  times->output += getTime();
}

void
NonLinDynamic::dynamOutput(GeomState* geomState, Vector& velocity,
                           Vector& vp, double time, int step, Vector& force, 
                           Vector &aeroF, Vector &acceleration, GeomState *refState) const
{
  times->output -= getTime();

  // update geomState bcx/vcx/acx for time dependent prescribed displacements and their time derivatives prior to output
  ControlLawInfo *claw = geoSource->getControlLaw();
  ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
  if(claw && claw->numUserDisp) {
    double *userDefineDisp = new double[claw->numUserDisp];
    double *userDefineVel  = new double[claw->numUserDisp];
    double *userDefineAcc  = new double[claw->numUserDisp];
    for(int i = 0; i < claw->numUserDisp; ++i) {
      userDefineVel[i] = 0;
      userDefineAcc[i] = 0;
    }
    userSupFunc->usd_disp(time, userDefineDisp, userDefineVel, userDefineAcc);
    DofSetArray *dsa = domain->getDSA();
    for(int i = 0; i < claw->numUserDisp; ++i) {
      int dof = dsa->locate(claw->userDisp[i].nnum,1 << claw->userDisp[i].dofnum);
      if(dof >= 0) {
        bcx[dof] = userDefineDisp[i];  // actually, prescribed displacements are output from the geomState, not bcx
        vcx[dof] = userDefineVel[i];
        acx[dof] = userDefineAcc[i];
      }
    }
    geomState->updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes(), userDefineVel, userDefineAcc);
    delete [] userDefineDisp; delete [] userDefineVel; delete [] userDefineAcc;
  }

  if(domain->reactionsReqd(time, step+1)) {
    domain->computeReactionForce(*reactions, geomState, allCorot, kelArray, time, refState, velocity,
                                 acceleration, vcx, acx, Cuc, Ccc, Muc, Mcc); 
  }

  domain->postProcessing(geomState, force, aeroF, time, (step+1), velocity.data(), vcx,
                         allCorot, acceleration.data(), acx, refState, reactions, M, C);
  times->output += getTime();
}

void
NonLinDynamic::updatePrescribedDisplacement(GeomState *geomState)
{
 if(domain->solInfo().initialTime == 0.0) {
   // Measure time necessary to update the prescribed displacments
   times->timePresc -= getTime();

   // note 2: "if both IDISP and IDISP6 are present in the input file, FEM selects IDISP6 to construct the geometric stiffness"
   if((domain->numInitDisp() > 0) && (domain->numInitDisp6() == 0))
     geomState->updatePrescribedDisplacement(domain->getInitDisp(), domain->numInitDisp(), domain->getNodes());
   
   if(domain->numInitDisp6() > 0) 
     geomState->updatePrescribedDisplacement(domain->getInitDisp6(), domain->numInitDisp6(), domain->getNodes());
   
   if(domain->nDirichlet() > 0)
     geomState->updatePrescribedDisplacement(domain->getDBC(), domain->nDirichlet(), domain->getNodes()); 

   times->timePresc += getTime();
 }
}

void
NonLinDynamic::setBC(double *userDefineDisplacement, double *userDefineVel, double *userDefineAcc)
{
  DofSetArray *dsa = domain->getDSA();
  for(int i = 0; i < claw->numUserDisp; ++i) {
    int dof = dsa->locate(claw->userDisp[i].nnum, 1 << claw->userDisp[i].dofnum);
    if(dof >= 0) {
      bcx[dof] = userDefineDisplacement[i];
      vcx[dof] = userDefineVel[i];
      acx[dof] = userDefineAcc[i];
    }
  }
}

void
NonLinDynamic::getInitialTime(int &initTimeIndex, double &initTime)
{
 initTimeIndex = domain->solInfo().initialTimeIndex;
 initTime      = domain->solInfo().initialTime;
}

void
NonLinDynamic::getNewmarkParameters(double &beta, double &gamma,
                                    double &alphaf, double &alpham)
{
 beta  = domain->solInfo().newmarkBeta;
 gamma = domain->solInfo().newmarkGamma;
 alphaf = domain->solInfo().newmarkAlphaF;
 alpham = domain->solInfo().newmarkAlphaM;
}

double
NonLinDynamic::getResidualNorm(const Vector &rhs, GeomState &geomState, double localDelta)
{
  double beta, gamma, alphaf, alpham, dt = 2*localDelta;
  getNewmarkParameters(beta, gamma, alphaf, alpham);
  Vector res(rhs);
  domain->applyResidualCorrection(geomState, allCorot, res, dt*dt*beta);
  return solver->getResidualNorm(res);
}

bool
NonLinDynamic::factorWhenBuilding() const { 
  return factor; //domain->solInfo().iacc_switch || domain->solInfo().mpcDirect != 0;
}

void
NonLinDynamic::initializeParameters(int step, GeomState *geomState)
{
  if(step == 1 || domain->solInfo().reinit_lm) {
    domain->initializeMultipliers(*geomState, allCorot);
  }
  domain->initializeParameters();
}

void
NonLinDynamic::updateParameters(GeomState *geomState)
{
  domain->updateMultipliers(*geomState, allCorot);
  domain->updateParameters();
}

bool
NonLinDynamic::checkConstraintViolation(double &err, GeomState *gs)
{
  err = domain->getError(allCorot, *gs);
  return (err <= domain->solInfo().penalty_tol);
}

LinesearchInfo&
NonLinDynamic::linesearch()
{
 return domain->solInfo().getNLInfo().linesearch;
}

bool
NonLinDynamic::getResizeFlag()
{
  return (domain->GetnContactSurfacePairs() > 0);
}

void
NonLinDynamic::postProcessSA(Vector &sol)
{
  domain->buildPostSensitivities<double>(solver, spm, K, *allSens, &sol, bcx, false);
}

void
NonLinDynamic::postProcessNLSA(GeomState *refState, GeomState *geomState)
{
  domain->buildNLPostSensitivities<double>(solver, *allSens, refState, geomState, allCorot, true); 
}

SensitivityInfo*
NonLinDynamic::getSensitivityInfo() { return domain->senInfo; }

int
NonLinDynamic::getNumSensitivities() { return domain->getNumSensitivities(); }

void
NonLinDynamic::sensitivityAnalysis(GeomState *geomState, GeomState *refState)
{
#ifdef USE_EIGEN3
  filePrint(stderr, " ... starting nonlinear sensitivity analysis\n");
  SensitivityInfo *senInfo = getSensitivityInfo();
  int numShapeVars = domain->getNumShapeVars();
  int numThicknessGroups = domain->getNumThicknessGroups();
  int numStructQuantTypes = domain->getNumSensitivityQuantityTypes();
  if(domain->solInfo().sensitivity) {
    postProcessNLSA(refState, geomState);
    AllSensitivities<double> *allSens = getAllSensitivities();
    for(int isen = 0; isen < getNumSensitivities(); ++isen) {
      switch (senInfo[isen].type) {

        case SensitivityInfo::DisplacementWRTshape:
 
          filePrint(stderr, " ... numShapeVars = %d              ...\n", numShapeVars);
          if( numShapeVars > 0) {
            if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
              if(!allSens->dispWRTshape) { 
                allSens->dispWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
                Vector *rhsSen;
                for(int ishap=0; ishap< numShapeVars; ++ishap) {
                  allSens->dispWRTshape[ishap] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(domain->numUncon(),1);
                  rhsSen = new Vector( solVecInfo() );
                  rhsSen->copy(allSens->linearstaticWRTshape[ishap]->data());
                  (*rhsSen) *= -1;
                  getSolver()->reSolve(*rhsSen);
                  *allSens->dispWRTshape[ishap] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                                  Eigen::Dynamic> >(rhsSen->data(),
                                                                                                    domain->numUncon(),1);
                  delete rhsSen;
                }
              }
            }
          }     
          break;

        case SensitivityInfo::DisplacementWRTthickness:

          filePrint(stderr," ... numThicknessGroups = %d         ...\n", numThicknessGroups);
          if( numThicknessGroups > 0 ) {  
            if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
              if(!allSens->dispWRTthick) { 
                allSens->dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
                Vector *rhsSen;
                for(int iparam=0; iparam< numThicknessGroups; ++iparam) {
                  allSens->dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(domain->numUncon(),1);
                  rhsSen = new Vector( solVecInfo() );
                  rhsSen->copy(allSens->linearstaticWRTthick[iparam]->data());
                  (*rhsSen) *= -1; 
                  getSolver()->reSolve(*rhsSen);
                  *allSens->dispWRTthick[iparam] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                                   Eigen::Dynamic> >(rhsSen->data(),
                                                                                                     domain->numUncon(),1);
                  delete rhsSen;
                }
              }
            } else if(domain->solInfo().sensitivityMethod == SolverInfo::Adjoint) { // Adjoint
              int numDispNodes = domain->getNumDispNodes();
              int numTotalDispDofs = domain->getTotalNumDispDofs();
              if(!allSens->lambdaDisp) 
                computeLambdaDisp(numDispNodes,numTotalDispDofs);
              allSens->dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
              std::vector<DispNode> *dispNodes = domain->getDispNodes();
              for(int iparam=0; iparam < numThicknessGroups; ++iparam) {
                allSens->dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numTotalDispDofs,1);
                allSens->dispWRTthick[iparam]->setZero();
                int dispDofIndex = 0;
                for(int inode=0; inode < numDispNodes; ++inode) {
                  int numDispDofs = (*dispNodes)[inode].numdofs;
                  for(int idof=0; idof < numDispDofs; ++idof) {
                     (*allSens->dispWRTthick[iparam])(dispDofIndex,0) -= 
                        allSens->lambdaDisp[dispDofIndex]->adjoint()*(allSens->linearstaticWRTthick[iparam]->col(0));
                     dispDofIndex++;
                  }
                }
              }
            }
          }
          if (allSens->residual !=0 && allSens->dwrDisp==0) {
            int numTotalDispDofs = domain->getTotalNumDispDofs();
            allSens->dwrDisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numTotalDispDofs);
            for (int i=0; i<numTotalDispDofs; ++i) {
              (*allSens->dwrDisp)[i] = allSens->lambdaDisp[i]->dot(*(allSens->residual));
            }
          }
          break;

        case SensitivityInfo::StressVMWRTthickness:

          if( numThicknessGroups > 0 ) {
            if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
              allSens->dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
              Vector *rhsSen;
              for(int iparam=0; iparam< numThicknessGroups; ++iparam) {
                allSens->dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(domain->numUncon(),1);
                rhsSen = new Vector( solVecInfo() );
                rhsSen->copy(allSens->linearstaticWRTthick[iparam]->data());
                (*rhsSen) *= -1;
                getSolver()->reSolve(*rhsSen);
                *allSens->dispWRTthick[iparam] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                                 Eigen::Dynamic> >(rhsSen->data(),
                                                                                                   domain->numUncon(),1);
                allSens->vonMisesWRTthick->col(iparam) += *allSens->vonMisesWRTdisp * (*allSens->dispWRTthick[iparam]);
                delete rhsSen;
              }
            }
            else if(domain->solInfo().sensitivityMethod == SolverInfo::Adjoint) {
              int numStressNodes = domain->getNumStressNodes();
              if(!allSens->lambdaStressVM) computeLambdaNLStressVM(numStressNodes);
              for(int inode = 0; inode < numStressNodes; ++inode)
                for(int iparam=0; iparam < numThicknessGroups; ++iparam)
                  (*allSens->vonMisesWRTthick)(inode,iparam) -= 
                      allSens->lambdaStressVM[inode]->adjoint()*(allSens->linearstaticWRTthick[iparam]->col(0));
            }
          }
          if (allSens->residual !=0 && allSens->dwrStressVM==0) {
            int numStressNodes = domain->getNumStressNodes(); 
            allSens->dwrStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(numStressNodes);
            for (int i=0; i<numStressNodes; i++) {
              (*allSens->dwrStressVM)[i] = allSens->lambdaStressVM[i]->dot(*(allSens->residual));
            }
          }
          break;

        case SensitivityInfo::AggregatedStressVMWRTthickness:
          if( numThicknessGroups > 0 ) {
            if(domain->solInfo().sensitivityMethod != SolverInfo::Adjoint) {
              filePrint(stderr, " *** Error: aggregated stress sensitivity is not supported by direct sensitivity, exiting\n");
              exit(-1);
            }
            if(!allSens->lambdaAggregatedStressVM) {
              computeLambdaAggregatedStress();
            }
            for(int iparam=0; iparam< numThicknessGroups; ++iparam) {
              (*allSens->aggregatedVonMisesWRTthick)(iparam) -= 
                  allSens->lambdaAggregatedStressVM->adjoint()*(allSens->linearstaticWRTthick[iparam]->col(0));
            }
          }
         if (allSens->residual !=0 && allSens->dwrAggregatedStressVM == 0) {
           allSens->dwrAggregatedStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(1);
           (*allSens->dwrAggregatedStressVM)[0] = allSens->lambdaAggregatedStressVM->dot(*(allSens->residual));
          }
          break;

        default:
          break; 
      }
    }
    Vector *fullDispBuffer = new Vector(solVecInfo());
    geomState->get_tot_displacement(*fullDispBuffer);
    domain->sensitivityPostProcessing(*allSens, fullDispBuffer, bcx, geomState, refState, allCorot); 
    delete fullDispBuffer;
  }
#endif
}

void
NonLinDynamic::computeLambdaNLStressVM(int numStressNodes) 
{
  #ifdef USE_EIGEN3
  AllSensitivities<double> *allSens = getAllSensitivities();
  std::vector<int> *stressNodes = domain->getStressNodes();
  allSens->lambdaStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>*[numStressNodes];
  Vector *rhsSen;
  setUpPODSolver(OutputInfo::VMstThic);
  for(int inode = 0; inode < numStressNodes; ++inode) {
    allSens->lambdaStressVM[inode] = new Eigen::Matrix<double, Eigen::Dynamic, 1>(domain->numUncon());
    rhsSen = new Vector( solVecInfo() );
    for(int i=0; i<domain->numUncon(); ++i) (*rhsSen)[i] = (*allSens->vonMisesWRTdisp)((*stressNodes)[inode],i);
    getSolver()->reSolve(*rhsSen);
    *allSens->lambdaStressVM[inode] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                      Eigen::Dynamic> >(rhsSen->data(),domain->numUncon(),1);
    delete rhsSen;
  }
  #endif
}

void 
NonLinDynamic::computeLambdaAggregatedStress()
{
  #ifdef USE_EIGEN3
  AllSensitivities<double> *allSens = getAllSensitivities();
  allSens->lambdaAggregatedStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(domain->numUncon());
  Vector *rhsSen = new Vector( solVecInfo() );
  setUpPODSolver(OutputInfo::AGstThic);

  for(int i=0; i<domain->numUncon(); ++i) { (*rhsSen)[i] = (*allSens->aggregatedVonMisesWRTdisp)(i); }
  getSolver()->reSolve(*rhsSen);
  *allSens->lambdaAggregatedStressVM = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                       Eigen::Dynamic> >(rhsSen->data(),
                                                                                         domain->numUncon(),1);
  delete rhsSen;
  #endif
}

void
NonLinDynamic::computeLambdaDisp(int numDispNodes, int numTotalDispDofs)
{
  #ifdef USE_EIGEN3
  AllSensitivities<double> *allSens = getAllSensitivities();
  std::vector<DispNode> *dispNodes = domain->getDispNodes();
  allSens->lambdaDisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>*[numTotalDispDofs];
  int dispDofIndex = 0;
  Vector *rhsSen;

  setUpPODSolver(OutputInfo::DispThic);
  for(int inode = 0; inode < numDispNodes; ++inode) {
    int node = (*dispNodes)[inode].nodeID, loc;
    int numDispDofs = (*dispNodes)[inode].numdofs;
    for(int idof = 0; idof<numDispDofs; ++idof) {
      allSens->lambdaDisp[dispDofIndex] = new Eigen::Matrix<double, Eigen::Dynamic, 1>(domain->numUncon());
      allSens->lambdaDisp[dispDofIndex]->setZero();
      rhsSen = new Vector( solVecInfo() );
      Vector lambdadisp(domain->numUncon(),0.0);
      int dof = (*dispNodes)[inode].dofs[idof];
      loc = domain->returnLocalDofNum(node, dof);
      if (loc >= 0) { rhsSen->zero(); (*rhsSen)[loc] = 1.0; }
      else continue;
      getSolver()->reSolve(*rhsSen);
      *allSens->lambdaDisp[dispDofIndex] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                           Eigen::Dynamic> >(rhsSen->data(),
                                                                                             domain->numUncon(),1);
      delete rhsSen;
      dispDofIndex++;
    }
  }
  #endif
}

void 
NonLinDynamic::setUpPODSolver(OutputInfo::Type type) {
  if(!domain->solInfo().readInAdjointROB.empty()) {
    Rom::PodProjectionSolver* podSolver = dynamic_cast<Rom::PodProjectionSolver*>(solver);
    if(podSolver) {
      std::map<OutputInfo::Type,int>::iterator it = domain->solInfo().adjointMap.find(type);
      if(it != domain->solInfo().adjointMap.end()) {
        int adjointBasisId = it->second;
        int blockCols = domain->solInfo().maxSizeAdjointBasis[adjointBasisId];
        int startCol = std::accumulate(domain->solInfo().maxSizeAdjointBasis.begin(), 
                                       domain->solInfo().maxSizeAdjointBasis.begin()+adjointBasisId, 0);
        podSolver->setLocalBasis(startCol, blockCols);
        podSolver->factor();
      }
      else { std::cerr << "ERROR: adjoint basis is not defined for " << type << " quantity of interest\n"; }
    }
  }
}
