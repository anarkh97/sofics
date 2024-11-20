#include <Paral.d/MDTemp.h>
#include <Driver.d/DecDomain.h>
#include <Dist.d/DistDom.h>
#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>
#include <Driver.d/SysState.h>
#include <Math.d/SparseMatrix.h>

MultiDomainTemp::MultiDomainTemp(Domain *d)
{
  domain = d;

#ifdef DISTRIBUTED
  decDomain = new GenDistrDomain<double>(domain);
#else
  decDomain = new GenDecDomain<double>(domain);
#endif

  kelArray = 0;
  allCorot = 0;
  geomState = 0;
  dynMat = 0;
}

MultiDomainTemp::~MultiDomainTemp()
{
  int nsub = decDomain->getNumSub();
  if(geomState) delete geomState;
  if(kelArray) { for(int i=0; i<nsub; ++i) delete [] kelArray[i]; delete [] kelArray; }
  if(allCorot) {
    for(int i=0; i<nsub; ++i) {
      if(allCorot[i]) {
        for(int iElem = 0; iElem < decDomain->getSubDomain(i)->numElements(); ++iElem) {
          if(allCorot[i][iElem] && (allCorot[i][iElem] != dynamic_cast<Corotator*>(decDomain->getSubDomain(i)->getElementSet()[iElem])))
            delete allCorot[i][iElem];
        }
        delete [] allCorot[i];
      }
    }
    delete [] allCorot;
  }
  delete decDomain;
}

void
MultiDomainTemp::temptrProject(DistrVector &f)
{
  filePrint(stderr, " *** WARNING: MultiDomainTemp::temptrProject is not implemented\n");
}

void
MultiDomainTemp::tempProject(DistrVector &v)
{
  filePrint(stderr, " *** WARNING: MultiDomainTemp::tempProject is not implemented\n");
}

DistrInfo &
MultiDomainTemp::solVecInfo()
{
  return decDomain->solVecInfo();
}

void
MultiDomainTemp::getTempTimes(double &dtemp, double &tmax)
{
  dtemp = domain->solInfo().getTimeStep();
  tmax = domain->solInfo().tmax;
}

void
MultiDomainTemp::getTempNeum(double &epsiln)
{
  epsiln = domain->solInfo().alphaTemp;
}

int
MultiDomainTemp::getTimeIntegration()
{
  return domain->solInfo().timeIntegration;
}

int
MultiDomainTemp::getAeroheatFlag()
{
  return domain->solInfo().aeroheatFlag;
}

int
MultiDomainTemp::getThermohFlag()
{
  return domain->solInfo().thermohFlag;
}

int
MultiDomainTemp::getHzemFlag()
{
  return domain->solInfo().hzemFilterFlag;
}

int
MultiDomainTemp::getZEMFlag()
{
  return domain->solInfo().hzemFlag;
}

void
MultiDomainTemp::getSteadyStateParam(int &steadyFlag, int &steadyMin,
                                         int &steadyMax, double &steadyTol)
{
  steadyFlag = domain->solInfo().steadyFlag;
  steadyMin  = domain->solInfo().steadyMin;
  steadyMax  = domain->solInfo().steadyMax;
  steadyTol  = domain->solInfo().steadyTol;
}

void
MultiDomainTemp::getQuasiStaticParameters(double &maxVel, double &qsbeta)
{
  maxVel = domain->solInfo().qsMaxvel;
  qsbeta = domain->solInfo().qsBeta;
}

void
MultiDomainTemp::tempInitState(TempState<DistrVector> &inState)
{
  execParal(decDomain->getNumSub(), this, &MultiDomainTemp::subTempInitState, inState);
}

void
MultiDomainTemp::subTempInitState(int isub, TempState<DistrVector> &inState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector d(inState.getDisp().subData(isub), inState.getDisp().subLen(isub));
  StackVector v(inState.getVeloc().subData(isub), inState.getVeloc().subLen(isub));
  StackVector v_p(inState.getPrevVeloc().subData(isub), inState.getPrevVeloc().subLen(isub));
  sd->initTempVector(d, v, v_p);
}

void
MultiDomainTemp::computeExtForce(DistrVector &ext_f, double t, int tIndex, DistrVector &prev_f)
{
  execParal(decDomain->getNumSub(), this, &MultiDomainTemp::subComputeExtForce, ext_f, t, tIndex, prev_f);

  // apply projector ONLY for Dynamics, hzemFilterFlag is set to zero
  // in the beginning for quasistatic
  int useFilter = domain->solInfo().hzemFilterFlag;
  if(useFilter) temptrProject(ext_f);
}

void
MultiDomainTemp::subComputeExtForce(int isub, DistrVector &ext_f, double t, int tIndex, DistrVector &prev_f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subext_f(ext_f.subData(isub), ext_f.subLen(isub));
  StackVector subprev_f(prev_f.subData(isub), prev_f.subLen(isub));
  SparseMatrix *kuc = (*(dynMat->Kuc))[isub];
  sd->computeExtForce(subext_f, t, tIndex, kuc, subprev_f);
}

void
MultiDomainTemp::preProcess()
{
  // Makes local renumbering, connectivities and dofsets
  decDomain->preProcess();

  // Make all element's dofs
  execParal(decDomain->getNumSub(), this, &MultiDomainTemp::makeAllDOFs);

  // Make the geomState (used for etemp and explicit nonlinear)
  if((domain->solInfo().gepsFlg == 1 && domain->numInitDisp6() > 0) || domain->solInfo().isNonLin()) {
    geomState = new DistrGeomState(decDomain);
  }

  // Update geomState with prescribed dirichlet boundary conditions (explicit nonlinear and contact)
  if(domain->solInfo().isNonLin())
    execParal(decDomain->getNumSub(), this, &MultiDomainTemp::initSubPrescribedTemperature);

  // Make corotators and kelArray (used for ETEMP and explicit nonlinear)
  if((domain->solInfo().gepsFlg == 1 && domain->numInitDisp6() > 0) || domain->solInfo().isNonLin()) {
    allCorot = new Corotator**[decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MultiDomainTemp::makeSubCorotators);

    kelArray = new FullSquareMatrix*[decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MultiDomainTemp::makeSubElementArrays);
  }
}

void
MultiDomainTemp::makeAllDOFs(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  sd->makeAllDOFs();
}

void
MultiDomainTemp::initSubPrescribedTemperature(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  if(sd->nDirichlet() > 0)
    (*geomState)[isub]->updatePrescribedDisplacement(sd->getDBC(), sd->nDirichlet(), sd->getNodes());
}

void
MultiDomainTemp::makeSubCorotators(int isub)
{
  SubDomain *sd  = decDomain->getSubDomain(isub);
  int numele     = sd->numElements();
  allCorot[isub] = new Corotator*[numele];
  sd->createCorotators(allCorot[isub]);
}

void
MultiDomainTemp::makeSubElementArrays(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);

  // allocate the element stiffness array
  sd->createKelArray(kelArray[isub]);

  // update geomState with ETEMP if specified (linear only)
  if((sd->numInitDisp6() > 0) && (domain->solInfo().gepsFlg == 1))
    (*geomState)[isub]->updatePrescribedDisplacement(sd->getInitDisp6(), sd->numInitDisp6(), sd->getNodes());

  // build the element stiffness matrices.
  Vector elementInternalForce(sd->maxNumDOF(), 0.0);
  Vector residual(sd->numUncon(), 0.0);
  sd->getStiffAndForce(*(*geomState)[isub], elementInternalForce, allCorot[isub], kelArray[isub], residual,
                       1.0, 0.0, (*geomState)[isub], (Vector*) NULL, NULL);
}

MDDynamMat
MultiDomainTemp::buildOps(double coeM, double coeC, double coeK)
{
  // Have each subdomain create their operators, then put the
  // dynamic matrices in the Feti Solver
  dynMat = new MDDynamMat;

  // XXX: check HZEM and HZEMFILTER
  decDomain->buildOps(*dynMat, coeM, coeC, coeK, (Rbm **) 0, kelArray);

  return *dynMat;
}

void
MultiDomainTemp::getInternalForce(DistrVector& d, DistrVector& f)
{
  if(domain->solInfo().isNonLin()) geomState->explicitUpdate(decDomain, d);
  execParal(decDomain->getNumSub(), this, &MultiDomainTemp::subGetInternalForce, d, f);
}

void
MultiDomainTemp::subGetInternalForce(int isub, DistrVector& d, DistrVector& f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subf(f.subData(isub), f.subLen(isub));
  if(domain->solInfo().isNonLin()) {
    Vector elementInternalForce(sd->maxNumDOF(), 0.0);
    Vector residual(sd->numUncon(), 0.0);
    sd->getInternalForce(*(*geomState)[isub], elementInternalForce, allCorot[isub], kelArray[isub], residual,
                         1.0, 1.0, (*geomState)[isub]);
    // Note: a dummy value of t is passed to getInternalForce above. This is ok because there are no time-dependent
    // follower forces for temperature dofs
    subf.linC(residual, -1.0);
  }
  else {
    StackVector subd(d.subData(isub), d.subLen(isub));
    subf.zero();
    sd->getKtimesU(subd, sd->getBcx(), subf, 1.0, ((kelArray) ? kelArray[isub] : NULL));
  }
}

void
MultiDomainTemp::aeroHeatPreProcess(DistrVector& d_n, DistrVector& v_n, DistrVector& v_p)
{
  filePrint(stderr, " *** WARNING: MultiDomainTemp::aeroHeatPreProcess is not implemented\n");
}

void
MultiDomainTemp::thermohPreProcess(DistrVector& d_n, DistrVector& v_n, DistrVector& v_p)
{
  filePrint(stderr, " *** WARNING: MultiDomainTemp::thermohPreProcess is not implemented\n");
}

void
MultiDomainTemp::getInitialTime(int &initTimeIndex, double &initTime)
{
  initTimeIndex = domain->solInfo().initialTimeIndex;
  initTime      = domain->solInfo().initialTime;
}

MDTempDynamPostProcessor *
MultiDomainTemp::getPostProcessor()
{
  return new MDTempDynamPostProcessor(decDomain, geomState, allCorot);
}

void
MDTempDynamPostProcessor
::tempdynamOutput(int tIndex, MDDynamMat& dMat, DistrVector& ext_f, TempState<DistrVector> &state)
{
  double t = tIndex*domain->solInfo().getTimeStep(); // XXX fixed timestep only
  SysState<DistrVector> distState(state.getDisp(), state.getVeloc(), state.getPrevVeloc());
  decDomain->postProcessing(state.getDisp(), ext_f, t, NULL, tIndex, &dMat, &distState);
}

int
MultiDomainTemp::getModeDecompFlag()
{
  return domain->solInfo().modeDecompFlag;
}

void
MultiDomainTemp::modeDecomp(double t, int tIndex, DistrVector& d_n)
{
  filePrint(stderr, " *** WARNING: MultiDomainTemp::modeDecomp is not implemented\n");
}

int
MultiDomainTemp::cmdComHeat(int cmdFlag)
{
  filePrint(stderr, " *** WARNING: MultiDomainTemp::cmdComHeat is not implemented\n");
  return 0;
}
