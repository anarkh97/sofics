#include <Utils.d/dbg_alloca.h>
#include <algorithm>
#include <cstdio>
#include <Utils.d/Memory.h>

#include <Driver.d/Domain.h>
#include <Problems.d/NonLinStatic.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/TemperatureState.h>
#include <Corotational.d/utilities.h>
#include <Solvers.d/Solver.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Timers.d/GetTime.h>
#include <Problems.d/StaticDescr.h>

extern int verboseFlag;

NonLinStatic::NonLinStatic(Domain *d)
{
  domain = d;
  kelArray = 0;
  allCorot = 0;
  bcx = 0;
  solver = 0;
  prec = 0;
  times = 0;
  reactions = 0;

  if(domain->GetnContactSurfacePairs()) {
    domain->InitializeStaticContactSearch(MortarHandler::CTC);
    updateCS = true;
  }
}

NonLinStatic::~NonLinStatic()
{
  clean();
  if(times) delete times;
  if(reactions) delete reactions;
}

int
NonLinStatic::solVecInfo()
{
  return domain->numUncon();
}

int
NonLinStatic::sysVecInfo()
{
  return 0;
}

void
NonLinStatic::clean()
{
  if(bcx)      { delete [] bcx; bcx = 0; }
  if(solver)   { delete solver; solver = 0; }
  if(prec)     { delete prec; prec = 0; }
  if(kelArray) { delete [] kelArray; kelArray = 0; }
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
NonLinStatic::updateContactSurfaces(GeomState& geomState)
{
  clean();
  domain->UpdateSurfaces(MortarHandler::CTC, &geomState);
  domain->PerformStaticContactSearch(MortarHandler::CTC);
  domain->deleteSomeLMPCs(mpc::ContactSurfaces);
  domain->ExpComputeMortarLMPC(MortarHandler::CTC);
  domain->UpdateContactSurfaceElements(&geomState);
  preProcess(false);
  geomState.resizeLocAndFlag(*domain->getCDSA());
}

void
NonLinStatic::updateStates(GeomState *refState, GeomState& geomState, double lambda)
{
  if(domain->solInfo().piecewise_contact) updateCS = true;
  domain->updateStates(refState, geomState, allCorot, lambda);
}

double
NonLinStatic::getStiffAndForce(GeomState& geomState, Vector& residual, Vector& elementInternalForce, 
                               Vector &, double lambda, GeomState *refState, bool forceOnly)
{
  times->buildStiffAndForce -= getTime();

  if(domain->GetnContactSurfacePairs()) {
    if(!domain->solInfo().piecewise_contact || updateCS) {
      updateContactSurfaces(geomState);
      updateCS = false;
    }
    residual.conservativeResize(domain->getCDSA()->size());
    elementInternalForce.resize(domain->maxNumDOF());
  }

  reactions->zero();
  if(forceOnly) {
    domain->getInternalForce(geomState, elementInternalForce, allCorot,
                             kelArray, residual, lambda, 0, refState, reactions);
  }
  else {
    domain->getStiffAndForce(geomState, elementInternalForce, allCorot, 
                             kelArray, residual, lambda, 0, refState, reactions);
  }

  times->buildStiffAndForce += getTime();

  return sqrt(residual*residual);
}

void
NonLinStatic::updatePrescribedDisplacement(GeomState *geomState, double)
{
 // Measure the time necessary to update the prescribed displacments

 times->timePresc -= getTime();

 int numDirichlet = domain->nDirichlet();

 BCond *dbc = domain->getDBC();

 double delta = domain->solInfo().getNLInfo().dlambda;

 geomState->updatePrescribedDisplacement(dbc, numDirichlet, delta);

 times->timePresc += getTime();
}

void
NonLinStatic::initializeParameters(int step, GeomState *geomState)
{
  if(step == 1 || domain->solInfo().reinit_lm) {
    domain->initializeMultipliers(*geomState, allCorot);
  }
  domain->initializeParameters();
}

void
NonLinStatic::updateParameters(GeomState *geomState)
{
  domain->updateMultipliers(*geomState, allCorot);
  domain->updateParameters();
}

bool
NonLinStatic::checkConstraintViolation(double &err, GeomState *gs)
{
  err = domain->getError(allCorot, *gs);
  return (err <= domain->solInfo().penalty_tol);
}

int
NonLinStatic::checkConvergence(int iter, double normDv, double normRes)
{
 // Measure time necessary to check for convergence of Newton algorithm

 times->timeCheck -= getTime();

 // Note when useTolInc is false, this function is called before normDv is calculated
 bool useTolInc = (domain->solInfo().getNLInfo().tolInc != std::numeric_limits<double>::infinity() ||
                   domain->solInfo().getNLInfo().absTolInc != std::numeric_limits<double>::infinity());

 if(iter == 0) {
   if(useTolInc) firstDv  = normDv;
   else { normDv = 0; firstDv = 1; }
   firstRes = normRes;
 }
 else if(iter == 1 && !useTolInc) {
  firstDv  = normDv;
 }

 double relativeDv  = normDv/firstDv;
 double relativeRes = normRes/firstRes;

 if(verboseFlag) {
   filePrint(stderr," ----------------------------------------------------\n");
   if(useTolInc || iter >= 1) {
     filePrint(stderr, " Newton Iter    #%d\tcurrent dv   = % e\n \t\t\t"
                       "first dv     = % e\n \t\t\trelative dv  = % e\n",
                       iter+1, normDv, firstDv, relativeDv);
     filePrint(stderr, "                \tcurrent Res  = % e\n \t\t\t"
                       "first Res    = % e\n \t\t\trelative Res = % e\n",
                       normRes, firstRes, relativeRes);
   }
   else {
     filePrint(stderr, " Newton Iter    #%d\tcurrent Res  = % e\n \t\t\t"
                       "first Res    = % e\n \t\t\trelative Res = % e\n",
                       iter+1, normRes, firstRes, relativeRes);
   }
   filePrint(stderr," ----------------------------------------------------\n");
 }

 int converged = 0;

 // Check convergence criteria
 if(iter > 0 && ((normRes <= tolerance*firstRes && normDv <= domain->solInfo().getNLInfo().tolInc*firstDv)
    || (normRes < domain->solInfo().getNLInfo().absTolRes && normDv < domain->solInfo().getNLInfo().absTolInc)))
   converged = 1;

 // Check divergence
 else if(iter > 0 && normRes > 10000000*firstRes)
   converged = -1;

 // Store residual norm and dv norm for output
 times->norms[iter].normDv      = normDv;
 times->norms[iter].relativeDv  = relativeDv;
 times->norms[iter].normRes     = normRes;
 times->norms[iter].relativeRes = relativeRes;

 times->timeCheck += getTime();

 return converged;
}

GeomState*
NonLinStatic::createGeomState()
{
 times->timeGeom -= getTime();

 GeomState *geomState;
 if(domain->solInfo().soltyp == 2) 
   geomState = (GeomState *) new TemperatureState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes());
 else
   geomState = new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes(), &domain->getElementSet(),
                             domain->getNodalTemperatures()); 

 times->timeGeom += getTime();

 return geomState;
}

int
NonLinStatic::reBuild(int iteration, int step, GeomState&)
{
 times->rebuild -= getTime();

 int rebuildFlag = 0;

 if(iteration % domain->solInfo().getNLInfo().updateK == 0 && (step-1) % domain->solInfo().getNLInfo().stepUpdateK == 0) {
   if(domain->solInfo().mpcDirect != 0) {
     if(solver) delete solver;
     if(prec) delete prec;
     preProcess();
   }
   else {
     spm->zeroAll();
     AllOps<double> ops;
     if(spp) { // rebuild preconditioner as well as the solver
       spp->zeroAll();
       ops.spp = spp;
     }
     domain->makeSparseOps<double>(ops, 1.0, 0.0, 0.0, spm, kelArray, (FullSquareMatrix *) NULL);
     domain->getTimers().factor -= getTime();
     solver->factor();
     if(prec) prec->factor();
     domain->getTimers().factor += getTime();
  }
  rebuildFlag = 1;
 }

 times->rebuild += getTime();

 return rebuildFlag;
}

int
NonLinStatic::elemVecInfo()
{
  return domain->maxNumDOF();
}

int
NonLinStatic::getMaxit()
{
  return domain->solInfo().getNLInfo().maxiter;
}

double
NonLinStatic::getScaleFactor()
{
 return domain->solInfo().getNLInfo().lfactor;
}

double
NonLinStatic::getDeltaLambda0()
{
 return domain->solInfo().getNLInfo().dlambda;
}

double
NonLinStatic::getMaxLambda()
{
 return domain->solInfo().getNLInfo().maxLambda;
}

LinesearchInfo&
NonLinStatic::linesearch()
{
 return domain->solInfo().getNLInfo().linesearch;
}

void
NonLinStatic::getRHS(Vector& rhs)
{
 // ... BUILD THE RHS FORCE (not including follower or internal forces)
 times->formRhs -= getTime();
 domain->computeConstantForce<double>(rhs);
 times->formRhs += getTime();
}

void
NonLinStatic::preProcess(bool factor)
{
	// Allocate space for the Static Timers
	if(!times) times = new StaticTimers;

	startTimerMemory(times->preProcess, times->memoryPreProcess);

	times->timePre -= getTime();
	domain->preProcessing();
	times->timePre += getTime();

	int numdof = domain->numdof();

	times->makeBCs -= getTime();
	std::vector<int> bc(numdof);
	if(!bcx) bcx = new double[numdof];

	// Make the boundary conditions info
	domain->make_bc(bc.data(), bcx);
	if(!reactions) {
		reactions = new Vector(domain->nDirichlet());
		reactions->zero();
	}

	times->makeBCs += getTime();

	// Now, call make_constrainedDSA(bc) to
	// built c_dsa that will incorporate all
	// the boundary conditions info

	times->makeDOFs -= getTime();
	domain->make_constrainedDSA(bc);
	domain->makeAllDOFs();
	times->makeDOFs += getTime();

	stopTimerMemory(times->preProcess, times->memoryPreProcess);
	AllOps<double> allOps;

	long buildMem = -memoryUsed();
	times->timeBuild -= getTime();

	if(domain->solInfo().rbmflg) {
		Rbm *rigidBodyModes = domain->constructRbm(); // new policy is to construct rbms if GRBM is requested in input file
		// but only use them when it is appropriate to do so. In nonlinear statics it is not
		// since the nullity of the tangent stiffness matrix may be less than the nullity
		// of the number of rigid body modes
		delete rigidBodyModes;
	}

	domain->buildOps<double>(allOps, 1.0, 0.0, 0.0, (Rbm *) NULL, kelArray,
	                         (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL, factor);
	times->timeBuild += getTime();
	buildMem += memoryUsed();

	solver = allOps.sysSolver;
	spm = allOps.spm;
	prec = allOps.prec;
	spp = allOps.spp;

	if(!allCorot) { // first time only
		// ... CREATE THE ARRAY OF ELEMENT COROTATORS
		startTimerMemory(times->preProcess, times->memoryPreProcess);
		times->corotatorTime -= getTime();
		allCorot = new Corotator *[domain->numElements()];
		domain->createCorotators(allCorot);
		times->corotatorTime += getTime();
	}

	if(!kelArray) { // first time only
		// ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
		times->kelArrayTime -= getTime();
		domain->createKelArray(kelArray);
		times->kelArrayTime += getTime();
		stopTimerMemory(times->preProcess, times->memoryPreProcess);
	}

	// Set the nonlinear tolerance used for convergence
	tolerance = domain->solInfo().getNLInfo().tolRes;
}

Solver *
NonLinStatic::getSolver()
{
  return solver;
}

SingleDomainPostProcessor<double, Vector, Solver> *
NonLinStatic::getPostProcessor()
{
  return new SingleDomainPostProcessor<double,Vector,Solver>(domain,bcx,times,solver);
}

void
NonLinStatic::printTimers()
{
  times->timeTimers -= getTime();

  long memoryUsed = solver->size();
  double solveTime  = solver->getSolutionTime();

  times->printStaticTimers(solveTime, memoryUsed, domain);

  times->timeTimers += getTime();

  if(domain->solInfo().massFlag) {
    double mass = domain->computeStructureMass();
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," ... Structure mass = %e  ...\n",mass);
    filePrint(stderr," --------------------------------------\n");
  }
}

double
NonLinStatic::getTolerance()
{
  return std::max(tolerance*firstRes, domain->solInfo().getNLInfo().absTolRes);
}

void
NonLinStatic::staticOutput(GeomState *geomState, double lambda, Vector& force,
                           Vector &, GeomState *refState)
{
  times->output -= getTime();
  Vector dummyForce(domain->numUncon(), 0.0);
  int step = (int)std::floor(lambda/domain->solInfo().getNLInfo().dlambda+0.5);
  domain->postProcessing(geomState, force, dummyForce, lambda, step, 0, 0, allCorot,
                         (double *) 0, (double *) 0, refState, reactions);
  times->output += getTime();
}

double
NonLinStatic::getEnergy(double lambda, Vector& force, GeomState* geomState)
{
  Vector sol(domain->numUncon(), 0.0);
  for(int i=0; i<domain->numNode(); ++i) {
    Node *node_i = domain->getNodes()[i];
    int xloc  = domain->getCDSA()->locate(i, DofSet::Xdisp);
    if(xloc >= 0 && node_i) {
      sol[xloc]  = ( (*geomState)[i].x - node_i->x);
    }
    int yloc  = domain->getCDSA()->locate(i, DofSet::Ydisp);
    if(yloc >= 0 && node_i)
      sol[yloc]  = ( (*geomState)[i].y - node_i->y);
    int zloc  = domain->getCDSA()->locate(i, DofSet::Zdisp);
    if(zloc >= 0 && node_i)
      sol[zloc]  = ( (*geomState)[i].z - node_i->z);
    double rot[3];
    mat_to_vec((*geomState)[i].R,rot);
    int xrot  = domain->getCDSA()->locate(i, DofSet::Xrot);
    if(xrot >= 0 && node_i)
      sol[xrot]  = rot[0];
    int yrot  = domain->getCDSA()->locate(i, DofSet::Yrot);
    if(yrot >= 0 && node_i)
      sol[yrot]  = rot[1];
    int zrot  = domain->getCDSA()->locate(i, DofSet::Zrot);
    if(zrot >= 0 && node_i)
      sol[zrot]  = rot[2];
  }

  // compute external energy not including follower forces
  double Wext = -lambda*force*sol;

  // compute internal energy and external energy due to follower forces
  // XXXX note: need to include pressure, gravity and thermal forces in getElementEnergy for this to work
  double Wela = 0.0;
  for(int i = 0; i < domain->numElements(); ++i)
     Wela += allCorot[i]->getElementEnergy(*geomState, domain->getNodes());

  // Total Energy = Wext + Wela
  return Wext + Wela;
}

double
NonLinStatic::getResidualNorm(Vector &rhs, GeomState &geomState)
{
  Vector res(rhs);
  domain->applyResidualCorrection(geomState, allCorot, res, 1.0);
  return solver->getResidualNorm(res);
}

bool
NonLinStatic::getResizeFlag()
{
  return (domain->GetnContactSurfacePairs() > 0);
}
