#include <Utils.d/dbg_alloca.h>
#include <iostream>
#include <cstdio>
#include <cmath>

#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>
#include <Solvers.d/Rbm.h>

#include <Math.d/FullMatrix.h>
#include <Sfem.d/Sfem.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/MathUtils.h>

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

extern Sfem *sfem;

template<class T, class VectorType, class SolverType>
int
SingleDomainStatic<T, VectorType, SolverType>::solVecInfo()
{
 int ret = domain->numUncon();
 if(domain->solInfo().inpc) ret *= sfem->getP();
 return ret;
}

template<class T, class VectorType, class SolverType>
int
SingleDomainStatic<T, VectorType, SolverType>::solVecInfo(int i)
{
 int ret = (domain->numUncon());
 ret *= i;
 return ret;
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::clean()
{
// RT
  solver->clean_up();
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getRHS(VectorType &rhs)
{
 if(domain->solInfo().loadcases.size() > 0)
   filePrint(stderr," ... Building the Force (Case %2d)   ...\n", domain->solInfo().loadcases.front());
 else
   filePrint(stderr," ... Building the Force             ...\n");

 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 startTimerMemory(times->formRhs, times->memoryRhs);
 domain->template buildRHSForce<T>(rhs, kuc);
 
 // rigid body mode projector (or eigen mode projector)
 bool useProjector = (domain->solInfo().filterFlags || domain->solInfo().modeFilterFlag);
 if(useProjector) 
   trProject(rhs);

 stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getRHSinpc(VectorType &rhs)
{
 filePrint(stderr," ... Building the Force   (inpc)    ...\n");

 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 startTimerMemory(times->formRhs, times->memoryRhs);
 rhs.zero();
 for (int i=0; i<(*allOps.rhs_inpc).size(); ++i) rhs[i] = (*allOps.rhs_inpc)[i]; 

 int useProjector=domain->solInfo().filterFlags;
 if(useProjector)
   trProject(rhs);

 stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::preProcessSA()
{
 domain->buildPreSensitivities<T>(allSens, bcx);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::postProcessSA(GenVector<T> &sol)
{
 domain->buildPostSensitivities<T>(allOps.sysSolver, allOps.K, allOps.spm, allSens, &sol, bcx);
 domain->sensitivityPostProcessing(allSens, &sol, bcx);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::preProcess()
{
	// Allocate space for the Static Timers
	times = new StaticTimers;

	startTimerMemory(times->preProcess, times->memoryPreProcess);

	// Makes renumbering, connectivities and dofsets
	domain->preProcessing();

	int numdof = domain->numdof();

	times->makeBCs -= getTime();
	std::vector<int> bc(numdof);
	bcx     = new T[numdof];

	// Make boundary conditions info
	if(domain->getImplicitFlag() || domain->nCDirichlet()) {
		if(domain->getImplicitFlag()) bcxC = new DComplex[numdof * domain->getNumWaveDirections()];
		else bcxC = new DComplex[numdof];
		((HData *)domain)->make_bc(domain, bc.data(), bcxC);
		for(int i=0; i<numdof; ++i) ScalarTypes::copy(bcx[i],bcxC[i]); // temp fix, needs to be done for every direction before post processing
	}
	else domain->make_bc(bc.data(),bcx);
	times->makeBCs += getTime();

	// Now, call make_constrainedDSA(bc) to  built c_dsa
	// that will incorporate all the boundary conditions info
	times->makeDOFs -= getTime();
	domain->make_constrainedDSA(bc);
	domain->makeAllDOFs();
	times->makeDOFs += getTime();

	// if we have initial displacements, we have to consider
	// the nonlinear tangent stiffness matrix instead of the
	// linear stiffness matrix. This is due to the prestress.

	kelArray  = 0;
	geomState = 0;
	allCorot  = 0;

	if(domain->numInitDisp6() > 0 && domain->solInfo().gepsFlg == 1) {
		FullSquareMatrix *geomKelArray=0, *dummy=0;
		domain->computeGeometricPreStress(allCorot, geomState, kelArray, times,
		                                  geomKelArray, dummy);
	}

	stopTimerMemory(times->preProcess, times->memoryPreProcess);

	int useProjector = domain->solInfo().filterFlags;
	int useHzemFilter = domain->solInfo().hzemFilterFlag;
	int useSlzemFilter = domain->solInfo().slzemFilterFlag;

	if(!rigidBodyModes) {
		// ... Construct geometric rigid body modes if necessary
		if(useProjector || domain->solInfo().rbmflg) {
			rigidBodyModes = domain->constructRbm();
			if(useProjector) std::cout << " ... RBM Filter Requested           ..." << std::endl;
		}
			// ... Construct "thermal rigid body mode" if necessary
		else if(useHzemFilter || domain->solInfo().hzemFlag) {
			rigidBodyModes = domain->constructHzem();
			if(useHzemFilter) std::cout << " ... HZEM Filter Requested          ..." << std::endl;
		}
			// ... Construct "sloshing rigid body mode" if necessary
		else if(useSlzemFilter || domain->solInfo().slzemFlag) {
			rigidBodyModes = domain->constructSlzem();
			if(useSlzemFilter) std::cout << " ... SLZEM Filter Requested         ..." << std::endl;
		}
	}

	if(domain->solInfo().rbmflg || domain->solInfo().hzemFlag || domain->solInfo().slzemFlag) {
		domain->template getSolverAndKuc<T>(allOps, kelArray, rigidBodyModes);
	}
	else {
		domain->template getSolverAndKuc<T>(allOps, kelArray, (Rbm*)NULL);
	}
	solver = allOps.sysSolver;
	kuc = allOps.Kuc;
	kcc = allOps.Kcc;

	if(useProjector || useHzemFilter || useSlzemFilter)
		projector_prep(rigidBodyModes, allOps.M);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::scaleDisp(VectorType &sol)
{
  domain->scaleDisp(sol.data());
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::scaleInvDisp(VectorType &sol)
{
  domain->scaleInvDisp(sol.data());
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::scaleDisp(VectorType &sol, double alpha)
{
  domain->scaleDisp(sol.data(), alpha);
}


template<class T, class VectorType, class SolverType>
SolverType *
SingleDomainStatic<T, VectorType, SolverType>::getSolver()
{
  return solver;
}

template<class T, class VectorType, class SolverType>
SingleDomainPostProcessor<T, VectorType, SolverType> *
SingleDomainStatic<T, VectorType, SolverType>::getPostProcessor()
{
 return new SingleDomainPostProcessor<T,VectorType,SolverType>(domain,bcx,times,solver,kuc,kcc);
}

// ... The next function is the only one needed for Statics and FAcoustics
template<class T, class VectorType, class SolverType>
void
SingleDomainPostProcessor<T, VectorType, SolverType>::staticOutput(VectorType &sol, 
                                                                   VectorType &force, 
                                                                   bool printTimers, int ndflag)
{
#ifdef DISTRIBUTED
 if(structCom->myID() != 0) return; // used for parallel mumps, only one process should write the output file
#endif
 startTimerMemory(times->output, times->memoryOutput);
 domain->template postProcessing<T>(sol,bcx,force,ndflag,0,0,0,kuc,kcc);
 stopTimerMemory(times->output, times->memoryOutput);

 long memoryUsed = 0;
 double solveTime = 0.0;

 memoryUsed = solver->size();
 solveTime  = solver->getSolutionTime();
 if(printTimers) {
   times->printStaticTimers(solveTime, memoryUsed, domain);
 }

 if(ndflag <= 1) filePrint(stderr," --------------------------------------\n");
}

template<class T, class VectorType, class SolverType>
void
SingleDomainPostProcessor<T, VectorType, SolverType>::staticOutput(GeomState &geomState, 
                                                                   double x)
{
  times->output -= getTime();
  domain->template postProcessing<T>(geomState,x);
  times->output += getTime();

  MatrixTimers mt   = domain->getTimers();
  double solveTime  = solver->getSolutionTime();
  double memoryUsed = solver->getMemory();
}

template<class T, class VectorType, class SolverType>
void
SingleDomainPostProcessor<T, VectorType, SolverType>::staticOutput(GeomState *geomState, 
                                                                   double x)
{
  startTimerMemory(times->output, times->memoryOutput);
  domain->template postProcessing<T>(*geomState,x);
  stopTimerMemory(times->output, times->memoryOutput);

  MatrixTimers mt  = domain->getTimers();
  double solveTime = solver->getSolutionTime();
}


//------------------------------------------------------------------------------

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getFreqSweepRHS(VectorType *rhs, VectorType **u, int k)
{ 
  startTimerMemory(times->formRhs, times->memoryRhs);
  double omega2 = geoSource->shiftVal();

  double omega = sqrt(omega2);

  VectorType *vec = new VectorType(solVecInfo());

  if (u==0) {
    for(int i=0; i<vec->size(); ++i)
      (*vec)[i] = 0;
  } else {
    for(int i=0; i<vec->size(); ++i)
      (*vec)[i] = double(k)*(double(k-1)*(*u[k-1])[i] + 2.0*omega*(*u[k])[i]);
    allOps.M->mult(vec->data(), rhs->data());

    if(allOps.C_deriv) {
      for(int j=0; j<=k-1; ++j) {
        if(allOps.C_deriv[k-j-1]) {
          double ckj = DCombination(k,j);
          for(int i=0; i<vec->size(); ++i) (*vec)[i] = -ckj*(*u[j+1])[i];
          allOps.C_deriv[k-j-1]->multAdd(vec->data(), rhs->data());
        }
      }
    }
    if(allOps.K_deriv) {
      for(int j=0; j<=k-1; ++j) {
        if(allOps.K_deriv[k-j]) {
          double ckj = DCombination(k,j);
          for(int i=0; i<vec->size(); ++i) (*vec)[i] = -ckj*(*u[j+1])[i];
          allOps.K_deriv[k-j]->multAdd(vec->data(), rhs->data());
        }
      }
    }
  }
  delete vec;
  domain->template buildFreqSweepRHSForce<T>(*rhs, allOps.Muc, allOps.Cuc_deriv,allOps.Kuc_deriv, k, omega);
  stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getRHS(VectorType &rhs, double omega, double deltaomega)
{ 
  VectorType *vec = new VectorType(solVecInfo());
  domain->template buildRHSForce<T>(rhs, *vec, kuc, allOps.Muc,
                                    allOps.Cuc_deriv, allOps.Kuc_deriv,
                                    allOps.Kuc_arubber_l,allOps.Kuc_arubber_m,
                                    omega, deltaomega);
  delete vec;
}

//------------------------------------------------------------------------------





template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getWCAWEFreqSweepRHS(VectorType *rhs, VectorType **wcawe_u, T* pU, T*pb, int maxRHS,  int iRHS)
{ 
  startTimerMemory(times->formRhs, times->memoryRhs);

  double omega2 = geoSource->shiftVal();
  double omega = sqrt(omega2);

  VectorType *vec = new VectorType(solVecInfo());
//  pu = -Z{2}*W(:,ii-1); Z{2} = (-2.0*freqn * M+i*C)/1.0;
  for(int i=0; i<vec->size(); ++i)
      (*vec)[i] = 2.0*omega*(*wcawe_u[iRHS-1])[i];
  allOps.M->mult(vec->data(), rhs->data());
  if(allOps.C_deriv) {
    if(allOps.C_deriv[0]) {
      for(int i=0; i<vec->size(); ++i)
        (*vec)[i] = -(*wcawe_u[iRHS-1])[i];
      allOps.C_deriv[0]->multAdd(vec->data(), rhs->data()); 
    }
  }
//  pu = pu + b{j+1}*PP(1,ii-j);
  double factorial = 1.0;
  for(int j=1;j<=iRHS;j++) {
    factorial *= double(j);
    for(int i=0; i<vec->size(); ++i)
      (*vec)[i] = 0;
    domain->template buildFreqSweepRHSForce<T>(*vec,
                           allOps.Muc, allOps.Cuc_deriv,allOps.Kuc_deriv, j, omega);
    for(int i=0; i<vec->size(); ++i)
      (*rhs)[i] += pb[j-1]* (*vec)[i]/factorial;
  }
//  pu = pu - Z{j+1}*W(:,1:(ii-j))*PP(:,ii-j);
//  assume j > 2 all 0 (not so for rubber), and for j=2: Z{3} = -M;
  if (iRHS>1) for(int j=2;j<=2;j++) {
    for(int i=0; i<vec->size(); ++i)
      (*vec)[i] = 0;
    for(int k=0;k<iRHS+1-j;k++) {
      for(int i=0; i<vec->size(); ++i)
        (*vec)[i] += pU[k+(j-1)*(iRHS+1)]* (*wcawe_u[k])[i];
    }
    allOps.M->multAdd(vec->data(), rhs->data());
  }

  delete vec;
  stopTimerMemory(times->formRhs, times->memoryRhs);
}

