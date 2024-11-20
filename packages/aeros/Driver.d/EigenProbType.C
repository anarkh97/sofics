#include <cstdio>
#include <iostream>
#include <cmath>
#include <list>
#include <algorithm>
#include <string>

#include <Driver.d/Dynam.h>
#include <Utils.d/SolverInfo.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>
#include <Math.d/DBSparseMatrix.h>
#include <Utils.d/DistHelper.h> 
#include <Utils.d/linkfc.h>
#include <Utils.d/ModeData.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/Communicator.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/SolverInfo.h>

#ifdef DISTRIBUTED
#include <mpi.h>
#endif
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

extern int verboseFlag;
extern GeoSource *geoSource;
extern SolverInfo &solInfo;

//------------------------------------------------------------------------------

// build eigensolver

/*
template <class EigOps, class VecType, class VecSet,
           class PostProcessor, class ProblemDescriptor>
  EigenSolver< EigOps, VecType, VecSet,
               PostProcessor, ProblemDescriptor>
* EigenSolver< EigOps, VecType, VecSet,
               PostProcessor, ProblemDescriptor>
  ::buildEigenSolver(ProblemDescriptor * pDesc)
{
   EigenSolver< EigOps, VecType, VecSet, 
	        PostProcessor, ProblemDescriptor> * newSolver;
   int solverType = pDesc->getEigenSolverType();
     
   switch (solverType)
   {
      case SolverInfo::Arpack:
          newSolver = new SymArpackSolver<EigOps, VecType, VecSet,
                                         PostProcessor, ProblemDescriptor>
                                         (pDesc);
#ifndef USE_ARPACK
          filePrint(stderr," *** ERROR: executable not linked with ARPACK. See flag USE_ARPACK in Makefile.\n");
	  exit(-1);
#endif
          break;
      case SolverInfo::SubSpace :
           newSolver = new SubSpaceSolver<EigOps, VecType, VecSet, 
                                         PostProcessor, ProblemDescriptor>
                                         (pDesc);
          break;
      case SolverInfo::LobPcg :
           newSolver = new LOBPCGSolver<EigOps, VecType, VecSet, 
                                        PostProcessor, ProblemDescriptor>
                                        (pDesc);
	  break;
      default:
          filePrint(stderr," *** ERROR: unknown eigensolver.\n");
	  exit(-1);
   }

   return newSolver;
}
*/

// Global variable for mode data
extern ModeData modeData;
    
//------------------------------------------------------------------------------

// ... Solving linear Eigenproblem (A - lambda B) P = 0
// ... for numEig eigenvalues and eigenvectors

template <class EigOps, class VecType, class VecSet, 
	  class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet,
		PostProcessor, ProblemDescriptor>
::setUp()
{
	// ... Construct connectivities, renumbering and dofset array
	probDesc->preProcess();

	// ... Get the eigen post processor
	postProcessor = probDesc->getPostProcessor();

	// ... Get number of eigenvalues and eigenvectors
	numEig = probDesc->getNumEigen();

	// ... Build matrices A and B
	eM = new EigOps;

	// build, factor matrix in buildEigOps
	probDesc->buildEigOps( *eM );
	if(solInfo.printMatLabExit) return;

	// ... get the number of rigid body modes
	// if GRBM is not specified, or if a shift is specified, or if the solver is skyline/sparse/feti then get rbms from solver
	if(domain->solInfo().rbmflg == 0 || geoSource->shiftVal() != 0.0 || isFeti(domain->solInfo().solvercntl->type) ||
	   (!isFeti(domain->solInfo().solvercntl->type) &&
	    (domain->solInfo().solvercntl->subtype == 0 || domain->solInfo().solvercntl->subtype == 1))) {
		nrmod = eM->dynMat->numRBM();
	}
		// otherwise, use the geometric modes
	else {
		nrmod = eM->rigidBodyModes->numRBM();
	}

	if(solInfo.test_ulrich) nrmod = 0;

	totalEig = numEig+nrmod;
	// ... this could be an issue for small problems
	//     so confirm that this->totalEig is NOT bigger than problem size
	if(totalEig > probDesc->solVecSize()) totalEig = probDesc->solVecSize();

	int nsmax;
	double tolEig, tolJac;
	bool explicitK;
	this->probDesc->getSubSpaceInfo(origSubSize, nsmax, tolEig, tolJac, explicitK);

	// The following code is used to output only rigid body
	// modes from a given structure. Selected by setting
	// the subspace size equal to zero in the input file.
	// The number output is selected by the number of eigen values
	// also in the input file.
	if(origSubSize == 0) {
		// Print those rigid body modes to the output file
		filePrint(stderr," ... Output Rigid Body Modes and exit ...\n");
		VecSet RBMs(nrmod, probDesc->solVecInfo());
		// if GRBM is not specified, or if a shift is specified, or if the solver is skyline/sparse/feti then get rbms from solver
		if(domain->solInfo().rbmflg == 0 || geoSource->shiftVal() != 0.0 || isFeti(domain->solInfo().solvercntl->type) ||
		   (!isFeti(domain->solInfo().solvercntl->type) && (domain->solInfo().solvercntl->subtype == 0 || domain->solInfo().solvercntl->subtype == 1))) {
			eM->dynMat->getRBMs(RBMs);
		}
			// otherwise, use the geometric modes
		else {
			eM->rigidBodyModes->getRBMs(RBMs);
		}
		Vector eValues(nrmod, 0.0);
		postProcessor->eigenOutput(eValues, RBMs);
		return;
	}

	// ... Construct vector and vectorset for eigenvalues and eigenvectors
	if (!solInfo.doEigSweep) {
		eigVal = new Vector(totalEig);
		eigVec = new VecSet(totalEig, probDesc->solVecInfo());

		eigVal->zero();
	}

/*
 int eigSolverType = probDesc->getEigenSolverType();
 switch(eigSolverType){
   case SolverInfo::SubSpace:
     initialize(); // ... Initialize eigensolver
     subsolve();   // ... Solve eigenproblem 
     break;
   case SolverInfo::LobPcg: 
     initializeLOBPCG();
     subsolveLOBPCG();
     break;
#ifdef USE_ARPACK
   case SolverInfo::Arpack: //HB
     solveArpack();
     break;
#endif
   default: //HB
     filePrint(stderr," *** WARNING: selected eigensolver is NOT supported.\n");
 }
 this->probDesc->printTimers(eM->dynMat);

 // ... Free memory
 cleanup();

 // ... Output results
 postProcessor->eigenOutput(*eigVal,*eigVec);
*/

}
//------------------------------------------------------------------------------

template <class EigOps, class VecType, class VecSet, 
	  class PostProcessor, class ProblemDescriptor>
void
LOBPCGSolver< EigOps, VecType, VecSet, 
	       PostProcessor, ProblemDescriptor>::initialize()
{
   this->probDesc->getSubSpaceInfo(subSpaceSize, nsmax, tolEig, tolJac, explicitK);
   if (explicitK)
     filePrint(stderr, " ... Using Explicit Stiffness Matrix...\n");
 
   // ... adjust subSpaceSize
   // subSpaceSize = std::max(subSpaceSize,2*this->totalEig);
   if(subSpaceSize < this->totalEig) subSpaceSize = std::min(2*(this->totalEig-this->nrmod),this->totalEig-this->nrmod+8)+this->nrmod; // PJSA
   // ... this could be an issue for small problems
   //     so confirm that subSpace is NOT bigger than problem size
   if (subSpaceSize > this->probDesc->solVecSize())
     subSpaceSize = this->probDesc->solVecSize();

   nsub = subSpaceSize - this->nrmod;  // number of flexible modes
   //int nProbSize = 3*nsub + this->nrmod;
   //filePrint(stderr," ... subSpaceSize = %d (%d)\n",subSpaceSize, nsub);
   // ... Declaration and Initialization of Vector sets Q and Z
   Q = new VecSet(subSpaceSize, this->probDesc->solVecInfo());
   Z = new VecSet(subSpaceSize, this->probDesc->solVecInfo());

   // ... Allocate memory for eigen values and initialize to zero;
   subVal = new Vector( 2*nsub+this->nrmod, 0.0);
}

//------------------------------------------------------------------------------

template <class EigOps, class VecType, class VecSet, 
	  class PostProcessor, class ProblemDescriptor>
void
SubSpaceSolver< EigOps, VecType, VecSet, 
	       PostProcessor, ProblemDescriptor>::initialize()
{
   this->probDesc->getSubSpaceInfo(subSpaceSize, nsmax, tolEig, tolJac, explicitK);
   if (explicitK)
     filePrint(stderr, " ... Using Explicit Stiffness Matrix...\n");
 
   // ... adjust subSpaceSize
   //subSpaceSize = std::max(subSpaceSize,2*this->totalEig);
   if(subSpaceSize < this->totalEig) subSpaceSize = std::min(2*(this->totalEig-this->nrmod),this->totalEig-this->nrmod+8)+this->nrmod; // PJSA
   // ... this could be an issue for small problems
   //     so confirm that subSpace is NOT bigger than problem size
   if (subSpaceSize > this->probDesc->solVecSize())
     subSpaceSize = this->probDesc->solVecSize();

   //filePrint(stderr,"subSpaceSize = %d\n",subSpaceSize);
   // ... Declaration and Initialization of Vector sets Q and Z
   Q = new VecSet(subSpaceSize, this->probDesc->solVecInfo());
   Z = new VecSet(subSpaceSize, this->probDesc->solVecInfo());

   // ... Allocate memory for eigen values and initialize to zero;
   subVal = new Vector(subSpaceSize, 0.0);
   subOld = new Vector(subSpaceSize, 0.0);

   // ... nsub = # of flexible modes
   nsub = subSpaceSize - this->nrmod;
}

//------------------------------------------------------------------------------

/*
template <class EigOps, class VecType, class VecSet, 
	  class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet, 
	       PostProcessor, ProblemDescriptor>::cleanup()  
{
  if(Q)      { delete Q;      Q      = 0; }
  if(Z)      { delete Z;      Z      = 0; }
  if(subVal) { delete subVal; subVal = 0; }
  if(subOld) { delete subOld; subOld = 0; }
}
*/

//------------------------------------------------------------------------------
template <class EigOps, class VecType, class VecSet, 
	  class PostProcessor, class ProblemDescriptor>
void
LOBPCGSolver< EigOps, VecType, VecSet, 
	      PostProcessor, ProblemDescriptor>::solve()
{
  //... Check if multiple-shift eigen-analysis is requested
  if(solInfo.doEigSweep) {
    filePrint(stderr," *** ERROR: Multiple-shift eigenvalue calculation not implemented for LOBPCGSolver. Please use ARPACK.\n");
    exit(-1);
  }

  //... Setup
  this->setUp();
  if(solInfo.printMatLabExit) return; // just compute and print M and K
  if(this->origSubSize == 0) return; // just compute and print the rigid body modes

  //... Initialize
  initialize();

  // number of flex modes*3 plus number of rbm
  //int nProbSize = 3*nsub+this->nrmod;
  // ... Check if nsub (# of flexible modes) is less than zero
  if(nsub < 0) {
    filePrint(stderr," *** ERROR in routine LOBPCGSolver::solve()    ***\n");
    filePrint(stderr," *** subspaceSize =  %d                        ***\n", subSpaceSize);
    filePrint(stderr," *** while # of Rigid Body Modes (this->nrmod) = %d ***\n", this->nrmod);
    filePrint(stderr," *** subspaceSize must be bigger than this->nrmod   ***\n");
    return;
  }

  // ... check that subspace is big enough
  if(subSpaceSize < this->nrmod) {
    this->probDesc->error(subSpaceSize,this->nrmod);
    return;
  }

  // ... Store the rbms in the first this->nrmod vectors of VecSet Z
  if (this->nrmod) {
    // if GRBM is not specified, or if a shift is specified, or if the solver is skyline/sparse/feti then get rbms from solver
    if(domain->solInfo().rbmflg == 0 || geoSource->shiftVal() != 0.0 || isFeti(domain->solInfo().solvercntl->type) ||
       (!isFeti(domain->solInfo().solvercntl->type) && (domain->solInfo().solvercntl->subtype == 0 || domain->solInfo().solvercntl->subtype == 1))) {
      this->eM->dynMat->getRBMs(*Z);
    }
    // otherwise, use the geometric modes
    else {
      this->eM->rigidBodyModes->getRBMs(*Z);
    }
  }

  // ... initialize the first set of vectors
  this->probDesc->initQ((*Z)+this->nrmod,nsub);

  // ... copy the intialized vectors to vector set Q
  int i, j;
  for(i=0; i<subSpaceSize; ++i) (*Q)[i] = (*Z)[i];

  int lengthNsub = (nsub*(nsub+1))/2;

  double *mu    = new double[lengthNsub];
  double *kappa = new double[lengthNsub];
  
  // Orthonormalize inital eigenvector guess to rbm
  // Multiply mass by Z:   [Q] = [M]*[Z];
  for (j = 0; j < this->nrmod; ++j) {
    for (i = 0; i < j; ++i) { double s = (*Q)[i]*(*Z)[j]; (*Z)[j].linAdd(-s, (*Z)[i]); } // MGS 
    this->eM->M->mult((*Z)[j],(*Q)[j]);
    double coef = 1.0/sqrt((*Z)[j]*(*Q)[j]);
    (*Z)[j] *= coef;
    (*Q)[j] *= coef;
  }

  // orthonormalize. Make R^t M Z = 0 then Z-> M Z so that R^t Z = 0 which
  // is necessary for the reSolve to have a meaning
  this->ortho((*Q), (*Z), nsub, this->nrmod);

  // compute Kappa, the coarse stiffness matrix 
  for (j = this->nrmod; j < subSpaceSize; ++j)
    this->eM->refK->mult((*Z)[j],(*Q)[j]);

  int index = 0;
  for (j = this->nrmod; j < subSpaceSize; ++j)
    for (i = this->nrmod; i <= j; ++i)
      kappa[index++] = (*Q)[i]*(*Z)[j];

  // initialize mu, the coarse mass matrix
  for (j = this->nrmod; j < subSpaceSize; ++j)
     this->eM->M->mult((*Z)[j],(*Q)[j]);
  
  index = 0;
  for (j = this->nrmod; j < subSpaceSize; ++j)
    for(i = this->nrmod; i <= j; ++i)
      mu[index++] = (*Q)[i]*(*Z)[j];
  
  FullSquareMatrix xx(nsub);
  this->getJacobi(kappa, mu, xx, (*subVal).data()+this->nrmod, nsmax, nsub, tolJac);

  for (j = this->nrmod; j < subSpaceSize; ++j)
    (*Q)[j] = (*Z)[j];

  for (j = 0; j < nsub; ++j)  {
    (*Z)[j+this->nrmod].zero();
    for (i = 0; i < nsub; ++i)
      (*Z)[j+this->nrmod].linAdd(xx[j][i],(*Q)[i+this->nrmod]);
  }
  
  // Now, Z is M-orthogonal and is the inital guess
/* 
  // Solve standard eigenvalue prob. to initialize residual
  for (j = this->nrmod; j < subSpaceSize; ++j)
    this->eM->refK->mult((*Z)[j],(*Q)[j]);
  index = 0;
  for (j = this->nrmod; j < subSpaceSize; ++j)
    for (i = this->nrmod; i <= j; ++i)
      kappa[index++] = (*Q)[j]*(*Z)[i];

  // initialize mu as the identity matrix
  for (j = 0; j < lengthNsub; ++j)
    mu[j] = 0.0;
  for (j = 1; j <= nsub; ++j)  {
    index = (j*(j+1))/2 - 1;
    mu[index] = 1.0;
  }
  getJacobi(kappa, mu, xx, (*subVal).data()+this->nrmod, nsmax, nsub, tolJac); 

  // expand the eigenvectors
  for (j = this->nrmod; j < subSpaceSize; ++j)
    (*Q)[j] = (*Z)[j];

  for (j = 0; j < nsub; ++j)  {
    (*Z)[j+this->nrmod].zero();
    for (i = 0; i < nsub; ++i)
      (*Z)[j+this->nrmod].linAdd(xx[j][i],(*Q)[i+this->nrmod]);
  }   
*/      
  // form residual
  double refNorm;
  VecSet *R = new VecSet(subSpaceSize, this->probDesc->solVecInfo());
  for (j = 0; j < this->nrmod; j++)
    (*R)[j] = (*Z)[j];

  for (j = this->nrmod; j < subSpaceSize; ++j)  {
    this->eM->refK->mult((*Z)[j], (*R)[j]);
    refNorm = (*R)[j].norm();
    this->eM->M->mult((*Z)[j],(*Q)[j]);
    //filePrint(stderr, "zmz(%d) = %e, zkz(%d) = %e\n", j, (*Z)[j]*(*Q)[j], j, (*Z)[j]*(*R)[j]);
    (*R)[j] -= (*subVal)[j] * (*Q)[j];
    if(verboseFlag) filePrint(stderr, "Initial Norm %d: %e (%e), lam = %e\n", j, (*R)[j].norm()/refNorm, refNorm, (*subVal)[j]);
  }
    
  // do 1st iteration 
  delete [] mu; 
  delete [] kappa;
  lengthNsub = (2*nsub*(2*nsub+1))/2;
  mu    = new double[lengthNsub];
  kappa = new double[lengthNsub];

  // solve preconditioned system and orthonormalize
  this->eM->dynMat->reSolve(nsub,(*R)+this->nrmod);
  this->ortho((*Q), (*R), nsub, this->nrmod);
 
  // form Kappa and Mu
  VecSet *tmpR = new VecSet(subSpaceSize, this->probDesc->solVecInfo());
  for (j = this->nrmod; j < subSpaceSize; ++j)  {
    this->eM->refK->mult((*Z)[j],(*Q)[j]);
    this->eM->refK->mult((*R)[j],(*tmpR)[j]);
  }

  index = 0;
  for (j = this->nrmod; j < subSpaceSize; ++j)
    for (i = this->nrmod; i <= j; ++i)
      kappa[index++] = (*Q)[j]*(*Z)[i];

  for (j = this->nrmod; j < subSpaceSize; ++j)  {
    for (i = this->nrmod; i < subSpaceSize; i++)
      kappa[index++] = (*tmpR)[j]*(*Z)[i];
    for (i = this->nrmod; i <= j; ++i)
      kappa[index++] = (*tmpR)[j]*(*R)[i];
  }

  // initialize mu, the coarse mass matrix
  for (j = this->nrmod; j < subSpaceSize; ++j)  {
     this->eM->M->mult((*Z)[j],(*Q)[j]);
     this->eM->M->mult((*R)[j],(*tmpR)[j]);
  }

  index = 0;
  for (j = this->nrmod; j < subSpaceSize; ++j)
    for(i = this->nrmod; i <= j; ++i)
      mu[index++] = (*Q)[j]*(*Z)[i];   

  for (j = this->nrmod; j < subSpaceSize; ++j)  {
    for(i = this->nrmod; i < subSpaceSize; ++i)
      mu[index++] = (*tmpR)[j]*(*Z)[i];   
    for(i = this->nrmod; i <= j; ++i)
      mu[index++] = (*tmpR)[j]*(*R)[i];   
  }
 
  xx.setSize(2*nsub);
  this->getJacobi(kappa, mu, xx, (*subVal).data()+this->nrmod, nsmax, 2*nsub, tolJac);
/*
  // form P_k, the search direction
  VectorSet *Pk = new VectorSet(subSpaceSize, this->probDesc->solVecInfo()); 
  for (j = 0; j < nsub; j++)  {
    (*Pk)[j+this->nrmod].zero();
    for (i = 0; i < nsub; i++)
      (*Pk)[j+this->nrmod].linAdd(xx[j][nsub+i], (*R)[i+this->nrmod]);
  }
*/ 
  // form new eigenvector estimate
  for (j = this->nrmod; j < subSpaceSize; ++j)
    (*Q)[j] = (*Z)[j];
    
  for (j = 0; j < nsub; ++j)  {
    (*Z)[j+this->nrmod].zero();
    for (i = 0; i < nsub; ++i)
      (*Z)[j+this->nrmod].linAdd(xx[j][i],(*Q)[i+this->nrmod]);
    for (i = 0; i < nsub; ++i)
      (*Z)[j+this->nrmod].linAdd(xx[j][nsub+i],(*R)[i+this->nrmod]);
    //(*Z)[j+this->nrmod] += (*Pk)[j+this->nrmod];
  } 
    
  // form residual
  for (j = this->nrmod; j < subSpaceSize; ++j)  {
    this->eM->refK->mult((*Z)[j], (*R)[j]);
    refNorm = (*R)[j].norm();
    this->eM->M->mult((*Z)[j],(*Q)[j]);
    //filePrint(stderr, "zmz(%d) = %e, zkz(%d) = %e\n", j, (*Z)[j]*(*Q)[j], j, (*Z)[j]*(*R)[j]);
    (*R)[j] -= (*subVal)[j] * (*Q)[j];
    if(verboseFlag) filePrint(stderr, "1st Iter, Norm %d: %e (%e), lam = %e\n", j, (*R)[j].norm()/refNorm, refNorm, (*subVal)[j]);
  }
/*
  delete [] mu; 
  delete [] kappa;
  lengthNsub = (3*nsub*(3*nsub+1))/2;
  mu    = new double[lengthNsub];
  kappa = new double[lengthNsub];    
  xx.setSize(3*nsub);
*/
  int iter;  
  //VecSet *tmpP = new VecSet(subSpaceSize, this->probDesc->solVecInfo());
  for (iter = 0; iter < nsmax; ++iter) {
    this->eM->dynMat->reSolve(nsub, (*R)+this->nrmod);
    this->ortho(*Q, *R, nsub, this->nrmod);
 
    // form Kappa and Mu
    for (j = this->nrmod; j < subSpaceSize; ++j)  {
      this->eM->refK->mult((*Z)[j],(*Q)[j]);
      this->eM->refK->mult((*R)[j],(*tmpR)[j]);
      //this->eM->refK->mult((*Pk)[j],(*tmpP)[j]);
    }

    index = 0;
    for (j = this->nrmod; j < subSpaceSize; ++j)
      for (i = this->nrmod; i <= j; ++i)
        kappa[index++] = (*Q)[j]*(*Z)[i];

    for (j = this->nrmod; j < subSpaceSize; ++j)  {
      for (i = this->nrmod; i < subSpaceSize; ++i)
        kappa[index++] = (*tmpR)[j] * (*Z)[i];
      for (i = this->nrmod; i <= j; ++i)
        kappa[index++] = (*tmpR)[j]*(*R)[i];
    }
/*
    for (j = this->nrmod; j < subSpaceSize; ++j)  {
      for (i = this->nrmod; i < subSpaceSize; ++i)
        kappa[index++] = (*tmpP)[j] * (*Z)[i];
      for (i = this->nrmod; i < subSpaceSize; ++i)
        kappa[index++] = (*tmpP)[j] * (*R)[i];
      for (i = this->nrmod; i <= j; ++i)
        kappa[index++] = (*tmpP)[j]*(*Pk)[i];
    }
*/
    // initialize mu, the coarse mass matrix
    for (j = this->nrmod; j < subSpaceSize; ++j)  {
      this->eM->M->mult((*Z)[j],(*Q)[j]);
      this->eM->M->mult((*R)[j],(*tmpR)[j]);
      //this->eM->M->mult((*Pk)[j],(*tmpP)[j]);
    }

    index = 0;
    for (j = this->nrmod; j < subSpaceSize; ++j)
      for(i = this->nrmod; i <= j; ++i)
        mu[index++] = (*Q)[j]*(*Z)[i];   
    for (j = this->nrmod; j < subSpaceSize; ++j)  {
      for (i = this->nrmod; i < subSpaceSize; ++i)
        mu[index++] = (*tmpR)[j] * (*Z)[i];
      for (i = this->nrmod; i <= j; ++i)
        mu[index++] = (*tmpR)[j]*(*R)[i];
    }
/*
    for (j = this->nrmod; j < subSpaceSize; ++j)  {
      for (i = this->nrmod; i < subSpaceSize; ++i)
        mu[index++] = (*tmpP)[j] * (*Z)[i];
      for (i = this->nrmod; i < subSpaceSize; ++i)
        mu[index++] = (*tmpP)[j] * (*R)[i];
      for (i = this->nrmod; i <= j; ++i)
        mu[index++] = (*tmpP)[j]*(*Pk)[i];
    }
*/	
    // returns sorted eigenvalues
    this->getJacobi(kappa, mu, xx, (*subVal).data()+this->nrmod, nsmax, 2*nsub, tolJac);
/*
    // form new search direction 
    for (j = this->nrmod; j < subSpaceSize; j++)
      (*tmpP)[j] = (*Pk)[j];

    for (j = 0; j < nsub; j++)  {
      (*Pk)[j+this->nrmod].zero();
      for (i = 0; i < nsub; i++)
        (*Pk)[j+this->nrmod].linAdd(xx[j][nsub+i], (*R)[this->nrmod+i]);
      for (i = 0; i < nsub; i++)
        (*Pk)[j+this->nrmod].linAdd(xx[j][2*nsub+i], (*tmpP)[this->nrmod+i]);
    }
*/
    // form new eigenvector estimate
    for (j = this->nrmod; j < subSpaceSize; ++j)
      (*Q)[j] = (*Z)[j];
    
    for (j = 0; j < nsub; ++j)  {
      (*Z)[j+this->nrmod].zero() ;  
      //(*Z)[j+this->nrmod] = (*Pk)[this->nrmod+j];
      for (i = 0; i < nsub; ++i)
        (*Z)[j+this->nrmod].linAdd(xx[j][i],(*Q)[i+this->nrmod]);
      for (i = 0; i < nsub; ++i)
        (*Z)[j+this->nrmod].linAdd(xx[j][nsub+i],(*R)[i+this->nrmod]);
    } 
   
    // form residual and test convergence
    int hasCon = 1;
    for (j = this->nrmod; j < subSpaceSize; ++j)  {
      this->eM->refK->mult((*Z)[j], (*R)[j]);
      this->eM->M->mult((*Z)[j],(*Q)[j]);
      //filePrint(stderr, "zmz(%d) = %e, zkz(%d) = %e\n", j, (*Z)[j]*(*Q)[j], j, (*Z)[j]*(*R)[j]);
      refNorm = (*R)[j].norm();
      (*R)[j] -= (*subVal)[j] * (*Q)[j];
      if(verboseFlag) filePrint(stderr, "Loop Iter Norm %d: %e (%e), lam = %e\n", j, (*R)[j].norm()/refNorm, refNorm, (*subVal)[j]);

      if ((*R)[j].sqNorm()/refNorm > tolEig)  {
//	filePrint(stderr, "Norm %d: %e (%e)\n", j, (*R)[j].sqNorm()/refNorm, refNorm);
        hasCon = 0;
      }

    }

    if (hasCon) break;
    if(verboseFlag) filePrint(stderr, "Done iteration %d\n", iter);
  } 
  
  // M-orthonormalize Z
  for (i = this->nrmod; i < subSpaceSize; ++i)
    this->eM->M->mult((*Z)[i], (*Q)[i]);
  
  for (j = this->nrmod; j < subSpaceSize; ++j) {
    for (i = this->nrmod; i < j; ++i) { double s = (*Q)[i]*(*Z)[j]; (*Z)[j].linAdd(-s, (*Z)[i]); } // MGS
    double coef = 1.0/sqrt((*Z)[j]*(*Q)[j]);
    (*Z)[j] *= coef;
  }

  // save converged eigenvalues and eigenvectors
  filePrint(stderr, "... Saving eigen solution for %d eigs ...\n", this->totalEig);
  for(i = 0; i < this->totalEig; ++i) {
    filePrint(stderr, "%d: %e   (%e)\n", i, (*subVal)[i], (*Z)[i].norm());
    (*this->eigVal)[i] = (*subVal)[i];
    (*this->eigVec)[i] = (*Z)[i];
  }

  // ... Output results
  if(this->probDesc->getFilter()) {
    VecSet rModeData(modeData.numModes,this->probDesc->solVecInfo());
    Vector rModeVal(modeData.numModes);
    this->pickMostCorrelatedModes(rModeVal,rModeData);
  }
  else this->postProcessor->eigenOutput(*this->eigVal,*this->eigVec);

  // ... Print timers
  this->probDesc->printTimers(this->eM->dynMat);

  // ... Free memory
  if(Q)      { delete Q;      Q      = 0; }
  if(Z)      { delete Z;      Z      = 0; }
  if(subVal) { delete subVal; subVal = 0; }

/*
  // ... Output results
  this->postProcessor->eigenOutput(*this->eigVal,*this->eigVec);
*/

}
//------------------------------------------------------------------------------

template <class EigOps, class VecType, class VecSet,
    class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet,
         PostProcessor, ProblemDescriptor>::performQR(Vector *eigVal, VecSet *eigVec, int totalEig)
{
   if( strcmp(solInfo.eigenvaluename, "") != 0) {
     const char* eigenvaluename = solInfo.eigenvaluename;
     std::ofstream eigenout(eigenvaluename, std::ios::out);
     if(!eigenout) { std::cerr << " ... Error: cannot open file " << eigenvaluename << std::endl;   exit(-1); } 
     std::ostringstream s;
     const double pi = 3.141592653589793;
     for(int i=0; i<totalEig; ++i) { eigenout << sqrt((*eigVal)[i])/(2.0*pi) << "\n"; }
     eigenout.close();
   }
#ifdef USE_EIGEN3
   if( strcmp(solInfo.xmatrixname, "") == 0 || strcmp(solInfo.qmatrixname, "") == 0 || strcmp(solInfo.rmatrixname,"") == 0) {
     filePrint(stderr, " *** ERROR: xmatrix, qmatrix and rmatrix keywords must be specified in input file\n");
   } else {
     Eigen::MatrixXd Xeigen(probDesc->solVecInfo(),totalEig);
     for(int j=0; j<totalEig; ++j) {
       VecType col = (*eigVec)[j];
       for(int i=0; i<probDesc->solVecInfo(); ++i) Xeigen(i,j) = col[i];
     }

     Eigen::HouseholderQR<Eigen::MatrixXd> qr(Xeigen);
     const typename Eigen::HouseholderQR<Eigen::MatrixXd>::HouseholderSequenceType &Q = qr.householderQ();

//     Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Xeigen);
//     const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> &P = qr.colsPermutation();
//    const typename Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::HouseholderSequenceType &Q = qr.householderQ();

     Eigen::MatrixXd Qupdated = Eigen::MatrixXd(Q).block(0,0,probDesc->solVecInfo(),totalEig);

     Eigen::MatrixXd Reigen = (Qupdated.transpose()*Xeigen).block(0,0,totalEig,totalEig);
/*
     Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0);
     std::cerr << Reigen.format(HeavyFmt) << std::endl;
     std::cerr << Xeigen.rows() << "x" << Xeigen.cols() << std::endl;
     std::cerr << Reigen.rows() << "x" << Reigen.cols() << std::endl;
*/
   // ... Output results
     this->postProcessor->eigenQROutput(Xeigen,Qupdated,Reigen);
   }
#else
   filePrint(stderr," *** ERROR: performQR() needs eigen3\n");
   exit(-1);
#endif
}
//------------------------------------------------------------------------------

template <class EigOps, class VecType, class VecSet, 
	  class PostProcessor, class ProblemDescriptor>
void
SubSpaceSolver< EigOps, VecType, VecSet, 
	       PostProcessor, ProblemDescriptor>::solve()
{
 //... Check if multiple-shift eigen-analysis is requested
 if(solInfo.doEigSweep) {
   filePrint(stderr," *** ERROR: Multiple-shift eigenvalue calculation not implemented for SubSpaceSolver. Please use ARPACK.\n");
   exit(-1);
 }

 // ... Set up
 this->setUp();
 if(solInfo.printMatLabExit) return; // just compute and print M and K
 if(this->origSubSize == 0) return; // just compute and print the rigid body modes

 // ... Initialize
 initialize();

 // ... Check if nsub (# of flexible modes) is less than zero
 if(nsub < 0) {
   filePrint(stderr," *** ERROR in routine SubSpaceSolver::solve()    ***\n");
   filePrint(stderr," *** subspaceSize =  %d                       ***\n", subSpaceSize);
   filePrint(stderr," *** while # of Rigid Body Modes (this->nrmod) = %d ***\n",this->nrmod);
   filePrint(stderr," *** subspaceSize must be bigger than this->nrmod   ***\n");
   return;
 }

 // ... Eigenvalues filtering
 bool filterEigen = solInfo.filtereig; // force "filtering" the eigenvalues: set this->eigVal[i] to 0.0
                                       // if eigVal[i]<0.0 or |eigVal[i]|<tol

 if(filterEigen) filePrint(stderr," ... Eigenvalue ''filtering'' enabled with tol = 1.E-6 ...\n");

 // ... check that subspace is big enough
 if(subSpaceSize < this->nrmod) {
   this->probDesc->error(subSpaceSize,this->nrmod);
   return;
 }

 // ... Store the rbms in the first this->nrmod vectors of VecSet Z
 if (this->nrmod) {
   // if GRBM is not specified, or if a shift is specified, or if the solver is skyline/sparse/feti then get rbms from solver
   if(domain->solInfo().rbmflg == 0 || geoSource->shiftVal() != 0.0 ||
      isFeti(domain->solInfo().solvercntl->type) ||
      ((domain->solInfo().solvercntl->subtype == 0 || domain->solInfo().solvercntl->subtype == 1))) {
     this->eM->dynMat->getRBMs(*Z);
   }
   // otherwise, use the geometric modes
   else {
     this->eM->rigidBodyModes->getRBMs(*Z);
   }
 }

 // ... initialize the first set of vectors
 this->probDesc->initQ((*Z)+this->nrmod,nsub);

 // ... copy the intialized vectors to vector set Q
 int i;
 for(i=0; i<subSpaceSize; ++i) (*Q)[i] = (*Z)[i];

 int lengthNsub = (nsub*(nsub+1))/2;

 double *mu    = new double[lengthNsub];
 double *kappa = new double[lengthNsub];

 // Multiply mass by Z:   [Q] = [M]*[Z];
 int j;
 for(j=0; j<this->nrmod; ++j) {
   for(i=0; i<j; ++i) { double s = (*Q)[i]*(*Z)[j]; (*Z)[j].linAdd(-s, (*Z)[i]); } // MGS
   this->eM->M->mult((*Z)[j],(*Q)[j]);
   double coef = 1.0/sqrt((*Z)[j]*(*Q)[j]);
   (*Z)[j] *= coef;
   (*Q)[j] *= coef;
 }

 // orthonormalize. Make R^t M Z = 0 then Z-> M Z so that R^t Z = 0 which
 // is necessary for the reSolve to have a meaning
 this->ortho((*Q), (*Z), nsub, this->nrmod);

 for(j=this->nrmod; j<subSpaceSize; ++j) {
   this->eM->M->mult((*Z)[j],(*Q)[j]);
   (*Z)[j] = (*Q)[j];
 }
 
 FullSquareMatrix xx(nsub); xx.zero();

 // Subspace Iteration Loop
 int iter;
 VecSet *tmpVec = new VecSet(1, this->probDesc->solVecInfo());
 VecSet *MZ = new VecSet(1, this->probDesc->solVecInfo());
 VecSet *residual = new VecSet(1, this->probDesc->solVecInfo());
 for(iter=0; iter<nsmax; ++iter) {

   for(i = this->nrmod; i<subSpaceSize; ++i)
     (*Q)[i] = (*Z)[i];

   this->eM->dynMat->reSolve(nsub,(*Z)+this->nrmod);
   this->ortho((*Q), (*Z), nsub, this->nrmod);

   if (explicitK)  {
     for(j=this->nrmod; j<subSpaceSize; ++j)
       this->eM->refK->mult((*Z)[j],(*Q)[j]);
   }

   int index = 0;
   for(j=this->nrmod; j<subSpaceSize; ++j)
     for(i = this->nrmod; i<=j; ++i)
       kappa[index++] = (*Q)[i]*(*Z)[j];

   for(j=this->nrmod; j<subSpaceSize; ++j)
     this->eM->M->mult((*Z)[j],(*Q)[j]);
     
   index = 0;
   for(j=this->nrmod; j<subSpaceSize; ++j)
     for(i = this->nrmod; i<= j; ++i)
       mu[index++] = (*Q)[i]*(*Z)[j];

   // returns sorted eigenvalues
   this->getJacobi(kappa,mu,xx,(*subVal).data()+this->nrmod,nsmax,nsub,tolJac);

   // test for convergence
   int hasCon = 1;
   double maxErr = 0;
   if (explicitK)  {
     // double refNorm;
     for (j=0; j<nsub; ++j)  {
       (*tmpVec)[0].zero();
       for(i = 0; i<nsub; ++i)
         (*tmpVec)[0].linAdd(xx[j][i],(*Z)[i+this->nrmod]);
       this->eM->refK->mult((*tmpVec)[0], (*residual)[0]);
       this->eM->M->mult((*tmpVec)[0], (*MZ)[0]);
       // refNorm = (*residual)[0].norm();
       (*residual)[0] -= (*subVal)[j+this->nrmod] * (*MZ)[0];
     }
   }
   else  {
     for (i = this->nrmod; i < this->totalEig; ++i) {
       double err = std::abs(((*subOld)[i] - (*subVal)[i])/(*subVal)[i]);
       if(err > tolEig) hasCon = 0;
       if(err > maxErr) maxErr = err;
       //if ( std::abs(((*subOld)[i] - (*subVal)[i])/(*subVal)[i]) > tolEig) hasCon = 0;
       (*subOld)[i] = (*subVal)[i];
     }
   }

   if((domain->solInfo().solvercntl->fetiInfo.numPrint() > 0) && (iter % domain->solInfo().solvercntl->fetiInfo.numPrint() == 0) && verboseFlag)
     filePrint(stderr,"Subspace Iteration %d (Max %d): Error = %e (Tol %e)\n", iter+1, nsmax, maxErr, tolEig);

   if (hasCon) break;

   if(iter < nsmax-1)
     for(j=0; j<nsub; ++j)  {
        (*Z)[j+this->nrmod].zero();
        for(i = 0; i<nsub; ++i)
           (*Z)[j+this->nrmod].linAdd(xx[j][i],(*Q)[i+this->nrmod]);
     }
 }

 for(i=0; i<this->nrmod; ++i) { (*Q)[i] = (*Z)[i]; }

 for(j=0; j<nsub; ++j) {
    (*Q)[j+this->nrmod].zero();
    for(i = 0; i<nsub;++i) (*Q)[j+this->nrmod].linAdd(xx[j][i],(*Z)[i+this->nrmod]);
 }

 // save converged eigenvalues and eigenvectors
 for(i=0; i<this->totalEig; ++i) {
    (*this->eigVal)[i] = (*subVal)[i]+geoSource->shiftVal(); // PJSA: modified for shifted eigen
    (*this->eigVec)[i] = (*Q)[i];
    if(filterEigen) { if((*this->eigVal)[i]<0 || fabs((*this->eigVal)[i])<1.E-6 ) (*this->eigVal)[i] = 0.0; }
 }


 // check this->eigVec/this->eigVal convergence
 if (explicitK)  {
   (*tmpVec)[0] = 0.0;
   double refNorm;
   for (i=0; i<this->totalEig; ++i)  {
     this->eM->M->mult( (*this->eigVec)[i], (*tmpVec)[0] );
     (*tmpVec)[0] *= (*this->eigVal)[i];
     this->eM->refK->mult( (*this->eigVec)[i], (*residual)[0] );
     refNorm = (*residual)[0].norm();
     (*residual)[0] -= (*tmpVec)[0];

     filePrint(stderr, "EigProb Norm %d: %e (%e)\n", i, (*residual)[0].norm()/refNorm, refNorm);
   }
 }

 if(solInfo.qrfactorization && !this->probDesc->getFilter()) this->performQR(this->eigVal,this->eigVec,this->totalEig);

 // ... Output results
 if(this->probDesc->getFilter()) {
   VecSet rModeData(modeData.numModes,this->probDesc->solVecInfo());
   Vector rModeVal(modeData.numModes);
   this->pickMostCorrelatedModes(rModeVal,rModeData);
   this->performQR(&rModeVal,&rModeData,modeData.numModes);
 } else this->postProcessor->eigenOutput(*this->eigVal,*this->eigVec);

 // ... Print timers
 this->probDesc->printTimers(this->eM->dynMat);

 // ... Free memory
  if(Q)      { delete Q;      Q      = 0; }
  if(Z)      { delete Z;      Z      = 0; }
  if(subVal) { delete subVal; subVal = 0; }
  if(subOld) { delete subOld; subOld = 0; }

/*
 // ... Output results
 this->postProcessor->eigenOutput(*this->eigVal,*this->eigVec);
*/

 delete [] mu; 
 delete [] kappa;
}

extern "C"      {
void _FORTRAN(cfjacobi)(double *,double *,double *, double *,int&,double&,int &);
void _FORTRAN(dspgv)(const int &, const char &, const char &,
                     const int &N, double *AP, double *BP, double *W,
                     double *Z, const int &LDZ, double *work, int &info);
void _FORTRAN(zhpgv)(const int &, const char &, const char &,
                     const int &N, complex<double> *AP, complex<double> *BP, double *W,
                     complex<double> *Z, const int &LDZ, complex<double> *work, double *work2, int &info);
void _FORTRAN(dggev)(const char &JOBVL, const char &JOBVR,
                     const int &N, double *AP, const int &LDA, double *BP, const int &LDB, 
                     double *ALPHAR, double *ALPHAI, double *BETA, 
                     double *VL, const int &LDVL, double *VR, const int &LDVR,
                     double *work, const int &LWORK, int &info);
}

template <class EigOps, class VecType, class VecSet,
          class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet,
             PostProcessor, ProblemDescriptor>
  ::getJacobi(double *kappa, double *mu, FullSquareMatrix &xx,
              double *eigVal, int nsmax, int subSpaceSize, double tolJac)
{
  int i,j;

  int info;
  int subType = solInfo.eigenSolverSubType;
  switch(subType) {
    case 0 : {
      double *work = new double[3*subSpaceSize];
      //_FORTRAN(cfjacobi)(kappa,mu,xx[0],eigVal,nsmax,tolJac,subSpaceSize);
      _FORTRAN(dspgv)(1, 'V', 'U', subSpaceSize, kappa, mu, eigVal, xx[0],
                      xx.dim(), work, info);
      if(info !=0 ) {
        //int N = subSpaceSize;
        //cerr << "N = " << N << std::endl;
        //cerr << "AP = "; for(int i=0; i<N*(N+1)/2; ++i) std::cerr << kappa[i] << " "; std::cerr << std::endl;
        //cerr << "BP = "; for(int i=0; i<N*(N+1)/2; ++i) std::cerr << mu[i] << " "; std::cerr << std::endl;
        //cerr << "LDZ = " << xx.dim() << std::endl;
        filePrint(stderr, "Error in dspgv: info = %d, N = %d\n", info, subSpaceSize);
      }
      delete [] work;
    }
    break;
    case 1 : {
      // PJSA: zhpgv for complex (is working)
      int N = subSpaceSize;
      complex<double> *work = new complex<double>[2*N-1];
      double *work2 = new double[3*N-2];
      complex<double> *a = new complex<double>[N*(N+1)/2];
      complex<double> *b = new complex<double>[N*(N+1)/2];
      GenFullSquareMatrix<complex<double> > zz(N);
      for(i=0; i<N*(N+1)/2; ++i) {
        a[i] = kappa[i];
        b[i] = mu[i];
      }
      _FORTRAN(zhpgv)(1, 'V', 'U', subSpaceSize, a, b, eigVal, zz[0],
                      zz.dim(), work, work2, info);
      if(info !=0)
        filePrint(stderr, "Error in zhpgv: %d\n", info);
      for(i=0; i<N; ++i)
        for(j=0; j<N; ++j) xx[i][j] = zz[i][j].real();
      delete [] a; delete [] b; delete [] work2; delete [] work;
    }
    break;
    case 2 : {
      // PJSA: dggev for indefinite (not working)
      int N = subSpaceSize;
      double *work = new double[8*N];
      double *a = new double[N*N];
      double *b = new double[N*N];
      double *alphai = new double[N]; 
      double *beta = new double[N];
      for(j=0; j<N; ++j)
        for(i=0; i<=j; ++i) {
          int I = i+1;
          int J = j+1;
          a[i+j*N] = a[j+i*N] = kappa[I+(J-1)*J/2-1];
          b[i+j*N] = b[j+i*N] = mu[I+(J-1)*J/2-1];
        }
      _FORTRAN(dggev)('N', 'V', N, a, N, b, N, eigVal, alphai, beta, xx[0], xx.dim(), xx[0], 
                      xx.dim(), work, 8*N, info);
      if(info !=0)
        filePrint(stderr, "Error in ddgev: %d\n", info);
      for(i=0; i<N; ++i) eigVal[i] /= beta[i];
      delete [] a; delete [] b; delete [] work; delete [] alphai; delete [] beta;
    }
    break;
  }
  // sort eigenvalues.
  int is = 1;
  while(is != 0) {
    is = 0;
    for(i=1; i<subSpaceSize; ++i) {
      if(eigVal[i] < eigVal[i-1] ) {
        is = 1;
        std::swap(eigVal[i-1], eigVal[i]);
        for(j=0; j<subSpaceSize; ++j) {
          std::swap( xx[i][j], xx[i-1][j] );
        }
      }
    }
  }

}

template <class EigOps, class VecType, class VecSet,
          class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet,
             PostProcessor, ProblemDescriptor>
::absoluteInnerproductNormalized(const VecType& v1, const VecType& v2, double &result)
{
  VecType normalizedV1(v1), normalizedV2(v2);
  double normv1 = normalizedV1.norm();
  double normv2 = normalizedV2.norm();
  normalizedV1 *= (1.0/normv1);
  normalizedV2 *= (1.0/normv2);
  result = std::abs(normalizedV1*normalizedV2);
}

template <class EigOps, class VecType, class VecSet,
          class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet,
             PostProcessor, ProblemDescriptor>
::pickMostCorrelatedModes(Vector &rModeVal, VecSet &rModeData)
{
  int i,j,k,selectedIndex[modeData.numModes];
  for(i=0; i<modeData.numModes; ++i) selectedIndex[i] = 0;
  double result, maxResult[modeData.numModes];
  bool isIncluded;
  VecSet vModeData(modeData.numModes,probDesc->solVecInfo());
  probDesc->convertModeDataToVecSet(vModeData);
  for(i=0; i<modeData.numModes; ++i) {
    maxResult[i] = 0.0;
    for(j=0; j<totalEig; ++j) {
      absoluteInnerproductNormalized(vModeData[i],(*eigVec)[j],result);
      if(result > maxResult[i]) {
//        isIncluded = false;
//        for(k=0; k<i; ++k) if(selectedIndex[k] == j) { isIncluded = true; }
//        if(!isIncluded) {
          maxResult[i] = result;
          selectedIndex[i] = j;
//        }
      }
    }
    filePrint(stderr," ... mode %d is paired with mode %d with normalized inner product of %e\n", i+1, selectedIndex[i]+1, maxResult[i]);
    rModeData[i] = (*eigVec)[selectedIndex[i]];
    rModeVal[i] = (*eigVal)[selectedIndex[i]];
  }
  this->postProcessor->eigenOutput(rModeVal, rModeData, modeData.numModes);
}

template <class EigOps, class VecType, class VecSet,
          class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet,
             PostProcessor, ProblemDescriptor>
::ortho(VecType *v1, VecType *vr, int nsub, int nrbm)
{
  int i,j;
  for(j=0; j<nsub; ++j) {
    for(i=0; i<nrbm; ++i) {
      double s = v1[i]*vr[j+nrbm];
      vr[j+nrbm] -= s*vr[i]; // Vector -= operation
    }
  }
}

template <class EigOps, class VecType, class VecSet,
          class PostProcessor, class ProblemDescriptor>
void
EigenSolver< EigOps, VecType, VecSet,
             PostProcessor, ProblemDescriptor>
::ortho(VecSet& v1, VecSet& vr, int nsub, int nrbm)
{
  int i,j;
  for(j=0; j<nsub; ++j) {
    for(i=0; i<nrbm; ++i) {
      double s = v1[i]*vr[j+nrbm];
      vr[j+nrbm] -= s*vr[i]; // Vector -= operation
    }
  }
}



#ifdef USE_ARPACK
extern "C" {
#ifdef DISTRIBUTED
void _FORTRAN(pdsaupd)(MPI_Fint *COMM, int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
                       double* RESID, int* NCV, double* V, int* LDV, int* IPARAM, int* INPTR, double* WORKD,
                       double* WORKL, int* LWORKL, int* INFO);

void  _FORTRAN(pdseupd)(MPI_Fint *COMM, const int& REC , char* HOWMMY , int* SELECT, double* D , double* Z,
                        int* LDZ, double* SIGMA, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL, double* RESID,
                        int* NCV, double* V, int* LDV, int* IPARAM, int* INPTR, double* WORKD, double* WORKL,
                        int* LWORKL, int* INFO);
#else
void _FORTRAN(dsaupd)(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL, double* RESID,
                      int* NCV, double* V, int* LDV, int* IPARAM, int* INPTR, double* WORKD, double* WORKL,
                      int* LWORKL, int* INFO);

void  _FORTRAN(dseupd)(const int& REC , char* HOWMMY , int* SELECT, double* D , double* Z, int* LDZ, double* SIGMA,
                       char* BMAT, int* N, char* WHICH, int* NEV, double* TOL, double* RESID,
                       int* NCV, double* V, int* LDV, int* IPARAM, int* INPTR, double* WORKD, double* WORKL,
                       int* LWORKL, int* INFO);
#endif
}

template <class EigOps, class VecType, class VecSet, 
	  class PostProcessor, class ProblemDescriptor>
void
SymArpackSolver< EigOps, VecType, VecSet, 
	         PostProcessor, ProblemDescriptor>::solve()
{

  // for multiple shift eigen analysis
  if (solInfo.doEigSweep) {
    if (solInfo.lbound == 0.0) geoSource->resetShift(solInfo.lbound);
    else geoSource->setShift(solInfo.lbound);
  }

  //... Set up
  this->setUp();
  if (solInfo.printMatLabExit) return; // just compute and print M and K
  if (this->origSubSize == 0) return; // just compute and print the rigid body modes

  filePrint(stderr, " ... Implicitly Restarted Arnoldi Method using Arpack ...\n");

  this->probDesc->getSubSpaceInfo(subSpaceSize,nsmax,tolEig,tolJac,explicitK);
  if (explicitK)
    filePrint(stderr, " ... Using Explicit Stiffness Matrix...\n");

  // ... Some declarations
  bool sconv = false;
  int nev    = this->totalEig - this->nrmod; // number of eigenvalues searched for
  int nmodes = nev;
  //int nmodes , nev;
  char which[4];
  double* resid; 
  double* workd;
  double* workl; 
  double* LanVects;
  double* RitzVects;

  bool reortho     = true; // force double M-orthogonalization (like in Salinas): CGS ~> ICGS
  bool printInfo   = true;
  bool filterEigen = solInfo.filtereig; // force "filtering" the eigenvalues: set this->eigVal[i] to 0.0
                                                  // if this->eigVal[i]<0.0 or |this->eigVal[i]|<tol
  double newShift = 0.0;
  const double pi = 3.141592653589793;
  int counter = 0, diffEig = 0, convEig = 0;
  int myTotalEig = 0;
  std::list<double> TotEigVal;
  int i,j;

  //--------------------Multiple shift eigen analysis------------------------

  if (solInfo.doEigSweep)
    filePrint(stderr, " ... Multiple-Shift Eigen Analysis  ... \n");

  while (!sconv) {

    if (solInfo.doEigSweep)  {
      if (solInfo.nshifts > 1) { // Apply multiple shifts
        counter++;
        sprintf(which,"LA"); //"LA" --> compute the NEV largest (algebraic) eigenvalues
        nmodes = (counter == solInfo.nshifts) ? nev-myTotalEig : nev/solInfo.nshifts+1; 
        if (diffEig != 0 && counter != solInfo.nshifts) { nmodes += diffEig; diffEig = 0; }
        subSpaceSize = 2*nmodes;
        if (newShift != 0.0) {
          rebuildSolver(newShift);
          this->nrmod = this->eM->dynMat->numRBM();
        }
        //if (geoSource->shiftVal() > 0.0) { nev += this->nrmod; this->nrmod = 0; }
/* PJSA: shifted matrix can be singular
        if (geoSource->shiftVal() > 0.0) { this->nrmod = 0; }
*/
      }
      else if (solInfo.lbound >= 0.0 && solInfo.ubound > 0.0) {
        sprintf(which,"LA");
        //nmodes = (solInfo.neigps) ? solInfo.neigps+1 : 51;
        nmodes = solInfo.neigps+1;
        //if (diffEig != 0) { nmodes += diffEig; diffEig = 0; }
        subSpaceSize = 2*nmodes;
        if (newShift != 0.0) {
          rebuildSolver(newShift);
          this->nrmod = this->eM->dynMat->numRBM();
        }
        //if (geoSource->shiftVal() > 0.0) { nev += this->nrmod; this->nrmod = 0; }
/* PJSA: shifted matrix can be singular
        if (geoSource->shiftVal() > 0.0) { this->nrmod = 0; }
*/
      } 
      else if (solInfo.nshifts <= 1) {
        filePrint(stderr," *** WARNING: Number of shifts = %d. Performing one single-shift eigen calculation.\n", solInfo.nshifts);
        sprintf(which,"LA");
        sconv = true;
      }
    
      // ... Construct vector and vectorset for eigenvalues and eigenvectors

      this->eigVal = new Vector(nmodes+this->nrmod);
      this->eigVec = new VecSet(nmodes+this->nrmod, this->probDesc->solVecInfo());

      this->eigVal->zero();

    }
    else { // one eigen calculation
      sconv = true;
      if(strcmp("",solInfo.which) == 0) {
        if(geoSource->shiftVal() <= 0.0) {
          sprintf(which,"LA");
          //cerr << "which is undefined, setting to LA\n";
        }
        else {
          sprintf(which,"BE");
          //cerr << "which is undefined, setting to BE\n";
        }
      } else
      sprintf(which,"%s",solInfo.which); // specifies which of the Ritz values of
                                         // OP to compute (see Arpack manual)
      // ... adjust subSpaceSize
      if(subSpaceSize <= this->totalEig) {
        subSpaceSize = std::min(2*(this->totalEig-this->nrmod),this->totalEig-this->nrmod+8)+this->nrmod;
      }
      if(subSpaceSize-this->nrmod > this->probDesc->solVecSize()) {
        subSpaceSize = this->probDesc->solVecSize()+this->nrmod;
      }
    }

    // ... Declaration and Initialization of Vector sets Q and Z
    LanVects = new double[subSpaceSize*this->probDesc->solVecSize()]; // Arpack assumes Lanczos & Ritz vectors
    RitzVects= new double[subSpaceSize*this->probDesc->solVecSize()]; // are stored continuously in memory

    Q = new VecSet(subSpaceSize, this->probDesc->solVecInfo(), LanVects, false);
    Z = new VecSet(subSpaceSize, this->probDesc->solVecInfo(), RitzVects, false);

    // ... nsub = # of flexible modes
    nsub = subSpaceSize - this->nrmod;

    // ... Check if nsub (# of flexible modes) is less than zero
    if(nsub < 0) {
      filePrint(stderr," *** ERROR in routine SymArpackSolver::solve() ***\n");
      filePrint(stderr," *** subspaceSize = %3d                        ***\n", subSpaceSize);
      filePrint(stderr," *** while # of Rigid Body Modes (nrmod) = %3d ***\n", this->nrmod);
      filePrint(stderr," *** subspaceSize must be bigger than nrmod    ***\n");
      return;
    }

    // ... check that subspace is big enough
    if(subSpaceSize < this->nrmod) {
      this->probDesc->error(subSpaceSize,this->nrmod);
      return;
    }

    // ... Store the rbms in the first this->nrmod vectors of VecSet Z
    if(this->nrmod) {
      // if GRBM is not specified, or if a shift is specified, or if the solver is skyline/sparse/feti then get rbms from solver
      if(domain->solInfo().rbmflg == 0 || geoSource->shiftVal() != 0.0 || domain->solInfo().solvercntl->type == SolverSelection::Feti ||
         (domain->solInfo().solvercntl->type != SolverSelection::Feti && (domain->solInfo().solvercntl->subtype == 0 || domain->solInfo().solvercntl->subtype == 1))) {
        this->eM->dynMat->getRBMs(*Z);
      }
      // otherwise, just use the geometric modes
      else {
        this->eM->rigidBodyModes->getRBMs(*Z);
      }
    }
 
    // ... copy the this->nrmod Z vectors to vector set Q
    for(i=0; i<this->nrmod; ++i) (*Q)[i] = (*Z)[i];

    // ... M-Orthonormalize the first this->nrmod Z. Q will store M.Z
    for(j=0; j<this->nrmod; ++j) {
      for(i=0; i<j; ++i) { double s = (*Q)[i]*(*Z)[j]; (*Z)[j].linAdd(-s, (*Z)[i]); } // MGS
      this->eM->M->mult((*Z)[j],(*Q)[j]);
      double coef = 1.0/sqrt((*Z)[j]*(*Q)[j]);
      (*Z)[j] *= coef;
      (*Q)[j] *= coef;
    }

    // Implicit Restarted Arnoldi Method (IRAM) loop using Arpack for solving A.x = lambda.M.x

    double shift = geoSource->shiftVal();

    if(printInfo) {
      if(shift != 0.0 && newShift == 0.0) filePrint(stderr," ... shift = %7.2e                ...\n",shift);
      if(reortho)     filePrint(stderr," ... Re-orthogonalization enabled (CGS -> ICGS) ...\n");
      if(filterEigen) filePrint(stderr," ... Eigenvalue ''filtering'' enabled with tol = 1.E-6 ...\n");
    }

    // using Arpack
    int nevold = this->nrmod;
    int ncv = nsub;  // number of Lanczos vectors generated at each iteration
    int nloc = this->probDesc->solVecSize();  // dimension of the eigenproblem
    int iparam[11], ipntr[11];
    int ido = 0, info = 0, iram = nevold*nloc;
    int lworkl = ncv * (ncv + 8); // this should be at least ncv^2 + 8*ncv (as explained in ARPACK manual) 
    resid = new double[nloc+1]; // residual vector
    workd = new double[3*nloc]; // distributed array to be used in the basic Arnoldi iteration
    workl = new double[lworkl]; // private array 
    char bmat[4], howmny[4];
    sprintf(bmat,"G");    // "G"  --> generalized eigenvalue problem A*x = labmda*B*x
    sprintf(howmny,"A");  // "A"  --> compute NEV Ritz vectors
    for(i=0; i<11; ++i) iparam[i] = 0;
    for(i=0; i<11; ++i) ipntr[i]  = 0;
    iparam[0] = 1;    // "exact" shift
    iparam[2] = (solInfo.maxArnItr) ? solInfo.maxArnItr : nsmax; // maxitr
    iparam[6] = solInfo.arpack_mode;  // Mode = 3 (default) shift-invert mode
                                      //          => A symmetric & M symmetric positive semi-definite
                                      // Mode = 4 buckling mode
    // note: for buckling mode, the geometric stiffness matrix KG takes the place of M and this can be indefinite.
    // for mode 4: OP = inv(K-sigma*KG)*K. The shift sigma must be non-zero

    for(i=0; i <= nloc ; i++) resid[i] = 0.0;
    for(i=0; i < 3*nloc; i++) workd[i] = 0.0;
    for(i=0; i < lworkl; i++) workl[i] = 0.0;

    VecType tmp(this->probDesc->solVecInfo()); // temporary vector used in M-orthogonalization
    tmp.zero();

    bool conv = false;
    while(!conv) {
#ifdef DISTRIBUTED
       MPI_Fint mpi_fint = MPI_Comm_c2f(*(structCom->getCommunicator()));
        _FORTRAN(pdsaupd)(&mpi_fint,
                         &ido, bmat, &nloc, which, &nmodes, &tolEig, resid,
                         &ncv, &LanVects[iram], &nloc, iparam, ipntr, workd, workl,
                         &lworkl, &info);
#else
        _FORTRAN(dsaupd)(&ido, bmat, &nloc, which, &nmodes, &tolEig, resid,
                         &ncv, &LanVects[iram], &nloc, iparam, ipntr, workd, workl,
                         &lworkl, &info);
#endif
      switch(ido) {
        case -1: // Performing y <- OP*M*x
        { 
#ifdef ARPACK_DEBUG
          std::cerr<<" perform y <- OP*M*x (1)\n"; 
#endif
          VecType x(this->probDesc->solVecInfo(),&workd[ipntr[0]-1],false);
          VecType y(this->probDesc->solVecInfo(),&workd[ipntr[1]-1],false);
          if(nevold) { // M-Orthogonalize x to old eigenvectors, then update z=Mx
            this->eM->M->mult(x,tmp); 
            for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; x.linAdd(-s,(*Z)[i]); } // CGS
            if(reortho) {
              this->eM->M->mult(x,tmp); 
              for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; x.linAdd(-s,(*Z)[i]); } // CGS
            }
          }
          this->eM->M->mult(x,y); 
          this->eM->dynMat->reSolve(y);
          //M-Orthogonalize y to old eigenvectors
          if(nevold) { // M-Orthogonalize y to old eigenvectors
            this->eM->M->mult(y,tmp);
            for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; y.linAdd(-s,(*Z)[i]); } // CGS 
            if(reortho) {
              this->eM->M->mult(y,tmp); 
              for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; y.linAdd(-s,(*Z)[i]); } // CGS
            }
          }
          break;
        }
        case 1:  // Performing y <- OP*M*x. z=M*x is already available
        {
#ifdef ARPACK_DEBUG
          std::cerr<<" perform y <- OP*M*x (2)\n";
#endif
          VecType x(this->probDesc->solVecInfo(),&workd[ipntr[0]-1],false);
          VecType y(this->probDesc->solVecInfo(),&workd[ipntr[1]-1],false);
          VecType z(this->probDesc->solVecInfo(),&workd[ipntr[2]-1],false);
          if(nevold) { // M-Orthogonalize x to old eigenvectors, then update z=Mx
            this->eM->M->mult(x,tmp);
            for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; x.linAdd(-s,(*Z)[i]); } // CGS
            if(reortho) {
              this->eM->M->mult(x,tmp); 
              for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; x.linAdd(-s,(*Z)[i]); } // CGS
            }
            this->eM->M->mult(x,z); 
          }
          this->eM->dynMat->solve(z,y);
          //M-Orthogonalize y to old eigenvectors
          if(nevold) { // M-Orthogonalize y to old eigenvectors
            this->eM->M->mult(y,tmp);
            for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; y.linAdd(-s,(*Z)[i]); } // CGS
            if(reortho) {
              this->eM->M->mult(y,tmp); 
              for(i=0; i<nevold; i++){ double s = (*Z)[i]*tmp; y.linAdd(-s,(*Z)[i]); } // CGS
            }
          }
          break;
        }  
        case 2: // Performing z <- M*x
        {  
#ifdef ARPACK_DEBUG
          std::cerr<<" perform z <- M*x\n";
#endif
          VecType x(this->probDesc->solVecInfo(),&workd[ipntr[0]-1],false);
          VecType z(this->probDesc->solVecInfo(),&workd[ipntr[1]-1],false);
          this->eM->M->mult(x,z);
          break;
        }
        default:
          conv = true;
      }
    }
  
    // Either Arnoldi iteration has converged or there is an error 
    if(info < 0){
      // PJSA: if this happens try increasing nsbspv
      filePrint(stderr," *** WARNING: Arpack error with _saupd, info = %d\n",info);
      filePrint(stderr," ... %d Arnoldi update iterations taken ...\n",iparam[2]);
      filePrint(stderr," ,,, current Arnoldi factorization size? %d ...\n",iparam[4]);
    } else if(nmodes > 0){
      //filePrint(stderr,"NLOC = %d, NEV = %d, TOL = %e, NCV = %d, INFO = %d\n", nloc, nmodes, tolEig, ncv, info);
      //filePrint(stderr,"IPARAM(3) = %d, IPARAM(5) = %d, IPARAM(9) = %d, IPARAM(10) = %d, IPARAM(11) = %d\n", iparam[2], iparam[4], iparam[8], iparam[9], iparam[11]);
      if(printInfo){
        filePrint(stderr," ... Number of Implicit Arnoldi update: %d ...\n",iparam[2]);
        filePrint(stderr," ... Number of linear (re-)solve: %d ...\n",iparam[8]);
      }
      int ierr;
      int* select = new int[ncv]; // for howmny = "A", select is used as a workspace for reordering the Ritz values
      for(i=0; i<ncv; i++) select[i] = 0;
      double* lambda = &(*this->eigVal)[nevold];

#ifdef DISTRIBUTED
        MPI_Fint mpi_fint = MPI_Comm_c2f(*(structCom->getCommunicator()));
        _FORTRAN(pdseupd)(&mpi_fint,
                         1, howmny, select, lambda, &RitzVects[iram], &nloc, &shift,
                         bmat, &nloc, which, &nmodes, &tolEig,
                         resid, &ncv, &LanVects[iram], &nloc, iparam, ipntr,
                         workd, workl, &lworkl, &ierr);
#else
        _FORTRAN(dseupd)(1, howmny, select, lambda, &RitzVects[iram], &nloc, &shift,
                         bmat, &nloc, which, &nmodes, &tolEig,
                         resid, &ncv, &LanVects[iram], &nloc, iparam, ipntr,
                         workd, workl, &lworkl, &ierr);
#endif

      if(ierr!=0)
        filePrint(stderr," *** WARNING: Arpack error with _seupd, info = %d\n",ierr);
      else {
        int nconv = iparam[4]; // number of converged Ritz values
        if(nconv<nmodes) {
          filePrint(stderr," *** WARNING: only %d/%d eigenvalues did converge.\n",nconv,nmodes);
          if (solInfo.doEigSweep) nmodes=nconv-1;
        }
        if (solInfo.doEigSweep) {
          double eigDiff;
          if (solInfo.nshifts > 1) {
            myTotalEig += (counter == solInfo.nshifts && solInfo.nshifts != 0) ? nmodes : (nmodes-1);
            if(myTotalEig == nev) { convEig = nmodes; sconv = true; }
          } else if (solInfo.lbound >= 0.0 && solInfo.ubound > 0.0) {
            myTotalEig += nmodes-1;
            for (i=0; i<nmodes; i++) {
              if((*this->eigVal)[i] < solInfo.lbound)
                filePrint(stderr,"*** WARNING: calculated eigenvalue (%f) below lower bound (%f)\n",(*this->eigVal)[i],solInfo.lbound);
              else if ((*this->eigVal)[i] > solInfo.ubound) {
                sconv = true;
                convEig = i;
                myTotalEig -= (nmodes-convEig-1);
                break;
              }
            }
          } else if (solInfo.nshifts <= 1)
            convEig = this->totalEig;

          // Calculate next shift
          if (!sconv) {
            for(i=nmodes-1+this->nrmod; i>=this->nrmod; i--){
              eigDiff = fabs((*this->eigVal)[i]-(*this->eigVal)[i-1])/(*this->eigVal)[i];
              if (eigDiff <= 1.0e-06) continue;
              else {
                diffEig = nmodes-1+this->nrmod-i;
                myTotalEig -= diffEig;
                convEig = i;
                newShift = ((*this->eigVal)[i]+(*this->eigVal)[i-1])/2.0;
                break;
              }
            }
          } 
        } 
        else 
          convEig = this->totalEig;

        int filtered = 0; // count number of "filtered" eigen frequencies

        for (i=0; i<convEig; i++) { // save converged eigenvalues and eigenvectors
          if(filterEigen) { if((*this->eigVal)[i]<0 || fabs((*this->eigVal)[i])<1.E-6 ) { (*this->eigVal)[i] = 0.0; filtered++; } }
          if (solInfo.doEigSweep) 
             TotEigVal.push_back((*this->eigVal)[i]);
          (*this->eigVec)[i] = (*Z)[i];
        }
      }
      if(select) { delete [] select; select = 0; }
    }
    if(this->probDesc->getFilter()) { 
      VecSet rModeData(modeData.numModes,this->probDesc->solVecInfo());
      Vector rModeVal(modeData.numModes);
      this->pickMostCorrelatedModes(rModeVal,rModeData);
    } else this->postProcessor->eigenOutput(*this->eigVal,*this->eigVec, convEig);


    // free locally allocated memory
    if(resid)    { delete [] resid    ; resid    = 0; }
    if(workd)    { delete [] workd    ; workd    = 0; }
    if(workl)    { delete [] workl    ; workl    = 0; }
    if(LanVects) { delete [] LanVects ; LanVects = 0; }
    if(RitzVects){ delete [] RitzVects; RitzVects= 0; }
  }

  if (solInfo.doEigSweep) {

    filePrint(stderr," --------------------------------------\n");
      
    // Print omega and frequency values to screen
    if(solInfo.buckling) {
      filePrint(stderr," Mode\tLambda\n");
      filePrint(stderr," --------------------------------------\n");
      int imode, maxmode = TotEigVal.size();
      for(imode=0; imode<maxmode; ++imode) {
        filePrint(stderr," %d\t%e\n",imode+1,TotEigVal.front());
        TotEigVal.pop_front();
      }
      filePrint(stderr," --------------------------------------\n");
    } else {
      filePrint(stderr," Mode\tOmega^2\t\tFrequency\n");
      filePrint(stderr," --------------------------------------\n");
      int imode, maxmode = TotEigVal.size();
      for(imode=0; imode<maxmode; ++imode) {
        filePrint(stderr," %d\t%e\t%e\n",imode+1,TotEigVal.front(),
                sqrt(TotEigVal.front())/(2.0*pi));
        TotEigVal.pop_front();
      }
      filePrint(stderr," --------------------------------------\n");
    }
  }

  // ... Print timers
  this->probDesc->printTimers(this->eM->dynMat);

  //... free memory
  if(Q) { delete Q;  Q = 0; }
  if(Z) { delete Z;  Z = 0; }
}

template <class EigOps, class VecType, class VecSet,
          class PostProcessor, class ProblemDescriptor>
void
SymArpackSolver< EigOps, VecType, VecSet,
                 PostProcessor, ProblemDescriptor>::rebuildSolver(double newShift)
{

  filePrint(stderr,"\n ... Rebuilding LHS with shift = %e ... \n", newShift);
  geoSource->setShift(newShift);  // set new shift

  if(this->eM) this->probDesc->reBuild( *this->eM );
  else this->probDesc->buildEigOps( *this->eM );

}


/* HB: to finish ...
template <class EigOps, class VecType, class VecSet,
          class PostProcessor, class ProblemDescriptor>
void
SubSpaceSolver< EigOps, VecType, VecSet,
               PostProcessor, ProblemDescriptor>::Morthogonalize(VecType& x, VecSet& V, int nV, int OrthoAlgo=0; VecType* tmp=0, VectSet* MV=0)
{
  if(nV) { // M-Orthogonalize x to old eigenvectors,
    switch(OrthoAlgo) {
     default: 
     case 0: // CGS & ICGS  
       if(!MV) { 
         VecType& z = (tmp) ? *tmp : *(new VectType(x)); 
         this->eM->M->mult(x,z);
         for(int i=0; i<nV; i++){ double s = V[i]*z; x.linAdd(-s,V[i]); } 
         if(reortho) {
           this->eM->M->mult(x,z);
           for(int i=0; i<nV; i++){ double s = Z[i]*z; x.linAdd(-s,V[i]); } 
         }
       } else {
         for(int i=0; i<nV; i++){ double s = MV[i]*x; x.linAdd(-s,V[i]); }  // actually this is CGS ...
         if(reortho) { for(int i=0; i<nV; i++){ double s = MV[i]*x; x.linAdd(-s,V[i]); } }
       } 
       break;
     case 1: // MGS & IMGS
       if(!MV) {
         VecType& z = (tmp) ? *tmp : *(new VectType(x));
         for(int i=0; i<nV; i++){ this->eM->M->mult(x,z); double s = V[i]*z; x.linAdd(-s,V[i]); }
         if(reortho) { for(int i=0; i<nV; i++){ this->eM->M->mult(x,z); double s = Z[i]*z; x.linAdd(-s,V[i]); } }
       } else {
         for(int i=0; i<nV; i++){ double s = MV[i]*x; x.linAdd(-s,V[i]); }
         if(reortho) { for(int i=0; i<nV; i++){ double s = MV[i]*x; x.linAdd(-s,V[i]); } }
       }
       break;
    }
  }      
}
*/
#endif
