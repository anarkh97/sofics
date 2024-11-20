#include <cstdio>
#include <algorithm>
#include <Utils.d/dbg_alloca.h>
#include <Paral.d/MDEigen.h>
#include <Paral.d/MDDynam.h>
#include <Driver.d/SubDomain.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Driver.d/Dynam.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/VectorSet.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Feti.d/DistrVectorSet.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/SolverInfo.h>
#include <Feti.d/Feti.h>
#include <Utils.d/DistHelper.h>

#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif

template<class Scalar>
GenMultiDomainEigen<Scalar>::GenMultiDomainEigen(Domain* d)
{
 domain = d;
#ifdef DISTRIBUTED
 decDomain = new GenDistrDomain<Scalar>(domain);
#else
 decDomain = new GenDecDomain<Scalar>(domain);
#endif
 times  = new StaticTimers;
}

template<class Scalar>
DistrInfo &
GenMultiDomainEigen<Scalar>::solVecInfo()
{
 return decDomain->solVecInfo();
}

template<class Scalar>
int
GenMultiDomainEigen<Scalar>::solVecSize()
{
 return decDomain->solVecInfo().len;
}


template<class Scalar>
void
GenMultiDomainEigen<Scalar>::preProcess()
{
 times->preProcess -= getTime();

 // Makes renumbering, connectivities and dofsets
 decDomain->preProcess();
 // Make all element's dofs
 MultiDomainOp mdop(&MultiDomainOp::makeAllDOFs, decDomain->getAllSubDomains());
#ifdef DISTRIBUTED
 execParal(decDomain->getNumSub(), &mdop, &MultiDomainOp::runFor);
#else
 threadManager->execParal(decDomain->getNumSub(), &mdop);
#endif
 times->preProcess += getTime();
}

template<class Scalar>
GenMultiDomainEigenPostProcessor<Scalar> *
GenMultiDomainEigen<Scalar>::getPostProcessor()
{
 return new GenMultiDomainEigenPostProcessor<Scalar>(decDomain, solver, times);
}

template<class Scalar>
void 
GenMultiDomainEigen<Scalar>::buildEigOps(MDDynamMat &dMat)
{
 times->getFetiSolverTime -= getTime(); // PJSA 5-25-05
 decDomain->buildOps(dMat, 0.0, 0.0, 1.0);
 solver = dMat.dynMat;
 times->getFetiSolverTime += getTime();
}

template<class Scalar>
void 
GenMultiDomainEigen<Scalar>::reBuild(MDDynamMat &dMat)
{
 times->getFetiSolverTime -= getTime();
 decDomain->rebuildOps(dMat, 0.0, 0.0, 1.0);
 times->getFetiSolverTime += getTime();
}

template<class Scalar>
void
GenMultiDomainEigen<Scalar>::error(int subspacesize, int numRbm)
{
 filePrint(stderr,"*** subspace size must be bigger than number of rbm   ***\n");
 filePrint(stderr,"subspace size = %d\n",subspacesize);
 filePrint(stderr,"number of rbm = %d\n",numRbm);
}

template <class Scalar>
void
GenMultiDomainEigen<Scalar>::printTimers(GenParallelSolver<double> *solver)
{
 fflush(stderr);
 times->numSubdomain = decDomain->getNumSub();

 // --- Memory computations should be done so that other problems using
 // --- FETI can use the same routines.

 // --- NOTE: the average, min, and max should be done per processor,
 // ---       not per subdomain.

 int i;
 long (*memory)=(long *) dbg_alloca(sizeof(long)*decDomain->getNumSub());
 long totMemPrec = 0, totMemK = 0;

 // Kii memory calculations
 for(i=0; i<decDomain->getNumSub(); ++i) memory[i] = 0;
 execParal(decDomain->getNumSub(), this,
           &GenMultiDomainEigen<Scalar>::getMemoryPrec, memory);
 for(i=0; i<decDomain->getNumSub(); ++i) totMemPrec += memory[i];

 // K memory calculations
 for(i=0; i<decDomain->getNumSub(); ++i) memory[i] = 0;
 execParal(decDomain->getNumSub(), this,
           &GenMultiDomainEigen<Scalar>::getMemoryK, memory);
 for(i=0; i<decDomain->getNumSub(); ++i) totMemK += memory[i];

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

 times->memoryPrecond = totMemPrec;
 times->memoryK       = totMemK;

 //filePrint(stderr," ... Print Timers                   ... \n");
 
   switch(domain->solInfo().solvercntl->fetiInfo.version) {
     default:
     case FetiInfo::feti1:
     case FetiInfo::feti2:
       times->printStaticTimers(domain->getTimers(),
                                solver->getSolutionTime(),
                                domain->solInfo() ,
                                solver->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;

     case FetiInfo::fetidp:
       times->printFetiDPtimers(domain->getTimers(),
                                solver->getSolutionTime(),
                                domain->solInfo() ,
                                solver->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;
   }
}

template<class Scalar>
void
GenMultiDomainEigen<Scalar>::getMemoryK(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryK();
}

template<class Scalar>
void
GenMultiDomainEigen<Scalar>::getMemoryPrec(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryPrec();
}

template<class Scalar>
int 
GenMultiDomainEigen<Scalar>::getNumEigen()
{
  return domain->solInfo().nEig;
}

template<class Scalar>
int 
GenMultiDomainEigen<Scalar>::getEigenSolverType()
{
  return domain->solInfo().eigenSolverType;
}

template<class Scalar>
void
GenMultiDomainEigen<Scalar>::getSubSpaceInfo(int& subspacesize, int& maxIter, 
                                  double& tolEig, double& tolJac, bool &explicitK)
{
 subspacesize = domain->solInfo().subspaceSize;
 tolEig       = domain->solInfo().tolEig;
 tolJac       = domain->solInfo().tolJac;

 int numIter1 = domain->solInfo().maxitEig;
 int numIter2 = 10*subspacesize;

 maxIter = std::max(numIter1, numIter2);

 if (numIter1 > 0) maxIter = numIter1;

 explicitK = domain->solInfo().explicitK;
}

template<class Scalar>
void
GenMultiDomainEigenPostProcessor<Scalar>::eigenOutput(GenVector<Scalar> &eigValues, GenDistrVectorSet<Scalar> &eigVectors, int convEig)
{
  if (!times)
    times = new StaticTimers;
  startTimerMemory(times->output, times->memoryOutput);

  const double pi = 3.141592653589793;

  // Create a dummy DistrVector
  GenDistrVector<double> *dummyVector = 0;

  // loop over all eigenvalues
  int imode;
  int maxmode = (convEig && convEig < eigValues.size()) ? convEig : eigValues.size();

  if(domain->solInfo().modeDecompFlag) {

    int cmode = 0;
    for (int imode=0; imode < maxmode; ++imode) {
      if(domain->solInfo().test_ulrich && (imode+1)%10 == 0) continue; // this is for testing Ulrich's method for finding missed eigenvalues
      cmode++;
    }

    filePrint(stderr," ... Outputting %d EigenVectors in Binary file requested ...\n", cmode);
#ifdef DISTRIBUTED
    char *filename = new char[40];
    sprintf(filename,"EIGENMODES%d",structCom->myID());
    BinFileHandler modefile(filename, "w");
#else
    BinFileHandler modefile("EIGENMODES", "w"); // XXXX need to fix this for eigensweep
#endif
    //modefile.write(&maxmode, 1);
    modefile.write(&cmode, 1);
    int veclength = eigVectors[1].size();
    modefile.write(&veclength, 1);

    // Write eigenmodes in file EIGENMODES
    for (int imode=0; imode < maxmode; ++imode) {
      if(domain->solInfo().test_ulrich && (imode+1)%10 == 0) { filePrint(stderr," ### MODE %d IS NOT BEING OUTPUT ###\n",imode+1); continue; } // this is for testing Ulrich's method for finding missed eigenvalues
      double *modes = eigVectors[imode].data();
      modefile.write(modes, veclength);
    }
  }

  for(imode=0; imode<maxmode; ++imode) {
    if(domain->solInfo().buckling)
      decDomain->postProcessing(eigVectors[imode],*dummyVector,eigValues[imode],dummyVector,imode);
    else
      decDomain->postProcessing(eigVectors[imode],*dummyVector,sqrt(eigValues[imode])/(2.0*pi),dummyVector,imode);
#ifdef DISTRIBUTED
    structCom->sync();
#endif
  }

 if(!domain->solInfo().doEigSweep) {
   // Print omega and frequency values to screen
   if(domain->solInfo().buckling) {
     filePrint(stderr," --------------------------------------\n");
     filePrint(stderr," Mode\tLambda\n");
     filePrint(stderr," --------------------------------------\n");
     //for(imode=0; imode<eigValues.size(); ++imode)
     for(imode=0; imode<maxmode; ++imode)
       filePrint(stderr," %d\t%e\n",imode+1,eigValues[imode]);
     filePrint(stderr," --------------------------------------\n");
   } else {
     filePrint(stderr," --------------------------------------\n");
     filePrint(stderr," Mode\tOmega^2\t\tFrequency\n");
     filePrint(stderr," --------------------------------------\n");
     //for(imode=0; imode<eigValues.size(); ++imode)
     for(imode=0; imode<maxmode; ++imode)
       filePrint(stderr," %d\t%e\t%e\n",imode+1,eigValues[imode],
               sqrt(eigValues[imode])/(2.0*pi));
     filePrint(stderr," --------------------------------------\n");
   }
 }
 stopTimerMemory(times->output, times->memoryOutput);
}
