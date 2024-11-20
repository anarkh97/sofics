#include <alloca.h>
#include <Threads.d/Paral.h>
#include <Threads.d/PHelper.h>
#include <Driver.d/Domain.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Solvers.d/PCGSolver.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Memory.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/DecDomain.h>
#include <Sfem.d/Sfem.h>

#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif

extern Sfem* sfem;

template<class Scalar>
GenMultiDomainInpcStatic<Scalar>::GenMultiDomainInpcStatic(Domain *d)
{
 domain = d;
/*
 switch(domain->solInfo().solvercntl->fetiInfo.version) {
  default:
  case FetiInfo::feti1:
    filePrint(stderr," ... FETI-1 is Selected             ...\n");
    break;
  case FetiInfo::feti2:
    filePrint(stderr," ... FETI-2 is Selected             ...\n");
    break;
  case FetiInfo::fetidp:
    if (!(domain->solInfo().solvercntl->fetiInfo.dph_flag)) 
      filePrint(stderr," ... FETI-Dual/Primal is Selected   ...\n");
    else 
      filePrint(stderr," ... FETI-DPH is Selected           ...\n");
    break;
  }
*/
#ifdef DISTRIBUTED
 decDomain = new GenDistrDomain<Scalar>(domain);
#else
 decDomain = new GenDecDomain<Scalar>(domain);
#endif

 times  = new StaticTimers;

 // debug: check static timers are initialized
 //MatrixTimers &mt = domain->getTimers();
}

template<class Scalar>
DistrBlockInfo &
GenMultiDomainInpcStatic<Scalar>::solVecInfo()
{
 return info;
}

template<class Scalar>
DistrBlockInfo &
GenMultiDomainInpcStatic<Scalar>::solVecInfo(int i)
{
 info.nblocks = i;
 return info;
}


template<class Scalar>
void
GenMultiDomainInpcStatic<Scalar>::getRHSinpc(DistrBlockVector<Scalar> &rhs)
{
 times->formRhs -= getTime();
 rhs.setBlock(0,*rhs_inpc);
 for(int i=1; i<info.nblocks; ++i) rhs.getBlock(i).zero();
 // copy rhs_inpc
 times->formRhs += getTime();
}

template<class Scalar>
void
GenMultiDomainInpcStatic<Scalar>::preProcess()
{
 // Makes renumbering, connectivities and dofsets
 startTimerMemory(times->preProcess, times->memoryPreProcess);
 decDomain->preProcess();
 // Make all element's dofs
 MultiDomainOp mdop(&MultiDomainOp::makeAllDOFs, decDomain->getAllSubDomains());
#ifdef DISTRIBUTED
 execParal(decDomain->getNumSub(), &mdop, &MultiDomainOp::runFor);
#else
 threadManager->execParal(decDomain->getNumSub(), &mdop);
#endif
 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 int L =  sfem->getL();
 int P = sfem->getP();
 int ndim = sfem->getndim();
 int output_order = sfem->getoutput_order();

 info.nblocks = P;
 info.nnzblkindex = new int[P];
 for (int i=0; i<P; ++i) info.nnzblkindex[i]=1;
 info.blockinfo=new DistrInfo[P];
 for(int i=0; i<P; ++i) info.blockinfo[i] = decDomain->solVecInfo();
 sfbm = new DistrSfemBlockMatrix<Scalar>(L,P,ndim,output_order,info); 

 for(int i=0; i< ndim+1; ++i) {
   startTimerMemory(times->sfemBuildOps, times->memorySfemBuildOps);
   decDomain->setNewProperties(i);
   MDDynamMat mdMat;
   bool make_feti = (i==0);
   decDomain->buildOps(mdMat, 0.0, 0.0, 1.0, (Rbm **) 0, (FullSquareMatrix **) 0, make_feti);
   sfbm->setKi(mdMat.K,i);
   stopTimerMemory(times->sfemBuildOps, times->memorySfemBuildOps);
   if(i==0) {
     feti_solver = (GenFetiSolver<Scalar> *) mdMat.dynMat; //decDomain->getFetiSolver();
     sfbm->setMeanSolver(feti_solver);
     times->formRhs -= getTime();
     rhs_inpc = new GenDistrVector<Scalar>(info.blockinfo[0]);
     //feti_solver->makeStaticLoad(*rhs_inpc);
     execParal(decDomain->getNumSub(), this, &GenMultiDomainInpcStatic<Scalar>::subGetRHS, *rhs_inpc, mdMat.Kuc);
     mdMat.K->getAssembler()->assemble(*rhs_inpc);
     times->formRhs += getTime();
   }
 }
 solver = new GenPCGSolver<Scalar, DistrBlockVector<Scalar>, DistrSfemBlockMatrix<Scalar> >(sfbm, domain->solInfo().solvercntl->precond, domain->solInfo().solvercntl->maxit,
                                                                                            domain->solInfo().solvercntl->tol, domain->solInfo().solvercntl->maxvecsize);
}

template<class Scalar>
void
GenMultiDomainInpcStatic<Scalar>::subGetRHS(int isub, GenDistrVector<Scalar>& rhs, GenSubDOp<Scalar> *Kuc)
{
 GenSubDomain<Scalar> *sd = decDomain->getSubDomain(isub);
 GenStackVector<Scalar> subrhs(rhs.subData(isub), rhs.subLen(isub));
 sd->buildRHSForce(subrhs, (*Kuc)[isub]);
}

template<class Scalar>
void
GenMultiDomainInpcStatic<Scalar>::rebuildSolver()
{
  std::cerr << "GenMultiDomainInpcStatic<Scalar>::rebuildSolver() is not implemented \n";
}

template<class Scalar>
void
GenMultiDomainInpcStatic<Scalar>::clean()
{
 decDomain->clean_up();
}

template<class Scalar>
GenSolver<Scalar> *
GenMultiDomainInpcStatic<Scalar>::getSolver()
{
 return solver;
}

template<class Scalar>
GenMultiDomainInpcPostProcessor<Scalar> *
GenMultiDomainInpcStatic<Scalar>::getPostProcessor()
{
 return new GenMultiDomainInpcPostProcessor<Scalar>(decDomain, solver, times);
}

template<class Scalar>
void
GenMultiDomainInpcPostProcessor<Scalar>::staticOutput(DistrBlockVector<Scalar> &sol, DistrBlockVector<Scalar> &force,
                                                      bool printTimers, int ndflag)
{
 startTimerMemory(times->output, times->memoryOutput);
 decDomain->postProcessing(sol.getBlock(0), force.getBlock(0), 0.0, 0, x, 0, 0, ndflag);
 if(ndflag == 3) x++;
// if(ndflag == 0) decDomain->postProcessing(sol.getBlock(0), force.getBlock(0));
// if(ndflag == 0) decDomain->postProcessing(sol.getBlock(0), force.getBlock(0), 0.0, 0, 0, 0, 0);
// else std::cerr << "XXXX\n";
 stopTimerMemory(times->output, times->memoryOutput);

/*
 times->numSubdomain = decDomain->getNumSub();

 // --- Memory computations should be done so that other problems using
 // --- FETI can use the same routines.

 // --- NOTE: the average, min, and max should be done per processor,
 // ---       not per subdomain.

 int i;
 long (*memory)=(long *) alloca(sizeof(long)*decDomain->getNumSub());
 long totMemPrec = 0, totMemK = 0;

 // Kii memory calculations
 for(i=0; i<decDomain->getNumSub(); ++i) memory[i] = 0;
 execParal(decDomain->getNumSub(), this,  
           &GenMultiDomainInpcPostProcessor<Scalar>::getMemoryPrec, memory);
 for(i=0; i<decDomain->getNumSub(); ++i) totMemPrec += memory[i];

 // K memory calculations
 for(i=0; i<decDomain->getNumSub(); ++i) memory[i] = 0;
 execParal(decDomain->getNumSub(), this, 
           &GenMultiDomainInpcPostProcessor<Scalar>::getMemoryK, memory);
 for(i=0; i<decDomain->getNumSub(); ++i) totMemK += memory[i];

 Timings &fetiTimers = feti_solver->getTimers();
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

 if(printTimers) {
   switch(domain->solInfo().solvercntl->fetiInfo.version) {
     default:
     case FetiInfo::feti1:
     case FetiInfo::feti2:
       times->printStaticTimers(domain->getTimers(), 
                                feti_solver->getSolutionTime(), 
                                domain->solInfo() , 
                                feti_solver->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;

     case FetiInfo::fetidp:
       times->printFetiDPtimers(domain->getTimers(),
                                feti_solver->getSolutionTime(),
                                domain->solInfo() ,
                                feti_solver->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;
   }
//   geoSource->closeOutputFiles();
 }
*/
 long memoryUsed = 0;
 double solveTime  = 0.0;

 memoryUsed = solver->size();
 solveTime  = solver->getSolutionTime();
 times->precond = solver->getTimers().precond;
 if(printTimers) {
   times->printStaticTimers(solveTime, memoryUsed, domain);
 }

 filePrint(stderr," --------------------------------------\n");
}

template<class Scalar>
void
GenMultiDomainInpcPostProcessor<Scalar>::getMemoryKii(int iSub, long *memory )
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryKii();
}

template<class Scalar>
void
GenMultiDomainInpcPostProcessor<Scalar>::getMemoryK(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryK();
}

template<class Scalar>
void
GenMultiDomainInpcPostProcessor<Scalar>::getMemoryPrec(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryPrec();
}

template<class Scalar>
void
GenMultiDomainInpcStatic<Scalar>::setIWaveDir(int _i)
{
  std::cerr << "GenMultiDomainInpcStatic<Scalar>::setIWaveDir(...) is not implemented \n";
}

