#include <Utils.d/dbg_alloca.h>
#include <Threads.d/Paral.h>
#include <Threads.d/PHelper.h>
#include <Driver.d/Domain.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Memory.h>
#include <Driver.d/GeoSource.h>
#include <Solvers.d/MultiDomainRbm.h>

#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif

template<class Scalar>
GenMultiDomainStatic<Scalar>::GenMultiDomainStatic(Domain *d)
 : MultiDomainBase(d->solInfo())
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
GenMultiDomainStatic<Scalar>::~GenMultiDomainStatic()
{
 delete decDomain;
 delete times;
 // solver deleted in StaticSolver
}

template<class Scalar>
DistrInfo &
GenMultiDomainStatic<Scalar>::solVecInfo()
{
 return decDomain->solVecInfo();
}

template<class Scalar>
DistrInfo &
GenMultiDomainStatic<Scalar>::solVecInfo(int i) 
{
 std::cerr << "Warning : GenMultiDomainStatic<Scalar>::solVecInfo(int i) should not be called" << std::endl;
 return decDomain->solVecInfo();
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getRHS(GenDistrVector<Scalar> &rhs)
{
 if(domain->solInfo().loadcases.size() > 0)
   filePrint(stderr," ... Building the Force (Case %2d)   ...\n", domain->solInfo().loadcases.front());
 else
   filePrint(stderr," ... Building the Force             ...\n");

 times->formRhs -= getTime();
 execParal(decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::subGetRHS, rhs);

 // rbm or eigen mode projector 
 if(domain->solInfo().filterFlags || domain->solInfo().modeFilterFlag)
   trProject(rhs);

 if(domain->solInfo().solvercntl->type == SolverSelection::Iterative) allOps.spMat->getAssembler()->assemble(rhs); // XXXX
 times->formRhs += getTime();
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::subGetRHS(int isub, GenDistrVector<Scalar>& rhs)
{
 GenSubDomain<Scalar> *sd = decDomain->getSubDomain(isub);
 GenStackVector<Scalar> subrhs(rhs.subData(isub), rhs.subLen(isub));
 sd->buildRHSForce(subrhs, (*allOps.Kuc)[isub]);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::preProcess()
{
 // Makes renumbering, connectivities and dofsets
 startTimerMemory(times->preProcess, times->memoryPreProcess);
 decDomain->preProcess();
 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 // Construct FETI solver and any other matrices required
 times->getFetiSolverTime -= getTime();
 for(int i=0; i<decDomain->getNumSub(); ++i) decDomain->getSubDomain(i)->makeAllDOFs();
 decDomain->buildOps(allOps, 0.0, 0.0, 1.0);
 solver = allOps.sysSolver;

 times->getFetiSolverTime += getTime();

 int useRbmFilter = domain->solInfo().filterFlags;
 if(useRbmFilter || domain->solInfo().rbmflg) {
   MultiDomainRbm<Scalar> *rigidBodyModes = decDomain->constructRbm();
   if(useRbmFilter) {
     filePrint(stderr," ... RBM Filter Requested           ...\n");
     projector_prep(rigidBodyModes, allOps.M);
   }
   delete rigidBodyModes;
 }

 int useModeFilter = domain->solInfo().modeFilterFlag;
 if(useModeFilter) {
   filePrint(stderr, " ... MODE Filter Requested          ...\n");
   eigmode_projector_prep(decDomain->solVecInfo());
 }
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::rebuildSolver()
{
  times->getFetiSolverTime -= getTime();
  GenMDDynamMat<Scalar> ops;
  ops.sysSolver = allOps.sysSolver;
  ops.K = allOps.K;
  ops.Kuc = allOps.Kuc;
  ops.M = allOps.M;
  ops.Muc = allOps.Muc;
  ops.C_deriv = allOps.C_deriv;
  ops.Cuc_deriv = allOps.Cuc_deriv;
  ops.K_deriv = allOps.K_deriv;
  ops.Kuc_deriv = allOps.Kuc_deriv;
  ops.num_K_deriv = allOps.num_K_deriv;
  ops.K_arubber_l = allOps.K_arubber_l;
  ops.K_arubber_m = allOps.K_arubber_m;
  ops.Kuc_arubber_l = allOps.Kuc_arubber_l;
  ops.Kuc_arubber_m = allOps.Kuc_arubber_m;
  ops.num_K_arubber = allOps.num_K_arubber;

  decDomain->rebuildOps(ops, 0.0, 0.0, 1.0);
  paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<Scalar>::setRebuildPade, true);
  times->getFetiSolverTime += getTime();
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::scaleDisp(GenDistrVector<Scalar> &u)
{
  decDomain->scaleDisp(u);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::scaleInvDisp(GenDistrVector<Scalar> &u)
{
  decDomain->scaleInvDisp(u);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::scaleDisp(GenDistrVector<Scalar> &u, double alpha)
{
  decDomain->scaleDisp(u, alpha);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::forceContinuity(GenDistrVector<Scalar> &u)
{
  decDomain->forceContinuity(u);
}


template<class Scalar>
void
GenMultiDomainStatic<Scalar>::forceAssemble(GenDistrVector<Scalar> &u)
{
  decDomain->forceAssemble(u);
}


template<class Scalar>
void
GenMultiDomainStatic<Scalar>::clean()
{
}

template<class Scalar>
GenParallelSolver<Scalar> *
GenMultiDomainStatic<Scalar>::getSolver()
{
  return solver;
}

template<class Scalar>
GenMultiDomainPostProcessor<Scalar> *
GenMultiDomainStatic<Scalar>::getPostProcessor()
{
  return new GenMultiDomainPostProcessor<Scalar>(decDomain, solver, times);
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::staticOutput(GenDistrVector<Scalar> &sol, GenDistrVector<Scalar> &force,
                                                  bool printTimers, int ndflag)
{
 startTimerMemory(times->output, times->memoryOutput);
 decDomain->postProcessing(sol, force);
 stopTimerMemory(times->output, times->memoryOutput);

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
           &GenMultiDomainPostProcessor<Scalar>::getMemoryPrec, memory);
 for(i=0; i<decDomain->getNumSub(); ++i) totMemPrec += memory[i];

 // K memory calculations
 for(i=0; i<decDomain->getNumSub(); ++i) memory[i] = 0;
 execParal(decDomain->getNumSub(), this, 
           &GenMultiDomainPostProcessor<Scalar>::getMemoryK, memory);
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

 if(printTimers) {
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
//   geoSource->closeOutputFiles();
 }

 filePrint(stderr," --------------------------------------\n");
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::getStressStrain(GenDistrVector<Scalar> &, int fileNumber,
                                                     int stressIndex, double time, int pflag)
{
  std::cerr << "GenMultiDomainPostProcessor::getStressStrain not implemented" << std::endl;
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::setsizeSfemStress(int fileNumber)
{
  std::cerr << "GenMultiDomainPostProcessor::setsizeSfemStress(int fileNumber) not implemented" << std::endl;
}

template<class Scalar>
int
GenMultiDomainPostProcessor<Scalar>::getsizeSfemStress()
{
  std::cerr << "GenMultiDomainPostProcessor::getsizeSfemStress() not implemented" << std::endl;
  return 0;
}

template<class Scalar>
Scalar*
GenMultiDomainPostProcessor<Scalar>::getSfemStress(int fileNumber)
{
  std::cerr << "GenMultiDomainPostProcessor::getSfemStress() not implemented" << std::endl;
  return 0;
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::updateSfemStress(Scalar* str, int fileNumber)
{
  std::cerr << "GenMultiDomainPostProcessor::updateSfemStress() not implemented" << std::endl;
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::getMemoryK(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryK();
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::getMemoryPrec(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryPrec();
}

// FETI-H
template<class Scalar>
void
GenMultiDomainStatic<Scalar>::setIWaveDir(int _i)
{
 int i;
 for(i=0;i<decDomain->getNumSub();i++)
   decDomain->getSubDomain(i)->iWaveDir = _i;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getFreqSweepRHS(GenDistrVector<Scalar> *rhs, 
                                              GenDistrVector<Scalar> **sol_prev, int iRHS)
{
  // TODO
  for(int i=0;i<decDomain->getNumSub();i++) {
    decDomain->getSubDomain(i)->M = (*allOps.M)[i];
    decDomain->getSubDomain(i)->Muc = (GenCuCSparse<Scalar> *)(*allOps.Muc)[i];
    if (allOps.C) decDomain->getSubDomain(i)->C = (*allOps.C)[i];
    if (allOps.C_deriv) {
      decDomain->getSubDomain(i)->C_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      (decDomain->getSubDomain(i)->C_deriv)[0] = (*(allOps.C_deriv[0]))[i];
      for(int j=1;j<iRHS+1;j++) (decDomain->getSubDomain(i)->C_deriv)[j] = 0;
    }
    if (allOps.Cuc_deriv) {
      decDomain->getSubDomain(i)->Cuc_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      (decDomain->getSubDomain(i)->Cuc_deriv)[0] = (*(allOps.Cuc_deriv[0]))[i];
      for(int j=1;j<iRHS+1;j++) (decDomain->getSubDomain(i)->Cuc_deriv)[j] = 0;
    }
    if (allOps.K_deriv) {
      decDomain->getSubDomain(i)->K_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      for(int j=0;j<iRHS+1;j++) if ((*(allOps.K_deriv[j]))[i]!=0) {
          (decDomain->getSubDomain(i)->K_deriv)[j] = (*(allOps.K_deriv[j]))[i];
      }
    }
    if (allOps.Kuc_deriv) {
      decDomain->getSubDomain(i)->Kuc_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      for(int j=0;j<iRHS+1;j++) if ((*(allOps.Kuc_deriv[j]))[i]!=0) 
        (decDomain->getSubDomain(i)->Kuc_deriv)[j] = (*(allOps.Kuc_deriv[j]))[i];
                                else
        (decDomain->getSubDomain(i)->Kuc_deriv)[j] = 0;
      
    }
  }


  Timings &timers = solver->getTimers();
  if(domain->solInfo().isCoupled && domain->solInfo().getFetiInfo().fsi_corner == 0) {
    timedParal(timers.buildRhs, decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::multMCoupled1, rhs, sol_prev, iRHS);
    decDomain->getWiCommPattern()->exchange();
    timedParal(timers.buildRhs, decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::multMCoupled2, rhs);
  }
  else timedParal(timers.buildRhs, decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::multM, rhs, sol_prev, iRHS);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::multM(int iSub, GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **u, int k)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(iSub);
  if (u==0) sd->multM(rhs->subData(iSub), 0, k);
  else {
    GenStackVector<Scalar> **sub_u = new GenStackVector<Scalar> * [k+1];
    for(int i=0; i<=k; ++i)
      sub_u[i]=
        new  GenStackVector<Scalar>(u[i]->subData(iSub), u[i]->subLen(iSub));
    sd->multM(rhs->subData(iSub), sub_u, k); // TODO
    delete [] sub_u;
  }
}


template<class Scalar>
void
GenMultiDomainStatic<Scalar>::multMCoupled1(int iSub, GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **u, int k)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(iSub);
  GenStackVector<Scalar> **sub_u = new GenStackVector<Scalar> * [k+1];
  for(int i=0; i<=k; ++i)
    sub_u[i]= new  GenStackVector<Scalar>(u[i]->subData(iSub), u[i]->subLen(iSub));
  sd->multMCoupled1(rhs->subData(iSub), sub_u, k, decDomain->getWiCommPattern()); // TODO
  delete [] sub_u;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::multMCoupled2(int iSub, GenDistrVector<Scalar> *rhs)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(iSub);
  sd->multMCoupled2(rhs->subData(sd->localSubNum()), decDomain->getWiCommPattern()); // TODO
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getWCAWEFreqSweepRHS(GenDistrVector<Scalar> *rhs, 
   GenDistrVector<Scalar> **wcawe_u, Scalar *pU, Scalar *pb,
   int maxRHS, int iRHS)
{
  for(int i=0;i<decDomain->getNumSub();i++) {
    decDomain->getSubDomain(i)->M = (*allOps.M)[i];
    decDomain->getSubDomain(i)->Muc = (GenCuCSparse<Scalar> *)(*allOps.Muc)[i];
    if (allOps.C) decDomain->getSubDomain(i)->C = (*allOps.C)[i];
    if (allOps.C_deriv) {
      decDomain->getSubDomain(i)->C_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      (decDomain->getSubDomain(i)->C_deriv)[0] = (*(allOps.C_deriv[0]))[i];
      for(int j=1;j<iRHS+1;j++) (decDomain->getSubDomain(i)->C_deriv)[j] = 0;
    }
    if (allOps.Cuc_deriv) {
      decDomain->getSubDomain(i)->Cuc_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      (decDomain->getSubDomain(i)->Cuc_deriv)[0] = (*(allOps.Cuc_deriv[0]))[i];
      for(int j=1;j<iRHS+1;j++) (decDomain->getSubDomain(i)->Cuc_deriv)[j] = 0;
    }
    if (allOps.K_deriv) {
      decDomain->getSubDomain(i)->K_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      for(int j=0;j<iRHS+1;j++) if ((*(allOps.K_deriv[j]))[i]!=0) {
          (decDomain->getSubDomain(i)->K_deriv)[j] = (*(allOps.K_deriv[j]))[i];
      }
    }
    if (allOps.Kuc_deriv) {
      decDomain->getSubDomain(i)->Kuc_deriv =
         new GenSparseMatrix<Scalar>*[iRHS+1];
      for(int j=0;j<iRHS+1;j++) if ((*(allOps.Kuc_deriv[j]))[i]!=0) 
        (decDomain->getSubDomain(i)->Kuc_deriv)[j] = (*(allOps.Kuc_deriv[j]))[i];
                                else
        (decDomain->getSubDomain(i)->Kuc_deriv)[j] = 0;
      
    }
  }
  Timings &timers = solver->getTimers();
  timedParal(timers.buildRhs, decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::multWCAWE, rhs, wcawe_u, pU, pb, maxRHS, iRHS);
 
}


template<class Scalar>
void
GenMultiDomainStatic<Scalar>::multWCAWE(int iSub, GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **u, Scalar *pU, Scalar *pb, int maxRHS, int iRHS)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(iSub);
  GenStackVector<Scalar> **sub_u = new GenStackVector<Scalar> * [iRHS+1];
  for(int i=0; i<=iRHS; ++i)
    sub_u[i]=
      new  GenStackVector<Scalar>(u[i]->subData(iSub), u[i]->subLen(iSub));
  sd->multWCAWE(rhs->subData(iSub), sub_u, pU,pb,maxRHS,iRHS); 
  delete [] sub_u;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getRHS(GenDistrVector<Scalar> &rhs, double omega,
                                     double deltaomega)
{
  filePrint(stderr," ... Building the Force             ...\n");
  rhs.zero();
  GenDistrVector<Scalar> *tmp = new GenDistrVector<Scalar>(solVecInfo());
  double o[2] = { omega, deltaomega };
// RT: 4/30/09 - should be timed
  execParal(decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::makeSubdomainStaticLoadGalPr, rhs, *tmp, o);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::makeSubdomainStaticLoadGalPr(int isub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &tmp, double *o)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(isub);
  GenStackVector<Scalar> subf(f.subData(isub), f.subLen(isub));
  GenStackVector<Scalar> subv(tmp.subData(isub), tmp.subLen(isub));

  // TODO subdomain shouldn't have Cuc_deriv pointer
  sd->buildRHSForce(subf, subv, (*allOps.Kuc)[isub], (*allOps.Muc)[isub],
                    sd->Cuc_deriv, sd->Kuc_deriv, sd->Kuc_arubber_l, sd->Kuc_arubber_m, o[0], o[1]);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::pade(GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **sol_prev, double *h, double x)
{
  execParal(decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::subPade, sol, sol_prev, h, x);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::subPade(int iSub, GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **u, double *h, double x)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(iSub);
  GenStackVector<Scalar> *sub_sol = new GenStackVector<Scalar>(sol->subData(iSub), sol->subLen(iSub));
  int usize = (domain->solInfo().getSweepParams()->nFreqSweepRHS+1)*domain->solInfo().getSweepParams()->padeN;
  GenStackVector<Scalar> **sub_u = new GenStackVector<Scalar> * [usize];
  for(int i=0; i<usize; ++i)
    sub_u[i]= new  GenStackVector<Scalar>(u[i]->subData(iSub), u[i]->subLen(iSub));
  sd->pade(sub_sol, sub_u, h, x); // TODO
  delete [] sub_u;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getRHSinpc(GenDistrVector<Scalar> &)
{
  std::cerr << "GenMultiDomainStatic::getRHSinpc not implemented" << std::endl;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::preProcessSA()
{
  std::cerr<< "GenMultiDomainStatic::PreProcessSA not implemented" << std::endl;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::postProcessSA(GenDistrVector<Scalar> &)
{
  std::cerr<< "GenMultiDomainStatic::PostProcessSA not implemented" << std::endl;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::assignRandMat()
{
  std::cerr << "GenMultiDomainStatic::assignRandMat() not implemented" << std::endl;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::retrieveElemset()
{
  std::cerr << "GenMultiDomainStatic::retrieveElemset() not implemented" << std::endl;
}

template<class Scalar>
AllSensitivities<Scalar> *
GenMultiDomainStatic<Scalar>::getAllSensitivities()
{
  std::cerr << "GenMultiDomainStatic::getAllSensitivities() not implemented" << std::endl;
  return 0;
}
