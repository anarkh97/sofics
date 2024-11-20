#ifndef _MD_INPC_STATIC_DESCR_H_
#define _MD_INPC_STATIC_DESCR_H_
#include <Driver.d/DecDomain.h>
#include <Sfem.d/DistrSfemBlockMatrix.h>
#include <Paral.d/MDOp.h>
#include <Paral.d/MDDynam.h>
#include <Paral.d/Assembler.h>
#include <iostream>

template <class Scalar> class GenFetiSolver;
class StaticTimers;
template <class Scalar> class GenDistrVector;
class DistrGeomState;
class DistrInfo;
template <class Scalar> class AllOps;
template <class Scalar> struct AllSensitivities;

template <class Scalar> 
class GenMultiDomainInpcPostProcessor 
{
    GenDecDomain<Scalar> *decDomain;
    GenSolver<Scalar> *solver;
    StaticTimers *times;
    int x;

  public:
    GenMultiDomainInpcPostProcessor(GenDecDomain<Scalar> *d, GenSolver<Scalar> *s, StaticTimers* _times=0)
      { decDomain = d; solver = s; times = _times; x = 0; }
    void staticOutput(DistrBlockVector<Scalar> &, DistrBlockVector<Scalar> &, bool printTimers = true, int ndflag=0);
    void staticOutput(DistrGeomState *u, double lambda = 1.0);
    void getMemoryK(int iSub, long *memory);
    void getMemoryPrec(int iSub, long *memory);
    void getMemoryKii(int iSub, long *memory);

    void getStressStrain(DistrBlockVector<Scalar> &sol, int fileNumber,
                                        int stressIndex, double time, int printFlag) 
          { decDomain->getStressStrain(sol.getBlock(0), fileNumber, stressIndex, time, printFlag); } 
              // YYY DG   (i) Assuming Findex in DecDomain & DistDom is stressIndex
              // YYY DG  (iv) For DistDom there is an extra argument "int x", also "iter", how to implement that ?
    void getStressStrain(GenDistrVector<Scalar> &sol, int fileNumber,
                                        int stressIndex, double time, int printFlag) 
          { decDomain->getStressStrain(sol, fileNumber, stressIndex, time, printFlag); } 
    void setsizeSfemStress(int fileNumber) { decDomain->setsizeSfemStress(fileNumber); }
    int getsizeSfemStress() { return decDomain->getsizeSfemStress(); }
    double* getSfemStress(int fileNumber) { return decDomain->getSfemStress(fileNumber);}
    void updateSfemStress(double* str, int fileNumber) { decDomain->updateSfemStress(str, fileNumber);}
    void setSolver(GenSolver<Scalar> *s) { solver = s; }
};


template<class Scalar>
class GenMultiDomainInpcStatic 
{
    Domain *domain;
    GenDecDomain<Scalar> *decDomain;
    GenSolver<Scalar> *solver;
    GenFetiSolver<Scalar> *feti_solver;
    DistrSfemBlockMatrix<Scalar> *sfbm;
    StaticTimers *times;
    GenDistrVector<Scalar> *rhs_inpc;
    DistrBlockInfo info;
    AllOps<Scalar> allOps;
 public:
    GenMultiDomainInpcStatic(Domain *d);
    ~GenMultiDomainInpcStatic() { delete decDomain; delete times; }  // solver deleted in StaticSolver

    DistrBlockInfo &solVecInfo();
    DistrBlockInfo &solVecInfo(int i);
    void getRHS(DistrBlockVector<Scalar> &) {std::cerr << "GenMultiDomainInpcStatic::getRHS not implemented" << std::endl;}
    void getRHSinpc(DistrBlockVector<Scalar>  &);
    void preProcessSA() {std::cerr << "GenMultiDomainInpcStatic::preProcessSA not implemented" << std::endl; }
    void postProcessSA(DistrBlockVector<Scalar> &) {std::cerr << "GenMultiDomainInpcStatic::postProcessSA not implemented" << std::endl; }
    void preProcess();
    void rebuildSolver();
    void scaleDisp(DistrBlockVector<Scalar> &) {std::cerr << "scaleDisp(DistrBlockVector not implemented" << std::endl;}
    void scaleInvDisp(DistrBlockVector<Scalar> &) {std::cerr << "scaleInvDisp(DistrBlockVector not implemented" << std::endl;}
    void scaleDisp(DistrBlockVector<Scalar> &, double alpha) {std::cerr << "scaleDisp(DistrBlockVector not implemented" << std::endl;}
    void forceContinuity(DistrBlockVector<Scalar> &) {std::cerr << "forceContinuity(DistrBlockVector not implemented" << std::endl;}
    void forceAssemble(DistrBlockVector<Scalar> &) {std::cerr << "forceAssemble(DistrBlockVector not implemented" << std::endl;}
    void clean();
    void setIWaveDir(int _i); // FETI-H
    void getFreqSweepRHS(DistrBlockVector<Scalar> *rhs, DistrBlockVector<Scalar> **sol_prev, int iRHS)
           {std::cerr << "getFreqSweepRHS(DistrBlockVector) not implemented" << std::endl; }	
    void getWCAWEFreqSweepRHS(DistrBlockVector<Scalar> *rhs, DistrBlockVector<Scalar> **sol_prev, Scalar *, Scalar *, int, int iRHS)
           {std::cerr << "getFreqSweepRHS(DistrBlockVector) not implemented" << std::endl; }	
    void buildDeltaK(double w0, double w) 
           {std::cerr << "buildDeltaK not implemented" << std::endl; }	
    void getRHS(DistrBlockVector<Scalar> &,double,double) {std::cerr << "GenMultiDomainInpcStatic::getRHS not implemented" << std::endl;}
    void pade(DistrBlockVector<Scalar> *sol, DistrBlockVector<Scalar> **sol_prev, double *h, double x)
           {std::cerr << "pade(DistrBlockVector) not implemented" << std::endl; }
    GenSolver<Scalar> *getSolver();
    AllOps<Scalar> *getAllOps() { return &allOps; }
    GenMultiDomainInpcPostProcessor<Scalar> *getPostProcessor();
    StaticTimers *getStaticTimers() { return times; }
    void assignRandMat() {decDomain->assignRandMat(); }     
    void retrieveElemset() {decDomain->retrieveElemset();}
    void project(DistrBlockVector<Scalar> &) {std::cerr << "project(DistrBlockVector) not implemented" << std::endl;}
    AllSensitivities<Scalar> *getAllSensitivities() { std::cerr << "GenMultiDomainInpcStatic::getAllSensitivities() not implemented" << std::endl; return 0; }
  private:
    void subGetRHS(int isub, GenDistrVector<Scalar>& rhs, GenSubDOp<Scalar> *Kuc);
};

#ifdef _TEMPLATE_FIX_
#include <Paral.d/MDInpcStatic.C>
#endif

#endif
