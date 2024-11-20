#ifndef _STATICPROBTYPE_H_
#define _STATICPROBTYPE_H_

#include <Utils.d/MathUtils.h>
#include <Sfem.d/SfemNonInpc.h>
#include <Sfem.d/SfemInpc.h>
#include <Sfem.d/MultInteg.h>

template<class Scalar> struct AllSensitivities;

template <class Scalar,
          class OpSolver,
          class VecType, 
          class PostProcessor, 
          class ProblemDescriptor,
          class ComplexVecType>
class StaticSolver 
{
 protected:
     ProblemDescriptor *probDesc;
     VecType *sol;
     VecType *rhs;
     VecType *adj;
//     OpSolver  *solver;
     OpSolver *allOps; // container with solver and any other operators required (eg. K, M, ...)
     PostProcessor *postProcessor;

     int ia, ib;
     VecType **a, **b;  // pade P, Q coefs
     VecType *P, *Q;    // pade P(x), Q(x)
     ComplexVecType **ca, **cb;  // fourier-pade P, Q coefs
     ComplexVecType *cP, *cQ;    // fourier-pade P(x), Q(x)

     SfemNonInpc<Scalar, VecType> *sfem_noninpc;
     SfemInpc<Scalar, VecType> *sfem_inpc;  

   public:
     StaticSolver() { sol = 0; rhs = 0; adj = 0; allOps = 0; probDesc = 0; postProcessor = 0;
                      a = 0; b = 0; P = 0; Q = 0; ca = 0; cb = 0; cP = 0; cQ = 0; }
     StaticSolver(ProblemDescriptor *PrbD) 
       { probDesc = PrbD; sol = 0; rhs = 0; adj = 0; allOps = 0; postProcessor = 0;
         a = 0; b = 0; P = 0; Q = 0; ca = 0; cb = 0; cP = 0; cQ = 0;
         if(domain->solInfo().noninpc) {
            sfem_noninpc = new SfemNonInpc<Scalar, VecType>();
            sfem_noninpc->init_XiPsi(); // needed only for pdf
            sfem->init_XiPsi(); // needed for other postprocessing
            sfem->genXiPsi(0);
//            sfem_noninpc->assignRandMat();
            probDesc->assignRandMat();
         }
         else if(domain->solInfo().inpc) {
            sfem_inpc = new SfemInpc<Scalar, VecType>();
            sfem_inpc->init_XiPsi(); // needed only for pdf
            sfem->init_XiPsi(); // needed for other postprocessing
         }
         else {
            sfem_inpc = 0;
            sfem_noninpc = 0;
         }
       }
     ~StaticSolver() { if(sol) delete sol; if(rhs) delete rhs; if(adj) delete adj; 
                       /*if(allOps) delete allOps;*/ if(postProcessor) delete postProcessor; 
                       if(a) { for(int i=0; i<ia; ++i) delete a[i]; delete [] a; }
                       if(b) { for(int i=0; i<ib; ++i) delete b[i]; delete [] b; }
                       if(P) delete P; if(Q) delete Q; 
                       if(ca) { for(int i=0; i<ia; ++i) delete ca[i]; delete [] ca; }
                       if(cb) { for(int i=0; i<ib; ++i) delete cb[i]; delete [] cb; }
                       if(cP) delete cP; if(cQ) delete cQ; }
     VecType * getpsol() {return sol;}
     void solve();
     void rebuildSolver(double frequency);
     void scaleDisp(VecType &u);
     void scaleInvDisp(VecType &u);
     void scaleDisp(VecType &u, double alpha);
     void forceContinuity(VecType &u);
     void forceAssemble(VecType &u);
     void pade1(VecType *sol, VecType **sol_prev, double x);
     void pade(VecType *sol, VecType **u, double *h, double x);
     double galProjection(bool,int,VecType *sol, VecType **u, VecType **v,
                        Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                        Scalar **&VhK_arubber_lV, Scalar **&VhK_arubber_mV,
                        double w, double deltaw);
     double krylovGalProjection(int,int,VecType *sol, VecType **u, VecType **v,
                        Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                        double w, double deltaw);
     void adaptGP(int dgpFlag, int minRHS, int maxRHS, int deltaRHS, int &nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  VecType **aa, VecType **bb, VecType **cc,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  Scalar **&VhK_arubber_lV, Scalar **&VhK_arubber_mV,
                  int ncheck, double *wcheck, double wc,
                  double alpha, double tol);
     double adaptGPSolRes(int dgpFlag, int nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  VecType **aa, VecType **bb, VecType **cc,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  Scalar **&VhK_arubber_lV, Scalar **&VhK_arubber_mV,
                  double w, double deltaw);
     void adaptWindowSweep();
     void adaptSweep();

     void fourier(VecType *sol, VecType **u, double *h, double x);
     void stochStress(VecType *sol);
 
     void PadeLanczos_BuildSubspace(int nRHS, VecType **sol_prev, VecType **OrthoVec, int numOrthoVec);
     void PadeLanczos_Evaluate(int nRHS, VecType **BasisVector, double* &VtKV, double* &Vtb, double w, VecType *sol, double shiftForK = 0.0);
};

struct FreqP {
 int findex;
 int pindex;
};

template <class Scalar,
          class VecType,
          class ProblemDescriptor>
struct FWindowSData { 
  int n0;
  int n;
  int imode;
  int irow;
  VecType **u; // globally orthogonalized / previous scaling space
  VecType **ups; // globally orthogonalized / previous scaling space 
  VecType **v; // localy orthogonalized /
               // scaling based on current frequency
  VecType **Ku;
  VecType **Mu;
  VecType **Cu;
  FWindowSData(int _n0, int _imode, ProblemDescriptor* probDesc) {
    n0=_n0; imode = _imode; irow = -1; n = 0;
    u = new VecType * [n0];
    for(int i = 0; i < n0; ++i)
      u[i] = new VecType(probDesc->solVecInfo());
    ups = new VecType * [n0];
    for(int i = 0; i < n0; ++i)
      ups[i] = new VecType(probDesc->solVecInfo());
    v = new VecType * [n0];
    for(int i = 0; i < n0; ++i)
      v[i] = new VecType(probDesc->solVecInfo());
    Ku = new VecType * [n0];
    for(int i = 0; i < n0; ++i)
      Ku[i] = new VecType(probDesc->solVecInfo());
    Mu = new VecType * [n0];
    for(int i = 0; i < n0; ++i)
      Mu[i] = new VecType(probDesc->solVecInfo());
    Cu = new VecType * [n0];
    for(int i = 0; i < n0; ++i) {
      Cu[i] = new VecType(probDesc->solVecInfo());
      Cu[i]->zero();
    }
  }
  ~FWindowSData() {
    for(int i=0;i<n0;i++) delete u[i];
    delete[] u;
    for(int i=0;i<n0;i++) delete ups[i];
    delete[] ups;
    for(int i=0;i<n0;i++) delete v[i];
    delete[] v;
    for(int i=0;i<n0;i++) delete Ku[i];
    delete[] Ku;
    for(int i=0;i<n0;i++) delete Mu[i];
    delete[] Mu;
    for(int i=0;i<n0;i++) delete Cu[i];
    delete[] Cu;
  }
};



#ifdef _TEMPLATE_FIX_
#include <Driver.d/StaticProbType.C>
#include <Driver.d/Pade.C>
#include <Driver.d/PadeLanczos.C>
#endif

#endif
