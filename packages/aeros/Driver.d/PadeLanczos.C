
// The following specialization is for the Pade-Lanczos
// UH 05/20/08 
// revised by PJSA 10/15/08

extern int verboseFlag;


#include <cassert>
#include <cstdlib>
#include <sys/time.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/DistHelper.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <valarray>

extern "C" {
void _FORTRAN(dgetrf) (const int &m, const int &n,
                       double *a, const int &lda, int *ipiv, int &info);
void _FORTRAN(dgetrs) (const char & trans, const int &n, const int &nrhs,
                       double *a, const int &lda, int *ipiv,
                       double *b, const int &ldb, int &info);
void _FORTRAN(zgetrf) (const int &m, const int &n,
                      complex<double> *a, const int &lda, int *ipiv, int &info);
void _FORTRAN(zgetrs) (const char & trans, const int &n, const int &nrhs,
                       complex<double> *a, const int &lda, int *ipiv,
                       complex<double> *b, const int &ldb, int &info);
}
inline void Tgetrf (const int &m, const int &n,
                    double *a, const int &lda, int *ipiv, int &info) {
 _FORTRAN(dgetrf)(m,n,a,lda,ipiv,info);
}
inline void Tgetrf (const int &m, const int &n,
                    complex<double> *a, const int &lda, int *ipiv, int &info) {
 _FORTRAN(zgetrf)(m,n,a,lda,ipiv,info);
}
inline void Tgetrs (const char & trans, const int &n, const int &nrhs,
                       double *a, const int &lda, int *ipiv,
                       double *b, const int &ldb, int &info) {
 _FORTRAN(dgetrs)(trans,n,nrhs,a,lda,ipiv,b,ldb,info);
}
inline void Tgetrs (const char & trans, const int &n, const int &nrhs,
                       complex<double> *a, const int &lda, int *ipiv,
                       complex<double> *b, const int &ldb, int &info) {
 _FORTRAN(zgetrs)(trans,n,nrhs,a,lda,ipiv,b,ldb,info);
}


/*
extern "C" {

  void _FORTRAN(dgemm)(const char &, const char &, const int &,const int &,
                       const int &, const double &, double *, const int &,
                       double *, const int &, const double &, double *,
                       const int &);

  void _FORTRAN(zgemm)(const char &, const char &, const int &,const int &,
                       const int &, const complex<double> &, complex<double> *,
                       const int &, complex<double> *, const int &, 
                       const complex<double> &, complex<double> *, const int &);
}

inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const double &f, double *g, const int &h,
        double *i, const int &j, const double &k, double *l, const int &m)
{
 _FORTRAN(dgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}

inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const complex<double> &f, complex<double> *g,
                  const int &h, complex<double> *i, const int &j,
                  const complex<double> &k, complex<double> *l, const int &m)
{
 _FORTRAN(zgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}
*/



template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::PadeLanczos_BuildSubspace(int nRHS, VecType **sol_prev, VecType **OrthoVec, int numOrthoVec)
{
//
// This routine builds a M-orthonormal basis for the Krylov subspace
// { (K - sigma M)^{-1} rhs, (K - sigma M)^{-1} M (K - sigma M)^{-1} rhs, ... , 
//  [(K - sigma M)^{-1} M]^n (K - sigma M)^{-1} rhs}
//
//--- UH (05/22/08) The random value should be modified for complex problems.
//
// UH (05/21/2008)
//

  if (allOps->C_deriv) {
    std::cerr << "\n !!! The Pade-Lanczos algorithm is not implemented yet";
    std::cerr << " for the damped case !!!\n\n";
    exit(-1);
  }

  srand((unsigned int) time(0));

//  startTimerMemory(times->formRhs, times->memoryRhs);

  //------
  VecType *u = new VecType(probDesc->solVecInfo());
  VecType *w = new VecType(probDesc->solVecInfo());
  //------
  for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
    filePrint(stderr," ... Compute Basis Vector #%3d        ... \n",iRHS);
    //---
    if (iRHS == 0) {
      probDesc->getRHS(*u);
      allOps->sysSolver->solve(*u, *w);
    }
    else {
      //--- Get the previous vector
      *w = *(sol_prev[iRHS-1]);
      //--- Compute M * V(:, iRHS-1) 
      allOps->M->mult(*w, *u);
      if(domain->solInfo().modeFilterFlag) probDesc->project(*u); // PJSA 10-15-08
      //--- Solve the linear system
      allOps->sysSolver->solve(*u, *w);
    }
    //---
    bool hasExpansionForV = false;
    while (hasExpansionForV == false) {
      double mNorm = 0.0;
      std::vector<Scalar> beta(iRHS + numOrthoVec);
      //--- Perform two steps of Gram-Schmidt
      for (int iijj = 0; iijj < 2; ++iijj) {
        allOps->M->mult(*w, *u);
        mNorm = sqrt(ScalarTypes::Real( (*w) * (*u) ));
        (*w) *= 1.0 / mNorm;
        if ((iRHS == 0) && (numOrthoVec == 0)) {
          hasExpansionForV = true;
          break;
        }
        double invMNorm = 1.0 / mNorm;
        //--------
        for (int jj = 0; jj < numOrthoVec; ++jj) {
          beta[jj] = ( (*u) * (*(OrthoVec[jj])) ) * invMNorm;
        }
        for (int jj = 0; jj < iRHS; ++jj) {
          beta[jj + numOrthoVec] = ( (*u) * (*(sol_prev[jj])) ) * invMNorm;
        }
        //--------
        for (int jj = 0; jj < numOrthoVec; ++jj) {
          *u = *(OrthoVec[jj]);
          *u *= -beta[jj];
          (*w) += *u;
        }
        for (int jj = 0; jj < iRHS; ++jj) {
          *u = *(sol_prev[jj]);
          *u *= -beta[jj + numOrthoVec];
          (*w) += *u;
        }
      } // for (int iijj = 0; iijj < 2; ++iijj)
      //---
      allOps->M->mult(*w, *u);
      mNorm = sqrt(ScalarTypes::Real( (*w) * (*u) ));
      if (mNorm < 1e-15) {
        //--- Reset w to random values
        for (int iirr = 0; iirr < w->size(); ++iirr)
          (*w)[iirr] = ((double) rand()) / RAND_MAX;
        continue;
      }
      else {
        hasExpansionForV = true;
        *w *= 1.0 / mNorm;
      }
    } // while (hasExpansionForV == false)

    *(sol_prev[iRHS]) = *w;

  } // for (int iRHS = 0; iRHS < nRHS; ++iRHS)
  //------
  
  delete u;
  delete w;

//  stopTimerMemory(times->formRhs, times->memoryRhs);

//--- Check the M-orthonormality ---  
/*
for (int iVec = 0; iVec < numOrthoVec; ++iVec) {
  VecType *Mv = new VecType(probDesc->solVecInfo());
  allOps->M->mult(*(OrthoVec[iVec]), *Mv);
  for (int jVec = 0; jVec < numOrthoVec; ++jVec) {
    double dot = ScalarTypes::Real( (*OrthoVec[jVec]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  for (int jRHS = 0; jRHS < nRHS; ++jRHS) {
    double dot = ScalarTypes::Real( (*sol_prev[jRHS]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  std::cout << std::endl;
  delete Mv;
}
//---
for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
  VecType *Mv = new VecType(probDesc->solVecInfo());
  allOps->M->mult(*(sol_prev[iRHS]), *Mv);
  for (int jVec = 0; jVec < numOrthoVec; ++jVec) {
    double dot = ScalarTypes::Real( (*OrthoVec[jVec]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  for (int jRHS = 0; jRHS < nRHS; ++jRHS) {
    double dot = ScalarTypes::Real( (*sol_prev[jRHS]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  std::cout << std::endl;
  delete Mv;
}
*/
//----------------------------------

}


extern "C" {
extern void _FORTRAN(dsysv)(const char &, const int &, const int &,
                 double *, const int &, int *, double *, const int &,
                 double *, int &, int &);
extern void _FORTRAN(zsysv)(const char &, const int &, const int &,
                 complex<double> *, const int &, int *, complex<double> *,
                 const int &, complex<double> *, int &, int &);
void _FORTRAN(dspev)(const char &JOBZ, const char &UPLO,
                     const int &N, double *AP, double *W,
                     double *Z, const int &LDZ, double *WORK, int &INFO);
}

#ifndef _TSYSV__
#define _TSYSV__
inline void Tsysv(const char &a, const int &b, const int &c,
                 double *d, const int &e, int *f, double *g, const int &h,
                 double *i, int &j, int &k) {
 _FORTRAN(dsysv)(a,b,c,d,e,f,g,h,i,j,k);
}
inline void Tsysv(const char &a, const int &b, const int &c,
                 complex<double> *d, const int &e, int *f, complex<double> *g,
                 const int &h,
                 complex<double> *i, int &j, int &k) {
 _FORTRAN(zsysv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif
/*
void _FORTRAN(dgesv)(const int &, const int &, double *, const int &, int *, double *,
                       const int &, int &);

  void _FORTRAN(zgesv)(const int &, const int &, complex<double> *, const int &, int *, complex<double> *,
                       const int &, int &);


inline void Tgesv(const int &a, const int &b, double *c, const int &d, int *e, double *f,
                  const int &g, int &h)
{
 _FORTRAN(dgesv)(a,b,c,d,e,f,g,h);
}

inline void Tgesv(const int &a, const int &b, complex<double> *c, const int &d, int *e, complex<double> *f,
                  const int &g, int &h)
{
 _FORTRAN(zgesv)(a,b,c,d,e,f,g,h);
}*/




template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::PadeLanczos_Evaluate(int nRHS, VecType **BasisVector, double* &VtKV, double* &Vtb, double w, VecType *sol, double shiftForK)
{
// 
// This routine builds a Pade approximation to H(w) = (K - w^2 M)^{-1} b
// We assume that the basis stored in V is M-orthonormal.
// (The vectors of V are stored in BasisVector).
// The Pade approximation is H_{Pade}(w) = V * (V^T * K * V - w^2 * V^T * M * V)^{-1} * V^T * b
// Note that V^T M V is not computed and set to the identity matrix.
//
// On the first call of this routine, the matrix V^T K V and the vector V^T b are computed.
//
// UH (05/21/2008)
//

  if (VtKV == 0) {
    VtKV = new double[nRHS * nRHS];
    //----
    VecType *u = new VecType(probDesc->solVecInfo());
    for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
      //--
      if (allOps->K == 0) {
        std::cerr << " !!! Stiffness matrix is not created !!! " << std::endl;
        assert(0 > 1);
      }
      allOps->K->mult(*(BasisVector[iRHS]), *u);
      //--
      for (int jRHS = iRHS; jRHS < nRHS; ++jRHS) {
        VtKV[jRHS + iRHS * nRHS] = ScalarTypes::Real( (*u) * (*(BasisVector[jRHS])) );
        VtKV[iRHS + jRHS * nRHS] = VtKV[jRHS + iRHS * nRHS];
      } // for (int jRHS = 0; jRHS < nRHS; ++jRHS)
      //--
    } // for (int iRHS = 0; iRHS < nRHS; ++iRHS)
    //--- 
    //--- If allOps->K is assembled with a shift proportional to M,
    //--- then the next lines cancel the shift.
    //--- 
    for (int jj = 0; jj < nRHS; ++jj)
      VtKV[jj * (1 + nRHS)] += shiftForK;
    //----
    Vtb = new double[nRHS];
    probDesc->getRHS(*u);
    for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
      Vtb[iRHS] = ScalarTypes::Real( (*u) * (*(BasisVector[iRHS])) );
    }
    //----

////--- Check VtKV ---  
//for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
//  for (int jRHS = 0; jRHS < nRHS; ++jRHS) {
//    std::cout.precision(2);
//    std::cout << VtKV[iRHS + jRHS * nRHS] << " ";
//  }
//  std::cout << std::endl;
//}
////----------------------------------

    // PJSA 10-14-08 solve the reduced eigenvalue problem to get the poles of the Pade approximant
    bool xxxx = false; // change to true if eigenvectors are required
                       // this idea needs a bit of work and testing. esp output freqsweep modes and eigenmodes go to same file
    if(domain->solInfo().getSweepParams()->pade_poles && nRHS == domain->solInfo().getSweepParams()->nFreqSweepRHS*domain->solInfo().getSweepParams()->padeN) {
      char jobz = (xxxx) ? 'V' : 'N';
      double *ap = new double[nRHS*(nRHS+1)/2];
      double *w = new double[nRHS];
      double *z = (xxxx) ? new double[nRHS*nRHS] : 0;
      double *work = new double[3*nRHS];
      int info;
      // fill the ap matrix
      int index = 0;
      for(int j = 0; j < nRHS; ++j)
        for(int i = 0; i <= j; ++i)
          ap[index++] = VtKV[i + j*nRHS];

      _FORTRAN(dspev)(jobz, 'U', nRHS, ap, w, z, nRHS, work, info);

      if(info < 0) filePrint(stderr, "Error in dspev: the %d-th argument had an illegal value.\n", -info);
      else if(info > 0) filePrint(stderr, "Warning in dspev: the algorithm failed to converge; %d\n"
                                          "                  off-diagonal elements of an intermediate tridiagonal"
                                          "                  form did not converge to zero.\n", info);
      filePrint(stderr," --------------------------------------\n");
      filePrint(stderr," Mode\tPoles of Pade Approximant\n");
      filePrint(stderr," --------------------------------------\n");
      int imode = 0;
      for(int i = 0; i < nRHS; ++i) {
        if(w[i] >= domain->solInfo().getSweepParams()->pade_poles_sigmaL && w[i] <= domain->solInfo().getSweepParams()->pade_poles_sigmaU) {
          filePrint(stderr, " %d\t%e\n", ++imode, w[i]);
          if(xxxx) {
            u->zero();
            for(int j = 0; j < nRHS; ++j) u->linAdd(z[nRHS*i+j], (*(BasisVector[j])));
            //XXXX domain->postProcessing<Scalar>(*u, (Scalar *)0, *u, 0, w[i]); 
          }
        }
      }
      filePrint(stderr," --------------------------------------\n");
      delete [] ap; delete [] w; delete [] work;
    }

    delete u;
  }

  //--- Shift the projected stiffness matrix
  std::vector<double> copyVtKV(nRHS * nRHS);
  for (int ii = 0; ii < nRHS * nRHS; ++ii)
    copyVtKV[ii] = VtKV[ii];
  for (int jj = 0; jj < nRHS; ++jj)
    copyVtKV[jj * (1 + nRHS)] += - w * w;

  std::vector<double> alpha(nRHS);
  for (int ii = 0; ii < nRHS; ++ii)
    alpha[ii] = Vtb[ii];

/*
      FILE *ff = fopen("pm","a");
      for(int ii=0;ii<nRHS;ii++) for(int jj=0;jj<nRHS;jj++) {
        fprintf(ff,"%d %d %.20e %.20e\n",ii+1,jj+1,
           ScalarTypes::Real(copyVtKV[ii*(nRHS)+jj]),
           ScalarTypes::Imag(copyVtKV[ii*(nRHS)+jj]));
      }
      fclose(ff);
*/


  //--- Solve the reduced linear system
  std::vector<int> ipiv(nRHS);
  int lwork = 3 * nRHS;
  std::vector<double> work(lwork);
  int info = 0;
  _FORTRAN(dsysv)('U', nRHS, 1, &copyVtKV[0], nRHS, &ipiv[0], &alpha[0], nRHS, &work[0], lwork, info);

//for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
//  std::cout << alpha[iRHS] << std::endl;
//}

  //--- Compute the Pade approximant vector
  for (int jj = 0; jj < sol->size(); ++jj) {
    Scalar value = 0.0;
    for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
      value += (*(BasisVector[iRHS]))[jj] * alpha[iRHS];
    } // for (int iRHS = 0; iRHS < nRHS; ++iRHS)
    (*sol)[jj] = value;
  } // for (int jj = 0; jj < sol->size(); ++jj)
    
}




template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
double
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::galProjection(bool gpReorthoFlag,
                  int nRHS, VecType *sol, VecType **u, VecType **v,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  Scalar **&VhK_arubber_lV, Scalar **&VhK_arubber_mV,
                  double w, double deltaw)
{
 if (verboseFlag) filePrint(stderr,"w deltaw size: %f %f %d %d\n",w,deltaw,sol->size(),int(gpReorthoFlag));

 VecType *f = new VecType(probDesc->solVecInfo());
 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());
 int num_arubber = geoSource->num_arubber;
 VecType *dl = (num_arubber>0)?new VecType(probDesc->solVecInfo()):0;
 VecType *dm = (num_arubber>0)?new VecType(probDesc->solVecInfo()):0;

 if (gpReorthoFlag) {

double time = 0.0;
time -= getTime();

   for(int i=0;i<nRHS;i++) {
     *u[i] = *v[i];
   }

   if(domain->solInfo().isCoupled) {
     for(int i=0;i<nRHS;i++) {
        scaleInvDisp(*u[i]);
     }
   }

   for(int i=0;i<nRHS;i++) {
     if (isFeti(domain->solInfo().solvercntl->type)) {
       forceContinuity(*u[i]);
        if (verboseFlag) filePrint(stderr,"Forcing continuity %d.\n",i);
     }
   }



 // Orthogonalize
//   int ngs = 5;
   int ngs = 2;
   for(int m=0;m<ngs;m++) { 
     for(int i=0;i<nRHS;i++) {
       Scalar nrm0 = *u[i] * *u[i];
       for(int j=0;j<i;j++) {
         Scalar dotp = *u[i] *  *u[j];
         (*u[i]).linAdd(-dotp,*u[j]);
       }
       Scalar nrm = *u[i] * *u[i];
       *u[i] *= 1.0/sqrt(ScalarTypes::Real(nrm));
       nrm = *u[i] * *u[i];
     }
   }

   if (VhMV!=0) delete[] VhMV; 
   if (VhKV!=0) delete[] VhKV; 
   if (VhCV!=0) delete[] VhCV; 
   if (VhK_arubber_lV!=0) {
     for(int ir=0;ir<num_arubber;ir++) {
       delete VhK_arubber_lV[ir];
       delete VhK_arubber_mV[ir];
     }
     delete VhK_arubber_lV;
     delete VhK_arubber_mV;
   }

   VhMV = new Scalar[(nRHS)*(nRHS)];
   VhKV = new Scalar[(nRHS)*(nRHS)];
   VhCV = new Scalar[(nRHS)*(nRHS)];
   if (num_arubber>0) {
     VhK_arubber_lV = new Scalar*[num_arubber];
     VhK_arubber_mV = new Scalar*[num_arubber];
     for(int ir=0;ir<num_arubber;ir++) {
       VhK_arubber_lV[ir] = new Scalar[(nRHS)*(nRHS)];
       VhK_arubber_mV[ir] = new Scalar[(nRHS)*(nRHS)];
     }
   } else {
     VhK_arubber_lV = VhK_arubber_mV = 0;
   }
   
   for(int i=0;i<nRHS;i++) {
     
     allOps->K->mult(*(u[i]), *a);
     allOps->M->mult(*(u[i]), *b);

     for(int k=0;k<sol->size();k++) (*c)[k] = 0;
     if (allOps->C_deriv) {
       if (allOps->C_deriv[0]) { allOps->C_deriv[0]->mult(*(u[i]), *c);
       }
     } else (*c).zero(); 
     //else if (allOps->C) allOps->C->mult(*(u[i]), *c);  
     

     for(int j=0;j<nRHS;j++) {
       VhMV[i*(nRHS)+j] = *b * *u[j];
       VhCV[i*(nRHS)+j] = *c * *u[j];
       VhKV[i*(nRHS)+j] = *a * *u[j];
// + (w-deltaw)*(w-deltaw) * VhMV[i*(nRHS)+j];
//       ScalarTypes::addComplex(VhKV[i*(nRHS)+j], -(w-deltaw)*VhCV[i*(nRHS)+j] );
     }

     for(int ir=0;ir<num_arubber;ir++) {
       allOps->K_arubber_l[ir]->mult(*(u[i]), *dl);
       allOps->K_arubber_m[ir]->mult(*(u[i]), *dm);
       for(int j=0;j<nRHS;j++) {
         VhK_arubber_lV[ir][i*(nRHS)+j] = *dl * *u[j];
         VhK_arubber_mV[ir][i*(nRHS)+j] = *dm * *u[j];
       }
     }

   }
   time += getTime();
   if (verboseFlag) filePrint(stderr,"Ortho + proj setup time: %e\n",time/1e3);
 }

//#define ARUBBER_O
#ifdef ARUBBER_O
 allOps->K_deriv[0]->zeroAll();
 probDesc->buildDeltaK(w-deltaw,w);
 Scalar* VdeltaKV = new Scalar[(nRHS)*(nRHS)];
 for(int i=0;i<(nRHS)*(nRHS);i++) VdeltaKV[i] = 0.0;
 if (allOps->K_deriv) if (allOps->K_deriv[0]) {
   for(int i=0;i<nRHS;i++) {
     allOps->K_deriv[0]->mult(*(u[i]), *a);
     for(int j=0;j<nRHS;j++)  
       VdeltaKV[i*(nRHS)+j] = *a * *u[j];
     }
 }
#endif


double time = 0.0;
time -= getTime();

 // Project 
 probDesc->getRHS(*f, w,deltaw);

 Scalar *z = new Scalar[nRHS];
 for(int i=0;i<nRHS;i++) {
   z[i] = *f * *u[i];
 }

 Scalar *zz = new Scalar[(nRHS)*(nRHS)];
 for(int i=0;i<nRHS;i++)
   for(int j=0;j<nRHS;j++) {
     zz[i*(nRHS)+j] = VhKV[i*(nRHS)+j] - w*w * VhMV[i*(nRHS)+j];
#ifdef ARUBBER_O
     zz[i*(nRHS)+j] += VdeltaKV[i*(nRHS)+j];
#endif
     ScalarTypes::addComplex(zz[i*(nRHS)+j], 
                      w * VhCV[i*(nRHS)+j]);
   }

 complex<double> *lambda=0, *mu=0, *deltalambda=0, *deltamu=0;
 if (num_arubber>0)  {
   lambda = new complex<double>[num_arubber];
   mu = new complex<double>[num_arubber];
   deltalambda = new complex<double>[num_arubber];
   deltamu = new complex<double>[num_arubber];
 }
 geoSource->getARubberLambdaMu(w,deltalambda,deltamu);
 geoSource->getARubberLambdaMu(w-deltaw,lambda,mu);
 for(int ir=0;ir<num_arubber;ir++) {
  deltalambda[ir] -= lambda[ir];
  deltamu[ir] -= mu[ir];
 }

 for(int ir=0;ir<num_arubber;ir++)      
 for(int i=0;i<nRHS;i++)
   for(int j=0;j<nRHS;j++) 
    ScalarTypes::addComplex(zz[i*(nRHS)+j],
       deltalambda[ir]*VhK_arubber_lV[ir][i*(nRHS)+j] +
       deltamu[ir]*VhK_arubber_mV[ir][i*(nRHS)+j]);

/*
 {
 int ir=0;
 for(int i=0;i<nRHS*nRHS;i++) fprintf(stderr,"VdeltaKV %d %e %e   %e %e  %e\n",
      i,ScalarTypes::Real(VdeltaKV[i]),ScalarTypes::Imag(VdeltaKV[i]),
       ScalarTypes::Real(deltalambda[ir]*VhK_arubber_lV[ir][i] +
       deltamu[ir]*VhK_arubber_mV[ir][i]),
       ScalarTypes::Imag(deltalambda[ir]*VhK_arubber_lV[ir][i] +
       deltamu[ir]*VhK_arubber_mV[ir][i]),
       ScalarTypes::Real(VhKV[i]));
 }
*/

//--- Solve the reduced linear system
 std::vector<int> ipiv(nRHS);
 int lwork = 3 * nRHS;
 std::vector<Scalar> work(lwork);
 int info = 0;
 Tgesv(nRHS, 1, &zz[0], nRHS, &ipiv[0], &z[0], nRHS, info);

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nRHS;i++) 
   (*sol).linAdd(z[i],*u[i]);

 time += getTime();
 filePrint(stderr,"Projection time: %e\n",time/1e3);

time = 0.0;
time -= getTime();
     Scalar nrmb = 0.0;
 allOps->K->mult(*sol, *a);
 allOps->M->mult(*sol, *b);
 if (allOps->C_deriv) {
    if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
 } else (*c).zero(); 
 (*a).linAdd(-w*w,*b);
 for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k], 
                      w * (*c)[k]);
#ifndef ARUBBER_O
 for(int ir=0;ir<num_arubber;ir++) {
   allOps->K_arubber_l[ir]->mult(*sol, *dl);
   allOps->K_arubber_m[ir]->mult(*sol, *dm);
   for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k],
                    deltalambda[ir]* (*dl)[k] + deltamu[ir]* (*dm)[k]);
 }
#else
 if (allOps->K_deriv) if (allOps->K_deriv[0]) { 
     allOps->K_deriv[0]->mult(*sol, *b);
     nrmb = *b * *b;
    (*a).linAdd(1.0,*b);
 }
#endif

// (*a).print();
 (*a).linAdd(-1.0,*f);
 forceAssemble(*a);
 forceAssemble(*f);
 Scalar nrma = *a * *a;
 Scalar nrmf = *f * *f;
  if (verboseFlag) filePrint(stderr,"residual: %e\n", sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf)));
  if (verboseFlag) filePrint(stderr,"deltaK res: %e\n", sqrt(ScalarTypes::Real(nrmb))/sqrt(ScalarTypes::Real(nrmf)));
// (*sol).print();
time += getTime();
  filePrint(stderr,"Residual compute time: %e\n",time/1e3);


 delete a;
 delete b;
 delete c;
 if (dl) delete dl;
 if (dm) delete dm;
 if (num_arubber>0)  {
   delete[] lambda;
   delete[] mu;
   delete[] deltalambda;
   delete[] deltamu;
 }


 delete f;

 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 delete[] z;
 delete[] zz;
#ifdef ARUBBER_O
 delete[] VdeltaKV;
#endif

 return sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
double
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::krylovGalProjection(int nRHS, int nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  double w, double deltaw)
{
 if(verboseFlag) filePrint(stderr,"w deltaw size: %f %f %d %d %d\n",w,deltaw,nRHS,nOrtho,sol->size());

 VecType *f = new VecType(probDesc->solVecInfo());

 if (nRHS>0) {
double time = 0.0;
double rhstime = 0.0;
time -= getTime();
   VecType *a = new VecType(probDesc->solVecInfo());

   for(int i=0;i<nOrtho;i++) *u[i] = *v[i];
   if(domain->solInfo().isCoupled) {
     for(int i=0;i<nOrtho;i++)
        scaleInvDisp(*u[i]);
     for(int i=0;i<nOrtho;i++) {
       *a = *u[i];
       int ngs = 2;
       for(int m=0;m<ngs;m++) { 
          for(int j=0;j<i;j++) {
            Scalar dotp = *a *  *u[j];
            (*a).linAdd(-dotp,*u[j]);
          }
          Scalar nrm = *a * *a;
          *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
       }
       *u[i] = *a;
     }
   }

   for(int i=0;i<nRHS;i++) {
     if (i == 0) {
       probDesc->getRHS(*f);
     }
     else {
       allOps->M->mult(*u[nOrtho+i-1], *f);
     }
     if(verboseFlag) filePrint(stderr,"\n ... Solving RHS #%3d               ...\n",i);
time += getTime();
rhstime -= getTime();
     allOps->sysSolver->solve(*f, *a);
rhstime += getTime();
time -= getTime();
     int ngs = 2;
     for(int m=0;m<ngs;m++) { 
        for(int j=0;j<nOrtho+i;j++) {
          Scalar dotp = *a *  *u[j];
          (*a).linAdd(-dotp,*u[j]);
        }
        Scalar nrm = *a * *a;
        *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
     }
     *u[nOrtho+i] = *a;
   }

   VecType *b = new VecType(probDesc->solVecInfo());
   VecType *c = new VecType(probDesc->solVecInfo());

   if (VhMV!=0) delete[] VhMV; 
   if (VhKV!=0) delete[] VhKV; 
   if (VhCV!=0) delete[] VhCV; 
   VhMV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
   VhKV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
   VhCV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];

   for(int i=0;i<nRHS+nOrtho;i++) {
     allOps->K->mult(*(u[i]), *a);
     allOps->M->mult(*(u[i]), *b);
     for(int k=0;k<sol->size();k++) (*c)[k] = 0;
     if (allOps->C_deriv) {
       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *c);
     } else (*c).zero(); 
     //else if (allOps->C) allOps->C->mult(*(u[i]), *c);  

     for(int j=0;j<nOrtho+nRHS;j++) {
       VhMV[i*(nOrtho+nRHS)+j] = *b * *u[j];
       VhCV[i*(nOrtho+nRHS)+j] = *c * *u[j];
       VhKV[i*(nOrtho+nRHS)+j] = *a * *u[j];
     }
   }

   delete b;
   delete c;
   delete a;

   for(int i=0;i<nOrtho+nRHS;i++) *v[i] = *u[i];
   if(domain->solInfo().isCoupled) {
     for(int i=0;i<nOrtho+nRHS;i++)
        scaleDisp(*v[i]);
   }
   time += getTime();
   if(verboseFlag) filePrint(stderr,"Krylov setup time: %e\n",time/1e3);
   if(verboseFlag) filePrint(stderr,"RHS solve time: %e\n",rhstime/1e3);
 }

double time = 0.0;
time -= getTime();
 // Project 
 probDesc->getRHS(*f, w,deltaw);

 Scalar *z = new Scalar[nOrtho+nRHS];
 for(int i=0;i<nOrtho+nRHS;i++) {
   z[i] = *f * *u[i];
 }

 Scalar *zz = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
 for(int i=0;i<nOrtho+nRHS;i++)
   for(int j=0;j<nOrtho+nRHS;j++) {
     zz[i*(nOrtho+nRHS)+j] = VhKV[i*(nOrtho+nRHS)+j] - w*w * VhMV[i*(nOrtho+nRHS)+j];
     ScalarTypes::addComplex(zz[i*(nOrtho+nRHS)+j], 
                      w * VhCV[i*(nOrtho+nRHS)+j]);
   }

//--- Solve the reduced linear system
 std::vector<int> ipiv(nOrtho+nRHS);
 int lwork = 3 * (nOrtho+nRHS);
 std::vector<Scalar> work(lwork);
 int info = 0;
 Tgesv(nOrtho+nRHS, 1, &zz[0], nOrtho+nRHS, &ipiv[0], &z[0], nOrtho+nRHS, info);

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nOrtho+nRHS;i++) 
   (*sol).linAdd(z[i],*u[i]);
 time += getTime();
  if (verboseFlag) filePrint(stderr,"Projection time: %e\n",time/1e3);

time = 0.0;
time -= getTime();
 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());
 allOps->K->mult(*sol, *a);
 allOps->M->mult(*sol, *b);
 if (allOps->C_deriv) {
    if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
 } else (*c).zero(); 
 (*a).linAdd(-w*w,*b);
 for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k], 
                      w * (*c)[k]);
 (*a).linAdd(-1.0,*f);
 forceAssemble(*a);
 forceAssemble(*f);
 Scalar nrma = *a * *a;
 Scalar nrmf = *f * *f;
  if (verboseFlag) filePrint(stderr,"residual: %e\n", sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf)));
 time += getTime();
  if (verboseFlag) filePrint(stderr,"Residual compute time: %e\n",time/1e3);

 delete a;
 delete b;
 delete c;

 delete f;


 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 delete[] z;
 delete[] zz;

 return sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void 
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::adaptGP(int dgpFlag, int minRHS, int maxRHS, int deltaRHS, int &nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  VecType **aa, VecType **bb, VecType **cc, 
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  Scalar **&VhK_arubber_lV, Scalar **&VhK_arubber_mV,
                  int ncheck, double *wcheck, double wc,
                  double alpha, double tol)
{
// filePrint(stderr,"minRHS,maxRHS,deltaRHS,nOrtho...: %d %d %d  %d %d\n",minRHS,maxRHS,deltaRHS,nOrtho,sol->size());

double time = 0.0;
double rhstime = 0.0;
double projmattime = 0.0;
double projmattime2 = 0.0;
double orthotime = 0.0;
double chrestime = 0.0;

time -= getTime();

 for(int i=0;i<nOrtho;i++) {
   *u[i] = *v[i];
   if(domain->solInfo().isCoupled) scaleDisp(*u[i],1.0/alpha);
 }

 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());
 VecType *f = new VecType(probDesc->solVecInfo());
 int num_arubber = geoSource->num_arubber;
 VecType *dl = (num_arubber>0)?new VecType(probDesc->solVecInfo()):0;
 VecType *dm = (num_arubber>0)?new VecType(probDesc->solVecInfo()):0;

 if (VhKV!=0) delete[] VhKV; 
 if (VhMV!=0) delete[] VhMV; 
 if (VhCV!=0) delete[] VhCV; 
 if (VhK_arubber_lV!=0) {
   for(int ir=0;ir<num_arubber;ir++) {
     delete VhK_arubber_lV[ir];
     delete VhK_arubber_mV[ir];
   }
   delete VhK_arubber_lV;
   delete VhK_arubber_mV;
 }
 Scalar *tmpVhKV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 Scalar *tmpVhMV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 Scalar *tmpVhCV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 Scalar **tmpVhK_arubber_lV = 0;
 Scalar **tmpVhK_arubber_mV = 0;
 if (num_arubber>0) {
   tmpVhK_arubber_lV = new Scalar*[num_arubber];
   tmpVhK_arubber_mV = new Scalar*[num_arubber];
   for(int ir=0;ir<num_arubber;ir++) {
     tmpVhK_arubber_lV[ir] = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
     tmpVhK_arubber_mV[ir] = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
   }
 } else {
   tmpVhK_arubber_lV = tmpVhK_arubber_mV = 0;
 }


 projmattime -= getTime();
 for(int i=0;i<nOrtho;i++) {
   allOps->K->mult(*(u[i]), *(aa[i]));
   allOps->M->mult(*(u[i]), *(bb[i]));
   (*c).zero();
   if (allOps->C_deriv) 
     if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *(cc[i]));
   for(int j=0;j<nOrtho;j++) {
     tmpVhKV[i*(nOrtho+maxRHS)+j] = *(aa[i]) * *(u[j]);
     tmpVhMV[i*(nOrtho+maxRHS)+j] = *(bb[i]) * *(u[j]);
     tmpVhCV[i*(nOrtho+maxRHS)+j] = *(cc[i]) * *(u[j]);

     for(int ir=0;ir<num_arubber;ir++) {
       allOps->K_arubber_l[ir]->mult(*(u[i]), *dl);
       allOps->K_arubber_m[ir]->mult(*(u[i]), *dm);
       for(int j=0;j<nOrtho;j++) {
         tmpVhK_arubber_lV[ir][i*(nOrtho+maxRHS)+j] = *dl * *u[j];
         tmpVhK_arubber_mV[ir][i*(nOrtho+maxRHS)+j] = *dm * *u[j];
       }
     }
   }
 }
 projmattime += getTime();

 Scalar *z = new Scalar[nOrtho+maxRHS];
 double *rescheck = new double[ncheck];
 double *oldrescheck = new double[ncheck];
 for(int i=0;i<ncheck;i++) oldrescheck[i] = 0.0;
 bool done = false;
 int lRHS = 0;
 int uRHS = minRHS;
 double maxres = 0.0;
 double oldmaxres = 0.0;

 Scalar *UU = 0;
 VecType **wcawe_u = 0 ;
 if (dgpFlag==2) {
   UU = new Scalar[maxRHS*maxRHS]; 
   wcawe_u = new VecType * [maxRHS];
   for(int i = 0; i < maxRHS; ++i) 
     wcawe_u[i] = new VecType(probDesc->solVecInfo());
   for(int i = 0; i < maxRHS*maxRHS; ++i) UU[i] = 0.0;
 }

 bool breakDown = false; 
 while (!done) {
 
// Add new vectors 
   orthotime -= getTime();
   int offset = nOrtho+maxRHS;
   for(int i=lRHS;i<uRHS;i++) {
     if (i == 0) probDesc->getRHS(*f);
     if (dgpFlag==1) {
       if (i == 0) {
         u[offset]->zero();   // ????
       } else {
         f->zero();
         probDesc->getFreqSweepRHS(f, u+offset, i);
       }
     } else if (dgpFlag==2) {
// WCAWE
       if (i>0) { 
         int ii = i+1;
      // for j=1:ii-1;
         std::vector<Scalar> pb(ii-1); 
         for(int j=1;j<ii;j++) {
//fprintf(stderr,"wcawe %d %d\n",ii,j);
           int iimj = ii-j;
      //     PP = diag(ones(ii-j,1));
           std::vector<Scalar> PP(iimj*iimj,0);
           for (int k=0; k<iimj; k++) PP[k+k*iimj] = 1.0;
      //   for k=1:j;
           for(int k=1;k<=j;k++) {
      //     PP = PP*inv(UU(k:(ii-j+k-1),k:(ii-j+k-1)));
      // Compute the inverse
            std::vector<Scalar> QQ(iimj*iimj,0);
            for (int m=k; m<iimj+k; m++) for (int n=k; n<iimj+k; n++)
              QQ[(m-k)+(n-k)*iimj] = UU[m-1+(n-1)*(maxRHS)];
            std::vector<Scalar> RR(iimj*iimj,0);
            for (int k=0; k<iimj; k++) RR[k+k*iimj] = 1.0;
            std::vector<int> ipiv(iimj);
            int info = 0;
            Tgesv(iimj, iimj, &QQ[0], iimj, &ipiv[0], &RR[0], iimj, info);
      // Multiply the matrices PP and RR and store in PP
            for (int kk=0; kk<iimj*iimj; kk++) QQ[kk] = PP[kk];
            char transa = 'N'; char transb = 'N';
            Tgemm(transa, transb, iimj, iimj, iimj,
                   1.0 , &QQ[0], iimj, &RR[0], iimj, 0.0, &PP[0], iimj);
      //    end
           }
      //   pu = pu + b{j+1}*PP(1,ii-j);
           pb[j-1] = PP[1-1+(iimj-1)*iimj];
      //    end
         } 
      
      //  end
      // Identical as above except j and k start from 2
         std::vector<Scalar> pU(ii*ii,0); 
      //  for j=2:ii-1;
         for(int j=2;j<ii;j++) {
           int iimj = ii-j;
      //    PP = diag(ones(ii-j,1));
           std::vector<Scalar> PP(iimj*iimj,0);
           for (int k=0; k<iimj; k++) PP[k+k*iimj] = 1.0;
      //    for k=2:j;
           for(int k=2;k<=j;k++) {
      //     PP = PP*inv(UU(k:(ii-j+k-1),k:(ii-j+k-1)));
      // Compute the inverse
            std::vector<Scalar> QQ(iimj*iimj,0);
            for (int m=k; m<iimj+k; m++) for (int n=k; n<iimj+k; n++)
              QQ[(m-k)+(n-k)*iimj] = UU[m-1+(n-1)*(maxRHS)];
            std::vector<Scalar> RR(iimj*iimj,0);
            for (int k=0; k<iimj; k++) RR[k+k*iimj] = 1.0;
            std::vector<int> ipiv(iimj);
            int info = 0;
            Tgesv(iimj, iimj, &QQ[0], iimj, &ipiv[0], &RR[0], iimj, info);
      // Multiply the matrices PP and RR and store in PP
            for (int kk=0; kk<iimj*iimj; kk++) QQ[kk] = PP[kk]; ;
            char transa = 'N'; char transb = 'N';
            Tgemm(transa, transb, iimj, iimj, iimj,
                   1.0 , &QQ[0], iimj, &RR[0], iimj, 0.0, &PP[0], iimj);
      //    end
           }
      //    pu = pu - Z{j+1}*W(:,1:(ii-j))*PP(:,ii-j);
           for(int k=0;k<iimj;k++)  {
             pU[k+(j-1)*ii] = PP[k+(iimj-1)*iimj];  
           }
      //    end
         } 
         probDesc->getWCAWEFreqSweepRHS(f, wcawe_u, &pU[0], &pb[0], maxRHS, i);
       }
     } else {
// KGP
       if (i>0) allOps->M->mult(*u[nOrtho+i-1], *f);
     }
      if (verboseFlag) filePrint(stderr,"\n ... Solving RHSa   #%3d               ...\n",i);
rhstime -= getTime();
     allOps->sysSolver->solve(*f, *a);
rhstime += getTime();
     if (dgpFlag) {
       if (isFeti(domain->solInfo().solvercntl->type)) {
//         filePrint(stderr,"Forcing continuity %d.\n",i);
         forceContinuity(*a);
       }
       *u[offset+i+1] = *a;
     }

     if (dgpFlag==2) {
// WCAWE orthogonalize
       *b = *a; 
        Scalar nrm0 = *a * *a;
       for(int j=0;j<i;j++) {
          Scalar dotp = *b *  *wcawe_u[j];
          UU[j+i*maxRHS] = dotp;
          (*b).linAdd(-dotp,*wcawe_u[j]);
       }
       UU[i*maxRHS+i] = sqrt(ScalarTypes::Real(*b * *b));
       *b *= 1.0/UU[i*maxRHS+i];
       *wcawe_u[i] = *b;
        double normRed = 
          ScalarTypes::Real(UU[i*maxRHS+i])/sqrt(ScalarTypes::Real(nrm0));
//filePrint(stderr,"norm reduction - wcawe space: %e\n",normRed);
        if (normRed<1e-17/tol) { breakDown = true; uRHS = i; break; }
     }

// GP space orthogonalize
     if(domain->solInfo().isCoupled) scaleDisp(*a,alpha);
     Scalar nrm0 = *a * *a;
     int ngs = 2;
     for(int m=0;m<ngs;m++) { 
        for(int j=0;j<nOrtho+i;j++) {
          Scalar dotp = *a *  *v[j];
          (*a).linAdd(-dotp,*v[j]);
        }
     }
     Scalar nrm = *a * *a;
     *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
      double normRed = 
        sqrt(ScalarTypes::Real(nrm))/sqrt(ScalarTypes::Real(nrm0));
//filePrint(stderr,"norm reduction - global space: %e\n",normRed);
        if (normRed<1e-20/tol && dgpFlag==2) { breakDown = true; uRHS = i; break; }
     *v[nOrtho+i] = *a;
     *u[nOrtho+i] = *a;
      if(domain->solInfo().isCoupled) scaleDisp(*u[nOrtho+i],1.0/alpha);
   }
   orthotime += getTime();


// Update matrices
   projmattime2 -= getTime();
   for(int i=nOrtho+lRHS;i<nOrtho+uRHS;i++) {
     allOps->K->mult(*(u[i]), *(aa[i]));
     allOps->M->mult(*(u[i]), *(bb[i]));
     (*c).zero();
     if (allOps->C_deriv) 
       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *(cc[i]));
     for(int j=0;j<nOrtho+uRHS;j++) {
       tmpVhKV[i*(nOrtho+maxRHS)+j] = *(aa[i]) * *(u[j]);
       tmpVhMV[i*(nOrtho+maxRHS)+j] = *(bb[i]) * *(u[j]);
       tmpVhCV[i*(nOrtho+maxRHS)+j] = *(cc[i]) * *(u[j]);
     }
     for(int ir=0;ir<num_arubber;ir++) {
       allOps->K_arubber_l[ir]->mult(*(u[i]), *dl);
       allOps->K_arubber_m[ir]->mult(*(u[i]), *dm);
       for(int j=0;j<nOrtho+uRHS;j++) {
         tmpVhK_arubber_lV[ir][i*(nOrtho+maxRHS)+j] = *dl * *u[j];
         tmpVhK_arubber_mV[ir][i*(nOrtho+maxRHS)+j] = *dm * *u[j];
       }
     }
   }
   for(int i=0;i<nOrtho+uRHS;i++) {
     for(int j=nOrtho+lRHS;j<nOrtho+uRHS;j++) {
       tmpVhKV[i*(nOrtho+maxRHS)+j] = *(aa[i]) * *(u[j]);
       tmpVhMV[i*(nOrtho+maxRHS)+j] = *(bb[i]) * *(u[j]);
       tmpVhCV[i*(nOrtho+maxRHS)+j] = *(cc[i]) * *(u[j]);
     }
     for(int ir=0;ir<num_arubber;ir++) {
       allOps->K_arubber_l[ir]->mult(*(u[i]), *dl);
       allOps->K_arubber_m[ir]->mult(*(u[i]), *dm);
       for(int j=nOrtho+lRHS;j<nOrtho+uRHS;j++) {
         tmpVhK_arubber_lV[ir][i*(nOrtho+maxRHS)+j] = *dl * *u[j];
         tmpVhK_arubber_mV[ir][i*(nOrtho+maxRHS)+j] = *dm * *u[j];
       }
     }
   }
   projmattime2 += getTime();

//Check residual at check points
   if (breakDown) { done = true; lRHS = uRHS; }
   else {
   chrestime -= getTime();
   Scalar *zz = new Scalar[(nOrtho+uRHS)*(nOrtho+uRHS)];
   for(int icheck=0;icheck<ncheck;icheck++) {
     probDesc->getRHS(*f, wcheck[icheck],wcheck[icheck]-wc);
     for(int i=0;i<nOrtho+uRHS;i++) {
       z[i] = *f * *u[i];
     }
     for(int i=0;i<nOrtho+uRHS;i++) {
       for(int j=0;j<nOrtho+uRHS;j++) {
         zz[i*(nOrtho+uRHS)+j] = tmpVhKV[i*(nOrtho+maxRHS)+j] -
                  wcheck[icheck]*wcheck[icheck] * tmpVhMV[i*(nOrtho+maxRHS)+j];
         ScalarTypes::addComplex(zz[i*(nOrtho+uRHS)+j], 
                    wcheck[icheck] * tmpVhCV[i*(nOrtho+maxRHS)+j]);
       }
     }
     complex<double> *lambda=0, *mu=0, *deltalambda=0, *deltamu=0;
     if (num_arubber>0)  {
       lambda = new complex<double>[num_arubber];
       mu = new complex<double>[num_arubber];
       deltalambda = new complex<double>[num_arubber];
       deltamu = new complex<double>[num_arubber];
       geoSource->getARubberLambdaMu(wcheck[icheck],deltalambda,deltamu);
       geoSource->getARubberLambdaMu(wc,lambda,mu);
     }
     for(int ir=0;ir<num_arubber;ir++) {
      deltalambda[ir] -= lambda[ir];
      deltamu[ir] -= mu[ir];
     }
    
     for(int ir=0;ir<num_arubber;ir++)      
     for(int i=0;i<nOrtho+uRHS;i++)
       for(int j=0;j<nOrtho+uRHS;j++) 
        ScalarTypes::addComplex(zz[i*(nOrtho+uRHS)+j],
           deltalambda[ir]*tmpVhK_arubber_lV[ir][i*(nOrtho+maxRHS)+j] +
           deltamu[ir]*tmpVhK_arubber_mV[ir][i*(nOrtho+maxRHS)+j]);


     std::vector<int> ipiv(nOrtho+uRHS);
     int lwork = 3 * (nOrtho+uRHS);
     std::vector<Scalar> work(lwork);
     int info = 0;
     Tgesv(nOrtho+uRHS, 1, &zz[0], nOrtho+uRHS, &ipiv[0], &z[0],
           nOrtho+uRHS, info);

     (*a).zero();
     for(int i=0;i<nOrtho+uRHS;i++) (*a).linAdd(z[i],*aa[i]);
     (*b).zero();
     for(int i=0;i<nOrtho+uRHS;i++) (*b).linAdd(z[i],*bb[i]);
     (*c).zero();
     for(int i=0;i<nOrtho+uRHS;i++) (*c).linAdd(z[i],*cc[i]);
/*
     (*sol).zero();
     for(int i=0;i<nOrtho+uRHS;i++) 
       (*sol).linAdd(z[i],*u[i]);
     allOps->K->mult(*sol, *a);
     allOps->M->mult(*sol, *b);
     if (allOps->C_deriv) {
        if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
     } else (*c).zero(); 
*/
     (*a).linAdd(-wcheck[icheck]*wcheck[icheck],*b);
     for(int k=0;k<sol->size();k++) 
       ScalarTypes::addComplex((*a)[k], wcheck[icheck] * (*c)[k]);

     (*sol).zero();
     for(int i=0;i<nOrtho+uRHS;i++) 
       (*sol).linAdd(z[i],*u[i]);
     for(int ir=0;ir<num_arubber;ir++) {
       allOps->K_arubber_l[ir]->mult(*sol, *dl);
       allOps->K_arubber_m[ir]->mult(*sol, *dm);
       for(int k=0;k<sol->size();k++) 
         ScalarTypes::addComplex((*a)[k],
                        deltalambda[ir]* (*dl)[k] + deltamu[ir]* (*dm)[k]);
     }

     (*a).linAdd(-1.0,*f);
     forceAssemble(*a);
     forceAssemble(*f);
     Scalar nrma = *a * *a;
     Scalar nrmf = *f * *f;
     rescheck[icheck] = 
      sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
     if (verboseFlag)
       filePrint(stderr,"rescheck[%d] = %e at %e\n",
                 icheck,rescheck[icheck],wcheck[icheck]);
   }
   delete[] zz;
   chrestime += getTime();

   maxres = 0.0;
   for(int i=0;i<ncheck;i++) if (rescheck[i]>maxres) maxres = rescheck[i];
   if (verboseFlag)
       filePrint(stderr,"maxres = %e, oldmaxres = %e  tol = %e \n",
             maxres,oldmaxres,tol);
   if (maxres<tol) done = true;
   int nsmaller=0;
//   for(int i=0;i<ncheck;i++) if (oldrescheck[i]!=0.0) if (rescheck[i]>tol) if (rescheck[i]<pow(0.8,deltaRHS)*oldrescheck[i]) nsmaller++; 
// For coupled
   for(int i=0;i<ncheck;i++) if (oldrescheck[i]!=0.0) if (rescheck[i]>tol) if (rescheck[i]<pow(0.9,deltaRHS)*oldrescheck[i]) nsmaller++; 
   if (oldmaxres!=0.0) if (nsmaller==0) done = true;
   for(int i=0;i<ncheck;i++) oldrescheck[i] = rescheck[i];
   oldmaxres = maxres;
   lRHS = uRHS;
   uRHS += deltaRHS;
   if (uRHS>maxRHS) done = true; 
   }
 }

 VhKV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 VhMV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 VhCV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 if (num_arubber>0) {
   VhK_arubber_lV = new Scalar*[num_arubber];
   VhK_arubber_mV = new Scalar*[num_arubber];
   for(int ir=0;ir<num_arubber;ir++) {
     VhK_arubber_lV[ir] = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
     VhK_arubber_mV[ir] = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
   }
 } else {
   VhK_arubber_lV = VhK_arubber_mV = 0;
 }


 for(int i=0;i<nOrtho+lRHS;i++) {
    for(int j=0;j<nOrtho+lRHS;j++) {
      VhKV[i*(nOrtho+lRHS)+j] = tmpVhKV[i*(nOrtho+maxRHS)+j];
      VhMV[i*(nOrtho+lRHS)+j] = tmpVhMV[i*(nOrtho+maxRHS)+j];
      VhCV[i*(nOrtho+lRHS)+j] = tmpVhCV[i*(nOrtho+maxRHS)+j];
    }
    for(int ir=0;ir<num_arubber;ir++) {
      for(int j=0;j<nOrtho+lRHS;j++) {
        VhK_arubber_lV[ir][i*(nOrtho+lRHS)+j] =
          tmpVhK_arubber_lV[ir][i*(nOrtho+maxRHS)+j];
        VhK_arubber_mV[ir][i*(nOrtho+lRHS)+j] =
          tmpVhK_arubber_mV[ir][i*(nOrtho+maxRHS)+j];
      }
    }
 }

 nOrtho += lRHS;
 if (verboseFlag) filePrint(stderr,"nOrtho = %d\n",nOrtho);
 
 delete[] tmpVhKV;
 delete[] tmpVhMV;
 delete[] tmpVhCV;
 if (num_arubber>0) {
   for(int ir=0;ir<num_arubber;ir++) {
     delete tmpVhK_arubber_lV[ir];
     delete tmpVhK_arubber_mV[ir];
   }
   delete tmpVhK_arubber_lV;
   delete tmpVhK_arubber_mV;
 }
 delete[] rescheck;
 delete[] oldrescheck;
 delete[] z;

 delete f;
 delete c;
 delete b;
 delete a;
 if (dl) delete dl;
 if (dm) delete dm;
 
 if (dgpFlag==2) {
   delete[] UU;
   for(int i=0;i<maxRHS;i++) delete wcawe_u[i];
   delete[] wcawe_u;
 }

 time += getTime();
 if (verboseFlag) {
   filePrint(stderr,"Total setup time: %e\n",time/1e3);
   filePrint(stderr,"Matrix setup time: %e\n",projmattime/1e3);
   filePrint(stderr,"Matrix setup time2: %e\n",projmattime2/1e3);
   filePrint(stderr,"Ortho+ time: %e\n",orthotime-rhstime/1e3);
   filePrint(stderr,"RHS solve time: %e\n",rhstime/1e3);
   filePrint(stderr,"Check res time: %e\n",chrestime/1e3);
 }
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
double
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::adaptGPSolRes(int dgpFlag, int nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  VecType **aa, VecType **bb, VecType **cc,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  Scalar **&VhK_arubber_lV, Scalar **&VhK_arubber_mV,
                  double w, double deltaw)
 //                  ,double alpha)
{
// filePrint(stderr,"w deltaw size: %f %f  %d %d\n",w,deltaw,nOrtho,sol->size());

double time = 0.0;
time -= getTime();

// This should be already done when this function is called
/*
 for(int i=0;i<nOrtho;i++) {
   *u[i] = *v[i];
   if(domain->solInfo().isCoupled) scaleDisp(*u[i],1.0/alpha);
 }
*/

 // Project 
 VecType *f = new VecType(probDesc->solVecInfo());
 probDesc->getRHS(*f, w,deltaw);

 Scalar *z = new Scalar[nOrtho];
 for(int i=0;i<nOrtho;i++) {
   z[i] = *f * *u[i];
 }

 Scalar *zz = new Scalar[(nOrtho)*(nOrtho)];
 for(int i=0;i<nOrtho;i++) {
   for(int j=0;j<nOrtho;j++) {
     zz[i*(nOrtho)+j] = VhKV[i*(nOrtho)+j] - w*w * VhMV[i*(nOrtho)+j];
     ScalarTypes::addComplex(zz[i*(nOrtho)+j], 
                      w * VhCV[i*(nOrtho)+j]);
   }
 }
 complex<double> *lambda=0, *mu=0, *deltalambda=0, *deltamu=0;
 int num_arubber = geoSource->num_arubber;
 VecType *dl = (num_arubber>0)?new VecType(probDesc->solVecInfo()):0;
 VecType *dm = (num_arubber>0)?new VecType(probDesc->solVecInfo()):0;
 if (num_arubber>0)  {
   lambda = new complex<double>[num_arubber];
   mu = new complex<double>[num_arubber];
   deltalambda = new complex<double>[num_arubber];
   deltamu = new complex<double>[num_arubber];
   geoSource->getARubberLambdaMu(w,deltalambda,deltamu);
   geoSource->getARubberLambdaMu(w-deltaw,lambda,mu);
 }
 for(int ir=0;ir<num_arubber;ir++) {
  deltalambda[ir] -= lambda[ir];
  deltamu[ir] -= mu[ir];
 }

 for(int ir=0;ir<num_arubber;ir++)      
 for(int i=0;i<nOrtho;i++)
   for(int j=0;j<nOrtho;j++) 
    ScalarTypes::addComplex(zz[i*(nOrtho)+j],
       deltalambda[ir]*VhK_arubber_lV[ir][i*(nOrtho)+j] +
       deltamu[ir]*VhK_arubber_mV[ir][i*(nOrtho)+j]);


//--- Solve the reduced linear system
 std::vector<int> ipiv(nOrtho);
 int lwork = 3 * (nOrtho);
 std::vector<Scalar> work(lwork);
 int info = 0;
 Tgesv(nOrtho, 1, &zz[0], nOrtho, &ipiv[0], &z[0], nOrtho, info);

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nOrtho;i++) 
   (*sol).linAdd(z[i],*u[i]);
  time += getTime();
  if (verboseFlag) filePrint(stderr,"Projection time: %e\n",time/1e3);

time = 0.0;
time -= getTime();
 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());

     (*a).zero();
     for(int i=0;i<nOrtho;i++) (*a).linAdd(z[i],*aa[i]);
     (*b).zero();
     for(int i=0;i<nOrtho;i++) (*b).linAdd(z[i],*bb[i]);
     (*c).zero();
     for(int i=0;i<nOrtho;i++) (*c).linAdd(z[i],*cc[i]);
/*
 allOps->K->mult(*sol, *a);
 allOps->M->mult(*sol, *b);
 if (allOps->C_deriv) {
    if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
 } else (*c).zero(); 
*/
 (*a).linAdd(-w*w,*b);
 for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k], w * (*c)[k]);

 for(int ir=0;ir<num_arubber;ir++) {
   allOps->K_arubber_l[ir]->mult(*sol, *dl);
   allOps->K_arubber_m[ir]->mult(*sol, *dm);
   for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k],
                    deltalambda[ir]* (*dl)[k] + deltamu[ir]* (*dm)[k]);
 }

 (*a).linAdd(-1.0,*f);
 forceAssemble(*a);
 forceAssemble(*f);
 Scalar nrma = *a * *a;
 Scalar nrmf = *f * *f;
 filePrint(stderr,"residual: %e\n", sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf)));
 time += getTime();
 if (verboseFlag) filePrint(stderr,"Residual compute time: %e\n",time/1e3);

 delete a;
 delete b;
 delete c;
 delete f;
 if (dl) delete dl;
 if (dm) delete dm;
 delete[] z;
 delete[] zz;

 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 return sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::adaptWindowSweep()
{
    // Adaptive window (of size 2 frequencies) sweep
    // with multiple right hand sides done locally or frequency-globally
    double time = 0.0;
    double rhstime = 0.0;
    double rhsctime1 = 0.0;
    double rhsctime2 = 0.0;
    double solutiontime = 0.0;
    double solvertime = 0.0;
    double outputtime = 0.0;
    double transitiontime= 0.0;
    double onefspacetime = 0.0;
    double fsweeptime = 0.0;
    double frsweeptime = 0.0;
    double Arsweeptime = 0.0;
    double ressweeptime = 0.0;
    double froneftime = 0.0;
    double Aroneftime = 0.0;
    double resoneftime = 0.0;
    int solvecount = 0;
    int rebuildscount = 0;


    domain->isCoarseGridSolve = false;

    time -= getTime();

    double w1 = domain->solInfo().getSweepParams()->adaptSweep.w1;
    double w2 = domain->solInfo().getSweepParams()->adaptSweep.w2;
    int numS = domain->solInfo().getSweepParams()->adaptSweep.numS;
    double dw = (w2-w1)/numS;

    int nDir = 1;
    bool isScattering = false;
    if (domain->numWaveDirections>=1 && domain->solInfo().loadcases.size()==0) {
        nDir = domain->numWaveDirections;
        isScattering = true;
        filePrint(stderr,"Solving %d scattering cases\n",nDir);
    }
    std::list<int> loadcases_backup;
    if (domain->solInfo().loadcases.size()>=1) {
        nDir = domain->solInfo().loadcases.size();
        loadcases_backup = domain->solInfo().loadcases;
        filePrint(stderr,"Solving %d loadcases\n",nDir);
    }
    if (nDir==0) nDir = 1;
    int maxLocalV = domain->solInfo().getSweepParams()->adaptSweep.minRHS;
    int dgpFlag = domain->solInfo().getSweepParams()->adaptSweep.dgp_flag;
    double tol = domain->solInfo().getSweepParams()->adaptSweep.atol;
    int localP = domain->solInfo().getSweepParams()->adaptSweep.maxP;
    // tolerance factor when looking for next coarse
    double ctolf = domain->solInfo().getSweepParams()->adaptSweep.ctolf;
    // tolerance factor when looking within one frequency
    double tol1f = domain->solInfo().getSweepParams()->adaptSweep.tol1f;
    if (ctolf==0.0) ctolf = 1e2;
    if (tol1f==0.0) tol1f = 1e-2;

    if (localP)
        filePrint(stderr,"Adaptive window sweep - local variant\n");
    else
        filePrint(stderr,"Adaptive window sweep\n");
    if (dgpFlag==1)
        filePrint(stderr,"DGP with %d local derivatives\n",maxLocalV);
    else if (dgpFlag==2)
        filePrint(stderr,"DGP with %d local WCAWE vectors\n",maxLocalV);
    else
        filePrint(stderr,"KGP with %d local vectors\n",maxLocalV);
    filePrint(stderr,"Parameter values: tol=%e ctolf=%e tol1f=%e\n",
              tol,ctolf,tol1f);

    // Temp storage for DGP and WCAWE
    VecType **temp_u=0;
    Scalar *UU = 0;
    if (dgpFlag) {
        UU = new Scalar[maxLocalV*maxLocalV];
        temp_u = new VecType * [maxLocalV+1];
        for(int i = 0; i < maxLocalV+1; ++i)
            temp_u[i] = new VecType(probDesc->solVecInfo());
        for(int i = 0; i < maxLocalV*maxLocalV; ++i) UU[i] = 0.0;
    }

    // Temp vectors
    VecType *f = new VecType(probDesc->solVecInfo());
    VecType *a = new VecType(probDesc->solVecInfo());
    VecType *fkgp = new VecType(probDesc->solVecInfo());
    VecType *akgp = new VecType(probDesc->solVecInfo());
    VecType *b = new VecType(probDesc->solVecInfo());
    VecType *c = new VecType(probDesc->solVecInfo());
    VecType *sol = new VecType(probDesc->solVecInfo());

    // Subspaces, 2 frequencies at a time
    const int nF = 2;
    FWindowSData<Scalar,VecType,ProblemDescriptor> **localFWSSpace =
        new FWindowSData<Scalar,VecType,ProblemDescriptor> *[nF*nDir];
    for(int i=0;i<nF*nDir;i++) {
        localFWSSpace[i] = 0;
    }

// vector of residuals, reduced matrix and vector, reduced matrices
    std::vector<double> res(nDir);
    std::vector<Scalar> fr(nF*nDir*maxLocalV);
    std::vector<Scalar> Ar(nF*nDir*maxLocalV*nF*nDir*maxLocalV);
    std::map<std::pair<int,int>,Scalar> Kr;
    std::map<std::pair<int,int>,Scalar> Mr;
    std::map<std::pair<int,int>,Scalar> Cr;

    int imode = 0;
    int icoarsew[nF];
    icoarsew[0] = -1;
    icoarsew[1] = 0;
    bool doneF = false;
    while (!doneF) {
        double w = w1 + icoarsew[1]*dw;
        // Scaling factor for coupled
        double alpha = 1.0; // ((w1+w2)/2.0)/w;
        if (icoarsew[0]!=-1) alpha = w/(w1 + icoarsew[0]*dw);
        // Setup solver for new frequency
        geoSource->setOmega(w);
        solvertime -= getTime();
        if (w>w1) {
            rebuildSolver(w);
            rebuildscount++;
        }
        solvertime += getTime();

        if (verboseFlag)
            filePrint(stderr, "icoarsew[0]=%d icoarsew[1]=%d alpha=%e isCoupled=%d\n",
                      icoarsew[0],icoarsew[1],alpha,int(domain->solInfo().isCoupled));

        // Erase old frequency spaces and
        // move new frequency spaces to old frequency spaces and update with
        // current scaling
        transitiontime -= getTime();
        for(int iDir=0;iDir<nDir;iDir++)  {
            if (localFWSSpace[iDir]!=0) delete localFWSSpace[iDir];
            localFWSSpace[iDir] = 0;
            if (localFWSSpace[nDir+iDir]!=0) {
                localFWSSpace[iDir] = localFWSSpace[nDir+iDir];
                localFWSSpace[nDir+iDir] = 0;
                // Update u's to current scaling and recompute Ku,Mu,Cu
                for(int i=0;i<localFWSSpace[iDir]->n;i++) {
                    *(localFWSSpace[iDir]->u[i]) = *(localFWSSpace[iDir]->v[i]);
                    *(localFWSSpace[iDir]->ups[i]) = *(localFWSSpace[iDir]->v[i]);
                    if(domain->solInfo().isCoupled)
                        scaleDisp(*(localFWSSpace[iDir]->u[i]),alpha);
                    allOps->K->mult(*(localFWSSpace[iDir]->u[i]),
                                    *(localFWSSpace[iDir]->Ku[i]));
                    allOps->M->mult(*(localFWSSpace[iDir]->u[i]),
                                    *(localFWSSpace[iDir]->Mu[i]));
                    if (allOps->C_deriv) if (allOps->C_deriv[0])
                            allOps->C_deriv[0]->mult(*(localFWSSpace[iDir]->u[i]),
                                                     *(localFWSSpace[iDir]->Cu[i]));
                }
            }
        }
        for(int iDir=0;iDir<nDir;iDir++)  {
            if (localFWSSpace[iDir]!=0) {
                // Recompute reduced matrices
                for(int i=0;i<localFWSSpace[iDir]->n;i++) {
                    // Compute new reduced matrix entries
                    int jDirL = (localP)?iDir:0;
                    int jDirU = (localP)?iDir+1:nDir;
                    for(int jDir=jDirL;jDir<jDirU;jDir++) {
                        if (localFWSSpace[jDir]!=0)
                            for(int j=0;j<localFWSSpace[jDir]->n;j++) {
                                std::pair<int,int> rowcol
                                    (localFWSSpace[jDir]->imode+j,localFWSSpace[iDir]->imode+i);
                                Kr[rowcol] = *(localFWSSpace[iDir]->Ku[i]) *
                                             *(localFWSSpace[jDir]->u[j]);
                                Mr[rowcol] = *(localFWSSpace[iDir]->Mu[i]) *
                                             *(localFWSSpace[jDir]->u[j]);
                                Cr[rowcol] = *(localFWSSpace[iDir]->Cu[i]) *
                                             *(localFWSSpace[jDir]->u[j]);
                                if (iDir!=jDir) {
                                    std::pair<int,int> rowcolt
                                        (localFWSSpace[iDir]->imode+i,
                                         localFWSSpace[jDir]->imode+j);
                                    Kr[rowcolt] = *(localFWSSpace[jDir]->Ku[j]) *
                                                  *(localFWSSpace[iDir]->u[i]);
                                    Mr[rowcolt] = *(localFWSSpace[jDir]->Mu[j]) *
                                                  *(localFWSSpace[iDir]->u[i]);
                                    Cr[rowcolt] = *(localFWSSpace[jDir]->Cu[j]) *
                                                  *(localFWSSpace[iDir]->u[i]);
                                }
                            }
                    }
                }
            }
        }
        transitiontime += getTime();

        // Count reduced variables and number them
        int nr = 0;
        for(int jDir=0;jDir<nF*nDir;jDir++) if (localFWSSpace[jDir]!=0) {
                localFWSSpace[jDir]->irow = nr;
                nr += localFWSSpace[jDir]->n;
            }

        if (localP)
            for(int jDir=0;jDir<nDir;jDir++) if (localFWSSpace[jDir]!=0)
                    localFWSSpace[jDir]->irow = 0;


        std::valarray<int> isDone(0,nDir);
        int newDir;
        if (icoarsew[0]==-1)
            newDir = localP?0:nDir/2;
        else {
            newDir = localP?0:std::max_element(res.begin(),res.end()) - res.begin();
        }


        onefspacetime -= getTime();
        if (!isScattering) if (domain->solInfo().loadcases.size()>0)
                domain->solInfo().loadcases = loadcases_backup;
        while (isDone.sum()<nDir) {
            if (verboseFlag)
                filePrint(stderr,
                          "nr=%d newDir=%d isDone.sum()=%d \n",nr,newDir,isDone.sum());

            // Build spaces at new frequency
            if (newDir!=-1) {
                if (isScattering)
                    probDesc->setIWaveDir(newDir);
                else {
                    if (newDir>0) domain->solInfo().loadcases.pop_front();
                }
                localFWSSpace[nDir+newDir] =
                    new FWindowSData<Scalar,VecType,ProblemDescriptor>
                        (maxLocalV,imode,probDesc);
                VecType **u = localFWSSpace[nDir+newDir]->u;
                VecType **v = localFWSSpace[nDir+newDir]->v;

                // Build one space
                int localDim = 0;
                for(int i=0;i<maxLocalV;i++) {
                    if (i == 0) {
                        rhsctime1 -= getTime();
                        probDesc->getRHS(*f);
                        rhsctime1 += getTime();
                        if (dgpFlag==1) temp_u[0]->zero();
                        if (dgpFlag==0) *fkgp = *f;
                    }
                    else if (dgpFlag==1) {
                        // DGP
                        f->zero();
                        rhsctime1 -= getTime();
                        probDesc->getFreqSweepRHS(f, temp_u, i);
                        rhsctime1 += getTime();
                    } else if (dgpFlag==2) {
                        // WCAWE
                        int ii = i+1;
                        std::vector<Scalar> pb(ii-1);
                        std::vector<Scalar> invU(i*i,0.0);
                        for (int k=0; k<i; k++) invU[k+k*i] = 1.0;
                        std::vector<Scalar> matrixU(i*i);
                        for(int k=0;k<i;k++) for(int l=0;l<i;l++)
                                matrixU[k+i*l] = UU[k+l*maxLocalV];
                        std::vector<int> ipiv(i);
                        int info = 0;
                        Tgesv(i, i, &matrixU[0], i, &ipiv[0], &invU[0], i, info);

                        for(int j=1;j<ii;j++) {
                            int iimj = ii-j;
                            std::vector<Scalar> e(iimj,0);
                            e[iimj-1] = 1.0;
                            for(int k=j;k>=1;k--) {
                                std::vector<Scalar> ee(iimj,0);
                                char transa = 'N';
                                Tgemv(transa, iimj, iimj, 1.0 , &invU[k-1+(k-1)*i], i,
                                      &e[0], 1, 0.0, &ee[0], 1);
                                for(int kk=0;kk<iimj;kk++) e[kk] = ee[kk];
                            }
                            pb[j-1] = e[0];
                        }

                        // Identical as above except j and k start from 2
                        std::vector<Scalar> pU(ii*ii,0);
                        for(int j=2;j<ii;j++) {
                            int iimj = ii-j;
                            std::vector<Scalar> e(iimj,0);
                            e[iimj-1] = 1.0;
                            for(int k=j;k>=2;k--) {
                                std::vector<Scalar> ee(iimj,0);
                                char transa = 'N';
                                Tgemv(transa, iimj, iimj, 1.0 , &invU[k-1+(k-1)*i], i,
                                      &e[0], 1, 0.0, &ee[0], 1);
                                for(int kk=0;kk<iimj;kk++) e[kk] = ee[kk];
                            }
                            for(int k=0;k<iimj;k++)
                                pU[k+(j-1)*ii] = e[k];
                        }
                        rhsctime1 -= getTime();
                        probDesc->getWCAWEFreqSweepRHS(f,
                                                       temp_u, &pU[0], &pb[0], maxLocalV, i);
                        rhsctime1 += getTime();
                    } else {
                        // KGP
                        allOps->M->mult(*u[i-1], *f);
                        if (dgpFlag==0) allOps->M->mult(*v[i-1], *fkgp);
                    }
                    if (verboseFlag)
                        filePrint(stderr," ... Solving RHS   #%3d               ...\n",i);
                    rhstime -= getTime();
                    allOps->sysSolver->solve(*f, *a);
                    solvecount++;
                    if (dgpFlag==0) {
                        allOps->sysSolver->solve(*fkgp, *akgp);
                        solvecount++;
                    }
                    rhstime += getTime();
                    if (dgpFlag) {
                        if (isFeti(domain->solInfo().solvercntl->type)) {
                            if (verboseFlag)
                                filePrint(stderr,"Forcing continuity %d.\n",i);
                            forceContinuity(*a);
                        }
                    }

                    if (dgpFlag==1)  *temp_u[i+1] = *a;
                    else if (dgpFlag==2) {
                        // WCAWE orthogonalize
                        Scalar nrm0 = *a * *a;
                        for(int j=0;j<i;j++) {
                            Scalar dotp = *a *  *temp_u[j];
                            UU[j+i*maxLocalV] = dotp;
                            (*a).linAdd(-dotp,*temp_u[j]);
                        }
                        UU[i*maxLocalV+i] = sqrt(ScalarTypes::Real(*a * *a));
                        *a *= 1.0/UU[i*maxLocalV+i];
                        *temp_u[i] = *a;
                        double normRedW =
                            ScalarTypes::Real(UU[i*maxLocalV+i])/sqrt(ScalarTypes::Real(nrm0));
                        if (normRedW<1e-17/tol) { break; }
                    }

                    localFWSSpace[nDir+newDir]->n = i;

                    *b = *a;
                    // Global space orthogonalize with previous coupled scaling
                    if(domain->solInfo().isCoupled) scaleDisp(*a,1.0/alpha);
                    Scalar nrm0 = *a * *a;
                    int ngs = 2;
                    for(int m=0;m<ngs;m++) {
                        for(int jDir=0;jDir<nF*nDir;jDir++)
                            if (localFWSSpace[jDir]!=0)
                                if (!localP || jDir==newDir || jDir==newDir+nDir)
                                    for(int j=0;j<localFWSSpace[jDir]->n;j++)
                                    {
                                        Scalar dotp = *a *  *(localFWSSpace[jDir]->ups[j]);
                                        (*a).linAdd(-dotp,*(localFWSSpace[jDir]->ups[j]));
                                    }
                    }
                    Scalar nrm = *a * *a;
                    double normRedG =
                        sqrt(ScalarTypes::Real(nrm))/sqrt(ScalarTypes::Real(nrm0));
                    *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
                    *(localFWSSpace[nDir+newDir]->ups[i]) = *a;
                    if(domain->solInfo().isCoupled) scaleDisp(*a,alpha);
                    *u[i] = *a;
                    // Local frequency space orthogonalize with current coupled scaling
                    if (dgpFlag==0) *a = *akgp;
                    else *a = *b;
                    nrm0 = *a * *a;
                    for(int m=0;m<ngs;m++) {
                        for(int jDir=nDir;jDir<nF*nDir;jDir++)
                            if (localFWSSpace[jDir]!=0)
                                if (!localP || jDir==newDir+nDir)
                                    for(int j=0;j<localFWSSpace[jDir]->n;j++)
                                    {
                                        Scalar dotp = *a *  *(localFWSSpace[jDir]->v[j]);
                                        (*a).linAdd(-dotp,*(localFWSSpace[jDir]->v[j]));
                                    }
                    }
                    nrm = *a * *a;
                    *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
                    *v[i] = *a;
                    double normRedL =
                        sqrt(ScalarTypes::Real(nrm))/sqrt(ScalarTypes::Real(nrm0));
                    if (normRedG<1e-20/tol) break;
                    localDim++;
                }

                localFWSSpace[nDir+newDir]->n = localDim;
                // Update matrices
                for(int i=0;i<localFWSSpace[nDir+newDir]->n;i++) {
                    allOps->K->mult(*(u[i]), *(localFWSSpace[nDir+newDir]->Ku[i]));
                    allOps->M->mult(*(u[i]), *(localFWSSpace[nDir+newDir]->Mu[i]));
                    if (allOps->C_deriv) if (allOps->C_deriv[0])
                            allOps->C_deriv[0]->mult(*(u[i]), *(localFWSSpace[nDir+newDir]->Cu[i]));
                }

                for(int i=0;i<localFWSSpace[nDir+newDir]->n;i++) {
                    // Compute new reduced matrix entries
                    for(int jDir=0;jDir<nF*nDir;jDir++) {
                        if (localFWSSpace[jDir]!=0)
                            if (!localP || jDir==newDir || jDir==newDir+nDir)
                                for(int j=0;j<localFWSSpace[jDir]->n;j++) {
                                    std::pair<int,int> rowcol
                                        (localFWSSpace[jDir]->imode+j,localFWSSpace[nDir+newDir]->imode+i);
                                    Kr[rowcol] = *(localFWSSpace[nDir+newDir]->Ku[i]) *
                                                 *(localFWSSpace[jDir]->u[j]);
                                    Mr[rowcol] = *(localFWSSpace[nDir+newDir]->Mu[i]) *
                                                 *(localFWSSpace[jDir]->u[j]);
                                    Cr[rowcol] = *(localFWSSpace[nDir+newDir]->Cu[i]) *
                                                 *(localFWSSpace[jDir]->u[j]);
                                    if (newDir+nDir!=jDir) {
                                        std::pair<int,int> rowcolt
                                            (localFWSSpace[nDir+newDir]->imode+i,
                                             localFWSSpace[jDir]->imode+j);
                                        Kr[rowcolt] = *(localFWSSpace[jDir]->Ku[j]) *
                                                      *(localFWSSpace[nDir+newDir]->u[i]);
                                        Mr[rowcolt] = *(localFWSSpace[jDir]->Mu[j]) *
                                                      *(localFWSSpace[nDir+newDir]->u[i]);
                                        Cr[rowcolt] = *(localFWSSpace[jDir]->Cu[j]) *
                                                      *(localFWSSpace[nDir+newDir]->u[i]);
                                    }
                                }
                    }
                }

                imode += localFWSSpace[nDir+newDir]->n;
                if (!localP) {
                    localFWSSpace[nDir+newDir]->irow = nr;
                    nr += localFWSSpace[nDir+newDir]->n;
                } else {
                    if (localFWSSpace[newDir]!=0)
                        localFWSSpace[nDir+newDir]->irow = localFWSSpace[newDir]->n;
                    else
                        localFWSSpace[nDir+newDir]->irow = 0;
                }
                isDone[newDir] = 1;
            }

            if (!localP) {
                // Sample residuals at midpoints of "done"
                std::valarray<int> iPick(0,nDir);
                int nPick = 0;
                int iLeft,iRight;
                for(iLeft=0;iLeft<nDir; iLeft++) if (isDone[iLeft]==1) break;
                if (iLeft/3!=iLeft) iPick[nPick++] = iLeft/3;
                bool last = false;
                while (!last) {
                    for(iRight=iLeft+1;iRight<nDir; iRight++) if (isDone[iRight]==1) break;
                    if (iRight==nDir) {
                        if ((iLeft+2*(nDir))/3!=iLeft && (iLeft+2*(nDir))/3<nDir)
                            iPick[nPick++] = (iLeft+2*(nDir))/3;
                        last = true;
                    } else {
                        if ((iLeft+iRight)/2!=iLeft) iPick[nPick++] = (iLeft+iRight)/2;
                    }
                    iLeft = iRight;
                }
                if (verboseFlag) {
                    filePrint(stderr,"iDir= ");
                    for(int ii=0;ii<nPick;ii++)
                        filePrint(stderr," %d",iPick[ii]);
                    filePrint(stderr,"\n");
                }
                // Evaluate residual at picked rhs's
                for(int iDir=0;iDir<nDir;iDir++) res[iDir] = 0;
                for(int ii=0;ii<nPick;ii++) {
                    int iDir = iPick[ii];
                    if (isScattering)
                        probDesc->setIWaveDir(iDir);
                    else
                        filePrint(stderr,
                                  "Error: multiple loadcases should be called with local only.\n");
                    rhsctime1 -= getTime();
                    probDesc->getRHS(*f);
                    rhsctime1 += getTime();
                    // Compute reduced vector
                    froneftime -= getTime();
                    for(int ir=0;ir<nr;ir++) fr[ir] = 0;
                    for(int jDir=0;jDir<nF*nDir;jDir++)
                        if (localFWSSpace[jDir]!=0)
                            for(int j=0;j<localFWSSpace[jDir]->n;j++)
                                fr[localFWSSpace[jDir]->irow+j] =
                                    *f * *(localFWSSpace[jDir]->u[j]);
                    froneftime += getTime();
                    // Compute reduced matrix
                    Aroneftime -= getTime();
                    for(int ir=0;ir<nr*nr;ir++) Ar[ir] = 0;
                    for(int jDir=0;jDir<nF*nDir;jDir++) if (localFWSSpace[jDir]!=0)
                            for(int j=0;j<localFWSSpace[jDir]->n;j++)  {
                                for(int kDir=0;kDir<nF*nDir;kDir++) if (localFWSSpace[kDir]!=0)
                                        for(int k=0;k<localFWSSpace[kDir]->n;k++)  {
                                            std::pair<int,int> rowcol
                                                (localFWSSpace[jDir]->imode+j,localFWSSpace[kDir]->imode+k);
                                            Ar[localFWSSpace[jDir]->irow+j+
                                               nr*(localFWSSpace[kDir]->irow+k)] =
                                                Kr[rowcol] - w*w*Mr[rowcol];
                                            ScalarTypes::addComplex(Ar[localFWSSpace[jDir]->irow+j+
                                                                       nr*(localFWSSpace[kDir]->irow+k)], w *Cr[rowcol]);
                                        }
                            }
                    Aroneftime += getTime();

                    // Solve the reduced linear system
                    std::vector<int> ipiv(nr);
                    int lwork = 3 * (nr);
                    std::vector<Scalar> work(lwork);
                    int info = 0;
                    Tgesv(nr, 1, &Ar[0], nr, &ipiv[0], &fr[0], nr, info);
                    // Evaluate the residual
                    resoneftime -= getTime();
                    (*a).zero();
                    for(int jDir=0;jDir<nF*nDir;jDir++)
                        if (localFWSSpace[jDir]!=0)
                            for(int j=0;j<localFWSSpace[jDir]->n;j++) {
                                int ir = localFWSSpace[jDir]->irow+j;
                                (*a).linAdd(fr[ir],*(localFWSSpace[jDir]->Ku[j]),
                                            -w*w*fr[ir],*(localFWSSpace[jDir]->Mu[j]));
                                Scalar dc = 0.0;
                                ScalarTypes::addComplex(dc,w*fr[ir]);
                                if (allOps->C_deriv) if (allOps->C_deriv[0])
                                        (*a).linAdd(dc,*(localFWSSpace[jDir]->Cu[j]));
                            }
                    (*a).linAdd(-1.0,*f);
                    forceAssemble(*a);
                    forceAssemble(*f);
                    Scalar nrma = *a * *a;
                    Scalar nrmf = *f * *f;
                    resoneftime += getTime();
                    res[iDir] = sqrt(ScalarTypes::Real(nrma)/ScalarTypes::Real(nrmf));
                    if (verboseFlag)
                        filePrint(stderr,"iDir=%d res=%e\n",iDir,res[iDir]);
                    if (res[iDir]<tol*tol1f) isDone[iDir] = 1;
                }
                double maxres = -1.0;
                for(int iDir=0;iDir<nDir;iDir++) if (res[iDir]>maxres) {
                        maxres = res[iDir];
                        newDir = iDir;
                        if (res[iDir]<tol*tol1f) newDir = -1;
                    }
                if (verboseFlag)
                    filePrint(stderr,"maxres=%e newDir=%d\n",maxres,newDir);
            } else {
                newDir++;
            }
        }
        onefspacetime += getTime();


        fsweeptime -= getTime();
        // Compute and output ROM solution and residual
        //  for interval icoarsew[0], icoarsew[1]
        //  and find new frequency above icoarsew[1]
        int iwr = (icoarsew[0]>=0)?icoarsew[0]-1:0;
        int deltaiwr;
        int liwr=-1, uiwr=-1;
        bool tolMet;
        bool seekMode = false;
        while (1) {
            // Find the increment to the next frequency
            if (verboseFlag)
                filePrint(stderr, "icoarsew[0]=%d icoarsew[1]=%d liwr=%d uiwr=%d deltaiwr=%d\n", icoarsew[0], icoarsew[1],
                          liwr, uiwr, deltaiwr);
            //   If done numS, we are done
            if (iwr == numS && !seekMode) {
                doneF = true;
                break;
            }
            //   If the first first time
            if (icoarsew[0] < 0 && liwr == -1) {
                seekMode = true;
                deltaiwr = numS;
                liwr = 0;
                if (verboseFlag)
                    filePrint(stderr,
                              "init liwr=%d uiwr=%d deltaiwr=%d\n", liwr, uiwr, deltaiwr);
            }
                // If generating solutions
            else if (iwr < icoarsew[1]) {
                deltaiwr = 1;
                if (verboseFlag)
                    filePrint(stderr,
                              "solution liwr=%d uiwr=%d deltaiwr=%d\n", liwr, uiwr, deltaiwr);
            }
                //   If just did the upper frequency, switch to seek mode
            else if (iwr == icoarsew[1]) {
                seekMode = true;
                deltaiwr = icoarsew[1] - icoarsew[0];
                if (iwr + deltaiwr > numS) deltaiwr = numS - iwr;
                liwr = icoarsew[1];
                if (verboseFlag)
                    filePrint(stderr,
                              "start seek liwr=%d uiwr=%d deltaiwr=%d\n", liwr, uiwr, deltaiwr);
            }
                // Seek mode
            else {
                if (tolMet) {
                    liwr = iwr;
                    if (uiwr != -1) deltaiwr = (uiwr - liwr) / 2;
                    if (iwr + deltaiwr > numS) deltaiwr = numS - iwr;
                    if (verboseFlag)
                        filePrint(stderr,
                                  "tolmet liwr=%d uiwr=%d deltaiwr=%d\n", liwr, uiwr, deltaiwr);
                } else {
                    uiwr = iwr;
                    deltaiwr = -(uiwr - liwr) / 2;
                    if (verboseFlag)
                        filePrint(stderr,
                                  "tolnotmet liwr=%d uiwr=%d deltaiwr=%d\n", liwr, uiwr, deltaiwr);
                }
                if (deltaiwr == 0) {
                    iwr = liwr;
                    icoarsew[0] = icoarsew[1];
                    icoarsew[1] = iwr;
                    if (verboseFlag)
                        filePrint(stderr,
                                  "stop liwr=%d uiwr=%d deltaiwr=%d\n", liwr, uiwr, deltaiwr);
                    break;
                }
            }

            iwr = iwr + deltaiwr;

            double wr = w1 + iwr * dw;
            if (verboseFlag)
                filePrint(stderr, "Computing residual at iwr= %d wr/2/pi=%e\n", iwr, wr / 2.0 / M_PI);
// RT uncommenting following line  12/8
            geoSource->setOmega(wr);
            std::vector<int> ipiv(nr);
            if (!localP) {
                // Compute reduced matrix
                Arsweeptime -= getTime();
                for (int ir = 0; ir < nr * nr; ir++) Ar[ir] = 0;
                for (int jDir = 0; jDir < nF * nDir; jDir++)
                    if (localFWSSpace[jDir] != 0)
                        for (int j = 0; j < localFWSSpace[jDir]->n; j++) {
                            for (int kDir = 0; kDir < nF * nDir; kDir++)
                                if (localFWSSpace[kDir] != 0)
                                    for (int k = 0; k < localFWSSpace[kDir]->n; k++) {
                                        std::pair<int, int> rowcol
                                            (localFWSSpace[jDir]->imode + j, localFWSSpace[kDir]->imode + k);
                                        //                 auto Kr_iter = Kr.find(rowcol);
                                        if (Kr.find(rowcol) == Kr.end())
                                            fprintf(stderr,
                                                    "rowcol %d %d not found\n",
                                                    localFWSSpace[jDir]->imode + j, localFWSSpace[kDir]->imode + k);
                                        Ar[localFWSSpace[jDir]->irow + j +
                                           nr * (localFWSSpace[kDir]->irow + k)] =
                                            Kr[rowcol] - wr * wr * Mr[rowcol];
                                        ScalarTypes::addComplex(Ar[localFWSSpace[jDir]->irow + j +
                                                                   nr * (localFWSSpace[kDir]->irow + k)], wr * Cr[rowcol]);
                                    }
                        }
                int info = 0;
                Tgetrf(nr, nr, &Ar[0], nr, &ipiv[0], info);
                Arsweeptime += getTime();
            }
            if (!isScattering)
                if (domain->solInfo().loadcases.size() > 0)
                    domain->solInfo().loadcases = loadcases_backup;
            for (int iDir = 0; iDir < nDir; iDir++) {
                // Get RHS
                if (isScattering)
                    probDesc->setIWaveDir(iDir);
                else {
                    if (iDir > 0) domain->solInfo().loadcases.pop_front();
                }
                rhsctime2 -= getTime();
                probDesc->getRHS(*f, wr, wr - (w1 + icoarsew[1] * dw));
                rhsctime2 += getTime();
                // Compute reduced vector
                if (localP) {
                    nr = 0;
                    if (localFWSSpace[iDir] != 0) nr += localFWSSpace[iDir]->n;
                    if (localFWSSpace[iDir + nDir] != 0) nr += localFWSSpace[iDir + nDir]->n;
                }
                frsweeptime -= getTime();
                for (int ir = 0; ir < nr; ir++) fr[ir] = 0;
                for (int jDir = 0; jDir < nF * nDir; jDir++)
                    if (localFWSSpace[jDir] != 0)
                        if (!localP || jDir == iDir || jDir == iDir + nDir)
                            for (int j = 0; j < localFWSSpace[jDir]->n; j++) {
                                fr[localFWSSpace[jDir]->irow + j] =
                                    *f * *(localFWSSpace[jDir]->u[j]);
                            }
                frsweeptime += getTime();
                if (localP) {
                    Arsweeptime -= getTime();
                    for (int ir = 0; ir < nr * nr; ir++) Ar[ir] = 0;
                    for (int jDir = 0; jDir < nF * nDir; jDir++)
                        if (localFWSSpace[jDir] != 0)
                            if (jDir == iDir || jDir == iDir + nDir)
                                for (int j = 0; j < localFWSSpace[jDir]->n; j++) {
                                    for (int kDir = 0; kDir < nF * nDir; kDir++)
                                        if (localFWSSpace[kDir] != 0)
                                            if (kDir == iDir || kDir == iDir + nDir)
                                                for (int k = 0; k < localFWSSpace[kDir]->n; k++) {
                                                    std::pair<int, int> rowcol
                                                        (localFWSSpace[jDir]->imode + j, localFWSSpace[kDir]->imode + k);
                                                    //                 auto Kr_iter = Kr.find(rowcol);
                                                    if (Kr.find(rowcol) == Kr.end())
                                                        fprintf(stderr,
                                                                "rowcol %d %d not found\n",
                                                                localFWSSpace[jDir]->imode + j,
                                                                localFWSSpace[kDir]->imode + k);
                                                    Ar[localFWSSpace[jDir]->irow + j +
                                                       nr * (localFWSSpace[kDir]->irow + k)] =
                                                        Kr[rowcol] - wr * wr * Mr[rowcol];
                                                    ScalarTypes::addComplex(Ar[localFWSSpace[jDir]->irow + j +
                                                                               nr * (localFWSSpace[kDir]->irow + k)],
                                                                            wr * Cr[rowcol]);
                                                }
                                }
                    Arsweeptime += getTime();
                    // Factor and solve the reduced linear system
                    std::vector<int> ipiv(nr);
                    int info = 0;
                    Tgesv(nr, 1, &Ar[0], nr, &ipiv[0], &fr[0], nr, info);
                } else {
                    // Solve the reduced linear system
                    int info = 0;
                    char trans = 'N';
                    Tgetrs(trans, nr, 1, &Ar[0], nr, &ipiv[0], &fr[0], nr, info);
                }

                // Evaluate the residual
                ressweeptime -= getTime();
                (*a).zero();
                for (int jDir = 0; jDir < nF * nDir; jDir++)
                    if (localFWSSpace[jDir] != 0)
                        if (!localP || jDir == iDir || jDir == iDir + nDir)
                            for (int j = 0; j < localFWSSpace[jDir]->n; j++) {
                                int ir = localFWSSpace[jDir]->irow + j;
                                (*a).linAdd(fr[ir], *(localFWSSpace[jDir]->Ku[j]),
                                            -wr * wr * fr[ir], *(localFWSSpace[jDir]->Mu[j]));
                                Scalar dc = 0.0;
                                ScalarTypes::addComplex(dc, wr * fr[ir]);
                                if (allOps->C_deriv)
                                    if (allOps->C_deriv[0])
                                        (*a).linAdd(dc, *(localFWSSpace[jDir]->Cu[j]));
                            }
/* 
       allOps->K->mult(*sol, *a);
       allOps->M->mult(*sol, *b);
       if (allOps->C_deriv) {
          if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
       } else (*c).zero(); 
       (*a).linAdd(-wr*wr,*b);
       for(int k=0;k<sol->size();k++) 
           ScalarTypes::addComplex((*a)[k], wr * (*c)[k]);
*/
                (*a).linAdd(-1.0, *f);
                forceAssemble(*a);
                forceAssemble(*f);
                Scalar nrma = *a * *a;
                Scalar nrmf = *f * *f;
                ressweeptime += getTime();
                res[iDir] = sqrt(ScalarTypes::Real(nrma) / ScalarTypes::Real(nrmf));
                if (iwr < icoarsew[1] || (iwr == numS && icoarsew[1] == numS)) {
                    solutiontime -= getTime();
                    if (verboseFlag)
                        filePrint(stderr, "fresidual %d: %e %d %d\n",
                                  iDir, res[iDir], iwr, icoarsew[1]);
                    // Compute the approximant vector
                    (*sol).zero();
                    for (int jDir = 0; jDir < nF * nDir; jDir++)
                        if (localFWSSpace[jDir] != 0)
                            if (!localP || jDir == iDir || jDir == iDir + nDir)
                                for (int j = 0; j < localFWSSpace[jDir]->n; j++)
                                    (*sol).linAdd(fr[localFWSSpace[jDir]->irow + j],
                                                  *(localFWSSpace[jDir]->u[j]));
                    if (domain->solInfo().isCoupled) scaleDisp(*sol);
                    solutiontime += getTime();
                    // Output the solution making sure the right frequency appears
                    domain->frequencies->push_front(wr);
                    if (domain->solInfo().isAcousticHelm()) {
                        SPropContainer &sProps = geoSource->getStructProps();
                        for (int iProp = 0; iProp < geoSource->getNumProps(); iProp++) {
                            if (sProps[iProp].kappaHelm != 0.0 ||
                                sProps[iProp].kappaHelmImag != 0.0) {
                                complex<double> k1 = wr / sProps[iProp].soundSpeed;
                                sProps[iProp].kappaHelm = real(k1);
                                sProps[iProp].kappaHelmImag = imag(k1);
                            }
                        }
                    }
                    outputtime -= getTime();
                    postProcessor->staticOutput(*sol, *rhs, (iwr == numS) && (iDir == nDir - 1));
                    outputtime += getTime();
                    domain->frequencies->pop_front();
                } else if (verboseFlag)
                    filePrint(stderr, "residual %d: %e %d %d\n",
                              iDir, res[iDir], iwr, icoarsew[1]);
            }
            double mres = *std::max_element(res.begin(), res.end());
            if (iwr < icoarsew[1]) {
                if (verboseFlag)
                    filePrint(stderr, "max fresidual: %e at %e \n", mres, wr / 2.0 / M_PI);
            } else {
                if (verboseFlag)
                    filePrint(stderr, "max residual: %e at %e \n", mres, wr / 2.0 / M_PI);
            }
            tolMet = mres<tol*ctolf;
        }
        fsweeptime += getTime();
    }

    time += getTime();
    if (verboseFlag) {
        filePrint(stderr,"Total time: %e\n",time/1e3);
        filePrint(stderr,"----\n");
        filePrint(stderr,"Transition time: %e\n",transitiontime/1e3);
        filePrint(stderr,"----\n");
        filePrint(stderr,"One frequency space construction time: %e\n",
                  onefspacetime/1e3);
        filePrint(stderr,"solver time: %e in %d solves\n",
                  solvertime/1e3,rebuildscount);
        filePrint(stderr,"rhs solve time: %e count: %d\n",rhstime/1e3,solvecount);
        filePrint(stderr,"fr time in one f: %e\n",froneftime/1e3);
        filePrint(stderr,"Ar time in one f: %e\n",Aroneftime/1e3);
        filePrint(stderr,"Residual computation time in one f: %e\n",resoneftime/1e3);
        filePrint(stderr,"RHS construction time in one f: %e\n",rhsctime1/1e3);
        filePrint(stderr,"----\n");
        filePrint(stderr,"Frequency sweep time: %e\n",fsweeptime/1e3);
        filePrint(stderr,"fr time in sweep: %e\n",frsweeptime/1e3);
        filePrint(stderr,"Ar time in sweep: %e\n",Arsweeptime/1e3);
        filePrint(stderr,"residual computation time in sweep: %e\n",ressweeptime/1e3);
        filePrint(stderr,"rhs construction time in sweep: %e\n",rhsctime2/1e3);
        filePrint(stderr,"solution construction time in sweep: %e\n",solutiontime/1e3);
        filePrint(stderr,"output time in sweep: %e\n",outputtime/1e3);
    }

    if (dgpFlag) {
        delete[] UU;
        for(int i=0;i<maxLocalV+1;i++) delete temp_u[i];
        delete[] temp_u;
    }
    delete f;
    delete fkgp;
    delete c;
    delete b;
    delete a;
    delete akgp;
    delete sol;

    for(int i=0;i<nF*nDir;i++) {
        if (localFWSSpace[i]!=0) delete localFWSSpace[i];
    }
    delete[] localFWSSpace;

}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::adaptSweep()
{
// Global adaptive sweep with multiple right hand sides done locally
 double time = 0.0;
 double romctime = 0.0;
 double transitiontime= 0.0;
 double onefspacetime = 0.0;
 double solvertime = 0.0;
 double rhstime = 0.0;
 double rhsctime1 = 0.0;
 double fsweeptime = 0.0;
 double frsweeptime = 0.0;
 double Arsweeptime = 0.0;
 double ressweeptime = 0.0;
 double solutiontime = 0.0;
 double outputtime = 0.0;
 double rhsctime2 = 0.0;
 int solvecount = 0;
 int rebuildscount = 0;

 domain->isCoarseGridSolve = false;

 time -= getTime();

 double w1 = domain->solInfo().getSweepParams()->adaptSweep.w1;
 double w2 = domain->solInfo().getSweepParams()->adaptSweep.w2;
 int numS = domain->solInfo().getSweepParams()->adaptSweep.numS;
 double dw = (w2-w1)/numS;

 int nDir = 1;
 bool isScattering = false;
 if (domain->numWaveDirections>=1 && domain->solInfo().loadcases.size()==0) {
   int nDir = domain->numWaveDirections;
   isScattering = true;
   filePrint(stderr,"Solving %d scattering cases\n",nDir);
 }
 std::list<int> loadcases_backup;
 if (domain->solInfo().loadcases.size()>=1) {
   nDir = domain->solInfo().loadcases.size();
   loadcases_backup = domain->solInfo().loadcases;
   filePrint(stderr,"Solving %d loadcases\n",nDir);
 }
 if (nDir==0) nDir = 1;
 int maxLocalV = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
 int maxRHS = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
 int minRHS = domain->solInfo().getSweepParams()->adaptSweep.minRHS;
 int deltaRHS = domain->solInfo().getSweepParams()->adaptSweep.deltaRHS;
 int dgpFlag = domain->solInfo().getSweepParams()->adaptSweep.dgp_flag;
 double tol = domain->solInfo().getSweepParams()->adaptSweep.atol;
 int maxP = domain->solInfo().getSweepParams()->adaptSweep.maxP;

 filePrint(stderr,"Adaptive sweep with multiple RHS - local variant\n");
 if (dgpFlag==1)
   filePrint(stderr,"DGP with %d,%d,%d local derivatives\n",
             minRHS,maxRHS,deltaRHS);
 else if (dgpFlag==2)
   filePrint(stderr,"DGP with %d,%d,%d local WCAWE vectors\n",
             minRHS,maxRHS,deltaRHS);
 else 
   filePrint(stderr,"KGP with %d,%d,%d local vectors\n",
             minRHS,maxRHS,deltaRHS);
 filePrint(stderr,"Parameter values: tol=%e, maxP=%d\n",
           tol,maxP);

 // Temp storage for DGP and WCAWE
 VecType **temp_u=0;
 Scalar *UU = 0;
 if (dgpFlag) {
   UU = new Scalar[maxLocalV*maxLocalV];
   temp_u = new VecType * [maxLocalV+1];
   for(int i = 0; i < maxLocalV+1; ++i)
     temp_u[i] = new VecType(probDesc->solVecInfo());
   for(int i = 0; i < maxLocalV*maxLocalV; ++i) UU[i] = 0.0;
 }

 // Temp vectors
 VecType *f = new VecType(probDesc->solVecInfo());
 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());
 VecType *sol = new VecType(probDesc->solVecInfo());

 // Subspaces, all maxP frequencies at a time
 FWindowSData<Scalar,VecType,ProblemDescriptor> **localFWSSpace = 
   new FWindowSData<Scalar,VecType,ProblemDescriptor> *[maxP*nDir];
 for(int i=0;i<maxP*nDir;i++) {
   localFWSSpace[i] = 0;
 }


// vector of residuals, reduced matrix and vector, reduced matrices
 const int maxRows = maxP*maxLocalV;
 std::vector<double> res(nDir);   
 std::vector<Scalar> fr(maxRows);
 std::vector<Scalar> Ar(maxRows*maxRows);
 std::map<std::pair<int,int>,Scalar> Kr; 
 std::map<std::pair<int,int>,Scalar> Mr; 
 std::map<std::pair<int,int>,Scalar> Cr; 
 
 std::vector<double> wc(maxP);
 std::vector<double> wcSorted(maxP);
 std::vector<double> newW(maxP);

 int numP = 0;
 int imode = 0;
 bool doneF = false;
 int numNewW = 2;
 newW[0] = w1;
 newW[1] = w2;
 int ncheck = 9;
 std::vector<double> wcheck;
 
 romctime -= getTime();
 while (!doneF) {
   for(int kP=0; kP<numNewW; kP++) {
     double w = newW[kP]; 
     wc[numP] = w;
     // Scaling factor for coupled
     double alpha = w/((w1+w2)/2.0);
     // Setup solver for new frequency
     geoSource->setOmega(w);
     solvertime -= getTime();
     if (w>w1) {
       rebuildSolver(w);
       rebuildscount++;
     }
     solvertime += getTime();

     // Update u's to current scaling and recompute Ku,Mu,Cu
     transitiontime -= getTime();
     for(int iP=0;iP<numP;iP++)  for(int iDir=0;iDir<nDir;iDir++)  {
       int iL = iP*nDir + iDir;    
       if (localFWSSpace[iL]!=0) {
         for(int i=0;i<localFWSSpace[iL]->n;i++) {
           *(localFWSSpace[iL]->u[i]) = *(localFWSSpace[iL]->v[i]);
//           *(localFWSSpace[iL]->ups[i]) = *(localFWSSpace[iL]->v[i]);
           if(domain->solInfo().isCoupled) 
             scaleDisp(*(localFWSSpace[iL]->u[i]),alpha);
           allOps->K->mult(*(localFWSSpace[iL]->u[i]),
                           *(localFWSSpace[iL]->Ku[i]));
           allOps->M->mult(*(localFWSSpace[iL]->u[i]),
                           *(localFWSSpace[iL]->Mu[i]));
           if (allOps->C_deriv) if (allOps->C_deriv[0])
             allOps->C_deriv[0]->mult(*(localFWSSpace[iL]->u[i]),
                                      *(localFWSSpace[iL]->Cu[i]));
         }
       }
     }
     // Recompute reduced matrices
     for(int iP=0;iP<numP;iP++)  for(int iDir=0;iDir<nDir;iDir++)  {
       int iL = iP*nDir + iDir;    
       if (localFWSSpace[iL]!=0) {
         for(int i=0;i<localFWSSpace[iL]->n;i++) {
           for(int jP=0;jP<numP;jP++) {
             int jL = jP*nDir + iDir;    
             if (localFWSSpace[jL]!=0) 
               for(int j=0;j<localFWSSpace[jL]->n;j++) {
                  std::pair<int,int> rowcol
                   (localFWSSpace[jL]->imode+j,localFWSSpace[iL]->imode+i); 
                  Kr[rowcol] = *(localFWSSpace[iL]->Ku[i]) *
                               *(localFWSSpace[jL]->u[j]);
                  Mr[rowcol] = *(localFWSSpace[iL]->Mu[i]) *
                               *(localFWSSpace[jL]->u[j]);
                  Cr[rowcol] = *(localFWSSpace[iL]->Cu[i]) *
                               *(localFWSSpace[jL]->u[j]);
                  if (iL!=jL) {
                    std::pair<int,int> rowcolt
                      (localFWSSpace[iL]->imode+i,
                       localFWSSpace[jL]->imode+j);
                    Kr[rowcolt] = *(localFWSSpace[jL]->Ku[j]) *
                                  *(localFWSSpace[iL]->u[i]);
                    Mr[rowcolt] = *(localFWSSpace[jL]->Mu[j]) *
                                  *(localFWSSpace[iL]->u[i]);
                    Cr[rowcolt] = *(localFWSSpace[jL]->Cu[j]) *
                                  *(localFWSSpace[iL]->u[i]);
                  }
               } 
           }
         }
       }
     }
     transitiontime += getTime();

     if (!isScattering) if (domain->solInfo().loadcases.size()>0)
       domain->solInfo().loadcases = loadcases_backup;
  
     onefspacetime -= getTime();
     for(int iDir=0;iDir<nDir;iDir++) {
       if (verboseFlag)
         filePrint(stderr, "numP=%d iDir=%d\n",numP,iDir);
  
     // Build spaces at new frequency
       if (isScattering)
         probDesc->setIWaveDir(iDir);
       else {
         if (iDir>0) domain->solInfo().loadcases.pop_front();
       }
       int newL = numP*nDir + iDir;
       localFWSSpace[newL] = 
         new FWindowSData<Scalar,VecType,ProblemDescriptor>
             (maxLocalV,imode,probDesc); 
       VecType **u = localFWSSpace[newL]->u;
       VecType **v = localFWSSpace[newL]->v;
    
       // Build one space 
       int localDim = 0;
       for(int i=0;i<maxLocalV;i++) {
         if (i == 0) { 
           rhsctime1 -= getTime();
           probDesc->getRHS(*f);
           rhsctime1 += getTime();
           if (dgpFlag==1) temp_u[0]->zero();
         } 
         else if (dgpFlag==1) {
           // DGP
           f->zero();
           rhsctime1 -= getTime();
           probDesc->getFreqSweepRHS(f, temp_u, i);
           rhsctime1 += getTime();
         } else if (dgpFlag==2) {
           // WCAWE
           int ii = i+1;
           std::vector<Scalar> pb(ii-1); 
           std::vector<Scalar> invU(i*i,0.0);
           for (int k=0; k<i; k++) invU[k+k*i] = 1.0;
           std::vector<Scalar> matrixU(i*i);
           for(int k=0;k<i;k++) for(int l=0;l<i;l++)
             matrixU[k+i*l] = UU[k+l*maxLocalV];
           std::vector<int> ipiv(i);
           int info = 0;
           Tgesv(i, i, &matrixU[0], i, &ipiv[0], &invU[0], i, info);

           for(int j=1;j<ii;j++) {
             int iimj = ii-j;
             std::vector<Scalar> e(iimj,0);
             e[iimj-1] = 1.0;
             for(int k=j;k>=1;k--) {
                std::vector<Scalar> ee(iimj,0);
                char transa = 'N';
                Tgemv(transa, iimj, iimj, 1.0 , &invU[k-1+(k-1)*i], i,
                      &e[0], 1, 0.0, &ee[0], 1);
                for(int kk=0;kk<iimj;kk++) e[kk] = ee[kk];
              }
              pb[j-1] = e[0];
           } 
        
           // Identical as above except j and k start from 2
           std::vector<Scalar> pU(ii*ii,0); 
           for(int j=2;j<ii;j++) {
             int iimj = ii-j;
             std::vector<Scalar> e(iimj,0);
             e[iimj-1] = 1.0;
             for(int k=j;k>=2;k--) {
                std::vector<Scalar> ee(iimj,0);
                char transa = 'N';
                Tgemv(transa, iimj, iimj, 1.0 , &invU[k-1+(k-1)*i], i,
                      &e[0], 1, 0.0, &ee[0], 1);
                for(int kk=0;kk<iimj;kk++) e[kk] = ee[kk];
              }
              for(int k=0;k<iimj;k++) 
               pU[k+(j-1)*ii] = e[k];  
           } 
           rhsctime1 -= getTime();
           probDesc->getWCAWEFreqSweepRHS(f,
                   temp_u, &pU[0], &pb[0], maxLocalV, i);
           rhsctime1 += getTime();
         } else {
           // KGP
           allOps->M->mult(*u[i-1], *f);
         }
         if (verboseFlag)
           filePrint(stderr," ... Solving RHS   #%3d               ...\n",i);
         rhstime -= getTime();
         allOps->sysSolver->solve(*f, *a);
         solvecount++;
         rhstime += getTime();
         if (dgpFlag) {
           if (isFeti(domain->solInfo().solvercntl->type)) {
             if (verboseFlag)
               filePrint(stderr,"Forcing continuity %d.\n",i);
             forceContinuity(*a);
           }
         }
    
         if (dgpFlag==1)  *temp_u[i+1] = *a;
         else if (dgpFlag==2) {
           // WCAWE orthogonalize
           Scalar nrm0 = *a * *a;
           for(int j=0;j<i;j++) {
              Scalar dotp = *a *  *temp_u[j];
              UU[j+i*maxLocalV] = dotp;
              (*a).linAdd(-dotp,*temp_u[j]);
           }
           UU[i*maxLocalV+i] = sqrt(ScalarTypes::Real(*a * *a));
           *a *= 1.0/UU[i*maxLocalV+i];
           *temp_u[i] = *a;
           double normRedW = 
             ScalarTypes::Real(UU[i*maxLocalV+i])/sqrt(ScalarTypes::Real(nrm0));
//           if (verboseFlag)
//             filePrint(stderr,"norm reduction - WCAWE: %e\n",normRedW);
           if (normRedW<1e-17/tol) { break; }
         }

         *b = *a;
         // Global space orthogonalize with mid-frequency coupled scaling
         if(domain->solInfo().isCoupled) scaleDisp(*a,1.0/alpha);
         Scalar nrm0 = *a * *a;
         int ngs = 2;
         for(int m=0;m<ngs;m++) { 
           for(int iP=0;iP<=numP;iP++) {
             int iL = iP*nDir + iDir;    
             if (localFWSSpace[iL]!=0) {
               for(int i=0;i<localFWSSpace[iL]->n;i++) {
                 Scalar dotp = *a *  *(localFWSSpace[iL]->v[i]);
                 (*a).linAdd(-dotp,*(localFWSSpace[iL]->v[i]));
               }
             }
           }
         }
         Scalar nrm = *a * *a;
         double normRedG = 
             sqrt(ScalarTypes::Real(nrm))/sqrt(ScalarTypes::Real(nrm0));
//         if (verboseFlag)
//           filePrint(stderr,"norm reduction - global space: %e\n",normRedG);
         *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
         *v[i] = *a;
         if(domain->solInfo().isCoupled) scaleDisp(*a,alpha);
         *u[i] = *a;
         if (normRedG<1e-20/tol) break;
         localDim++;
         localFWSSpace[newL]->n = localDim; 
       }

       // Compute new Ku,Mu,Cu
       for(int i=0;i<localFWSSpace[newL]->n;i++) {
         allOps->K->mult(*(u[i]), *(localFWSSpace[newL]->Ku[i]));
         allOps->M->mult(*(u[i]), *(localFWSSpace[newL]->Mu[i]));
         if (allOps->C_deriv) if (allOps->C_deriv[0])
           allOps->C_deriv[0]->mult(*(u[i]), *(localFWSSpace[newL]->Cu[i]));
       }
  
       // Compute new reduced matrix entries
       for(int i=0;i<localFWSSpace[newL]->n;i++) {
         for(int jP=0;jP<=numP;jP++) {
           int jL = jP*nDir + iDir;    
           if (localFWSSpace[jL]!=0) 
             for(int j=0;j<localFWSSpace[jL]->n;j++) {
               std::pair<int,int> rowcol
               (localFWSSpace[jL]->imode+j,localFWSSpace[newL]->imode+i); 
               Kr[rowcol] = *(localFWSSpace[newL]->Ku[i]) *
                            *(localFWSSpace[jL]->u[j]);
               Mr[rowcol] = *(localFWSSpace[newL]->Mu[i]) *
                            *(localFWSSpace[jL]->u[j]);
               Cr[rowcol] = *(localFWSSpace[newL]->Cu[i]) *
                            *(localFWSSpace[jL]->u[j]);
               if (newL!=jL) {
                 std::pair<int,int> rowcolt
                   (localFWSSpace[newL]->imode+i,
                    localFWSSpace[jL]->imode+j); 
                 Kr[rowcolt] = *(localFWSSpace[jL]->Ku[j]) *
                               *(localFWSSpace[newL]->u[i]);
                 Mr[rowcolt] = *(localFWSSpace[jL]->Mu[j]) *
                               *(localFWSSpace[newL]->u[i]);
                 Cr[rowcolt] = *(localFWSSpace[jL]->Cu[j]) *
                               *(localFWSSpace[newL]->u[i]);
               }
             } 
         }
       }

       imode += localFWSSpace[newL]->n;
filePrint(stderr,"imode=%d\n",imode);

       int prevL = (numP-1)*nDir + iDir;
       if (numP>0) 
         localFWSSpace[newL]->irow = localFWSSpace[prevL]->irow +
                                      localFWSSpace[prevL]->n;
       else localFWSSpace[newL]->irow = 0;
     }
     numP++; 
     onefspacetime += getTime();
     if (numP==maxP) break;
   }

   if (numP==maxP) break;
 
   // Sample intervals and find new coarse points
   numNewW = 0; 
   for(int nP=0;nP<numP;nP++) wcSorted[nP] = wc[nP];
   sort(wcSorted.begin(),wcSorted.begin()+numP);
   doneF = true;
   for(int nP=0;nP<numP-1;nP++) {
     double maxres = -1.0;
     int maxDir, maxIS;
     for(int iS=0;iS<ncheck;iS++) {
       double wS = wcSorted[nP] +
                   double(iS+1)/double(ncheck+1)*(wcSorted[nP+1]-wcSorted[nP]);
       // Check the residual at wS
       if (!isScattering) if (domain->solInfo().loadcases.size()>0)
          domain->solInfo().loadcases = loadcases_backup;
       for(int iDir=0;iDir<nDir;iDir++) {
         // Compute RHS
         if (isScattering)
           probDesc->setIWaveDir(iDir);
         else {
           if (iDir>0) domain->solInfo().loadcases.pop_front();
         }
         rhsctime2 -= getTime();
         probDesc->getRHS(*f, wS,wS-wc[numP-1]);
         rhsctime2 += getTime();

         // Compute reduced vector
         frsweeptime -= getTime();
         int nr = 0;
         for(int iP=0;iP<numP;iP++) {
           int iL = iP*nDir + iDir;
           if (localFWSSpace[iL]!=0) {
             nr += localFWSSpace[iL]->n;
             for(int j=0;j<localFWSSpace[iL]->n;j++) {
               fr[localFWSSpace[iL]->irow+j] = *f * *(localFWSSpace[iL]->u[j]);
             }
           }
         }
         frsweeptime += getTime();
    
         // Compute reduced matrix
         Arsweeptime -= getTime();
         for(int ir=0;ir<nr*nr;ir++) Ar[ir] = 0;
         for(int jP=0;jP<numP;jP++) {
           int jL = jP*nDir + iDir;
           if (localFWSSpace[jL]!=0) for(int j=0;j<localFWSSpace[jL]->n;j++)  {
             for(int kP=0;kP<numP;kP++) {
               int kL = kP*nDir + iDir;
               if (localFWSSpace[kL]!=0)
                 for(int k=0;k<localFWSSpace[kL]->n;k++) {
                 std::pair<int,int> rowcol
                   (localFWSSpace[jL]->imode+j,localFWSSpace[kL]->imode+k); 
                 if (Kr.find(rowcol)==Kr.end()) fprintf(stderr,
                             "rowcol %d %d not found\n",
                    localFWSSpace[jL]->imode+j,localFWSSpace[kL]->imode+k);
                 Ar[localFWSSpace[jL]->irow+j+
                 nr*(localFWSSpace[kL]->irow+k)] = 
                     Kr[rowcol] - wS*wS*Mr[rowcol];
                   ScalarTypes::addComplex(Ar[localFWSSpace[jL]->irow+j+
                     nr*(localFWSSpace[kL]->irow+k)], wS *Cr[rowcol]);
                   }
                 }
             }
         }
         Arsweeptime += getTime();
         // Factor and solve the reduced linear system
         std::vector<int> ipiv(nr);
         int info = 0;
         Tgesv(nr, 1, &Ar[0], nr, &ipiv[0], &fr[0], nr, info);
    
         // Evaluate the residual
         ressweeptime -= getTime();
         (*a).zero();
         for(int jP=0;jP<numP;jP++) {
           int jL = jP*nDir + iDir;
           if (localFWSSpace[jL]!=0) for(int j=0;j<localFWSSpace[jL]->n;j++)  {
             int ir = localFWSSpace[jL]->irow+j;
             (*a).linAdd(fr[ir],*(localFWSSpace[jL]->Ku[j]),
                         -wS*wS*fr[ir],*(localFWSSpace[jL]->Mu[j]));
             Scalar dc = 0.0;
             ScalarTypes::addComplex(dc,wS*fr[ir]);
             if (allOps->C_deriv) if (allOps->C_deriv[0])
               (*a).linAdd(dc,*(localFWSSpace[jL]->Cu[j]));
           }
         }
         (*a).linAdd(-1.0,*f);
         forceAssemble(*a);
         forceAssemble(*f);
         Scalar nrma = *a * *a;
         Scalar nrmf = *f * *f;
         ressweeptime += getTime();
         res[iDir] = sqrt(ScalarTypes::Real(nrma)/ScalarTypes::Real(nrmf));
    
         solutiontime -= getTime();
         if (verboseFlag) filePrint(stderr,"residual %d %d: %e\n",
                   numP,iDir,res[iDir]);

       }
       for(int iDir=0;iDir<nDir;iDir++) if (res[iDir]>maxres) {
         maxres = res[iDir];
         maxDir = iDir;
         maxIS = iS; 
       }
     }
     if (maxres > tol) {
       newW[numNewW] = wcSorted[nP] + double(maxIS+1)/double(ncheck+1)*
                          (wcSorted[nP+1]-wcSorted[nP]);
       if (verboseFlag) filePrint(stderr,"Selecting new frequency %e\n",
                                  newW[numNewW]);
                   
       numNewW++;
       doneF = false;
     }
   }
 }
 romctime += getTime();

 // Compute the residual and output the  ROM solution 
 //  for interval w1,w2
 fsweeptime -= getTime();
 for(int iwr=0;iwr<=numS;iwr++) {
   double wr = w1 + iwr*dw;
   filePrint(stderr,"Solution for iwr= %d wr/2/pi=%e\n",iwr,wr/2.0/M_PI);
// RT uncommenting following line  12/8
   geoSource->setOmega(wr);
   if (!isScattering) if (domain->solInfo().loadcases.size()>0)
      domain->solInfo().loadcases = loadcases_backup;
   for(int iDir=0;iDir<nDir;iDir++) {
     // Compute RHS
     if (isScattering)
       probDesc->setIWaveDir(iDir);
     else {
       if (iDir>0) domain->solInfo().loadcases.pop_front();
     }
     rhsctime2 -= getTime();
     probDesc->getRHS(*f, wr,wr-wc[numP-1]);
     rhsctime2 += getTime();
     // Compute reduced vector
     frsweeptime -= getTime();
     int nr = 0;
     for(int iP=0;iP<numP;iP++) {
       int iL = iP*nDir + iDir;
       if (localFWSSpace[iL]!=0) {
         nr += localFWSSpace[iL]->n;
         for(int j=0;j<localFWSSpace[iL]->n;j++) 
           fr[localFWSSpace[iL]->irow+j] = *f * *(localFWSSpace[iL]->u[j]);
       }
     }
     frsweeptime += getTime();

     // Compute reduced matrix
     Arsweeptime -= getTime();
     for(int ir=0;ir<nr*nr;ir++) Ar[ir] = 0;
     for(int jP=0;jP<numP;jP++) {
       int jL = jP*nDir + iDir;
       if (localFWSSpace[jL]!=0) for(int j=0;j<localFWSSpace[jL]->n;j++)  {
         for(int kP=0;kP<numP;kP++) {
           int kL = kP*nDir + iDir;
           if (localFWSSpace[kL]!=0) for(int k=0;k<localFWSSpace[kL]->n;k++) {
             std::pair<int,int> rowcol
               (localFWSSpace[jL]->imode+j,localFWSSpace[kL]->imode+k); 
             if (Kr.find(rowcol)==Kr.end()) fprintf(stderr,
                         "rowcol %d %d not found\n",
                localFWSSpace[jL]->imode+j,localFWSSpace[kL]->imode+k);
             Ar[localFWSSpace[jL]->irow+j+
             nr*(localFWSSpace[kL]->irow+k)] = 
                 Kr[rowcol] - wr*wr*Mr[rowcol];
               ScalarTypes::addComplex(Ar[localFWSSpace[jL]->irow+j+
                 nr*(localFWSSpace[kL]->irow+k)], wr *Cr[rowcol]);
               }
             }
         }
     }
     Arsweeptime += getTime();
     // Factor and solve the reduced linear system
     std::vector<int> ipiv(nr);
     int info = 0;
     Tgesv(nr, 1, &Ar[0], nr, &ipiv[0], &fr[0], nr, info);

     // Evaluate the residual
     ressweeptime -= getTime();
     (*a).zero();
     for(int jP=0;jP<numP;jP++) {
       int jL = jP*nDir + iDir;
       if (localFWSSpace[jL]!=0) for(int j=0;j<localFWSSpace[jL]->n;j++)  {
         int ir = localFWSSpace[jL]->irow+j;
         (*a).linAdd(fr[ir],*(localFWSSpace[jL]->Ku[j]),
                     -wr*wr*fr[ir],*(localFWSSpace[jL]->Mu[j]));
         Scalar dc = 0.0;
         ScalarTypes::addComplex(dc,wr*fr[ir]);
         if (allOps->C_deriv) if (allOps->C_deriv[0])
           (*a).linAdd(dc,*(localFWSSpace[jL]->Cu[j]));
       }
     }
/* 
     allOps->K->mult(*sol, *a);
     allOps->M->mult(*sol, *b);
     if (allOps->C_deriv) {
        if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
     } else (*c).zero(); 
     (*a).linAdd(-wr*wr,*b);
     for(int k=0;k<sol->size();k++) 
         ScalarTypes::addComplex((*a)[k], wr * (*c)[k]);
*/
     (*a).linAdd(-1.0,*f);
     forceAssemble(*a);
     forceAssemble(*f);
     Scalar nrma = *a * *a;
     Scalar nrmf = *f * *f;
     res[iDir] = sqrt(ScalarTypes::Real(nrma)/ScalarTypes::Real(nrmf));
     ressweeptime += getTime();

     solutiontime -= getTime();
     if (verboseFlag) filePrint(stderr,"fresidual %d: %e\n",
               iDir,res[iDir]);
     // Compute the approximant vector
     (*sol).zero();
     for(int jP=0;jP<numP;jP++) {
       int jL = jP*nDir + iDir;
       if (localFWSSpace[jL]!=0) for(int j=0;j<localFWSSpace[jL]->n;j++)  {
         (*sol).linAdd(fr[localFWSSpace[jL]->irow+j],
                          *(localFWSSpace[jL]->u[j]));
       }
     }
     if(domain->solInfo().isCoupled) scaleDisp(*sol);
     solutiontime += getTime();
     // Output the solution making sure the right frequency appears
     domain->frequencies->push_front(wr);
     if(domain->solInfo().isAcousticHelm()) {  
       SPropContainer& sProps = geoSource->getStructProps();
       for(int iProp=0;iProp<geoSource->getNumProps();iProp++) {
         if(sProps[iProp].kappaHelm!=0.0 ||
            sProps[iProp].kappaHelmImag!=0.0) {
           complex<double> k1 = wr/sProps[iProp].soundSpeed;
           sProps[iProp].kappaHelm = real(k1);
           sProps[iProp].kappaHelmImag = imag(k1);
         } 
       }
     }
     outputtime -= getTime();
     postProcessor->staticOutput(*sol, *rhs,(iwr==numS)&&(iDir==nDir-1));
     outputtime += getTime();
     domain->frequencies->pop_front();
   }
 }
 fsweeptime += getTime();

 time += getTime();

 if (verboseFlag) {
   filePrint(stderr,"Total time: %e\n",time/1e3);
   filePrint(stderr,"----\n");
   filePrint(stderr,"ROM construction time: %e\n",romctime/1e3);
   filePrint(stderr,"transition time: %e\n",transitiontime/1e3);
   filePrint(stderr,"solver time: %e in %d solves\n",
             solvertime/1e3,rebuildscount);
   filePrint(stderr,"one frequency space construction time: %e\n",
             onefspacetime/1e3);
   filePrint(stderr,"rhs construction time in one f: %e\n",rhsctime1/1e3);
   filePrint(stderr,"rhs solve time: %e count: %d\n",rhstime/1e3,solvecount);
   filePrint(stderr,"----\n");
   filePrint(stderr,"Final frequency sweep time: %e\n",fsweeptime/1e3);
   filePrint(stderr,"fr time in sweep: %e\n",frsweeptime/1e3);
   filePrint(stderr,"Ar time in sweep: %e\n",Arsweeptime/1e3);
   filePrint(stderr,"residual computation time in sweep: %e\n",ressweeptime/1e3);
   filePrint(stderr,"rhs construction time in sweep: %e\n",rhsctime2/1e3);
   filePrint(stderr,"solution construction time in sweep: %e\n",solutiontime/1e3);
   filePrint(stderr,"output time: %e\n",outputtime/1e3);
 }

 if (dgpFlag) {
   delete[] UU;
   for(int i=0;i<maxLocalV+1;i++) delete temp_u[i];
   delete[] temp_u;
 }
 delete f;
 delete c;
 delete b;
 delete a;
 delete sol;

 for(int iP=0;iP<numP;iP++) for(int iDir=0;iDir<nDir;iDir++) {
    int iL = iP*nDir + iDir;
    if (localFWSSpace[iL]!=0) delete localFWSSpace[iL];
 }
 delete[] localFWSSpace;

}
