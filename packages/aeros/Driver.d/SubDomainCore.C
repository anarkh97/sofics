#include <cstdio>
#include <algorithm>
#include <Driver.d/SubDomain.h>
#include <Driver.d/CornerMaker.h>

extern int salinasFlag;

extern "C"      {
void _FORTRAN(cfjacobi)(double *,double *,double *, double *,int&,double&,int &);
}

void
getJacobi(double *kappa, double * mu, FullSquareMatrix &xx,
          double *eigVal,int nsmax, int subSpaceSize, double tolJac)
{
  int i,j;

  _FORTRAN(cfjacobi)(kappa,mu,xx[0],eigVal,nsmax,tolJac,subSpaceSize);

  // sort eigenvalues.

  int is = 1;
  while(is != 0) {
    is = 0;
    for(i=1; i<subSpaceSize; ++i) {
      if(eigVal[i] < eigVal[i-1] ) {
        is = 1;
        std::swap( eigVal[i-1], eigVal[i] );
        for(j=0; j<subSpaceSize; ++j) {
          std::swap( xx[i][j], xx[i-1][j] );
        }
      }
    }
  }
}

void
getJacobi(DComplex *kappa, DComplex *mu, GenFullSquareMatrix<DComplex> &xx,
          DComplex *eigVal, int nsmax, int subSpaceSize, double tolJac)
{
  fprintf(stderr, " *** WARNING: getJacobi(...) not implemented for DComplex type \n");
}

template<>
void
GenSubDomain<DComplex>::getSRMult(const DComplex *lvec, const DComplex *interfvec,
                                  int nRBM, const double *locRBMs, DComplex *alpha) const
{
  fprintf(stderr, " *** WARNING: GenSubDomain<DComplex>::getSRMult(...) not implemented \n");
}

template<>
void
GenSubDomain<double>::getSRMult(const double *lvec, const double *interfvec,
                                int nRBM, const double *locRBMs,
                                double *alpha) const
{
 double *localvec = (double *) dbg_alloca(sizeof(double)*localLen());

 // Add the interface vector (interfvec) contribution to localvec
 for(int iDof = 0; iDof <localLen(); ++iDof)
   localvec[iDof] = lvec[iDof];

 for(int iDof = 0; iDof < scomm->lenT(SComm::std); ++iDof)
   localvec[scomm->boundDofT(SComm::std,iDof)] += interfvec[iDof];

 Tgemv('T',localLen(), nRBM, -1.0, locRBMs,
       localLen(),localvec, 1, 0.0, alpha, 1);
}


template<> double GenSubDomain<double>::Bcx(int i) { return (bcx) ? bcx[i] : bcxC[i].real(); }
template<> DComplex GenSubDomain<DComplex>::Bcx(int i) { return (bcxC) ? bcxC[i] : DComplex(bcx[i],0.0); }

