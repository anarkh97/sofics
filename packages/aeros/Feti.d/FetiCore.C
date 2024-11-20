#include <Feti.d/Feti.h>
#include <Math.d/BLAS.h>
#include "FetiOp.h"
#include <Driver.d/SubDomain.h>

template<>
void
GenFetiSolver<DComplex>::getRMult(int iSub, GenDistrVector<DComplex> *localvec, DComplex *alpha) const
{
  fprintf(stderr, "WARNING: GenFetiSolver<DComplex>::getRMult(...) not implemented \n");
}

template<>
void
GenFetiSolver<double>::getRMult(int iSub, GenDistrVector<double> *localvec, double *alpha) const
{
 double *lvec = localvec->subData(sd[iSub]->localSubNum());
 double *lAlpha = alpha + fetiOps[iSub]->alphaOffset[0];

 int numRBM = fetiOps[iSub]->numRBM;
 const double *locRBMs = fetiOps[iSub]->locRBMs.data();

 if(numRBM > 0)
 Tgemv('T',sd[iSub]->localLen(), numRBM, 1.0, locRBMs,
       sd[iSub]->localLen(),lvec, 1, 0.0, lAlpha, 1);
}


template<>
void
GenFetiSolver<DComplex>::addRP(int iSub, GenDistrVector<DComplex> *vec1, DComplex *vec2) const
{
  fprintf(stderr, "WARNING: GenFetiSolver<DComplex>::addRP(...) not implemented \n");
}

template<>
void
GenFetiSolver<double>::addRP(int iSub, GenDistrVector<double> *vec1, double *vec2) const
{
 double *localvec = vec1->subData(sd[iSub]->localSubNum());
 double *alpha    = vec2 + fetiOps[iSub]->alphaOffset[0];

 if(fetiOps[iSub]->numRBM > 0 && sd[iSub]->localLen() > 0)
   Tgemv('N', sd[iSub]->localLen(),fetiOps[iSub]->numRBM, 
         -1.0, fetiOps[iSub]->locRBMs.data(), sd[iSub]->localLen(),
         alpha, 1, 1.0, localvec, 1);
}

template<>
void
GenFetiSolver<DComplex>::getRBMs(DComplex *globRBM)
{
  fprintf(stderr, "WARNING: GenFetiSolver<DComplex>::getRBMs(...) not implemented \n");
}

template<>
void
GenFetiSolver<double>::getRBMs(double *globRBM)
{
 int nRBM = GtGsolver->numRBM();
 double *alphas = new double[nRBM*numrbms];
 GtGsolver->getRBMs(alphas);
 int iRBM;
 for(iRBM = 0; iRBM < nRBM; ++iRBM) {
    GenStackVector<double> iAlpha(alphas+iRBM*numrbms, numrbms);
    GenStackDistVector<double> R(internalDI, globRBM+iRBM*(internalDI.len));
    R.zero();
    addR(R, iAlpha);
 }

//CBM--Bug??
/*
 for(iRBM = 0; iRBM < nRBM; ++iRBM) {
    GenStackVector<double> iAlpha(alphas+iRBM*numrbms, numrbms);
    GenStackDistVector<double> R(internalDI, globRBM+iRBM*(internalDI.len));
    R.zero();
    addR(R, iAlpha);
 }
*/
}

#ifndef SALINAS
//CBM
template<>
void
GenFetiSolver<DComplex>::getRBMs(GenDistrVectorSet<DComplex> &globRBM)
{
  fprintf(stderr, "WARNING: GenFetiSolver<DComplex>::getRBMs(...) not implemented \n");
}

template<>
void
GenFetiSolver<double>::getRBMs(GenDistrVectorSet<double> &globRBM)
{
 int nRBM = GtGsolver->numRBM();
 double *alphas = new double[nRBM*numrbms];
 GtGsolver->getRBMs(alphas);
 int iRBM;
 for(iRBM = 0; iRBM < nRBM; ++iRBM) {
    GenStackVector<double> iAlpha(alphas+iRBM*numrbms, numrbms);
    addR(globRBM[iRBM], iAlpha);
 }
}
#endif

template<>
void
GenFetiSolver<DComplex>::addRS(int iSub, GenDistrVector<DComplex> *vec1, DComplex *vec2) const
{
 fprintf(stderr, "WARNING: GenFetiSolver<DComplex>::addRS(...) not implemented \n");
}

template<>
void
GenFetiSolver<double>::addRS(int iSub, GenDistrVector<double> *vec1, double *vec2) const
{
 double *localvec = vec1->subData(iSub);
 double *alpha    = vec2 + eqNums->firstdof(sd[iSub]->subNum() + gOffset);
 int m = sd[iSub]->localLen();
 int n = fetiOps[iSub]->numRBM;
 // localvec = R*alpha
 if(n > 0)
 Tgemv('N', m, n, 1.0, fetiOps[iSub]->locRBMs.data(), m,
                 alpha, 1, 1.0, localvec, 1);
}
