#include <Sfem.d/cijk.h>

template <class Scalar, class VecType >
void SfemNonInpc<Scalar, VecType>::update_Psi_u(VecType* u, VecType* psi_u) 
{
  psi = sfem->getpsi();
  for (int i=0;i<P;i++)  psi_u->updateBlock(i,psi[i],(*u));
}


template <class Scalar, class VecType >
void SfemNonInpc<Scalar, VecType>::compute_Coefs(VecType* psi_u)
{
  sfem->build_psisqr();
  psisqr = sfem->getPsiSqr();
  for (int i=0;i<P;i++) {
    Scalar f =  (domain->solInfo().nsample)*psisqr[i];
    psi_u->scaleBlock(i,f); 
  }
}


template <class Scalar, class VecType >
void SfemNonInpc<Scalar, VecType>::computeMean(VecType* psi_u, VecType *mean_u)
{
  mean_u->copyBlock((*psi_u),0); // 0 means the first block
}


template <class Scalar, class VecType >
void SfemNonInpc<Scalar, VecType>::computeStdDev(VecType* psi_u, VecType* sdev_u)
{
  sdev_u->zero();
  for (int i=1;i<P;i++) {
    sdev_u->addBlockSqr(i,psisqr[i],(*psi_u));
  }
  sdev_u->computeSqrt();
}


template <class Scalar, class VecType >
void SfemNonInpc<Scalar, VecType>::computePdf(int seed, VecType* psi_u, VecType *realz_u)
{
  realz_u->zero();
  genXiPsi(seed);
  for(int i=0;i<P;i++)   realz_u->computeRealz(i,psi[i],(*psi_u));  
}

/*
template <class Scalar, class VecType >
Scalar SfemNonInpc<Scalar, VecType>::computePdf(int seed, int dofno)
{
 int i;
 realz_u_1=0;
 genXiPsi(seed);
 for(i=0;i<P;i++) {
   realz_u_1=realz_u_1+(*(psi_u[i]))[dofno]*psi[i];
 }
 return realz_u_1; 
}*/
