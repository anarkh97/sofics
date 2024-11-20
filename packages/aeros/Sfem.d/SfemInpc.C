#include <Sfem.d/cijk.h>

template<class Scalar, class VecType>
void SfemInpc<Scalar, VecType>::computeMean(VecType* psi_u, VecType* mean_u)
{
  mean_u->copyBlock((*psi_u),0); // 0 means the first block
}

template<class Scalar, class VecType>
void SfemInpc<Scalar, VecType>::computeStdDev(VecType* psi_u, VecType* sdev_u)
{
//  build_psisqr(); 
  sfem->build_psisqr(); // YYY DG isn't it computed already ? check
  psisqr = sfem->getPsiSqr();
  sdev_u->zero();

  if (sfem->isreduced==true) {
   int* nnzindx=sfem->getNnzBlocks();

   for (int i=1;i<P;i++) {
     if (nnzindx[i] ==1) sdev_u->addBlockSqr(i,psisqr[i],(*psi_u)); // YYY DG
   }
  }
  else {
   for (int i=1;i<P;i++) {
     sdev_u->addBlockSqr(i,psisqr[i],(*psi_u));
   }
  }
  sdev_u->computeSqrt();
}

template <class Scalar, class VecType >
void SfemInpc<Scalar, VecType>::computePdf(int seed, VecType* psi_u, VecType* realz_u)
{
  realz_u->zero();
  genXiPsi(seed);

  for(int i=0;i<P;i++)   realz_u->computeRealz(i,psi[i],(*psi_u));
}

/*
template <class Scalar, class VecType >
Scalar SfemInpc<Scalar, VecType>::computePdf(int seed, int dofno)
{
 int i;
 realz_u_1=0;
 genXiPsi(seed);
 for(i=0;i<P;i++) {
   realz_u_1=realz_u_1+(*(psi_u[i]))[dofno]*psi[i];
 }
 return realz_u_1;
}*/

