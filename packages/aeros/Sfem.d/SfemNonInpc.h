#ifndef _SFEMNONINPC_H
#define _SFEMNONINPC_H

#include <Sfem.d/Sfem.h>

extern Sfem *sfem;

template <class Scalar, class VecType>
class SfemNonInpc : public Sfem {
  int n;
 public:
  SfemNonInpc() { P=sfem->getP(); L=sfem->getL(); ndim=sfem->getndim(); output_order=sfem->getoutput_order(); makealpha();};
  ~SfemNonInpc() {};
  void setn(int n1) {n=n1;};
  void update_Psi_u(VecType* u, VecType* psi_u); 
  void compute_Coefs(VecType* psi_u);
  void printCoefs(VecType* psi_u) {}; // print the coefficients in a file
  void computeMean(VecType* psi_u, VecType *mean_u);
  void computeStdDev(VecType* psi_u, VecType *sdev_u);
  void computePdf(int seed, VecType* psi_u, VecType *realz_u);
//  Scalar  computePdf(int seed, int dofno);
};

#ifdef _TEMPLATE_FIX_
#include <Sfem.d/SfemNonInpc.C>
#endif

#endif
