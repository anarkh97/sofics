#ifndef _SFEMINPC_H
#define _SFEMINPC_H

#include <Sfem.d/Sfem.h>

extern Sfem *sfem;
	
template <class Scalar, class VecType>
class SfemInpc : public Sfem {
  int n;
//  double ** E_elem;
 public:
  SfemInpc() { P=sfem->getP(); L= sfem->getL(); ndim=sfem->getndim(); output_order = sfem->getoutput_order(); makealpha(); };
  ~SfemInpc() {};
  void setn(int n1) {n=n1;};
  void printCoefs(VecType* psi_u) {}; // print the coefficients in a file
  void computeMean(VecType* psi_u, VecType *mean_u);
  void computeStdDev(VecType* psi_u, VecType *sdev_u);
  void computePdf(int seed, VecType* psi_u, VecType *realz_u);
//  Scalar  computePdf(int seed, int dofno);
};

#ifdef _TEMPLATE_FIX_
#include <Sfem.d/SfemInpc.C>
#endif

#endif
