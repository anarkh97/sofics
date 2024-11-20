#ifndef CIJK_H
#define CIJK_H
/*****************************************************************************/
/*          A Class to to Calculate the Expectated Value                     */
/*              of a product of Hermite polynomials                          */
/*****************************************************************************/

#include <Sfem.d/chaos.h>

class Cijk {

 private:
  int ndim;         /* number of dimensions in expansion */
  int order;        /* order of expansion */
  int maxorder, maxnterms;
  int i, j, k;
  Chaos PC1;
 
 public:
  Cijk(int ndim, int order, int P);
  double expectation(int i, int j, int k);
};

#endif
