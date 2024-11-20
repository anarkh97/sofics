#ifndef CHAOS_H
#define CHAOS_H

/*****************************************************************************/
/*             A Class to Manipulate Hermite Polynomial                      */
/*                    Written by George Saad                                 */
/*****************************************************************************/

template <class Scalar> class GenVector;
template <class Scalar> class GenFullM;

int nterm(int x, int y);        /* function to determine the number of terms corresponding to each specific order */

class Chaos {

private:
  int ndim;                    /* number of dimensions in the Chaos Expansionn */
  int order;                   /* order of the expansion */
                   
 public:
  Chaos(int dim, int order);
  int tnterms();
  void HermiteToPowers(int &maxorder, GenFullM<double> &HP);
  void ScHermiteToPowers(int &maxorder, double &scale, GenFullM<double> &HP);
  void mult1dHermite(GenVector<double> &A, GenVector<double> &B, GenVector<double> &C);
  void OrthoOrder(GenFullM<double> &D);
};

#endif
