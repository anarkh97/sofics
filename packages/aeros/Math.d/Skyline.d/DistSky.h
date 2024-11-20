#ifndef _DIST_SKY_H_
#define _DIST_SKY_H_

#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Utils.d/MyComplex.h>

template<class Scalar>
class GenDistSky : public GenSkyMatrix<Scalar> 
{
  int firstRow;
  int numRows;
  Scalar *nlines;

 public:
   GenDistSky(const Connectivity *cn, const EqNumberer *dsa, double trbm,
              int fRow, int numRows);
   virtual ~GenDistSky();

   void parallelFactor();
   void reSolve(Scalar *rhs);
   void reSolve(GenVector<Scalar> &rhs);
   void reSolve(int nRHS, Scalar **RHS);
   void reSolve(int nRHS, GenVector<Scalar> *RHS);
   void reSolve(int nRHS, Scalar *RHS);
   void zeroAll();
}; 

typedef GenDistSky<double> DistSky;
typedef GenDistSky<DComplex> ComplexDistSky;

#ifdef _TEMPLATE_FIX_
#include <Math.d/Skyline.d/DistSky.C>
#endif

#endif

