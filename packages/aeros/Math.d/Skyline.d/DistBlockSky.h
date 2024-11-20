#ifndef _DIST_BLOCK_SKY_H_
#define _DIST_BLOCK_SKY_H_

#include <Math.d/Skyline.d/BlockSky.h>
#include <Utils.d/MyComplex.h>

template<class Scalar> class GenVector;

template<class Scalar>
class GenDistBlockSky : public GenBlockSky<Scalar> 
{
  int firstRow;
  int numRows;
  Scalar *nlines;

 public:
   GenDistBlockSky(const Connectivity *cn, EqNumberer *dsa, double trbm,
                   int fRow, int numRows);
   virtual ~GenDistBlockSky();
   void parallelFactor();
   void reSolve(Scalar *rhs);
   void reSolve(GenVector<Scalar> &rhs);
}; 

typedef GenDistBlockSky<double> DistBlockSky;
typedef GenDistBlockSky<DComplex> ComplexDistBlockSky;

#ifdef _TEMPLATE_FIX_
  #include <Math.d/Skyline.d/DistBlockSky.C>
#endif

#endif

