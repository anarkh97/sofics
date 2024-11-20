#ifndef _DIST_SPARSE_H_
#define _DIST_SPARSE_H_

#include <Math.d/BLKSparseMatrix.h>
#include <Utils.d/MyComplex.h>

template<class Scalar>
class GenDistBLKSparse : public GenBLKSparseMatrix<Scalar> 
{
  int firstRow;
  int numRows;
  Scalar *nlines;

 public:
   GenDistBLKSparse(const Connectivity *cn, const EqNumberer *dsa, double trbm, const SolverCntl &scntl,
                    int fRow, int numRows);
   virtual ~GenDistBLKSparse();

   void factor();
   void reSolve(Scalar *rhs);
   void reSolve(GenVector<Scalar> &rhs);

}; 

typedef GenDistBLKSparse<double> DistBLKSparse;
typedef GenDistBLKSparse<DComplex> ComplexDistBLKSparse;

#ifdef _TEMPLATE_FIX_
  #include <Math.d/DistBLKSparse.C>
#endif

#endif

