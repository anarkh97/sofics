#ifndef _SPARSE_SET_
#define _SPARSE_SET_

#include <memory>
#include <Utils.d/MyComplex.h>
#include <Utils.d/resize_array.h>

template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenFullM;
typedef GenFullM<DComplex> FullMC;

template<class Scalar>
class GenSparseSet 
{
   ResizeArray<std::shared_ptr<GenSparseMatrix<Scalar>> > sm;
   int numSM;
 public:
   GenSparseSet(int n=1);
   ~GenSparseSet();
   int numRow();
   int numCol();
   long size() { return 0;}
   int num() { return numSM; }
   void setSparseMatrix(int i, std::shared_ptr<GenSparseMatrix<Scalar>> _sm);
   int addSparseMatrix(std::shared_ptr<GenSparseMatrix<Scalar>> _sm);
 
   void mult(Scalar *rhs, Scalar *result);
   void multAdd(Scalar *rhs, Scalar *result);
   void multIdentity(Scalar *result);
   void multIdentity(Scalar **result);
   void transposeMult(Scalar *rhs, Scalar *result);
   void transposeMultAdd(Scalar *rhs, Scalar *result);
   void transposeMultSubtract(Scalar *rhs, Scalar *result);
   void multSub(int nRHS, Scalar **rhs, Scalar **result);

   void multIdentity(DComplex **result, int col1, int col2);
   void transposeMultSubtract(int nRHS, DComplex **rhs, DComplex **result);
   void transposeMultSubtract(int nRHS, DComplex **rhs, FullMC *result, 
                              int col1=0);
   void transposeMultSubtract(int nRHS, DComplex **rhs, DComplex *result, 
                              int col1=0);
   void multSub(DComplex *rhs, DComplex *result);
};

typedef GenSparseSet<double> SparseSet;
typedef GenSparseSet<DComplex> ComplexSparseSet;

#ifdef _TEMPLATE_FIX_
#include <Math.d/SparseSet.C>
#endif

#endif
