#ifndef _BIGMATRIX_H_
#define _BIGMATRIX_H_

#if defined(sgi) && ! defined(_OPENMP)
#include <ulocks.h>
#endif

#include <Math.d/matrix.h>

template<class Scalar>
class GenBigMatrix : public GenFullM<Scalar> {
   int size;
   int numBlocks;
   int *blockIndex;
   int **kpvt;

#if defined(sgi) && ! defined(_OPENMP)
   barrier_t *new_barrier;

   void subFactor(int iThread, int numThreads, barrier_t *barrier);
#else
   void subFactor(int iThread, int numThreads);
#endif
   void factorDiagonal(int iThread, int numThreads, int iBlock);
   void lineUpdate(int iThread, int numThreads, int iBlock);
   void lineToColumnCopy(int iThread,int numThreads,int iBlock);
   void rankUpdate(int iThread, int numThreads, int iBlock);
   void symCopy(Scalar *, Scalar *, int, int);
	void subReSolve(int iThread, int numThreads, Scalar *rhs) const;
	void subForward(int iBlock, int, int, Scalar *rhs) const;
	void subBackward(int iBlock, int, int, Scalar *rhs) const;
	void diagSolve(int, int, int, Scalar *rhs) const;
	Scalar *block(int,int);
	const Scalar *block(int,int) const;
	int blockSize(int) const;
  
  public:
   GenBigMatrix(int size);
   virtual ~GenBigMatrix();
   Scalar *data() { return this->v; }
   void parallelFactor();
   void reSolve(Scalar *rhs) const;
   Scalar *operator[](int i);
};

template<class Scalar>
inline Scalar *
GenBigMatrix<Scalar>::operator[](int i)
{
 return this->v+i*size;
}

template<class Scalar>
inline Scalar *
GenBigMatrix<Scalar>::block(int i, int j)
{
	return this->v+blockIndex[i]*size + blockIndex[j];
}

template<class Scalar>
inline const Scalar *
GenBigMatrix<Scalar>::block(int i, int j) const
{
	return this->v+blockIndex[i]*size + blockIndex[j];
}

template<class Scalar>
inline int
GenBigMatrix<Scalar>::blockSize(int i) const
{
 return blockIndex[i+1]-blockIndex[i];
}

typedef GenBigMatrix<double> BigMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/BigMatrix.C>
#endif

#endif
