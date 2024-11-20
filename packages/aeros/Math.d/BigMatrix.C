#if defined(sgi) && ! defined(_OPENMP)
#include <ulocks.h>
#endif
#include <sys/types.h>

#ifndef WINDOWS
 #include <sys/mman.h>
#endif

#include <Threads.d/PHelper.h>
#include <Utils.d/linkfc.h>

// BLAS level three real Matrix Product

extern "C" {
  void _FORTRAN(dsisl)(const double *, const int &, const int &, int *, double *);
  void _FORTRAN(dsifa)(const double *, const int &, const int &, int *, int &);
}


extern int zeroFd;

template<class Scalar>
GenBigMatrix<Scalar>::~GenBigMatrix()
{
  if(blockIndex) { delete blockIndex; blockIndex = 0; }
  if(kpvt) { delete kpvt; kpvt = 0; } 
}

template<class Scalar>
GenBigMatrix<Scalar>::GenBigMatrix(int _size): GenFullM<Scalar>(_size) {
 size = _size;

 // compute number of blocks
 int sizeOfBlock = 96;
 numBlocks = (size+sizeOfBlock-1)/sizeOfBlock;

 // do the block splitting
 blockIndex = new int[numBlocks+1];
 int *pvts  = new int[size];
 kpvt = new int*[numBlocks];
 int iBlock;
 for(iBlock = 0; iBlock < numBlocks; ++iBlock) {
   blockIndex[iBlock] = sizeOfBlock*iBlock;
   kpvt[iBlock] = pvts+blockIndex[iBlock];
 }
 blockIndex[numBlocks] = size;
}

template<class Scalar>
void
GenBigMatrix<Scalar>::symCopy(Scalar *source, Scalar *dest, int nRow, int ncol)
{
 int i,j;
 for(i = 0; i+7 < ncol; i+= 8)
   for(j = 0; j <nRow; ++j) {
      dest[i*size + j] = source[j*size + i];
      dest[(i+1)*size + j] = source[j*size + (i+1)];
      dest[(i+2)*size + j] = source[j*size + (i+2)];
      dest[(i+3)*size + j] = source[j*size + (i+3)];
      dest[(i+4)*size + j] = source[j*size + (i+4)];
      dest[(i+5)*size + j] = source[j*size + (i+5)];
      dest[(i+6)*size + j] = source[j*size + (i+6)];
      dest[(i+7)*size + j] = source[j*size + (i+7)];
  }
 for(; i < ncol; ++i)
   for(j = 0; j <nRow; ++j)
      dest[i*size + j] = source[j*size + i];
}

// Operation of dgemm
// C := alpha*op( A )*op( B ) + beta*C

// Arguments of dgemm
// transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc

template<class Scalar>
void
GenBigMatrix<Scalar>::rankUpdate(int iThread, int numThreads, int iBlock)
{
 int rem = iBlock%numThreads;
 int myBlock = iBlock-rem+iThread;
 if(myBlock <= iBlock) myBlock += numThreads;

 Scalar alpha = -1.0;
 Scalar beta  =  1.0;

 for(;myBlock < numBlocks; myBlock += numThreads) {
 // the matrix a[myBlock:myBlock+1][iBlock+1:size]
 // needs to be updated by -a[iBlock+1:size][iBlock:iBlock+1]*
 //           a[iBlock:iBlock+1][myBlock:myBlock+1]

    Scalar *a = block(iBlock,myBlock);

    Scalar *b = block(myBlock,iBlock);

    Scalar *c = block(myBlock,myBlock);

   Tgemm('N','N',blockSize(myBlock),size-blockIndex[myBlock],
                             blockSize(iBlock),
                             alpha,a,size,b,size,beta,c,size);
 }

}

#if defined(sgi) && ! defined(_OPENMP)
template<class Scalar>
void
GenBigMatrix<Scalar>::subFactor(int iThread, int numThreads, barrier_t *bar)
#else
template<class Scalar>
void
GenBigMatrix<Scalar>::subFactor(int iThread, int numThreads)
#endif
{
 int iBlock;
 for(iBlock = 0; iBlock < numBlocks; ++iBlock) {
   // factor the diagonal block
   factorDiagonal(iThread, numThreads, iBlock);
#if defined(sgi) && ! defined(_OPENMP)
   barrier(bar, numThreads);
#endif
   // LineUpdate
   lineUpdate(iThread, numThreads, iBlock);
#if defined(sgi) && ! defined(_OPENMP)
   barrier(bar, numThreads); 
#endif

   // rank update
   rankUpdate(iThread, numThreads, iBlock);
#if defined(sgi) && ! defined(_OPENMP)
   barrier(bar, numThreads);
#endif
 }
}

template<class Scalar>
void
GenBigMatrix<Scalar>::parallelFactor()
{
 // create barrier using threadManager
#if defined(sgi) && ! defined(_OPENMP)
 new_barrier = threadManager->getBarrier();

 // call subFactor in parallel using numThreads
 int numThreads = threadManager->numThr();
 execParal(numThreads, this, &GenBigMatrix<Scalar>::subFactor, numThreads, new_barrier );
#else
 int numThreads = threadManager->numThr();
 execParal(numThreads, this, &GenBigMatrix<Scalar>::subFactor, numThreads);
#endif

}

template<class Scalar>
void
GenBigMatrix<Scalar>::subForward(int iBlock, int iThread, int numThreads, Scalar *rhs) const
{
 int rem = iBlock%numThreads;
 int myBlock = iBlock-rem+iThread;
 if(myBlock <= iBlock) myBlock += numThreads;

 for(; myBlock < numBlocks; myBlock += numThreads)
   Tgemv('T', blockSize(iBlock),blockSize(myBlock),
             -1.0, block(myBlock,iBlock), size, rhs+blockIndex[iBlock],
          1, 1.0, rhs+blockIndex[myBlock], 1);
}

template<class Scalar>
void
GenBigMatrix<Scalar>::subBackward(int iBlock, int iThread, int numThreads, Scalar *rhs) const
{
 int rem = iBlock%numThreads;
 int myBlock = iBlock-rem+iThread;
 if(myBlock >= iBlock) myBlock -= numThreads;
 for(; myBlock >= 0; myBlock -= numThreads)
    Tgemv('T', blockSize(iBlock),blockSize(myBlock),
             -1.0, block(myBlock,iBlock), size, rhs+blockIndex[iBlock],
          1, 1.0, rhs+blockIndex[myBlock], 1);
}

template<class Scalar>
void
GenBigMatrix<Scalar>::subReSolve(int iThread, int numThreads, Scalar *rhs) const
{
 int iBlock;
 for(iBlock = 0; iBlock < numBlocks; ++iBlock) {
   // Forward substitute of block iBlock
   subForward(iBlock, iThread, numThreads, rhs);
#if defined(sgi) && ! defined(_OPENMP)
   barrier(new_barrier , numThreads);
#endif
 }
 for(iBlock = numBlocks; iBlock--; ) {
#if defined(sgi) && ! defined(_OPENMP)
   barrier(new_barrier , numThreads);
#endif
   diagSolve(iBlock, iThread, numThreads, rhs);
#if defined(sgi) && ! defined(_OPENMP)
   barrier(new_barrier , numThreads);
#endif
   subBackward(iBlock, iThread, numThreads, rhs);
 }
}


template<class Scalar>
void
GenBigMatrix<Scalar>::reSolve(Scalar *rhs) const
{
 int numThreads = threadManager->numThr();
 execParal(numThreads, this, &GenBigMatrix<Scalar>::subReSolve, numThreads, rhs);
}
