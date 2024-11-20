#include <Math.d/Skyline.d/BlockSky.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/matrix.h>
#include <Math.d/SymFullMatrix.h>
#include <Math.d/Skyline.d/BlockSky.h>
#include <Math.d/FullMatrix.h>
#include <Math.d/Vector.h>
#include <Threads.d/Paral.h>
#include <Threads.d/PHelper.h>
#include <Driver.d/Communicator.h>
/*****************************************************************************
 *                   Copyright (C) 1999 CMSoft                               *
 *                                                                           *
 *  These lines of code and declarations contain unpublished proprietary     *
 *  information of CMSoft. They may not be copied or duplicated in whole     *
 *  or part without prior authorization from CMSoft.                         *
 *                                                                           *
 *****************************************************************************/


#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _TGEMM__
#define _TGEMM__
inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const double &f, double *g, const int &h,
                  double *i, const int &j, const double &k, double *l,
                  const int &m)
{
 _FORTRAN(dgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}

inline void Tgemm(const char &a, const char &b, const int &c,const int &d,
                  const int &e, const complex<double> &f, complex<double> *g, const int &h,
                  complex<double> *i, const int &j, const complex<double> &k, complex<double> *l,
                  const int &m)
{
 _FORTRAN(zgemm)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}
#endif

#ifndef _TGEMV__
#define _TGEMV__
inline void Tgemv(const char &a, const int &b, const int &c,
                  const double &d, double *e, const int &f,
                  double *g, const int &h, const double &i, double *j, const int &k)
{
 _FORTRAN(dgemv)(a,b,c,d,e,f,g,h,i,j,k);
}

inline void Tgemv(const char &a, const int &b, const int &c,
                  const complex<double> &d, complex<double> *e, const int &f,
                  complex<double> *g, const int &h, const complex<double> &i, complex<double> *j, const int &k)
{
 _FORTRAN(zgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif

template<class Scalar>
GenBlockSky<Scalar>::~GenBlockSky()
{
 if(dlp) { delete [] dlp; dlp=0; }
 if(perm) { delete [] perm; perm=0; }
 if(blWeight) { delete [] blWeight; blWeight=0; }
 if(blHeight) { delete [] blHeight; blHeight=0; }
 if(blTop) { delete [] blTop; blTop=0; }
 if(firstDof) { delete [] firstDof; firstDof=0; }
 if(lastCol) { delete [] lastCol; lastCol=0; }
 if(invDiag) { delete [] invDiag; invDiag=0; }
 if(skyA) { delete [] skyA; skyA=0; }
 if(sing) { delete [] sing; sing=0; }
 if(myRCN && rowColNum) { delete [] rowColNum; rowColNum = 0; }
}

template<class Scalar>
void
GenBlockSky<Scalar>::initialize() 
{
  dlp=0; perm=0; blWeight=0; blHeight=0; blTop=0; firstDof=0;
  lastCol=0; invDiag=0; skyA=0; sing=0; myRCN=0; rowColNum=0;
}

//
// Constructor for Coarse Problems
//
template<class Scalar>
GenBlockSky<Scalar>::GenBlockSky(const Connectivity *nodeToNode, const EqNumberer *eqnums,
                                 double _t)
{
 initialize();

 tol = _t;
 neq = eqnums->size();
 dlp = new int[neq];
 perm = new int[neq];

 int i, j;
 nBlocks = nodeToNode->csize();

 // get inverse renumbering
 auto invRenum = eqnums->renumPtr();
 int *renum = new int[nBlocks];
 int *blNum = new int[nBlocks];

 for(i=0; i<nBlocks; ++i)
   renum[i] = -1;
 // get the direct renumbering
 for(i=0; i<nBlocks; ++i) {
   renum[invRenum[i]] = i; 
 } 

 blWeight = new int[nBlocks];
 blHeight = new int[nBlocks];
 blTop    = new int[nBlocks+1];
 firstDof = new int[nBlocks];
 lastCol  = new int[nBlocks];
 int blN = 0;
 for(i = 0; i < nBlocks; ++i) {
   int nodeI = renum[i];
   if(nodeI < 0) break;
   blWeight[blN] = eqnums->weight( nodeI );
   if(blWeight[blN] <= 0) continue;
   firstDof[blN] = eqnums->firstdof( nodeI );
   
   int minDof = firstDof[blN];
   lastCol[blN] = blN;
   for(j = 0; j < nodeToNode->num(nodeI); ++j) {
     int jNode = (*nodeToNode)[nodeI][j];
     int jDof = eqnums->firstdof( jNode );
     if(jDof < 0) continue;
     if(jDof < minDof)
       minDof = jDof;
     if(invRenum[jNode] < i)
       lastCol[blNum[invRenum[jNode]]] = blN;
   }
   blHeight[blN] = firstDof[blN] - minDof + blWeight[blN];
   blNum[i] = blN++;
 }

 nBlocks = blN;
 blTop[0] = 0;
 int eq  =0;
 maxBlockSize = 0;
 for(i = 0; i < nBlocks; ++i) {
   blTop[i+1] = blTop[i] + blWeight[i]*blHeight[i];
   if(blWeight[i]*blHeight[i] > maxBlockSize)
       maxBlockSize = blWeight[i]*blHeight[i];
   for(j=0; j<blWeight[i]; ++j) {
     dlp[eq] = blTop[i] + j + j*blHeight[i] + blHeight[i] - blWeight[i]; 
     eq++;
   }
 }

 // Finally update lastcol, because it's currently only based on connectivity
 int lastSeen = 0;
 for(i = 0; i < nBlocks; ++i)
   if(lastSeen > lastCol[i]) lastCol[i] = lastSeen;
   else lastSeen = lastCol[i];

 skyA = new Scalar[blTop[nBlocks]];
#ifdef DISTRIBUTED
 for(i = 0; i < blTop[nBlocks]; ++i)
   skyA[i] = 0.0;
#else
 #ifdef _OPENMP
 #pragma omp parallel private(i)
 #endif
 {
   int myid = 0, numThreads = 1;
   #ifdef _OPENMP
   myid = omp_get_thread_num();
   numThreads = omp_get_num_threads();
   #endif
   for(i = myid; i < nBlocks; i+=numThreads) {
     for(int j = 0; j < blWeight[i]*blHeight[i]; j++)
         skyA[blTop[i]+j] = 0.0;
   }
 }
#endif

}

//
// Constructor for local problems
//
template<class Scalar>
GenBlockSky<Scalar>::GenBlockSky(const Connectivity *nodeToNode, const DofSetArray *eqnums,
 double _t)
{
 initialize(); 
 tol = _t;
 neq = eqnums->size();
 dlp = new int[neq];
 perm = new int[neq];
 // get unconstrained numbering
 rowColNum = eqnums->getUnconstrNum().data();

 int i, j;
 nBlocks = nodeToNode->csize();

 // get inverse renumbering
 auto invRenum = eqnums->renumPtr();
 int *renum = new int[nBlocks];
 int *blNum = new int[nBlocks];

 for(i=0; i<nBlocks; ++i)
   renum[i] = -1;
 // get the direct renumbering
 for(i=0; i<nBlocks; ++i) {
   renum[invRenum[i]] = i; 
 } 
 
 blWeight = new int[nBlocks];
 blHeight = new int[nBlocks];
 blTop    = new int[nBlocks+1];
 firstDof = new int[nBlocks];
 lastCol  = new int[nBlocks];
 int blN = 0;
 for(i = 0; i < nBlocks; ++i) {
   int nodeI = renum[i];
   if(nodeI < 0) break;
   blWeight[blN] = eqnums->weight( nodeI );
   if(blWeight[blN] <= 0) continue;
   firstDof[blN] = eqnums->firstdof( nodeI );
   
   int minDof = firstDof[blN];
   lastCol[blN] = blN;
   for(j = 0; j < nodeToNode->num(nodeI); ++j) {
     int jNode = (*nodeToNode)[nodeI][j];
     int jDof = eqnums->firstdof( jNode );
     if(jDof < 0) continue;
     if(jDof < minDof)
       minDof = jDof;
     if(invRenum[jNode] < i)
       lastCol[blNum[invRenum[jNode]]] = blN;
   }
   blHeight[blN] = firstDof[blN] - minDof + blWeight[blN];
   blNum[i] = blN++;
 }

 nBlocks = blN;
 blTop[0] = 0;
 int eq  =0;
 maxBlockSize = 0;
 for(i = 0; i < nBlocks; ++i) {
   blTop[i+1] = blTop[i] + blWeight[i]*blHeight[i];
   if(blWeight[i]*blHeight[i] > maxBlockSize)
       maxBlockSize = blWeight[i]*blHeight[i];
   for(j=0; j<blWeight[i]; ++j) {
     dlp[eq] = blTop[i] + j + j*blHeight[i] + blHeight[i] - blWeight[i]; 
     eq++;
   }
 }

 // Finally update lastcol, because it's currently only based on connectivity
 int lastSeen = 0;
 for(i = 0; i < nBlocks; ++i)
   if(lastSeen > lastCol[i]) lastCol[i] = lastSeen;
   else lastSeen = lastCol[i];

 skyA = new Scalar[blTop[nBlocks]];
#ifdef DISTRIBUTED
 for(i = 0; i < blTop[nBlocks]; ++i)
   skyA[i] = 0.0;
#else
 #ifdef _OPENMP
 #pragma omp parallel private(i)
 #endif
 {
   int myid = 0, numThreads = 1;
   #ifdef _OPENMP
   myid = omp_get_thread_num();
   numThreads = omp_get_num_threads();
   #endif
   for(i = myid; i < nBlocks; i+=numThreads) {
     for(int j = 0; j < blWeight[i]*blHeight[i]; j++)
         skyA[blTop[i]+j] = 0.0;
   }
 }
#endif
 delete [] renum;
 delete [] blNum;
}

//
// constructor for Dirichlet preconditioner
//
template<class Scalar>
GenBlockSky<Scalar>::GenBlockSky(const Connectivity *nodeToNode, const DofSetArray *dsa, double _t,
                                 int *rCN)
{
  initialize();
  tol = _t;

  int i, j, n, k;
  int ndof = dsa->size();

  nBlocks = dsa->numNodes();
  if( nodeToNode->csize() < nBlocks ) nBlocks = nodeToNode->csize();

// Count & renumber unconstrained dofs & mark constrained dofs

  neq  = 0;
  myRCN = 1;
  int *build_rowColNum = new int[neq];
  for(i = 0; i < neq; ++i) {
      build_rowColNum[i] = rCN[i];
    if(build_rowColNum[i] >= 0) neq++;
  }
  rowColNum = build_rowColNum;

// Build the dofTON in unconstrained numbering
  int *dofToN = new int[neq]; 

  for(i = 0; i < nBlocks; ++i) {
    int fdof = dsa->firstdof(i);
    int ndof = dsa->weight(i);
    for(j = 0; j < ndof; ++j) {
      if(rowColNum[j+fdof] >= 0) dofToN[rowColNum[j + fdof]] = i;
    }
  }

  int *uncFirstDof = new int[nBlocks];

  for(i=0; i < nBlocks; ++i)
    uncFirstDof[i] = -1;

  for(i=0; i < neq; ++i) {
    int nodeNum = dofToN[i];
    if(uncFirstDof[ nodeNum ] < 0)
       uncFirstDof[ nodeNum ] = i;
   }

// Allocate memory for diagonal location pointer
  dlp  = new int[neq];
  perm = new int[neq];

// Initialize first element of diagonal location pointers to one.
  dlp[0] = 1;

// Loop over Equations to find minimum Dof connection.
  for(n = 1; n < neq; ++n) {
    int thisNode = dofToN[n];
    int minDof   = n;
    for(i = 0; i < nodeToNode->num(thisNode); ++i) {
      k = uncFirstDof[(*nodeToNode)[thisNode][i] ];
      if(k >= 0 && k < minDof) minDof = k;
    }
    minDof = minDof - (minDof % 4);
    dlp[n] = dlp[n-1] + n - minDof + 1;
  }

  blWeight = new int[nBlocks];
  blHeight = new int[nBlocks];
  blTop    = new int[nBlocks+1];
  firstDof = new int[nBlocks];
  lastCol  = new int[nBlocks];

  skyA = new Scalar[blTop[nBlocks]];

  for(i = 0; i < blTop[nBlocks]; ++i)
    skyA[i] = 0.0;

  delete [] dofToN;
  delete [] uncFirstDof;
}

template<class Scalar>
void
GenBlockSky<Scalar>::factor()
{
 invDiag = new Scalar[neq];
 Scalar *ltmp = new Scalar[maxBlockSize];

 int i,j;

 int iTop, jTop;

 nzem = 0;
 Scalar *origDiag = new Scalar [neq];
 sing = new int [neq]; // PJSA
 for(i = 0; i < neq; ++i) {
   origDiag[i] = skyA[dlp[i]];
   sing[i] = 0; // PJSA
 }

 for(int iBlock = 0; iBlock < nBlocks; ++iBlock) {
  // First dof of the iBlock
  int fIDof = firstDof[iBlock];
  // Dimension of the iBlock column, not including the diagonal block
  int ldl = blHeight[iBlock] - blWeight[iBlock];
  // Degree of freedom corresponding to the top of iBlock
  iTop = fIDof - ldl;
  // Build the row of L corresponding to this block
  for(i = 0; i < ldl; ++i)
    for(j = 0; j <  blWeight[iBlock]; ++j) 
      ltmp[j*ldl+i] = 
          invDiag[iTop+i] * skyA[blTop[iBlock]+j*blHeight[iBlock]+i];

  // factor the diagonal block
  GenStackFSFullMatrix<Scalar> diagBlock(blWeight[iBlock], blWeight[iBlock], 
     blHeight[iBlock], skyA+blTop[iBlock]+ldl);
  if(ldl >0)
      Tgemm('T', 'N', blWeight[iBlock], blWeight[iBlock], ldl,
        -1.0, ltmp, ldl, skyA+blTop[iBlock],
         blHeight[iBlock], 1.0, skyA+blTop[iBlock]+ldl, blHeight[iBlock]);
  nzem += diagBlock.symLuFactor(perm+fIDof, tol, origDiag+fIDof, sing+fIDof);  // PJSA

  for(i = 0; i < blWeight[iBlock]; ++i) {
#ifdef USE_PERM
    invDiag[firstDof[iBlock] + perm[fIDof+i]] = diagBlock[i][i];
#else
    invDiag[firstDof[iBlock] + i] = diagBlock[i][i];
#endif
  }

  // Now update the iBlock row
  for(int jBlock = iBlock+1; jBlock <= lastCol[iBlock]; ++jBlock) {
    // jTop is the top degree of freedom of the j supernode
    jTop = firstDof[jBlock] - blHeight[jBlock] + blWeight[jBlock];

    if(jTop > firstDof[iBlock] + blWeight[iBlock]) continue;
    int iOffset, jOffset;
    int commonHeight;

    if(iTop >= jTop) {
     commonHeight = ldl;
     iOffset = 0;
     jOffset = iTop-jTop;
    } else {
     commonHeight = firstDof[iBlock]-jTop;
     jOffset = 0;
     iOffset = jTop-iTop;
    }
    
    // Final update of the aBlock
    if(commonHeight > 0)
     Tgemm('T', 'N', blWeight[iBlock], blWeight[jBlock], commonHeight,
        -1.0, ltmp+iOffset, ldl, skyA+blTop[jBlock]+jOffset,
         blHeight[jBlock], 1.0, skyA+blTop[jBlock]+
               firstDof[iBlock]-jTop, blHeight[jBlock]);
    if(jTop <= firstDof[iBlock])
      diagBlock.Lm1Mult(skyA+blTop[jBlock]+firstDof[iBlock]-jTop,
                     blWeight[jBlock], blHeight[jBlock], perm+fIDof);
  }
 }
 delete [] ltmp;
 delete [] origDiag;
 // if(nzem > 0) fprintf(stderr, " ... Found %d singularities\n", nzem);
}


template<class Scalar>
void
GenBlockSky<Scalar>::parallelFactor()
{
#if defined(sgi) && !defined(_OPENMP)
  int avg = size() / neq;
  int maxCPU = (avg*avg)/20000;
//if(maxCPU > 1 ) {
  if(maxCPU > 1 && threadManager->numThr() > 1) {
    barrier_t *barrier = threadManager->getBarrier();
    int npr = threadManager->numThr();
    if(npr > maxCPU) npr = maxCPU;
    filePrint(stderr," ... BlockSky parallel factor: neq %d avg. band. %d #CPUs %d\n",
              neq, avg, npr);
    // fprintf(stderr, " ... BlockSky parallel factor with %d CPUs (avg. bandwidth: %d) \n", maxCPU, avg);
    invDiag = new Scalar[neq];
    Scalar *ltmp = new Scalar[maxBlockSize];
    Scalar *origDiag = new Scalar[neq];
    execParal(npr, this, GenBlockSky<Scalar>::pFactor, npr, barrier, ltmp, invDiag, origDiag);
  } 
  else {
    filePrint(stderr, " ... BlockSky sequential factor (avg. bandwidth: %d) \n", avg);
    factor();
  }
#else
 invDiag = new Scalar[neq];
 Scalar *ltmp = new Scalar[maxBlockSize];

 int i,j;
 int iTop, jTop;

 nzem = 0;
 Scalar *origDiag = new Scalar [neq];
 for(i = 0; i < neq; ++i)
   origDiag[i] = skyA[dlp[i]];
 int iBlock;

 #ifdef _OPENMP
 #pragma omp parallel private(i,iBlock,j,iTop,jTop)
 #endif
 {
   // double *ltmp = bigLtmp;
   // double *ltmp2 = bigLtmp + maxBlockSize;
   int myid = 0, numThreads = 1;
   #ifdef _OPENMP
   myid = omp_get_thread_num();
   numThreads = omp_get_num_threads();
   fprintf(stderr, "I am %d out of %d\n", myid, numThreads);
   #endif
   int myFirstBlock = myid;
   for(iBlock = 0; iBlock < nBlocks; ++iBlock) {
    // First dof of the iBlock
    int fIDof = firstDof[iBlock];
    // Dimention of the iBlock column, not including the diagonal block
    int ldl = blHeight[iBlock] - blWeight[iBlock];
    // Degree of freedom corresponding to the top of iBlock
    iTop = fIDof - ldl;
    GenStackFSFullMatrix<Scalar> diagBlock(blWeight[iBlock], blWeight[iBlock], 
       blHeight[iBlock], skyA+blTop[iBlock]+ldl);
 #ifdef _OPENMP
 #pragma omp barrier
 #endif
    // Build the row of L corresponding to this block
   if(iBlock == myFirstBlock)
    {
      for(i = 0; i < ldl; ++i)
        for(j = 0; j <  blWeight[iBlock]; ++j) 
          ltmp[j*ldl+i] = 
            invDiag[iTop+i] * skyA[blTop[iBlock]+j*blHeight[iBlock]+i];

      // Update the diagonal block
      if(ldl >0)
        Tgemm('T', 'N', blWeight[iBlock], blWeight[iBlock], ldl,
          -1.0, ltmp, ldl, skyA+blTop[iBlock],
           blHeight[iBlock], 1.0, skyA+blTop[iBlock]+ldl, blHeight[iBlock]);
      // factor the diagonal block
      nzem += diagBlock.symLuFactor(perm+fIDof, tol, origDiag+fIDof);

      for(i = 0; i < blWeight[iBlock]; ++i) {
#ifdef USE_PERM
        invDiag[firstDof[iBlock] + perm[fIDof+i]] = diagBlock[i][i];
#else
        invDiag[firstDof[iBlock] + i] = diagBlock[i][i];
#endif
      }
      myFirstBlock += numThreads;
    }

   #ifdef _OPENMP
   #pragma omp barrier
   #endif
   // Now update the iBlock row
   for(int jBlock = myFirstBlock; 
		   jBlock <= lastCol[iBlock]; jBlock+=numThreads) {
     // jTop is the top degree of freedom of the j supernode
     jTop = firstDof[jBlock] - blHeight[jBlock] + blWeight[jBlock];

     if(jTop > firstDof[iBlock] + blWeight[iBlock]) continue;
     int iOffset, jOffset;
     int commonHeight;

     if(iTop >= jTop) {
       commonHeight = ldl;
       iOffset = 0;
       jOffset = iTop-jTop;
     } else {
       commonHeight = firstDof[iBlock]-jTop;
       jOffset = 0;
       iOffset = jTop-iTop;
     }
    
     // Final update of the aBlock
     if(commonHeight > 0)
      Tgemm('T', 'N', blWeight[iBlock], blWeight[jBlock], commonHeight,
         -1.0, ltmp+iOffset, ldl, skyA+blTop[jBlock]+jOffset,
          blHeight[jBlock], 1.0, skyA+blTop[jBlock]+
               firstDof[iBlock]-jTop, blHeight[jBlock]);
      if(jTop <= firstDof[iBlock])
        diagBlock.Lm1Mult(skyA+blTop[jBlock]+firstDof[iBlock]-jTop,
                     blWeight[jBlock], blHeight[jBlock], perm+fIDof);
    }
//   double *lxx = ltmp;
//   ltmp = ltmp2;
//   ltmp2 = lxx;
  }
 }
 delete [] origDiag;
 delete [] ltmp;

// fprintf(stderr, "Found %d singularities\n", nzem);
#endif
}


#ifdef sgi
template<class Scalar>
void
GenBlockSky<Scalar>::pFactor(int myid, int numThreads, 
    barrier_t *b, Scalar *ltmp, Scalar *invDiag, Scalar *origDiag)
{
 int i,j;

 int iTop, jTop;

 if(myid == 0) {
   nzem = 0;
   for(i = 0; i < neq; ++i)
     origDiag[i] = skyA[dlp[i]];
 }
 int iBlock;
 //fprintf(stderr, "I am %d out of %d\n", myid, numThreads);
 int myFirstBlock = myid;
 for(iBlock = 0; iBlock < nBlocks; ++iBlock) {
   // First dof of the iBlock
   int fIDof = firstDof[iBlock];
   // Dimention of the iBlock column, not including the diagonal block
   int ldl = blHeight[iBlock] - blWeight[iBlock];
   // Degree of freedom corresponding to the top of iBlock
   iTop = fIDof - ldl;
   GenStackFSFullMatrix<Scalar> diagBlock(blWeight[iBlock], blWeight[iBlock],
       blHeight[iBlock], skyA+blTop[iBlock]+ldl);
   // Barrier
   barrier(b, numThreads);
   // Build the row of L corresponding to this block
   if(iBlock == myFirstBlock)
    {
      for(i = 0; i < ldl; ++i)
        for(j = 0; j <  blWeight[iBlock]; ++j)
          ltmp[j*ldl+i] =
            invDiag[iTop+i] * skyA[blTop[iBlock]+j*blHeight[iBlock]+i];

      // Update the diagonal block
      if(ldl >0)
        Tgemm('T', 'N', blWeight[iBlock], blWeight[iBlock], ldl,
          -1.0, ltmp, ldl, skyA+blTop[iBlock],
           blHeight[iBlock], 1.0, skyA+blTop[iBlock]+ldl, blHeight[iBlock]);
      // factor the diagonal block
      nzem += diagBlock.symLuFactor(perm+fIDof, tol, origDiag+fIDof);

      for(i = 0; i < blWeight[iBlock]; ++i) {
#ifdef USE_PERM
        invDiag[firstDof[iBlock] + perm[fIDof+i]] = diagBlock[i][i];
#else
        invDiag[firstDof[iBlock] + i] = diagBlock[i][i];
#endif
      }
      myFirstBlock += numThreads;
    }

   barrier(b, numThreads);
   // Now update the iBlock row
   for(int jBlock = myFirstBlock;
                   jBlock <= lastCol[iBlock]; jBlock+=numThreads) {
     // jTop is the top degree of freedom of the j supernode
     jTop = firstDof[jBlock] - blHeight[jBlock] + blWeight[jBlock];

     if(jTop > firstDof[iBlock] + blWeight[iBlock]) continue;
     int iOffset, jOffset;
     int commonHeight;

     if(iTop >= jTop) {
       commonHeight = ldl;
       iOffset = 0;
       jOffset = iTop-jTop;
     } else {
       commonHeight = firstDof[iBlock]-jTop;
       jOffset = 0;
       iOffset = jTop-iTop;
     }

     // Final update of the aBlock
     if(commonHeight > 0)
      Tgemm('T', 'N', blWeight[iBlock], blWeight[jBlock], commonHeight,
         -1.0, ltmp+iOffset, ldl, skyA+blTop[jBlock]+jOffset,
          blHeight[jBlock], 1.0, skyA+blTop[jBlock]+
               firstDof[iBlock]-jTop, blHeight[jBlock]);
      if(jTop <= firstDof[iBlock])
        diagBlock.Lm1Mult(skyA+blTop[jBlock]+firstDof[iBlock]-jTop,
                     blWeight[jBlock], blHeight[jBlock], perm+fIDof);
    }
  }

}
#endif

template<class Scalar>
void
GenBlockSky<Scalar>::add(const FullSquareMatrix &kel, const int *dofs)
{
 // Construct stiffness matrix K (skyA)

 int i, j, ri, rj;
 int kndof = kel.dim();                  // Element stiffness dimension
 for(i = 0; i < kndof; ++i) {            // Loop over rows.
   if((ri = rowColNum[dofs[i]]) == -1) continue;// Skip constrained dofs
   for(j = 0; j < kndof; ++j) {          // Loop over columns.
     if(dofs[i] > dofs[j]) continue;     // Work with upper symmetric half.
     if((rj = rowColNum[dofs[j]]) == -1) continue; // Skip constrained dofs
     skyA[dlp[rj] - rj + ri ] += kel[i][j];
   }
 }
}


template<class Scalar>
void
GenBlockSky<Scalar>::add(int row_dof, int col_dof, Scalar s) {
 // Construct stiffness matrix K (skyA)

 int /*i, j,*/ ri, rj;
 if(row_dof > col_dof) {
   std::cerr << "WARNING: row_dof > col_dof in GenBlockSky<Scalar>::add() \n";
   return;  // Work with upper symmetric half.
 }

 if((ri = rowColNum[row_dof]) == -1) return;// Skip constrained dofs
 if((rj = rowColNum[col_dof]) == -1) return; // Skip constrained dofs
     skyA[dlp[rj] - rj + ri ] += s;
}



template<class Scalar>
void
GenBlockSky<Scalar>::add(const AssembledFullM &kel, const int *dofs)
{
 int i, j, ri, rj;
 int kndof = kel.dim();                   // Element stiffness dimension
 for(i = 0; i < kndof; ++i) {             // Loop over rows.
    if((ri = dofs[i]) == -1) continue;    // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {          // Loop over columns.
       if(dofs[i] > dofs[j]) continue;    // Work with upper symmetric half.
       if((rj = dofs[j]) == -1) continue; // Skip constrained dofs
       skyA[dlp[rj] - rj + ri ] += kel[i][j];
    }
 }
}

template<class Scalar>
void
GenBlockSky<Scalar>::add(const GenAssembledFullM<complex<double> > &kel, const int *dofs)
{
 int i, j, ri, rj;
 int kndof = kel.dim();                     // Element stiffness dimension
 for(i = 0; i < kndof; ++i) {               // Loop over rows.
    if((ri = dofs[i]) == -1) continue;      // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {            // Loop over columns.
       if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.
       if( (rj = dofs[j]) == -1 ) continue; // Skip constrained dofs
       skyA[dlp[rj] - rj + ri ] += kel[i][j];
    }
 }
}

template<class Scalar>
void
GenBlockSky<Scalar>::add(const FullM &knd, int fRow, int fCol)
{
  int nrow = knd.numRow(); // number of rows
  int ncol = knd.numCol(); // number of columns

  int iCol, iRow;
  for(iCol = 0; iCol < ncol; ++iCol) {
    int fk = dlp[fCol+iCol] + fRow - fCol - iCol ;
    int rowStop = (nrow < fCol+iCol-fRow+1) ? nrow : fCol+iCol-fRow+1;
    for(iRow = 0; iRow < rowStop; ++iRow)
      skyA[fk+iRow] += knd[iRow][iCol];
  }
}

// Assembly from a Boeing Sparse data structure in fortran
template<class Scalar>
void
GenBlockSky<Scalar>::addBoeing(int nlines, const int *Kai, const int *Kaj,
                               const double *nz, const int *map, Scalar multiplier)
{
 int i, j;
 for(i = 0; i < nlines; ++i) {
   int dofI = rowColNum[map[i]];
   if(dofI < 0) continue;
   for(j = Kai[i]; j < Kai[i+1]; ++j) {
     int dofJ = rowColNum[map[Kaj[j-1]-1]];
     if(dofJ < 0) continue;
     if(dofI < dofJ)
        skyA[dlp[dofJ] - dofJ + dofI] += (nz[j-1]*multiplier);
     else
        skyA[dlp[dofI] - dofI + dofJ] += (nz[j-1]*multiplier);
   }
 }
}


template<class Scalar>
void
GenBlockSky<Scalar>::solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution)
{
   solution = rhs;

   reSolve(solution.data());
}

template<class Scalar>
void
GenBlockSky<Scalar>::solve(const Scalar *rhs, Scalar *solution)
{
   int i;
   for(i=0; i<dim(); ++i)
     solution[i] = rhs[i];

   reSolve(solution);
}


template<class Scalar>
void
GenBlockSky<Scalar>::reSolve(Scalar *rhs)
{
 int iBlock, i;
 // Forward substitution
 for(iBlock = 0; iBlock < nBlocks; ++iBlock) {
   int ldl = blHeight[iBlock] - blWeight[iBlock];
   GenStackFSFullMatrix<Scalar> diagBlock(blWeight[iBlock], blWeight[iBlock], 
                               blHeight[iBlock], skyA+blTop[iBlock]+ldl);
   Tgemm('T', 'N', blWeight[iBlock], 1, ldl,
        -1.0, skyA+blTop[iBlock], blHeight[iBlock], rhs+firstDof[iBlock]-ldl,
         neq, 1.0, rhs+firstDof[iBlock], neq);   
   diagBlock.Um1TMult(rhs+firstDof[iBlock], 1, neq, perm+firstDof[iBlock]);
 }
 // scaling
 for(i = 0; i < neq; ++i)
  if(invDiag[i] != 0.0)
   rhs[i] /= invDiag[i];
  else rhs[i] = 0.0;

 // backward substitute
 for(iBlock = nBlocks; iBlock--; ) {
   int ldl = blHeight[iBlock] - blWeight[iBlock];
   GenStackFSFullMatrix<Scalar> diagBlock(blWeight[iBlock], blWeight[iBlock], 
     blHeight[iBlock], skyA+blTop[iBlock]+ldl);
   diagBlock.Um1Mult(rhs+firstDof[iBlock], 1, neq, perm+firstDof[iBlock]);
   Tgemm( 'N', 'N', ldl, 1, blWeight[iBlock], -1.0,
         skyA+blTop[iBlock], blHeight[iBlock], rhs+firstDof[iBlock], neq,
         1.0, rhs+firstDof[iBlock]-ldl, neq);
 }

}

template<class Scalar>
void
GenBlockSky<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
 reSolve(rhs.data());
}


template<class Scalar>
void
GenBlockSky<Scalar>::reSolve(int numRHS, Scalar **rhs)
{
 int iBlock, i;
 // Forward substitution
 for(iBlock = 0; iBlock < nBlocks; ++iBlock) {
   int ldl = blHeight[iBlock] - blWeight[iBlock];
   GenStackFSFullMatrix<Scalar> diagBlock(blWeight[iBlock], blWeight[iBlock],
     blHeight[iBlock], skyA+blTop[iBlock]+ldl);
   Tgemm('T', 'N', blWeight[iBlock], numRHS, ldl,
        -1.0, skyA+blTop[iBlock], blHeight[iBlock], rhs[0]+firstDof[iBlock]-ldl,
         neq, 1.0, rhs[0]+firstDof[iBlock], neq);
   diagBlock.Um1TMult(rhs[0]+firstDof[iBlock], numRHS, neq, perm+firstDof[iBlock]);
 }
 // scaling
 int n;
 for(n=0; n<numRHS; ++n)
   for(i = 0; i < neq; ++i)
    if(invDiag[i] != 0.0)
     rhs[n][i] /= invDiag[i];
    else rhs[n][i] = 0.0;

 // backward substitute
 for(iBlock = nBlocks; iBlock--; ) {
   int ldl = blHeight[iBlock] - blWeight[iBlock];
   GenStackFSFullMatrix<Scalar> diagBlock(blWeight[iBlock], blWeight[iBlock],
     blHeight[iBlock], skyA+blTop[iBlock]+ldl);
   diagBlock.Um1Mult(rhs[0]+firstDof[iBlock], numRHS, neq, perm+firstDof[iBlock]);
   Tgemm( 'N', 'N', ldl, numRHS, blWeight[iBlock], -1.0,
         skyA+blTop[iBlock], blHeight[iBlock], rhs[0]+firstDof[iBlock], neq,
         1.0, rhs[0]+firstDof[iBlock]-ldl, neq);
 }

}

template<class Scalar>
void
GenBlockSky<Scalar>::print(FILE *f)
{
 for(int i =0; i < blTop[nBlocks]; ++i)
   fprintf(f, "K[%d] = %e\n",i, ScalarTypes::Real(skyA[i]));
}

template<class Scalar>
void
GenBlockSky<Scalar>::zeroAll()
{
  int i;
#ifdef DISTRIBUTED
  for(i = 0; i < blTop[nBlocks]; ++i) 
    skyA[i] = 0.0;
#else
  for(i = 0; i <= dlp[neq-1]; ++i)
    skyA[i] = 0.0;
#endif
}

template<class Scalar>
void
GenBlockSky<Scalar>::clean_up()
{
 if(skyA) {
   delete [] skyA;
   skyA = 0;
 }
}

template<class Scalar>
void
GenBlockSky<Scalar>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
 communicator->globalSum(blTop[nBlocks], skyA);
#endif

/*
 FILE *file;
 if(syscom->myID() == 0)
   file = fopen("K.1","w");
 else if(syscom->myID() == 1)
   file = fopen("K.2","w");
 else if(syscom->myID() == 2)
   file = fopen("K.3","w");
 else
   file = fopen("K.4","w");
 print(file);
*/
}

template<class Scalar>
void
GenBlockSky<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  // int this function dof is in unconstrained dof numbering
  if(dof < 0) return;
  int ri = rowColNum[dof];
  if(ri < 0 ) return;  // Skip constrained dofs
  skyA[dlp[ri]] += dmass;
}

