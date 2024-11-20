#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/SubDomain.h>
#include <Feti.d/OrthoSet.h>
#include <Utils.d/linkfc.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/matrix.h>
#include <Math.d/IntFullM.h>
#include <Math.d/SymFullMatrix.h>
#include <Solvers.d/Rbm.h>
#include "CoarseSet.h"

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
void
GenCoarseSet<Scalar>::clean_up()
{
 if(neighbGs) {
   delete [] neighbGs;
   neighbGs=0;
 }
 if(leadingDimGs) {
   delete [] leadingDimGs;
   leadingDimGs = 0;
 }
}

template<class Scalar>
int
GenCoarseSet<Scalar>::offset(int subnum)
{
 int i;
 for(i=0; i < numNeighb; ++i)
   if(subnum == neighbs[i]) return subOffset[i];
 return 0;
}

template<class Scalar>
int
GenCoarseSet<Scalar>::edgeSize(int i)
{
 return subOffset[i+1] - subOffset[i];
}

// alpha (gSize x 1) = C' (gSize x numGs) * vec (numGs x 1)

template<class Scalar>
void
GenCoarseSet<Scalar>::getGtMult(const Scalar *vec, Scalar *alpha) const
{
 if(gSize > 0)
 Tgemv('T',gSize, numGs, 1.0, locGs, gSize, vec, 1, 0.0, alpha, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::getGtQMult(const Scalar *vec, Scalar *alpha) const
{
 if(gSize > 0)
 Tgemv('T',gSize, numGs, 1.0, locQGs, gSize, vec, 1, 0.0, alpha, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subAlphaGtQ(int iSub, const Scalar *vec, Scalar *alpha) const
{
 if(gSize > 0)
 Tgemv('T',gSize, neighbNumRBMs[iSub], -1.0, neighbQGs[iSub], gSize,
       vec, 1, 1.0, alpha, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::getGtQMult(int iSub, const Scalar *vec, Scalar *alpha) const
{
 if(gSize > 0)
 Tgemv('T',gSize, neighbNumRBMs[iSub], 1.0, neighbQGs[iSub], gSize, 
       vec, 1, 1.0, alpha, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::getGtFMult(const Scalar *vec, Scalar *alpha) const
{
 if(gSize > 0)
 Tgemv('T',gSize, numGs, 1.0, locFGs, gSize, vec, 1, 0.0, alpha, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::getGtFMult(int iSub, const Scalar *vec, Scalar *alpha) const
{
 if(gSize > 0)
 Tgemv('T',gSize, neighbNumRBMs[iSub], +1.0, neighbQGs[iSub], gSize,
       vec, 1, 1.0, alpha, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subLocAlphaG(Scalar *vec, Scalar *alpha)
{
 if(gSize > 0)
 Tgemv('N', gSize, numGs, -1.0, locGs, gSize,alpha, 1, 1.0, vec, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::addLocAlphaG(Scalar *vec, Scalar *alpha)
{
 if(gSize > 0)
 Tgemv('N', gSize, numGs, 1.0, locGs, gSize,alpha, 1, 1.0, vec, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subLocNuC(Scalar *vec, Scalar *nu)
{
 int i;
 for (i=0; i< numBCs; i++) {
   int index = (*locBCs)[0][i];
   vec[index] += nu[i];
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subLocNuFC(Scalar *vec, Scalar *nu)
{
 // Sign convention ????
 if(gSize > 0)
 Tgemv('N', gSize, numBCs, +1.0, locFBCs, gSize,
	nu, 1, 1.0, vec, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subLocAlphaQG(Scalar *vec, Scalar *alpha)
{
 if(gSize > 0)
 Tgemv('N', gSize,numGs, -1.0, locQGs, gSize,
	alpha, 1, 1.0, vec, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subLocAlphaFG(Scalar *vec, Scalar *alpha)
{
 if(gSize > 0)
 Tgemv('N', gSize,numGs, -1.0, locFGs, gSize,
       alpha, 1, 1.0, vec, 1);
}


template<class Scalar>
void
GenCoarseSet<Scalar>::subAlphaNeighbG(int sub, Scalar *vec, Scalar *neighbG,
                                      int lda, Scalar *alpha)
{
 if(lda > 0)
 Tgemv('N', subSize(sub), neighbNumRBMs[sub], 1.0, neighbG,
       lda, alpha, 1, 1.0, vec + subOffset[sub], 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::addAlphaNeighbG(int sub, Scalar *vec, Scalar *neighbG,
                                      int lda, Scalar *alpha)
{
if(lda > 0)
 Tgemv('N', subSize(sub), neighbNumRBMs[sub], -1.0, neighbG,
       lda, alpha, 1, 1.0, vec + subOffset[sub], 1);
}


template<class Scalar>
void
GenCoarseSet<Scalar>::getNeighbGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd,
                                  FSCommPattern<Scalar> *rbmPat)
{
 neighbGs     = new Scalar * [numNeighb];
 leadingDimGs = new int[numNeighb];

 int iSub;
 for(iSub = 0; iSub < numNeighb; ++iSub) {
#ifdef DISTRIBUTED
   leadingDimGs[iSub] = sd->getSComm()->sharedDOFsCount(iSub);
   neighbGs[iSub] = new Scalar[neighbNumRBMs[iSub]*leadingDimGs[iSub]]; // PJSA
   sd->recvInterfRBMs(iSub, neighbNumRBMs[iSub], neighbGs[iSub], rbmPat);  // PJSA
#else
   neighbGs[iSub]     = csets[neighbs[iSub]].subG(myNum);
   leadingDimGs[iSub] = csets[neighbs[iSub]].gSize;
#endif
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::getNeighbQGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd, int projType, 
                                   FSCommPattern<Scalar> *rbmPat, GenSolver<Scalar> *solver)
{
 int iSub, j, len; 

 neighbQGs    = new Scalar *[numNeighb];
 neighbGs     = new Scalar *[numNeighb];
 leadingDimGs = new int[numNeighb];

 int totNumQG = 0;
 for(iSub = 0; iSub < numNeighb; ++iSub)
   totNumQG += neighbNumRBMs[iSub];

 pool1 = new Scalar[totNumQG*gSize];

 for(iSub = 0; iSub < numNeighb; ++iSub) {
#ifdef DISTRIBUTED
   len = sd->getSComm()->sharedDOFsCount(iSub);
   Scalar *srcG = new Scalar[neighbNumRBMs[iSub]*len]; // PJSA
   sd->recvInterfRBMs(iSub, neighbNumRBMs[iSub], srcG, rbmPat);  // PJSA
#else
   Scalar *srcG = csets[neighbs[iSub]].subG(myNum);
   len = csets[neighbs[iSub]].gSize;
#endif
   neighbGs[iSub] = srcG;
   leadingDimGs[iSub] = len;
   neighbQGs[iSub] = pool1; pool1 += gSize*neighbNumRBMs[iSub];
   int i;
   for(i = 0; i < neighbNumRBMs[iSub]; ++i) {
     for(j = 0; j < gSize; ++j)
       neighbQGs[iSub][i*gSize + j] = 0.0;
     for(j = subOffset[iSub]; j <  subOffset[iSub+1] && j <gSize; ++j)
       neighbQGs[iSub][i*gSize + j] = -srcG[i*len+j-subOffset[iSub]];

     // KHP projType = 1 -> basic, Q=I
     //     projType = 2 -> Q=preconditioner 
     //                     (either Kbb, or Kbb - Kib^T Kii^-1 Kib)
     //     projType = 3 -> Q = kc
     //     projType = 4 -> Q = Diag(Kbb)
     if(!solver) {
      if(projType != 4) {
        sd->multKbb(neighbQGs[iSub]+i*gSize, neighbQGs[iSub]+i*gSize);
      } else {
        sd->multDiagKbb(neighbQGs[iSub]+i*gSize, neighbQGs[iSub]+i*gSize);
      }
    }
   }

  // multiple rhs version of multFi
  if(solver)
    sd->multMFi(solver, neighbQGs[iSub], neighbQGs[iSub],neighbNumRBMs[iSub]);
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subAlphaNeighbQG(int sub, Scalar *vec, Scalar *alpha)
{
 if(gSize > 0)
 Tgemv('N', gSize, neighbNumRBMs[sub], -1.0, neighbQGs[sub],
       gSize, alpha, 1, 1.0, vec, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subAlphaNeighbFG(int sub, Scalar *vec, Scalar *alpha)
{
 if(gSize > 0)
 Tgemv('N', gSize, neighbNumRBMs[sub], -1.0, neighbQGs[sub],
       gSize, alpha, 1, 1.0, vec, 1);
}

// beta = C'F * vec

template<class Scalar>
void
GenCoarseSet<Scalar>::getCtFMult(const Scalar *vec, Scalar *beta) const
{
 if(numBCs > 0 && gSize > 0)
 Tgemv('T', gSize, numBCs, -1, locFBCs, gSize, vec, 1, 0, beta, 1);
}

// beta (gSize x 1) = C' (gSize x numBCs) * vec (numBCs x 1)  

template<class Scalar>
void
GenCoarseSet<Scalar>::getCtMult(const Scalar *vec, Scalar *beta) const
{
 int i;
 for (i=0; i< numBCs; i++) {
   int index = (*locBCs)[0][i];
   beta[i] = -vec[index];
 }

}

template<class Scalar>
void
GenCoarseSet<Scalar>::getCtFMult(int iSub, const Scalar *vec, Scalar *beta) const
{
 int nCright = neighbCSubIndex[iSub+1] - neighbCSubIndex[iSub];
 if(nCright > 0 && gSize > 0)
 Tgemv('T', gSize, nCright, +1.0, neighbFBCs2[iSub], gSize,
       vec, 1, 1.0, beta + neighbCs[neighbCSubIndex[iSub]][1], 1);
}


////////////////////////////////////////////////////////////////////////////////

/***************************************************************
 * buildBNeighbCs                                              *
 *   build the neighbCSubIndex and neighbCs vectors            *
 *                                                             *
 ***************************************************************/
template<class Scalar>
void
GenCoarseSet<Scalar>::buildBNeighbCs(GenCoarseSet<Scalar> *csets)
{
 int iSub, iC;
 neighbCSubIndex = new int[numNeighb+1];
 neighbCSubIndex[0] = 0; 
 for(iSub = 0; iSub < numNeighb; ++iSub) {
#ifdef DISTRIBUTED
   // KENDALL 9/1/99
   // NEED A FIX FOR THIS!
   // neighbCSubIndex[iSub+1] = neighbCSubIndex[iSub] + ???;
#else
   neighbCSubIndex[iSub+1] = neighbCSubIndex[iSub] +
                             csets[neighbs[iSub]].getNumCs(myNum);    
#endif
 }
 neighbCs = new int[neighbCSubIndex[numNeighb]][3];

 for(iSub = 0; iSub < numNeighb; ++iSub) {
   csets[neighbs[iSub]].getCOffset(myNum, neighbCs+neighbCSubIndex[iSub]);
   for(iC = neighbCSubIndex[iSub]; iC < neighbCSubIndex[iSub+1]; ++iC) {
      neighbCs[iC][0] += subOffset[iSub];
      neighbCs[iC][2] = iSub;
   }
 }

}

template<class Scalar>
int
GenCoarseSet<Scalar>::getNumCs(int subNum)
{
 int r = 0;
 int iC;
 for(iC = 0; iC < numBCs ; ++iC)
   if(neighbs[(*locBCs)[3][iC]] == subNum)
     r++;
 return r;
}

template<class Scalar>
void
GenCoarseSet<Scalar>::getCOffset(int subNum, int (*offsets)[3]) 
{
 int index = 0;
 int iC;
 for(iC = 0; iC < numBCs ; ++iC)
   if(neighbs[(*locBCs)[3][iC]] == subNum) {
     offsets[index][0] = (*locBCs)[0][iC] - subOffset[(*locBCs)[3][iC]];
     offsets[index][1] = iC;
     index++;
   }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subNuNeighbC(Scalar *vec, Scalar *nu, int *offsets)
{
 int iSub, iC;
 // ML WARNING !!! CAREFUL ABOUT THE SIGN HERE! NOT SAME CONVENTION AS G
 for(iSub = 0; iSub < numNeighb; ++iSub) 
   for(iC = neighbCSubIndex[iSub]; iC < neighbCSubIndex[iSub+1]; ++iC) {
      vec[neighbCs[iC][0]] -= nu[offsets[iSub] + neighbCs[iC][1]];
   }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::subNuNeighbFC(Scalar *vec, Scalar *nu, int *offsets)
{
 int iSub, iC;
 Scalar *collectedNu = (Scalar *)
       dbg_alloca(sizeof(Scalar)*neighbCSubIndex[numNeighb]);
 for(iSub = 0; iSub < numNeighb; ++iSub)
   for(iC = neighbCSubIndex[iSub]; iC < neighbCSubIndex[iSub+1]; ++iC) 
     collectedNu[iC] = nu[offsets[iSub] + neighbCs[iC][1]];

// use the fact that all FCs were allocated consecutively
 Tgemv('N', gSize, neighbCSubIndex[numNeighb], -1.0, neighbFBCs2[0],
       gSize, collectedNu, 1, 1.0, vec, 1);
}

template<class Scalar>
void
GenCoarseSet<Scalar>::computeNeighbFBCs(GenSubDomain<Scalar> *sd, GenSolver<Scalar> *solver)
{
 int iSub, j;
 neighbFBCs2 = new Scalar *[numNeighb];
 for(iSub = 0; iSub < numNeighb; ++iSub)
   neighbFBCs2[iSub]=0;
 int totNumFBC = neighbCSubIndex[numNeighb];

 pool2 = new Scalar[totNumFBC*gSize];

 for(iSub = 0; iSub < numNeighb; ++iSub) {

   neighbFBCs2[iSub] = pool2; 
   int neighbNumCs = neighbCSubIndex[iSub+1]-neighbCSubIndex[iSub];
   pool2 += gSize*neighbNumCs;

   int i;
   for(i = 0; i < neighbNumCs; ++i) {
     for(j = 0; j < gSize; ++j)
        neighbFBCs2[iSub][i*gSize + j] = 0.0;

     int iC = i+ neighbCSubIndex[iSub];
     neighbFBCs2[iSub][i*gSize + neighbCs[iC][0]] = +1;
   }
   // muliple rhs version of multFi
   sd->multMFi(solver, neighbFBCs2[iSub], neighbFBCs2[iSub],neighbNumCs);
 }

}

template<class Scalar>
void
GenCoarseSet<Scalar>::addCtFCContrib(GenFullM<Scalar> *CtFC, int iSub, EqNumberer* eqNumber)
{
 Scalar * FC;
 int Ccol;
 int nCright;
 Scalar cornerSign;
 if(iSub == 0) { // then it is me
   nCright = numBCs;
   if(nCright <= 0) return;
   FC = locFBCs;
   Ccol = eqNumber->firstdof(myNum);
   cornerSign = 1;
 }else {
   iSub -= 1;
   nCright = neighbCSubIndex[iSub+1] - neighbCSubIndex[iSub];
   if(nCright <= 0) return;
   Ccol = eqNumber->firstdof(neighbs[iSub])+neighbCs[neighbCSubIndex[iSub]][1];
   FC = neighbFBCs2[iSub];
   cornerSign = -1;
 }

 // Find the maximum matrix allocation size
 int maxLeft = numBCs;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   if((neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub]) > maxLeft)
     maxLeft = (neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub]) ;
 }
 Scalar *stackMem = (Scalar *) dbg_alloca(sizeof(Scalar)*(maxLeft*nCright));

 if(numBCs > 0 && nCright > 0) {
     GenStackFullM<Scalar> ctfc(numBCs, nCright,stackMem);
     int i,j;
     for(i = 0; i < numBCs; ++i)
       for(j = 0; j < nCright; ++j) {
           ctfc[i][j] = cornerSign * FC[j*gSize + (*locBCs)[0][i] ];
       }
     CtFC->add(ctfc,eqNumber->firstdof(myNum),Ccol);
 }

 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int leftCcol;
   int nCleft = neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub];
   if(nCleft <= 0) continue;
   leftCcol = eqNumber->firstdof(neighbs[jSub]) +
           neighbCs[neighbCSubIndex[jSub]][1];
   int (*leftCs)[3] = neighbCs + neighbCSubIndex[jSub];

   GenStackFullM<Scalar> ctfc(nCleft, nCright,stackMem);
       int i,j;
       for(i = 0; i < nCleft; ++i)
         for(j = 0; j < nCright; ++j)
             ctfc[i][j] = -cornerSign * FC[j*gSize + leftCs[i][0] ];
    CtFC->add(ctfc,leftCcol,Ccol);
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::addCtGContrib(GenFullM<Scalar> *CtG, EqNumberer *cNums,
                                    EqNumberer *gNums)
{
 int i,j;

 // Find the maximum matrix allocation size
 int maxGs = numGs;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   if(gNums->weight(neighbs[jSub]) > maxGs)
     maxGs = gNums->weight(neighbs[jSub]);
 }
 Scalar *stackMem = (Scalar *) dbg_alloca(sizeof(Scalar)*(numBCs*maxGs));

 GenStackFullM<Scalar> ctg(numBCs, numGs,stackMem);
 ctg.zero();
 int myRow = cNums->firstdof(myNum);

 if (numGs>0 && numBCs >0)  {
   for (i=0; i< numGs; i++)
     for (j=0; j< numBCs; j++) {
        int GoffSet = i * gSize;
        ctg[j][i] = -locGs[GoffSet + (*locBCs)[0][j]];
       }
 }
 CtG->add(ctg,myRow,gNums->firstdof(myNum));
 // Now loop on the neighbors
 int iC = 0;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int neighbNumGs = gNums->weight(neighbs[jSub]);
   GenStackFullM<Scalar> ctg(numBCs, neighbNumGs,stackMem);
   ctg.zero();
   while(iC < numBCs && (*locBCs)[3][iC] == jSub) {
     int i;
     for(i =0; i < neighbNumGs; ++i)
       ctg[iC][i] = +neighbGs[jSub][i*leadingDimGs[jSub]+(*locBCs)[0][iC]-subOffset[jSub]];
     iC = iC+1;
   }
  CtG->add(ctg,myRow, gNums->firstdof(neighbs[jSub]));
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::addGtFCContrib(GenFullM<Scalar> *GtFC,int iSub, EqNumberer *cNums,
                                     EqNumberer *gNums)
{
 Scalar * QC;
 int Ccol;
 int nCright;
 Scalar cornerSign;

 if(iSub == 0) { // then it is me
   QC = locFBCs;
   Ccol = cNums->firstdof(myNum);
   nCright = numBCs;
   cornerSign = 1.0;
   iSub -= 1;
 } else {
   iSub -= 1;
   cornerSign = -1.0;
   QC = neighbFBCs2[iSub];
   nCright = neighbCSubIndex[iSub+1] - neighbCSubIndex[iSub];
   if(nCright > 0)
      Ccol = cNums->firstdof(neighbs[iSub]) +neighbCs[neighbCSubIndex[iSub]][1];
   else
      return;
 }

 // Find the maximum matrix allocation size
 int maxGs = numGs;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   if(neighbNumRBMs[jSub] > maxGs)
     maxGs = neighbNumRBMs[jSub];
 }
 Scalar *stackMem = (Scalar *) dbg_alloca(sizeof(Scalar)*(nCright*maxGs));

    if(numGs > 0 && nCright > 0) {
     GenStackFullM<Scalar> gtqc(numGs, nCright,stackMem);
     Tgemm('T','N', nCright, numGs, gSize, -cornerSign,
            QC, gSize, locGs, gSize, 0.0, gtqc.data(), nCright);
     GtFC->add(gtqc, gNums->firstdof(myNum), Ccol);
   }
 
 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
     int leftGcol = gNums->firstdof(neighbs[jSub]);
     if(neighbNumRBMs[jSub] > 0 && nCright > 0) {
       GenStackFullM<Scalar> gtqc(neighbNumRBMs[jSub], nCright,stackMem);
       Tgemm('T', 'N', nCright, neighbNumRBMs[jSub], subSize(jSub), +cornerSign,
             QC+subOffset[jSub], gSize, neighbGs[jSub], leadingDimGs[jSub], 0.0, gtqc.data(), nCright);
       GtFC->add(gtqc, leftGcol, Ccol);
     }
 }
}

/*****************************************************************
 * This subroutine adds all the contributions of (GC)t Fme       *
 *   (GC)t Qme (G_i C_i)                                         *
 * in input iSub is 0 if it is myself or 1+neighbor              *
 *                                                               *
 *****************************************************************/

template<class Scalar>
void
GenCoarseSet<Scalar>::addQContrib(GenSparseMatrix<Scalar> *GtQG, int iSub, EqNumberer* eqNumber, 
                                  int hasC)
{
 Scalar * QG;
 Scalar * QC;
 int Gcol, Ccol;
 int nGright, nCright;
 Scalar cornerSign;

 if(iSub == 0) { // then it is me
   QG = locQGs;
   QC = locFBCs;
   Gcol = eqNumber->firstdof(myNum);
   Ccol = Gcol + numGs;
   nGright = numGs;
   nCright = numBCs;
   cornerSign = 1.0;
   iSub -= 1;
 } else {
   iSub -= 1;
   QG = neighbQGs[iSub];
   Gcol = eqNumber->firstdof(neighbs[iSub]);
   nGright = neighbNumRBMs[iSub];
   cornerSign = -1.0;
   if(hasC) {
     QC = neighbFBCs2[iSub];
     nCright = neighbCSubIndex[iSub+1] - neighbCSubIndex[iSub];
     if(nCright > 0)
       Ccol = Gcol + neighbNumRBMs[iSub] + neighbCs[neighbCSubIndex[iSub]][1];
   }
 }
 // Find the maximum matrix allocation size
 int maxLeft  = (hasC && (numBCs > numGs)) ? numBCs : numGs;
 int maxRight = (hasC && (nCright > nGright)) ? nCright : nGright;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
     int nGleft = neighbNumRBMs[jSub];
     if(nGleft > maxLeft) maxLeft = nGleft;
     if(hasC) {
       int nCleft = neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub];
       if(nCleft > maxLeft) maxLeft = nCleft;
     }
 }
 Scalar *stackMem = (Scalar *) dbg_alloca(sizeof(Scalar)*(maxLeft*maxRight));

 // now that we have QG and QC, lets multiply them with all the others Gs and Cs

 if((numGs > 0) && (nGright > 0) && (gSize > 0)) {
    GenStackFullM<Scalar> gtqg(numGs, nGright, stackMem);
    Tgemm('T', 'N', nGright, numGs, gSize, 1.0,
          QG, gSize, locGs, gSize, 0.0, gtqg.data(), nGright);
    GtQG->add(gtqg, eqNumber->firstdof(myNum), Gcol);
 } 

 if(hasC) {
   if(numGs > 0 && nCright > 0) {
     GenStackFullM<Scalar> gtqc(numGs, nCright,stackMem);
      Tgemm('T', 'N', nCright, numGs, gSize, -cornerSign, QC, gSize,
            locGs, gSize, 0.0, gtqc.data(), nCright);
     GtQG->add(gtqc, eqNumber->firstdof(myNum), Ccol);
   }
 
   if(numBCs > 0 && nGright > 0) {
     GenStackFullM<Scalar> ctqg(numBCs, nGright, stackMem);
     int i,j;
     for(i = 0; i < numBCs; ++i)
       for(j = 0; j < nGright; ++j)
           ctqg[i][j] = -QG[j*gSize + (*locBCs)[0][i] ];
     GtQG->add(ctqg,eqNumber->firstdof(myNum)+numGs, Gcol);
   } 

   if(numBCs > 0 && nCright > 0) {
     GenStackFullM<Scalar> ctqc(numBCs, nCright, stackMem);
     int i,j;
     for(i = 0; i < numBCs; ++i)
       for(j = 0; j < nCright; ++j) {
           ctqc[i][j] = cornerSign*QC[j*gSize + (*locBCs)[0][i] ];
       }
     GtQG->add(ctqc,eqNumber->firstdof(myNum)+numGs, Ccol);
   }
 }

 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int leftGcol = eqNumber->firstdof(neighbs[jSub]),
       leftCcol;
   int nGleft = neighbNumRBMs[jSub];

   if(neighbNumRBMs[jSub] > 0 && nGright > 0) {
     GenStackFullM<Scalar> gtqg(nGleft, nGright,stackMem);
     Tgemm('T', 'N', nGright, nGleft, subSize(jSub), -1.0,
           QG+subOffset[jSub], gSize, neighbGs[jSub], leadingDimGs[jSub],
           0.0, gtqg.data(), nGright);
     GtQG->add(gtqg, leftGcol, Gcol);
   }
   if(hasC) {
     int nCleft = neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub];
     if(nCleft > 0)
       leftCcol = leftGcol + neighbNumRBMs[jSub] + 
                    neighbCs[neighbCSubIndex[jSub]][1];

     if(neighbNumRBMs[jSub] > 0 && nCright > 0) {
       GenStackFullM<Scalar> gtqc(neighbNumRBMs[jSub], nCright, stackMem);
        Tgemm('T','N', nCright, neighbNumRBMs[jSub], subSize(jSub), +cornerSign,
                        QC+subOffset[jSub], gSize,
                        neighbGs[jSub], leadingDimGs[jSub],
                        0.0, gtqc.data(), nCright);
       GtQG->add(gtqc, leftGcol, Ccol);
     }

     if(nCleft == 0) continue;

     int (*leftCs)[3] = neighbCs + neighbCSubIndex[jSub];
     if(nCleft > 0 && nGright > 0) {
       GenStackFullM<Scalar> ctqg(nCleft, nGright,stackMem);
       int i,j;
       for(i = 0; i < nCleft; ++i)
         for(j = 0; j < nGright; ++j)
             ctqg[i][j] = +QG[j*gSize + leftCs[i][0] ];
       GtQG->add(ctqg,leftCcol,Gcol);
     }

     if(nCleft > 0 && nCright > 0) {
       GenStackFullM<Scalar> ctqc(nCleft, nCright,stackMem);
       int i,j;
       for(i = 0; i < nCleft; ++i)
         for(j = 0; j < nCright; ++j)
             ctqc[i][j] = -cornerSign*QC[j*gSize + leftCs[i][0] ];
       GtQG->add(ctqc,leftCcol,Ccol);
     }
   }
 }
}

/*****************************************************************
 * This subroutine adds all the contributions of (GC)t Fme       *
 *   (G)t Fme (G_i)                                              *
 * in input iSub is 0 if it is myself or 1+neighbor              *
 *                                                               *
 *****************************************************************/
template<class Scalar>
void
GenCoarseSet<Scalar>::addFContrib(GenSymFullMatrix<Scalar> *GtFG, int iSub, EqNumberer* eqNumber)
{
 Scalar * QG;
 int Gcol;
 int nGright;

 if(iSub == 0) { // then it is me
   QG = locFGs;
   Gcol = eqNumber->firstdof(myNum);
   nGright = numGs;
 } else {
   iSub -= 1;
   QG = neighbQGs[iSub];
   Gcol = eqNumber->firstdof(neighbs[iSub]);
   nGright = neighbNumRBMs[iSub];
 }
 // Find the maximum matrix allocation size
 int maxGs = numGs;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   if(neighbNumRBMs[jSub] > maxGs)
     maxGs = neighbNumRBMs[jSub];
 }

 // now that we have QG and QC, lets multiply them with all the others Gs and Cs

 if(numGs > 0 && nGright > 0) {
    GenFullM<Scalar> gtqg(numGs, nGright);
    Tgemm('T', 'N', nGright, numGs, gSize, 1.0,
          QG, gSize, locGs, gSize,
          0.0, gtqg.data(), nGright);
    GtFG->add(gtqg, eqNumber->firstdof(myNum), Gcol);
 } 

 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   if(neighbNumRBMs[jSub] > 0 && nGright > 0) {
     GenFullM<Scalar> gtqg(neighbNumRBMs[jSub], nGright);
     Tgemm('T', 'N', nGright, neighbNumRBMs[jSub], subSize(jSub), -1.0,
           QG+subOffset[jSub], gSize, neighbGs[jSub], leadingDimGs[jSub],
           0.0, gtqg.data(), nGright);
     GtFG->add(gtqg, eqNumber->firstdof(neighbs[jSub]), Gcol);
   }
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::addGContrib(GenSparseMatrix<Scalar> *coarseMat, int iSub, 
	                          EqNumberer* eqNumber, int gOffset)
{
 Scalar *FG, *G;
 int Gcol, FGcol;
 int nGright;

 if(iSub == 0) { // then it is me
   FG = locFGs;
   G  = locGs;
   FGcol = eqNumber->firstdof(myNum);
   Gcol  = eqNumber->firstdof(myNum+gOffset);
   nGright = numGs;
   iSub -= 1;
 } else {
   iSub -= 1;
   FG = neighbQGs[iSub];
   G  = neighbGs[iSub];
   FGcol = eqNumber->firstdof(neighbs[iSub]);
   Gcol  = eqNumber->firstdof(neighbs[iSub]+gOffset);
   nGright = neighbNumRBMs[iSub];
 }
 // Find the maximum matrix allocation size
 int maxLeft  =  numGs;
 int maxRight =  nGright;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
     int nGleft = neighbNumRBMs[jSub];
     if(nGleft > maxLeft) maxLeft = nGleft;
 }
 Scalar *stackMem = (Scalar *) dbg_alloca(sizeof(Scalar)*(maxLeft*maxRight));

 if((numGs > 0) && (nGright > 0) && (gSize > 0) && (isDynamic==0)) {
    GenStackFullM<Scalar> gtg(numGs, nGright,stackMem);
    Tgemm('T', 'N', nGright, numGs, gSize, 1.0,
          FG, gSize, locGs, gSize, 0.0, gtg.data(), nGright);
    coarseMat->add(gtg, eqNumber->firstdof(myNum), FGcol);
    
    if(iSub == -1) {
      Tgemm('T', 'N', nGright, numGs, gSize, 1.0,
            G, gSize, locGs, gSize, 0.0, gtg.data(), nGright);
      // the gi^Tgj appears in two places
      coarseMat->add(gtg, eqNumber->firstdof(myNum), Gcol);
      coarseMat->add(gtg, eqNumber->firstdof(myNum+gOffset), FGcol);
    }
 } 

 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int leftFGcol = eqNumber->firstdof(neighbs[jSub]);
   int leftGcol = eqNumber->firstdof(neighbs[jSub]+gOffset);
   int nGleft = neighbNumRBMs[jSub];

   if(neighbNumRBMs[jSub] > 0 && nGright > 0) {
     GenStackFullM<Scalar> gtqg(nGleft, nGright,stackMem);
     Tgemm('T', 'N', nGright, nGleft, subSize(jSub), -1.0,
           FG+subOffset[jSub], gSize, neighbGs[jSub], leadingDimGs[jSub],
           0.0, gtqg.data(), nGright);

     coarseMat->add(gtqg, leftFGcol, FGcol);
     if((iSub == -1) && (isDynamic==0)) {
        Tgemm('T', 'N', nGright, nGleft, subSize(jSub), -1.0,
               G+subOffset[jSub], gSize, neighbGs[jSub], leadingDimGs[jSub],
               0.0, gtqg.data(), nGright);
        coarseMat->add(gtqg,leftFGcol, Gcol);
        coarseMat->add(gtqg, leftGcol, FGcol);
     }
   }
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::addCContrib(GenSparseMatrix<Scalar> *coarseMat, int iSub,
                  EqNumberer* eqNumber, int cOffset, int gOffset)
{
 Scalar *FG, *FC;
 int FGcol, FCcol;
 int nGright, nCright;
 Scalar cornerSign;
 int i,j;

 if(iSub == 0) { // then it is me
   FG = locFGs;
   FC = locFBCs;
   FGcol = eqNumber->firstdof(myNum);
   FCcol = eqNumber->firstdof(myNum+cOffset);
   nGright = numGs;
   nCright = numBCs;
   cornerSign = 1.0;
   iSub -= 1;
 } else {
   iSub -= 1;
   FG = neighbQGs[iSub];
   FGcol = eqNumber->firstdof(neighbs[iSub]);
   cornerSign = -1.0;
   nGright = neighbNumRBMs[iSub];
   FC = neighbFBCs2[iSub];
   nCright = neighbCSubIndex[iSub+1] - neighbCSubIndex[iSub];
   int cSubIndex = (neighbCSubIndex[numNeighb] > 0) ?
                        neighbCs[neighbCSubIndex[iSub]][1]  : 0;
   FCcol = eqNumber->firstdof(neighbs[iSub]+cOffset) + cSubIndex; 
 }
 // Find the maximum matrix allocation size
 int maxLeft  =  (numBCs > numGs) ? numBCs : numGs;
 int maxRight =  (nCright > nGright) ? nCright : nGright;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
     int nGleft = neighbNumRBMs[jSub];
     if(nGleft > maxLeft) maxLeft = nGleft;
     int nCleft = neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub];
     if(nCleft > maxLeft) maxLeft = nCleft;
 }
 Scalar *stackMem = (Scalar *) dbg_alloca(sizeof(Scalar)*(maxLeft*maxRight));

 if(numGs > 0 && nCright > 0) {
   // This can be cheaper by 'picking'
   GenStackFullM<Scalar> gtFc(numGs, nCright,stackMem);
   //MLX SIGN adjusted
   Tgemm('T','N', nCright, numGs, gSize, cornerSign,
                      FC, gSize,
                      locGs, gSize,
                      0.0, gtFc.data(), nCright);
   coarseMat->add(gtFc, eqNumber->firstdof(myNum), FCcol);
 } 
 if(numBCs > 0 && nGright > 0) {
    GenStackFullM<Scalar> ctFg(numBCs, nGright,stackMem);
    // MLX SIGN adjusted
    for(i = 0; i < numBCs; ++i)
       for(j = 0; j < nGright; ++j)
           ctFg[i][j] = FG[j*gSize + (*locBCs)[0][i] ];
     coarseMat->add(ctFg,eqNumber->firstdof(myNum+cOffset),FGcol);
 }

 if(numBCs > 0 && nCright > 0) {
     GenStackFullM<Scalar> ctFc(numBCs, nCright,stackMem);
     int i,j;
     for(i = 0; i < numBCs; ++i)
       for(j = 0; j < nCright; ++j) {
           ctFc[i][j] = cornerSign*FC[j*gSize + (*locBCs)[0][i] ];
       }
     coarseMat->add(ctFc,eqNumber->firstdof(myNum+cOffset),FCcol);
 }


 // Now loop on the neighbors
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int leftFGcol = eqNumber->firstdof(neighbs[jSub]),
       leftFCcol;

   int nCleft = neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub];
   if(nCleft > 0)
       leftFCcol = eqNumber->firstdof(neighbs[jSub]+cOffset) +
                         neighbCs[neighbCSubIndex[jSub]][1];

     if(neighbNumRBMs[jSub] > 0 && nCright > 0) {
     // MLX SIGN adjusted
       GenStackFullM<Scalar> gtFc(neighbNumRBMs[jSub], nCright, stackMem);
       Tgemm('T', 'N', nCright, neighbNumRBMs[jSub], subSize(jSub), -cornerSign,
             FC+subOffset[jSub], gSize, neighbGs[jSub], leadingDimGs[jSub],
             0.0, gtFc.data(), nCright);
       coarseMat->add(gtFc, leftFGcol, FCcol);
     }

     if(nCleft == 0) continue;

     int (*leftCs)[3] = neighbCs + neighbCSubIndex[jSub];
     if(nCleft > 0 && nGright > 0) {
       GenStackFullM<Scalar> ctFg(nCleft, nGright,stackMem);
       int i,j;
       // MLX SIGN adjusted
       for(i = 0; i < nCleft; ++i)
         for(j = 0; j < nGright; ++j)
             ctFg[i][j] = -FG[j*gSize + leftCs[i][0] ];
       coarseMat->add(ctFg,leftFCcol,FGcol);
     }

     if(nCleft > 0 && nCright > 0) {
       GenStackFullM<Scalar> ctFc(nCleft, nCright,stackMem);
       int i,j;
       for(i = 0; i < nCleft; ++i)
         for(j = 0; j < nCright; ++j)
             ctFc[i][j] = -cornerSign*FC[j*gSize + leftCs[i][0] ];
       coarseMat->add(ctFc,leftFCcol,FCcol);
     }
 }

}

template<class Scalar>
void
GenCoarseSet<Scalar>::addCGContrib(GenSparseMatrix<Scalar> *coarseMat,
                                   EqNumberer* eqNumber, int cOffset, int gOffset)
{
 if(numBCs == 0 || isDynamic ) return;
 int i,j;

 // Find the maximum matrix allocation size
 int maxGs = numGs;
 int jSub;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   if(eqNumber->weight(neighbs[jSub]+gOffset) > maxGs)
     maxGs = eqNumber->weight(neighbs[jSub]+gOffset);
 }
 Scalar *stackMem  = (Scalar *) dbg_alloca(sizeof(Scalar)*(numBCs*maxGs));
 Scalar *stackMem2 = (Scalar *) dbg_alloca(sizeof(Scalar)*(numBCs*maxGs));

 GenStackFullM<Scalar> ctg(numBCs, numGs,stackMem);
 GenStackFullM<Scalar> gtc(numGs, numBCs,stackMem2);
 ctg.zero();
 gtc.zero();
 int myFCRow = eqNumber->firstdof(myNum+cOffset);
 int myGRow  = eqNumber->firstdof(myNum+gOffset);

 if ((numGs>0) && (numBCs >0) )  {
   for (i=0; i< numGs; i++)
     for (j=0; j< numBCs; j++) {
        int GoffSet = i * gSize;
     // MLX SIGN adjusted
        gtc[i][j] = ctg[j][i] = locGs[GoffSet + (*locBCs)[0][j]];
       }
 }
 coarseMat->add(ctg, myFCRow, myGRow);
 coarseMat->add(gtc, myGRow, myFCRow);

  // Now loop on the neighbors
 int iC = 0;
 for(jSub =0; jSub < numNeighb; ++jSub) {
   int neighbNumGs = eqNumber->weight(neighbs[jSub]+gOffset);
   GenStackFullM<Scalar> ctg(numBCs, neighbNumGs,stackMem);
   GenStackFullM<Scalar> gtc(neighbNumGs, numBCs, stackMem2);
   ctg.zero();
   gtc.zero();
   while(iC < numBCs && (*locBCs)[3][iC] == jSub) {
     int i;
     for(i =0; i < neighbNumGs; ++i)
       // MLX SIGN adjusted
       gtc[i][iC] = ctg[iC][i] = 
         -neighbGs[jSub][i*leadingDimGs[jSub]+(*locBCs)[0][iC]-subOffset[jSub]];
     iC = iC+1;
   }
  coarseMat->add(ctg,myFCRow, eqNumber->firstdof(neighbs[jSub]+gOffset));
  coarseMat->add(gtc,eqNumber->firstdof(neighbs[jSub]+gOffset), myFCRow);
 }
}


template<class Scalar>
void
GenCoarseSet<Scalar>::addNuC(Scalar *vec, Scalar *sv, EqNumberer *eqNums, int cOffset)
{
 int i;
 Scalar *nu;
 nu = sv + eqNums->firstdof(myNum+cOffset);
 for (i=0; i< numBCs; i++) {
   int index = (*locBCs)[0][i];
   vec[index] += nu[i];
 }
 int iSub, iC;
 // ML WARNING !!! CAREFUL ABOUT THE SIGN HERE! NOT SAME CONVENTION AS G
 for(iSub = 0; iSub < numNeighb; ++iSub)
   for(iC = neighbCSubIndex[iSub]; iC < neighbCSubIndex[iSub+1]; ++iC) {
      Scalar *nu = sv + eqNums->firstdof(neighbs[iSub]+cOffset);
      vec[neighbCs[iC][0]] -= nu[neighbCs[iC][1]];
   }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::reGetNeighbGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *)
{
 int iSub;
 for(iSub = 0; iSub < numNeighb; ++iSub) {
   neighbGs[iSub]     = csets[neighbs[iSub]].subG(myNum);
   leadingDimGs[iSub] = csets[neighbs[iSub]].gSize;
 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::reGetNeighbQGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd, int, 
                                     GenSolver<Scalar> *solver)
{
 int iSub, j;

 int totNumQG = 0;
 for(iSub = 0; iSub < numNeighb; ++iSub)
   totNumQG += neighbNumRBMs[iSub];

 for(iSub = 0; iSub < numNeighb; ++iSub) {
   Scalar *srcG = csets[neighbs[iSub]].subG(myNum);
   neighbGs[iSub] = srcG;
   int len = csets[neighbs[iSub]].gSize;
   leadingDimGs[iSub] = len;

   int i;
   for(i = 0; i < neighbNumRBMs[iSub]; ++i) {
     for(j = 0; j < gSize; ++j)
       neighbQGs[iSub][i*gSize + j] = 0.0;
     for(j = subOffset[iSub]; j <  subOffset[iSub+1] && j <gSize; ++j)
        neighbQGs[iSub][i*gSize + j] = -srcG[i*len+j-subOffset[iSub]];
     if(!solver)
       sd->multKbb(neighbQGs[iSub]+i*gSize, neighbQGs[iSub]+i*gSize);
   }

  if(solver)
    sd->multMFi(solver, neighbQGs[iSub], neighbQGs[iSub],neighbNumRBMs[iSub]);

 }
}

template<class Scalar>
void
GenCoarseSet<Scalar>::reComputeNeighbFBCs(GenSubDomain<Scalar> *sd, GenSolver<Scalar> *solver)
{

 int iSub, j;

 for(iSub = 0; iSub < numNeighb; ++iSub) {

   int neighbNumCs = neighbCSubIndex[iSub+1]-neighbCSubIndex[iSub];

   int i;
   for(i = 0; i < neighbNumCs; ++i) {
     for(j = 0; j < gSize; ++j)
        neighbFBCs2[iSub][i*gSize + j] = 0.0;

     int iC = i+ neighbCSubIndex[iSub];
     neighbFBCs2[iSub][i*gSize + neighbCs[iC][0]] = +1;
   }
   sd->multMFi(solver, neighbFBCs2[iSub], neighbFBCs2[iSub],neighbNumCs);
 }
}

template class GenCoarseSet<double>;
template class GenCoarseSet<std::complex<double>>;