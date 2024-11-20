#ifndef _SKYMATRIX_C_
#define _SKYMATRIX_C_

#include <cstdio>
#include <iostream>

#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/Skyline.d/utility.h>
#include <Threads.d/PHelper.h>
#include <Timers.d/GetTime.h>
#ifdef USE_MPI
#include <Driver.d/Communicator.h>
#endif


extern "C"      {
void _FORTRAN(slacol)( int*, int*, int&, int&, int& );


// *** REAL FUNCTIONS ***
void _FORTRAN(svbu4rb)(double*,int*,int*,const double*,double*,
                       int*,double&,int&,int&,int&,int&,double*);

// parallel skyline factoring routine
void _FORTRAN(pfact)(double*,int*,int*,const double*,double*,
                     int*,double&,int&,int&,int&,int&,double*, int&, int&,
                     void *, void *);

void _FORTRAN(svbu4gmb)(double*,int*,int*,const double*,double*,double*,double*,
                        int*,double&,int&,int&,double&,int&,
		          double*,int&,double*,double*,int&,int*);

// forward/backward routine for 1 to 8 rhs vectors (nonsingular skymatrix)
void _FORTRAN(forbackr1nscf)(const double*,const int*,double*, const int &);
void _FORTRAN(forbackr2ns)(const double*,const int*,double*,double*,const int &);
void _FORTRAN(forbackr3ns)(const double*,const int*,double*,double*,double*,const int &);
void _FORTRAN(forbackr4ns)(const double*,const int*,double*,double*,double*,double*,const int &);
void _FORTRAN(forbackr5ns)(const double*,const int*,double*,double*,double*,
                                        double*,double*,const int &);
void _FORTRAN(forbackr6ns)(const double*,const int*,double*,double*,double*,
                                        double*,double*,double*,const int &);
void _FORTRAN(forbackr7ns)(const double*,const int*,double*,double*,double*,
                                        double*,double*,double*,
                                        double*,const int &);
void _FORTRAN(forbackr8ns)(const double*,const int*,double*,double*,double*,
                                        double*,double*,double*,
                                        double*,double*,const int &);

// forward/backward routine for 1 to 8 rhs vectors
void _FORTRAN(forbackr1cf)(const double*, const int*,double*,const int *, const int &);
void _FORTRAN(forbackr2)(const double*,const int*,double*,double*,const int *, const int &);
void _FORTRAN(forbackr3)(const double*,const int*,double*,double*,double*, const int *, const int &);
void _FORTRAN(forbackr) (const double*,const int*,double*,double*,double*,double*,
                        const int *, const int &);
void _FORTRAN(forbackr5)(const double*,const int*,double*,double*,double*,double*,
                         double*, const int *, const int &);
void _FORTRAN(forbackr6)(const double*,const int*,double*,double*,double*,double*,
                         double*, double*,const int *, const int &);
void _FORTRAN(forbackr7)(const double*,const int*,double*,double*,double*,double*,
                         double*, double*, double*,const int *, const int &);
void _FORTRAN(forbackr8)(const double*,const int*,double*,double*,double*,double*,
                         double*, double*, double*, double*,const int *, const int &);

// forward routine for 1 rhs vector (nonsingular skymatrix)
void _FORTRAN(forr1ns)(double*,int*,double*,int &);

// backward routine for 1 rhs vector (nonsingular skymatrix)
void _FORTRAN(backr1ns)(double*,int*,double*,int &);

// forward routine for 1 rhs vector (nonsingular skymatrix)
void _FORTRAN(forr1)(double*,int*,double*,int*,int &);

// backward routine for 1 rhs vector (nonsingular skymatrix)
void _FORTRAN(backr1)(double*,int*,double*,int*,int &);


// *** COMPLEX FUNCTIONS ***
void _FORTRAN(svbu4cb)(DComplex*,int*,int*,const DComplex*,
                       DComplex*,int*,
                       double&,int&,int&,int&,int&,DComplex*);

// parallel skyline factoring routine
void _FORTRAN(pfactc)(DComplex*,int*,int*,const DComplex*,DComplex*,int*,
                      double&,int&,int&,int&,int&,DComplex*, int&, int&,
                      void *, void *);

// forward/backward routine for 1 to 4 rhs vectors
void _FORTRAN(forbackc1)(const DComplex*,const int*,DComplex*,const int *, const int &);
void _FORTRAN(forbackc2)(const DComplex*,const int*,DComplex*,DComplex*,const int *,const int &);
void _FORTRAN(forbackc3)(const DComplex*,const int*,DComplex*,DComplex*,DComplex*,
                        const int *, const int &);
void _FORTRAN(forbackc4)(const DComplex*,const int*,DComplex*,DComplex*,DComplex*,
                         DComplex*,const int *, const int &);
}

// *** TEMPLATE FUNCTIONS ***
inline void Tsvbu4b(double *a, int *b, int *c, const double *d, double *e,
                    int *f, double &g, int &h, int &i, int &j, int &k, double *l)
{ _FORTRAN(svbu4rb)(a,b,c,d,e,f,g,h,i,j,k,l); }
inline void Tsvbu4b(DComplex *a, int *b, int *c, const DComplex *d, DComplex *e,
                    int *f, double &g, int &h, int &i, int &j, int &k, DComplex *l)
{ _FORTRAN(svbu4cb)(a,b,c,d,e,f,g,h,i,j,k,l); }

inline void Tpfact(double *a, int *b, int *c, const double *d, double *e,
                   int *f, double &g, int &h, int &i, int &j, int &k, double *l,
                   int &m, int &n, void *o, void *p)
{ _FORTRAN(pfact)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p); }
inline void Tpfact(DComplex *a, int *b, int *c, const DComplex *d, DComplex *e,
                   int *f, double &g, int &h, int &i, int &j, int &k, DComplex *l,
                   int &m, int &n, void *o, void *p)
{ _FORTRAN(pfactc)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p); }

inline void Tforback1(const double *a, const int *b, double *c1, const int *d, const int &e, int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr1cf)(a,b,c1,d,e);
  else _FORTRAN(forbackr1nscf)(a,b,c1,e);
}
inline void Tforback1(DComplex *a, const int *b, DComplex *c1, const int *d, const int &e, int nzem = 0)
{ _FORTRAN(forbackc1)(a,b,c1,d,e); }

inline void Tforback2(double *a, const int *b, double *c1, double *c2, const int *d, const int &e, const int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr2)(a,b,c1,c2,d,e);
  else _FORTRAN(forbackr2ns)(a,b,c1,c2,e);
}
inline void Tforback2(DComplex *a, const int *b, DComplex *c1, DComplex *c2, const int *d, const int &e, const int nzem = 0)
{ _FORTRAN(forbackc2)(a,b,c1,c2,d,e); }

inline void Tforback3(double *a, const int *b, double *c1, double *c2, double *c3,
                      const int *d, const int &e, const int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr3)(a,b,c1,c2,c3,d,e);
  else _FORTRAN(forbackr3ns)(a,b,c1,c2,c3,e);
}
inline void Tforback3(DComplex *a, const int *b, DComplex *c1, DComplex *c2, DComplex *c3,
                      const int *d, const int &e, const int nzem = 0)
{ _FORTRAN(forbackc3)(a,b,c1,c2,c3,d,e); }

inline void Tforback4(double *a, const int *b, double *c1, double *c2, double *c3, double *c4,
                      const int *d, const int &e, const int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr)(a,b,c1,c2,c3,c4,d,e);
  else _FORTRAN(forbackr4ns)(a,b,c1,c2,c3,c4,e);
}
inline void Tforback4(DComplex *a, const int *b, DComplex *c1, DComplex *c2, DComplex *c3, DComplex *c4,
                      const int *d, const int &e, const int nzem = 0)
{ _FORTRAN(forbackc4)(a,b,c1,c2,c3,c4,d,e); }

inline void Tforback5(double *a, const int *b, double *c1, double *c2, double *c3, double *c4,
                      double *c5, const int *d, const int &e, const int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr5)(a,b,c1,c2,c3,c4,c5,d,e);
  else _FORTRAN(forbackr5ns)(a,b,c1,c2,c3,c4,c5,e);
}
inline void Tforback5(DComplex *a, const int *b, DComplex *c1, DComplex *c2, DComplex *c3, DComplex *c4,
                      DComplex *c5, const int *d, const int &e, const int nzem = 0)
{ fprintf(stderr, "forbackc5 not implemented"); }

inline void Tforback6(double *a, const int *b, double *c1, double *c2, double *c3, double *c4,
                      double *c5, double *c6, const int *d, const int &e, const int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr6)(a,b,c1,c2,c3,c4,c5,c6,d,e);
  else _FORTRAN(forbackr6ns)(a,b,c1,c2,c3,c4,c5,c6,e);
}
inline void Tforback6(DComplex *a, const int *b, DComplex *c1, DComplex *c2, DComplex *c3, DComplex *c4,
                      DComplex *c5, DComplex *c6, const int *d, const int &e, const int nzem = 0)
{ fprintf(stderr, "forbackc6 not implemented"); }

inline void Tforback7(double *a, const int *b, double *c1, double *c2, double *c3, double *c4,
                      double *c5, double *c6, double *c7, const int *d, const int &e, const int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr7)(a,b,c1,c2,c3,c4,c5,c6,c7,d,e);
  else _FORTRAN(forbackr7ns)(a,b,c1,c2,c3,c4,c5,c6,c7,e);
}
inline void Tforback7(DComplex *a, const int *b, DComplex *c1, DComplex *c2, DComplex *c3, DComplex *c4,
                      DComplex *c5, DComplex *c6, DComplex *c7, const int *d, const int &e, const int nzem = 0)
{ fprintf(stderr, "forbackc7 not implemented"); }

inline void Tforback8(double *a, const int *b, double *c1, double *c2, double *c3, double *c4,
                      double *c5, double *c6, double *c7, double *c8, const int *d, const int &e, const int nzem = 0)
{
  if(nzem) _FORTRAN(forbackr8)(a,b,c1,c2,c3,c4,c5,c6,c7,c8,d,e);
  else _FORTRAN(forbackr8ns)(a,b,c1,c2,c3,c4,c5,c6,c7,c8,e);
}
inline void Tforback8(DComplex *a, const int *b, DComplex *c1, DComplex *c2, DComplex *c3, DComplex *c4,
                      DComplex *c5, DComplex *c6, DComplex *c7, DComplex *c8, const int *d, const int &e,
                      const int nzem = 0)
{ fprintf(stderr, "forbackc8 not implemented"); }

inline void Tfor1(double *a, int *b, double *c1, int *d, int &e, int nzem = 0)
{
  if(nzem) _FORTRAN(forr1)(a,b,c1,d,e);
  else _FORTRAN(forr1ns)(a,b,c1,e);
}
inline void Tfor1(DComplex *a, int *b, DComplex *c1, int *d, int &e, int nzem = 0)
{
  std::cerr << "WARNING ...";
}

inline void Tback1(double *a, int *b, double *c1, int *d, int &e, int nzem = 0)
{
  if(nzem) _FORTRAN(backr1)(a,b,c1,d,e);
  else _FORTRAN(backr1ns)(a,b,c1,e);
}
inline void Tback1(DComplex *a, int *b, DComplex *c1, int *d, int &e, int nzem = 0)
{
  std::cerr << "WARNING ...";
}

template<class Scalar>
GenSkyMatrix<Scalar>::~GenSkyMatrix()
{
 if(skyA) { delete [] skyA; skyA = 0; }
 if(scale) { delete [] scale; scale=0; }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::clean_up()
{
 numUncon = 0;
 if(skyA) {
   delete [] skyA;
   skyA=0;
 }
 if(pivot) {
   delete [] pivot;
   pivot=0;
 }
 if(dlp) {
   delete [] dlp;
   dlp=0;
 }
 if(rbm && myRbm) { rbm->clean_up(); }
 if(scale) { delete [] scale; scale=0; }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::zeroAll()
{
  if(numUncon == 0) return;
  int i;
  for(i = 0; i < dlp[numUncon-1]; ++i)
    skyA[i] = 0.0;
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::allMult(Scalar x)
{
 int i;
 for(i = 0; i < dlp[numUncon-1]; ++i)
   skyA[i] *= x;
}

template<class Scalar>
GenSkyMatrix<Scalar>::GenSkyMatrix(const Connectivity *cn, const DofSetArray *c_dsa,
                                   double trbm, Rbm *rigid) :
SkyData(cn,c_dsa,trbm,rigid)
{
  if(numUncon==0) { skyA=0; nzem=0; return; }
  constructTime = -getTime();

  // ALLOCATE MEMORY FOR SKYA
  skyA = new Scalar[ dlp[numUncon-1] ];

  // initialize to zero
  zeroAll();
  solveTime = 0.0;
  constructTime += getTime();
  if(rigid) isTRBM=0;
  else isTRBM=1;

  isScaled = 0;
  scale = 0;
  wasScaled = false;
}

template<class Scalar>
GenSkyMatrix<Scalar>::GenSkyMatrix(const Connectivity *cn, const EqNumberer *_dsa,
                                   double trbm, const int *bc) :
SkyData(cn,_dsa,trbm,bc)
{
  constructTime = -getTime();
  if(numUncon==0) { skyA=0; nzem=0; return; }

  // ALLOCATE MEMORY FOR SKYA
  skyA = new Scalar[ dlp[numUncon-1] ];

  // INITIALIZE SKY MATRIX TO ZERO
  zeroAll();
  solveTime = 0.0;
  constructTime += getTime();
  isTRBM = 1;

  isScaled = 0;
  scale = 0;
  wasScaled = false;
}

template<class Scalar>
GenSkyMatrix<Scalar>::GenSkyMatrix(const Connectivity *cn, const EqNumberer *_dsa,
                                   double trbm, int _isScaled) :
SkyData(cn,_dsa,trbm)
{
  constructTime = -getTime();
  if(numUncon == 0) { skyA=0; nzem=0; return; }

  // ALLOCATE MEMORY FOR SKYA
  skyA = new Scalar[ dlp[numUncon-1] ];

  // INITIALIZE SKY MATRIX TO ZERO
  zeroAll();
  solveTime = 0.0;
  constructTime += getTime();
  isTRBM = 1;

  isScaled = _isScaled;
  scale = 0;
  wasScaled = false;
}

template<class Scalar>
GenSkyMatrix<Scalar>::GenSkyMatrix(const Connectivity *cn, const EqNumberer *_dsa,
                                   double trbm, const int *rCN, int) :
SkyData( _dsa, cn, trbm, rCN )
{
  constructTime = -getTime();

  // ALLOCATE MEMORY FOR SKYA
  if(numUncon == 0) { skyA=0; nzem=0; return; }
  skyA = new Scalar[ dlp[numUncon-1] ];

  // INITIALIZE SKYA TO ZERO
  zeroAll();
  solveTime = 0.0;
  constructTime += getTime();
  isTRBM = 1;

  isScaled = 0;
  scale = 0;
  wasScaled = false;
}

template<class Scalar>
GenSkyMatrix<Scalar>::GenSkyMatrix(GenFullM<Scalar> *mat, double tolerance) :
SkyData(mat->numRow(), tolerance)
{
  constructTime = -getTime();

  // WARNING: DEFAULT value of tolerance is 1.0E-4
  if(numUncon==0) { skyA=0; nzem=0; return; }

  if(mat->numRow() != mat->numCol()) {
    skyA = 0;
    fprintf(stderr, " ERROR: SkyMatrix::SkyMatrix(FullM *mat) constructor\n");
  }
  else {
    skyA = new Scalar[ dlp[numUncon-1] ];
    int offset = 0;
    int i,j;
    for(i = 0; i < numUncon; ++i)
      for(j = 0; j <= i; ++j) {
        skyA[offset++] = (*mat)[i][j];
      }
  }
  constructTime += getTime();
  isTRBM = 1;

  isScaled = 0;
  scale = 0;
  wasScaled = false;
}

template<class Scalar>
GenSkyMatrix<Scalar>::GenSkyMatrix(const Connectivity *cn, const EqNumberer *,
                                   const ConstrainedDSA *c_dsa, double trbm) :
  SkyData(cn, c_dsa, trbm)
{
  // ALLOCATE MEMORY FOR SKYA
  skyA = new Scalar[ dlp[numUncon-1] ];

  // INITIALIZE TO ZERO
  zeroAll();

  isScaled = 0;
  scale = 0;
  wasScaled = false;
}

template<class Scalar>
GenSkyMatrix<Scalar>::GenSkyMatrix(int n, double tolerance) :
  SkyData(n, tolerance)
{
   // WARNING: DEFAULT value of tolerance is 1.0E-4
  if(numUncon==0) { skyA=0; nzem=0; return; }

  skyA = new Scalar[ dlp[numUncon-1] ];

  auto build_rcn = new int[numUncon];
  for(int i=0; i<numUncon; ++i)
      build_rcn[i] = i;
  rowColNum = build_rcn;

  zeroAll();
  isTRBM = 1;

  isScaled = 0;
  scale = 0;
  wasScaled = false;
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::addPoint(Scalar value, int i, int j) {

 int ri, rj;

 ri = rowColNum[i];
 rj = rowColNum[j];

 if (i>j)
   return;

 if ((ri>-1) && (rj>-1)) {
   skyA[dlp[rj] - rj + ri - 1] += value;
 }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::printMemory()
{
   fprintf(stderr,"Memory necessary for Skyline array (Mb): %10.3f\n",
          8.0*dlp[numUncon-1]/(1024.0*1024.0));
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::printConstructTime()
{
  fprintf(stderr,"SkyMatrix Construction Time = %14.5f\n",constructTime/1000.0);
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::mult(const GenVector<Scalar> &, GenVector<Scalar> &) const
{
  fprintf(stderr,"This shouldn't be called--SkyMatrix::mult\n");

 // For this routine to work, fortranDlp has to be created from dlp in
 // the constructor, for details of this, look at an old version of the
 // Skyline constructor. But a sky matrix - vector multiplication should
 // never be done since it is an expensive computation versus a
 // sparse - matrix vector multiply.

 // _FORTRAN(skymul)(skyA,numUncon,fortranDlp,rhs.data(),result.data(),iop,rhs.data(),rhs.data());
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result) const
{
  // local objects
  int i, j, istr, iend, index;

  // initialize
  for( i = 0; i < numUncon; i++)
    result[ i] = 0.0;

  // result = SkyMatrix * rhs
  for( j = 0; j < numUncon; j++) {
    istr = dlp[ j];
    if ( j > 0) istr -= dlp[ j - 1];
    istr = j - istr + 1;
    iend = j + 1;
    for( i = istr; i < iend; i++) {
      index = dlp[ j] - j + i - 1;
      result[ i] += skyA[ index] * rhs[ j]; // upper half (incl. diagonal)
      if ( i != j) // part of skymatrix below diagonal
        result[ j] += skyA[ index] * rhs[ i];
    }
  }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::factor()
{
  if(numUncon==0) return;

  if((isScaled)&&(!wasScaled))
    symmetricScaling();

 // Prints number of unconstrained dof and average column height
 // fprintf(stderr, "Stats %d %f\n",numUncon,float(dlp[numUncon-1])/numUncon);

  if(isTRBM) {
    //std::cerr << " ... Using TRBM method to factor Skyline Matrix ...\n";
    Factor();    // tolerance method
  }
  else {
    //std::cerr << " ... Using GRBM method to factor Skyline Matrix ...\n";
    Factor(rbm); // geometric method
  }
}

#include <unistd.h>
template<class Scalar>
void
GenSkyMatrix<Scalar>::Factor()
{
   //filePrint(stderr," ... Skyline factor: neq %d size %d avg. band. %d \n",
   //          neqs(), size(), size() / neqs());

   //Scalar *w = (Scalar*) dbg_alloca(5*numUncon*sizeof(Scalar));
   Scalar *w = new Scalar[5*numUncon];

   Scalar *dummyZEM = w + 4*numUncon;
   Scalar *dummyRHS = NULL; // send a dummy argument for rhs during FACTORING

// Factor matrix, sending a dummy rhs & dummy ZEM.
   int flag = 1, nops = 0;
   int numrbm = 0;

   Tsvbu4b(skyA, dlp, lacol, dummyRHS, w, pivot, TOLERANCE, numUncon,
           flag, nops, numrbm, dummyZEM);

   if(numrbm > 0 && this->print_nullity)
     std::cerr << " ... Matrix is singular: size = " << numUncon << ", rank = " << numUncon-numrbm << ", nullity = " << numrbm << " ...\n";

   // set number of zero energy modes
   nzem = numrbm;

   GenVector<Scalar> *zem = 0;
   if(nzem > 0) zem = getNullSpace();
   if(rbm && myRbm) delete rbm;
   myRbm = 1;
   rbm = new Rbm(zem, nzem, numUncon, myRbm);
   delete [] w;
}

#if defined(sgi) && !defined(_OPENMP)
extern ulock_t allocLock;
template<class Scalar>
void
GenSkyMatrix<Scalar>::pfact(int me, int nprocs, barrier_t *bar, Scalar *w)
#else
template<class Scalar>
void
GenSkyMatrix<Scalar>::pfact(int me, int nprocs, Scalar *w)
#endif
{
   Scalar *dummyZEM = w + 4*numUncon;
   Scalar *dummyRHS = NULL; // send a dummy argument for rhs during FACTORING

// Factor the matrix, sending a dummy rhs & dummy ZEM.
   int flag = 1, nops = 0;
   int numrbm = 0;
#if defined(sgi) && !defined(_OPENMP)
   Tpfact(skyA, dlp, lacol, dummyRHS, w, pivot, TOLERANCE, numUncon, flag,
          nops, numrbm, dummyZEM, me, nprocs, &bar, &allocLock);
   if(me == 0) nzem = numrbm;
#else
// Sequential version of skyline factoring
   Tsvbu4b(skyA, dlp, lacol, dummyRHS, w, pivot, TOLERANCE, numUncon,
           flag, nops, numrbm, dummyZEM);
   if(me == 0) nzem = numrbm;
#endif
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::parallelFactor()
{
  if((isScaled)&&(!wasScaled))
    symmetricScaling();

  //Scalar *w  = (Scalar *) dbg_alloca(5*numUncon*sizeof(Scalar));
  Scalar *w  = new Scalar[5*numUncon];

  Scalar *dummyZEM = w + 4*numUncon;
  Scalar *dummyRHS = NULL; // send a dummy argument for rhs during FACTORING

// Factor the matrix, sending a dummy rhs & dummy ZEM.
  int flag = 1, nops = 0;
  int numrbm = 0;

#if defined(sgi) && !defined(_OPENMP)
  //int avg = size() / neqs();
  int avg1= int(rmsBandwidth()); //HB: better estimate of the avg. bandwidth
  //int maxCPU = (avg*avg)/2500;
  int maxCPU = (avg1*avg1)/2500;
  if(maxCPU > 1 && threadManager->numThr() > 1) {
     barrier_t *barrier = threadManager->getBarrier();
     int npr = threadManager->numThr();
     if(npr > maxCPU) npr = maxCPU;
     // filePrint(stderr," ... Skyline parallel factor: neq %d size %d avg. band. %d rms band. %d #CPUs %d\n",
     //           neqs(), size(), avg, avg1, npr);
     execParal(npr, this, GenSkyMatrix::pfact, npr, barrier, w);
  } else {
#endif
     // filePrint(stderr," ... Skyline parallel factor: neq %d size %d avg. band. %d rms band. %d #CPUs %d\n",
     //           neqs(), size(), avg, avg1, 1);
     Tsvbu4b(skyA, dlp, lacol, dummyRHS, w, pivot,
             TOLERANCE, numUncon, flag, nops, numrbm, dummyZEM);
     nzem = numrbm;
#if defined(sgi) && !defined(_OPENMP)
   }
#endif
  if(numrbm > 0 && this->print_nullity)
     std::cerr << " ... Matrix is singular: size = " << numUncon << ", rank = " << numUncon-numrbm << ", nullity = " << numrbm << " ...\n";

  // set number of zero energy modes
  GenVector<Scalar> *zem = getNullSpace();
  if(rbm && myRbm) delete rbm;
  myRbm = 1;
  rbm = new Rbm(zem, nzem, numUncon, myRbm);
  delete [] w;
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution)
{
   solution = rhs;
   reSolve(solution.data());
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::solve(const Scalar *rhs, Scalar *solution)
{
   int i;
   for(i=0; i<dim(); ++i)
     solution[i] = rhs[i];

   reSolve(solution);
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
   if(numUncon==0) return;
   solveTime -= getTime();
   applyScaling(rhs.data());
   Tforback1(skyA, dlp, rhs.data(), pivot, numUncon, nzem);
   applyScaling(rhs.data());
   solveTime += getTime();
}
template<class Scalar>
void
GenSkyMatrix<Scalar>::forward(GenVector<Scalar> &rhs)
{
   if(numUncon==0) return;
   solveTime -= getTime();
   //applyScaling(rhs.data());
   Tfor1(skyA, dlp, rhs.data(), pivot, numUncon, nzem);
   //applyScaling(rhs.data());
   solveTime += getTime();
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::backward(GenVector<Scalar> &rhs)
{
   if(numUncon==0) return;
   solveTime -= getTime();
   //applyScaling(rhs.data());
   Tback1(skyA, dlp, rhs.data(), pivot, numUncon, nzem);
   //applyScaling(rhs.data());
   solveTime += getTime();
}


template<class Scalar>
void
GenSkyMatrix<Scalar>::reSolve(Scalar *rhs)
{
   if(numUncon==0) return;
   solveTime -= getTime();
   applyScaling(rhs);
   Tforback1(skyA, dlp, rhs, pivot, numUncon, nzem);
   applyScaling(rhs);
   solveTime += getTime();
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::reSolve(int nRHS, Scalar **rhs)
{
  if(numUncon==0) return;

  solveTime -= getTime();
  int n;
  for(n=0; n<nRHS; ++n)
    applyScaling(rhs[n]);

  int i = 0;
  int multiple = 4;
  if(nRHS <= 4)  // PJSA: changed to 4 (was 6) since > 4 not implemented for complex yet
    multiple = nRHS;

  switch(multiple) {
    default:
    case 8:
      {
        for( ; i < nRHS-7; i += 8)
          Tforback8(skyA,dlp,rhs[i],rhs[i+1],rhs[i+2],rhs[i+3],
                    rhs[i+4],rhs[i+5],rhs[i+6],rhs[i+7],pivot,numUncon,nzem);
      }
    case 7:
      {
        for( ; i < nRHS-6; i += 7)
          Tforback7(skyA,dlp,rhs[i],rhs[i+1],rhs[i+2],rhs[i+3],
                    rhs[i+4],rhs[i+5],rhs[i+6],pivot,numUncon,nzem);
      }
    case 6:
      {
        for( ; i < nRHS-5; i += 6)
          Tforback6(skyA,dlp,rhs[i],rhs[i+1],rhs[i+2],rhs[i+3],
                    rhs[i+4],rhs[i+5],pivot,numUncon,nzem);
      }
    case 5:
      {
        for( ; i < nRHS-4; i += 5)
          Tforback5(skyA,dlp,rhs[i],rhs[i+1],rhs[i+2],rhs[i+3],
                    rhs[i+4],pivot,numUncon,nzem);
      }
    case 4:
      {
        for( ; i < nRHS-3; i += 4)
          Tforback4(skyA,dlp,rhs[i],rhs[i+1],rhs[i+2],rhs[i+3],
                    pivot,numUncon,nzem);
      }
    case 3:
      {
        for( ; i < nRHS-2; i += 3)
          Tforback3(skyA,dlp,rhs[i],rhs[i+1],rhs[i+2],
                    pivot,numUncon,nzem);
      }
    case 2:
      {
        for( ; i < nRHS-1; i += 2)
          Tforback2(skyA,dlp,rhs[i],rhs[i+1],pivot,numUncon,nzem);
      }
    case 1:
      {
        for( ; i < nRHS; ++i)
          Tforback1(skyA,dlp,rhs[i],pivot,numUncon,nzem);
      }
      break;
  }

  for(n=0; n<nRHS; ++n)
    applyScaling(rhs[n]);

  solveTime += getTime();
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::reSolve(int nRHS, GenVector<Scalar> *rhs)
{
  if(numUncon==0) return;

  solveTime -= getTime();
  int n;
  for(n=0; n<nRHS; ++n)
    applyScaling(rhs[n].data());

  int i=0;
  int multiple = 4;
  if(nRHS <= 4) // PJSA: changed to 4 (was 6) since > 4 not implemented for complex yet
    multiple = nRHS;

  switch(multiple) {
    default:
    case 8:
      {
        for( ; i < nRHS-7; i += 8)
          Tforback8(skyA,dlp, rhs[i].data(),rhs[i+1].data(),rhs[i+2].data(),
                    rhs[i+3].data(),rhs[i+4].data(),rhs[i+5].data(),rhs[i+6].data(),
                    rhs[i+7].data(),pivot,numUncon,nzem);
      }
    case 7:
      {
        for( ; i < nRHS-6; i += 7)
         Tforback7(skyA,dlp,rhs[i].data(),rhs[i+1].data(),rhs[i+2].data(),
                   rhs[i+3].data(),rhs[i+4].data(),rhs[i+5].data(),rhs[i+6].data(),
                   pivot,numUncon,nzem);
      }
    case 6:
      {
        for( ; i < nRHS-5; i += 6)
          Tforback6(skyA,dlp,
                    rhs[i].data(),rhs[i+1].data(),rhs[i+2].data(),rhs[i+3].data(),
                    rhs[i+4].data(),rhs[i+5].data(),pivot,numUncon,nzem);
      }
    case 5:
      {
        for( ; i < nRHS-4; i += 5)
          Tforback5(skyA,dlp,
                    rhs[i].data(),rhs[i+1].data(),rhs[i+2].data(),rhs[i+3].data(),
                    rhs[i+4].data(),pivot,numUncon,nzem);
      }
    case 4:
      {
        for( ; i < nRHS-3; i += 4)
          Tforback4(skyA,dlp,rhs[i].data(),
                    rhs[i+1].data(),rhs[i+2].data(),rhs[i+3].data(),
                    pivot,numUncon,nzem);
      }
    case 3:
      {
        for( ; i < nRHS-2; i += 3)
          Tforback3(skyA,dlp,rhs[i].data(),rhs[i+1].data(),rhs[i+2].data(),
                    pivot,numUncon,nzem);
      }
    case 2:
      {
        for( ; i < nRHS-1; i += 2)
          Tforback2(skyA,dlp,rhs[i].data(),rhs[i+1].data(),
                    pivot,numUncon,nzem);
      }
    case 1:
      {
        for( ; i < nRHS; ++i)
          Tforback1(skyA,dlp,rhs[i].data(),pivot,numUncon,nzem);
      }
      break;
  }
  for(n=0; n<nRHS; ++n)
    applyScaling(rhs[n].data());

  solveTime += getTime();
}


template<class Scalar>
void
GenSkyMatrix<Scalar>::reSolve(int nRHS, Scalar *rhs)
{
  if(numUncon==0) return;

  solveTime -= getTime();

  int n;
  for(n=0; n<nRHS; ++n)
    applyScaling(rhs+n*numUncon);

  int i = 0;
  int multiple = 4;

  switch(multiple) {
    default:
    case 4:
      {
        for( ; i < nRHS-3; i += 4)
          Tforback4(skyA,dlp,rhs+i*numUncon,rhs+(i+1)*numUncon,
                    rhs+(i+2)*numUncon, rhs+(i+3)*numUncon,pivot,numUncon,nzem);
      }
    case 3:
      {
        for( ; i < nRHS-2; i += 3)
          Tforback3(skyA,dlp,rhs+i*numUncon,rhs+(i+1)*numUncon,
                    rhs+(i+2)*numUncon,pivot,numUncon,nzem);
      }
    case 2:
      {
        for( ; i < nRHS-1; i += 2)
          Tforback2(skyA,dlp,rhs+i*numUncon,rhs+(i+1)*numUncon,
                    pivot,numUncon,nzem);
      }
    case 1:
      {
        for( ; i < nRHS; ++i)
          Tforback1(skyA,dlp,rhs+i*numUncon,pivot,numUncon,nzem);
      }
      break;
  }

  for(n=0; n<nRHS; ++n)
    applyScaling(rhs+n*numUncon);

  solveTime += getTime();
}


template<class Scalar>
void
GenSkyMatrix<Scalar>::reSolve(GenFullM<Scalar> *mat)
{
  solveTime -= getTime();

  int i, j;
  int nRHS = mat->numCol();
  int length = mat->numRow();

  Scalar **tmp = new Scalar*[nRHS];

  // 1. Copy the matrix into a set of vectors
  for (j=0; j<nRHS; ++j)
     tmp[j] = mat->Column(j);

  // 2. Solve the system
  reSolve(nRHS,tmp);

  // 3. Overwrite into the matrix
  for (i=0; i<length; ++i)
    for (j=0; j<nRHS; ++j)
      (*mat)[i][j] = tmp[j][i];

  solveTime += getTime();

  for (i=0; i<nRHS; ++i)
     delete [] tmp[i];

  delete [] tmp;
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::getNullSpace(Scalar *rbmv)
{
  Scalar *dummyRHS = 0; // dummy rhs vector
  Scalar *w = 0;
  int nops, flag = 2;
  Tsvbu4b(skyA, dlp, lacol, dummyRHS, w, pivot,
          TOLERANCE, numUncon, flag, nops, nzem, rbmv);
}

template<class Scalar>
GenVector<Scalar> *
GenSkyMatrix<Scalar>::getNullSpace()
{
  int i, m;

// Declare & Initialize ZEM:
  GenVector<Scalar> *ZEM = new GenVector<Scalar>[nzem];
  Scalar *rbmv = new Scalar[numUncon*nzem];
  getNullSpace(rbmv);

// Copy rigid body modes (rbm) into the return vector ZEM:
  GenVector<Scalar> v(numUncon,0.0);
  for(m=0; m<nzem; ++m) {
    ZEM[m] = v;
    for(i=0; i<numUncon; ++i)
      ZEM[m][i] = rbmv[i+m*numUncon];
  }

  delete [] rbmv;
  return ZEM;
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::getRBMs(double *rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::getRBMs(Vector *rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::getRBMs(VectorSet& rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
Scalar &
GenSkyMatrix<Scalar>::diag(int dof)
{
  return skyA[dlp[dof] - 1];
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::add(const FullSquareMatrix &kel, const int *dofs)
{
// Construct stiffness matrix K (skyA)

  int i, j, ri, rj;
  int kndof = kel.dim();                	// Element stiffness dimension

  for( i = 0; i < kndof; ++i ) {           // Loop over rows.
    if( (ri = rowColNum[dofs[i]]) == -1 ) continue; // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {          // Loop over columns.
       if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.
       if( (rj = rowColNum[dofs[j]]) == -1 ) continue; // Skip constrained dofs
     /*  // MLX DEBUG
       if(rj > 0 && dlp[rj] - rj + ri - 1 < dlp[rj-1])
         std:cerr << "Adding outside of skyline!!! " << kel[i][j] << " "
         << dofs[i] << " " << dofs[j] << " "
            << ri << " " << rj << std::endl;*/
       skyA[dlp[rj] - rj + ri - 1] += kel[i][j];
    }
  }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::addone(Scalar d, int dofi, int dofj)
{
// Construct stiffness matrix K (skyA)
  int ri, rj;
  if(dofi > dofj) {
    std::cerr << "WARNING: dofi > dofj in SkyMatrix::addone() \n";
    return;  // Work with upper symmetric half.
  }
  if((ri = rowColNum[dofi]) == -1 ) return;  // Skip constrained dofs
  if((rj = rowColNum[dofj]) == -1 ) return;  // Skip constrained dofs
  // HB & PJSA: to be checked for CCt
  //if((ri = dofi) == -1 ) return;  // Skip constrained dofs
  //if((rj = dofj) == -1 ) return;  // Skip constrained dofs

  skyA[dlp[rj] - rj + ri - 1] += d;
}

//HB: add an array of data given by row-col indices & associated values
// NOT VALIDATED
/*template<class Scalar>
void
GenSkyMatrix<Scalar>::addData(int nData, int* dofi, int* dofj, Scalar* d)
{
  int ri, rj, i;
  bool skip;
  for(i=0; i<nData; i++){
    if(dofi[i] > dofj[i]) {
      std::cerr << "WARNING: dofi > dofj in SkyMatrix::addData() \n";
      return;  // Work with upper symmetric half.
    }
    skip = ((ri = rowColNum[dofi]) == -1 )? true : false ; // Skip constrained dofs
    skip = ((rj = rowColNum[dofj]) == -1 )? true : false ; // Skip constrained dofs
    // HB & PJSA: to be checked for CCt
    //if((ri = dofi) == -1 ) return;  // Skip constrained dofs
    //if((rj = dofj) == -1 ) return;  // Skip constrained dofs
    if(!skip) skyA[dlp[rj] - rj + ri - 1] += d[i];
  }
}*/

template<class Scalar>
void
GenSkyMatrix<Scalar>::add(const GenAssembledFullM<Scalar> &kel, const int *dofs)
{
 int i, j, ri, rj;
 int kndof = kel.dim();                  // Element stiffness dimension
 for( i = 0; i < kndof; ++i ) {          // Loop over rows.
    if( (ri = dofs[i]) == -1 ) continue; // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {          // Loop over columns.
       if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.
       if( (rj = dofs[j]) == -1 ) continue; // Skip constrained dofs
       skyA[dlp[rj] - rj + ri - 1] += kel[i][j];
    }
 }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::add(const GenFullM<Scalar> &knd, int fRow, int fCol)
{
  int nrow = knd.numRow(); // number of rows
  int ncol = knd.numCol(); // number of columns

  int iCol, iRow;
  for(iCol = 0; iCol < ncol; ++iCol) {
    int fk = dlp[fCol+iCol] + fRow - fCol - iCol - 1;
    int rowStop = (nrow < fCol+iCol-fRow+1) ? nrow : fCol+iCol-fRow+1;
    for(iRow = 0; iRow < rowStop; ++iRow) {
      skyA[fk+iRow] += knd[iRow][iCol];
    }
  }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::add(Scalar *_skyA)
{
  for(int i=0; i<size(); ++i) skyA[i] += _skyA[i];
}

template<class Scalar>
Scalar
GenSkyMatrix<Scalar>::getone(int fRow, int fCol)
{
  // PJSA: used to assemble CCt
  if(fRow > fCol) {
    int tmp = fRow; fRow = fCol; fCol = tmp;
  }
  Scalar d;
  int iRow;
  int fk = dlp[fCol] + fRow - fCol - 1;
  int rowStop = (1 < fCol-fRow+1) ? 1 : fCol-fRow+1;
  for(iRow = 0; iRow < rowStop; ++iRow) {
    d = skyA[fk+iRow];
  }
  return d;
}

/*
template<class Scalar>
void
GenSkyMatrix<Scalar>::addDiscreteMass(int cdof, Scalar dmass)
{
  // in this function cdof is in constrained dof numbering
  skyA[dlp[cdof]-1] += dmass;
}
*/

template<class Scalar>
void
GenSkyMatrix<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  // int this function dof is in unconstrained dof numbering
  if(dof < 0) return;
  int cdof = rowColNum[dof];
  if(cdof < 0) return;  // Skip constrained dofs
  skyA[dlp[cdof]-1] += dmass;
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
 communicator->globalSum(dlp[numUncon-1], skyA);
#endif
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::addBoeing(int nlines, const int *Kai, const int *Kaj,
                                const double *nz, const int *map, Scalar multiplier)
{
 if(numUncon==0) return;
 int i, j;
 for(i = 0; i < nlines; ++i) {
   if(map[i] == -1) continue;
   int dofI = rowColNum[map[i]];
   if(dofI < 0) continue;
   for(j = Kai[i]; j < Kai[i+1]; ++j) {
     if(map[Kaj[j-1]-1] == -1) continue;
     int dofJ = rowColNum[map[Kaj[j-1]-1]];
     if(dofJ < 0) continue;
     if(dofI < dofJ)
        skyA[dlp[dofJ] - dofJ + dofI - 1] += (nz[j-1]*multiplier);
     else
        skyA[dlp[dofI] - dofI + dofJ - 1] += (nz[j-1]*multiplier);
   }
 }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::symmetricScaling()
{
  if(!wasScaled){
  if(numUncon == 0) return;

  scale = new Scalar[numUncon];
  int i;
  for(i=0; i<numUncon; ++i) {
    if(diag(i) == 0.0) std::cerr << " *** WARNING: in GenSkyMatrix<Scalar>::symmetricScaling(), diag(i) = 0.0\n";
    scale[i] = Scalar(1.0) / ScalarTypes::sqrt(diag(i));
  }

  skyA[0] *= scale[0]*scale[0];
  for(i=1; i<numUncon; ++i) {
    int numEntries = dlp[i] - dlp[i-1];
    int j;
    for(j=0; j<numEntries; ++j)
      skyA[dlp[i]-1-j] *= scale[i-j]*scale[i];
 }
 wasScaled = true;
 }
}

template<class Scalar>
void
GenSkyMatrix<Scalar>::applyScaling(Scalar *vector) const
{
 if(isScaled) {
   int i;
   for(i=0; i<numUncon; ++i)
     vector[i] *= scale[i];
 }
}

//HB: compute the root-mean-square bandwidth
template<class Scalar>
double
GenSkyMatrix<Scalar>::rmsBandwidth() const
{
  if(numUncon == 0) return 0.0;

  double rmsBandwith = 1.0;
  for(int i=1; i<numUncon; ++i) {
    int numEntries = dlp[i] - dlp[i-1];
    rmsBandwith += numEntries*numEntries;
  }
  rmsBandwith = sqrt(rmsBandwith/numUncon);
  return rmsBandwith;
}

#endif
