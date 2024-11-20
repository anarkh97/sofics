#include <cstdio>
#include <Utils.d/linkfc.h>
#include <Math.d/Vector.h>

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

template<class Scalar, class AnyVector>
KrylovProjector<Scalar,AnyVector>::KrylovProjector(int _numDofs, int _maxVecStorage)
{
 numDofs   = _numDofs;
 numKrylov = 1;
 numSystem = 0;

 maxVecStorage = _maxVecStorage;

 kindex    = new int[maxVecStorage+1];

 wKw       = new Scalar[maxVecStorage];

 int totalSize = maxVecStorage*numDofs;

 wMat      = new Scalar[totalSize];
 kWMat     = new Scalar[totalSize];

 kindex[0] = 0;
 kindex[1] = 0;

 y	   = new Scalar[maxVecStorage];

}

template<class Scalar, class AnyVector>
void
KrylovProjector<Scalar,AnyVector>::newSystem()
{
 lastIndex = kindex[numKrylov];
 numSystem += 1;
}

template<class Scalar, class AnyVector>
int
KrylovProjector<Scalar,AnyVector>::addDirection(const AnyVector &w, const AnyVector &Kw, Scalar dotProd)
{
 // w   - direction
 // Kw  - projection
 // wKw - dot product

 int res = 0;

 if(kindex[numKrylov] == maxVecStorage) {
   res = -1;     // To notify the caller that we are restarting
   zeroKrylov(); // To zero the Krylov Space
   fprintf(stderr,"... Restarting the Krylov Projection.\n");
 }

 int curIndex  = kindex[numKrylov];
 wKw[curIndex] = dotProd;

 int index = curIndex*numDofs;
 int i, ii;
 for(i=0; i < numDofs; ++i) {
   ii = index + i;
   wMat[ii] = w[i];
   kWMat[ii] = Kw[i];
 }

 kindex[numKrylov] = kindex[numKrylov] + 1;
 return res;
}

template<class Scalar, class AnyVector>
void
KrylovProjector<Scalar,AnyVector>::zeroKrylov()
{
 numKrylov = 1;
 numSystem = 1;
 kindex[0] = 0;
 kindex[1] = 0;
 lastIndex = 0;
}

template<class Scalar, class AnyVector>
void
KrylovProjector<Scalar,AnyVector>::newKrylov() 
{
  if(kindex[numKrylov] == kindex[numKrylov-1]) return; // Empty previous Krylov
  numKrylov++;
  kindex[numKrylov] = kindex[numKrylov-1];
}

template<class Scalar, class AnyVector>
void
KrylovProjector<Scalar,AnyVector>::project(const AnyVector *resid, AnyVector *precResid)
{
 // For the first system, there are no stored directions
 if(numSystem == 1) return;

 int i, iKrylov;

 int maxit = numKrylov-1;

 if(numSystem > 1 && numKrylov == 1) maxit = 1;

 for(iKrylov = 0; iKrylov < maxit; ++iKrylov) { 

   int kSize = kindex[iKrylov+1] - kindex[iKrylov];

   // Compute y = W^t resid
   Tgemv('T', numDofs, kSize, 1.0, wMat + kindex[iKrylov]*numDofs,
         numDofs, resid->data(), 1, 0.0, y, 1);

   // Compute y = y - W^t K P r 
   Tgemv('T', numDofs, kSize, -1.0, kWMat + kindex[iKrylov]*numDofs,
         numDofs, precResid->data(), 1, 1.0, y, 1);

   // Compute y = y / wKw
   for(i = 0; i < kSize; ++i)
      y[i] /= wKw[kindex[iKrylov]+i];

   // Compute precResid = precResid + W y
   Tgemv('N', numDofs, kSize, 1.0, wMat + kindex[iKrylov]*numDofs,
         numDofs, y, 1, 1.0, precResid->data(), 1);
 }
}

template<class Scalar, class AnyVector>
int
KrylovProjector<Scalar,AnyVector>::initialization(AnyVector *rhs, AnyVector *x0)
{
 // ... Kx = b where b = rhs 
 // ... INITILIZATION BY K ORTHOGONAL PROJECTION 
 // ... WITH PREVIOUS KRYLOV SPACES

 if(numSystem == 1) return 0;

 int i, iKrylov;

 int maxit = numKrylov-1;

 if(numSystem > 1 && numKrylov == 1) maxit = 1;

 for(iKrylov = 0; iKrylov < maxit; ++iKrylov) {
   int kSize =  kindex[iKrylov+1]-kindex[iKrylov];

   // ... COMPUTE y = W^t rhs
   Tgemv('T', numDofs, kSize, 1.0, wMat + kindex[iKrylov]*numDofs,
         numDofs, rhs->data(), 1, 0.0, y, 1);

   for(i = 0; i < kSize; ++i)
      y[i] /= wKw[kindex[iKrylov]+i];

   // ... COMPUTE x0(new) = W y + x0
   Tgemv('N', numDofs, kSize, 1.0, wMat + kindex[iKrylov]*numDofs,
       numDofs, y, 1, 1.0, x0->data(), 1);

 }

 return 1;
}

template<class Scalar, class AnyVector>
void
KrylovProjector<Scalar,AnyVector>::ortho(AnyVector *direction, int)
{
 int kSize = kindex[numKrylov] - lastIndex;

 // ... Compute y = y - W^t K direction
 Tgemv('T', numDofs, kSize, -1.0, kWMat + lastIndex*numDofs,
       numDofs, direction->data(), 1, 0.0, y, 1);

 int i;
 for(i = 0; i < kSize; ++i)
   y[i] /= wKw[lastIndex + i];

 // ... Compute direction = direction + W y
 Tgemv('N', numDofs, kSize, 1.0, wMat + lastIndex*numDofs,
       numDofs, y, 1, 1.0, direction->data(), 1);
}
