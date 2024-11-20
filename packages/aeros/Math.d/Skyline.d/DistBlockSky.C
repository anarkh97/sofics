#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <Math.d/Vector.h>

#include <Comm.d/Communicator.h>
extern SysCom *scom;

#include <Utils.d/linkfc.h>

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
GenDistBlockSky<Scalar>::~GenDistBlockSky()
{
 if(nlines) { delete [] nlines; nlines=0; }
}


template<class Scalar>
GenDistBlockSky<Scalar>::GenDistBlockSky(const Connectivity *cn, EqNumberer *dsa,
                                         double trbm, int fRow, int nRow) :
 GenBlockSky<Scalar>(cn, dsa, trbm)
{
 firstRow = fRow;
 numRows  = nRow;

 // allocate memory for the number of rows we wish to store
 nlines = new Scalar[numRows*this->neqs()];
}

template<class Scalar>
void
GenDistBlockSky<Scalar>::parallelFactor()
{
 GenBlockSky<Scalar>::parallelFactor();

 // zero the n-lines
 int i;
 for(i=0; i<numRows*this->neqs(); ++i)
   nlines[i]=0.0;   

 Scalar **rows = (Scalar **) dbg_alloca(sizeof(Scalar*)*numRows);

 for(i=0; i<numRows; ++i) {
   nlines[i*this->neqs()+firstRow+i] = 1.0;
   rows[i] = nlines+i*this->neqs();
 }

 GenBlockSky<Scalar>::reSolve(numRows, rows);

 // delete entire Skyline matrix 
 delete [] this->skyA;
 this->skyA = 0;
}

template<class Scalar>
void
GenDistBlockSky<Scalar>::reSolve(Scalar *rhs)
{
 Scalar *partialSum = (Scalar *)dbg_alloca(sizeof(Scalar)*numRows);

 Tgemv('T', this->neqs(), numRows, 1.0, nlines, this->neqs(),
       rhs, 1, 0.0, partialSum, 1);

 // zero the rhs
 int i;
 for(i=0; i<this->neqs(); ++i)
   rhs[i] = 0.0;

 for(i=0; i<numRows; ++i)
   rhs[firstRow+i] = partialSum[i];

#ifdef DISTRIBUTED
  structCom->globalSum(this->neqs(), rhs);
#endif
}

template<class Scalar>
void
GenDistBlockSky<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
 Scalar *partialSum = (Scalar *)dbg_alloca(sizeof(Scalar)*numRows);

 Tgemv('T', this->neqs(), numRows, 1.0, nlines, this->neqs(),
       rhs.data(), 1, 0.0, partialSum, 1);

 // zero the rhs
 int i;
 for(i=0; i<this->neqs(); ++i)
   rhs[i] = 0.0;

 for(i=0; i<numRows; ++i)
   rhs[firstRow+i] = partialSum[i];

#ifdef DISTRIBUTED
  structCom->globalSum(this->neqs(), rhs.data());
#endif
}


