#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>
#include <Feti.d/GCROrthoSet.h>
#include <Timers.d/DistTimer.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <iostream>

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
GenGCROrthoSet<Scalar>::GenGCROrthoSet(int _len, int maxSize, FSCommunicator *_fetiCom)
{
  this->len = _len;
  this->maxP = maxSize;
  this->fetiCom = _fetiCom;
  this->numP = 0;
  allFPFiP = new Scalar[this->maxP];

  this->numTasks = threadManager->numThr();
  this->oos = new TaskDescr*[this->numTasks];
  int locLen, idx=0;
  int avl = this->len/this->numTasks;
  int mod = this->len % this->numTasks;
  for(int i = 0; i < this->numTasks; ++i) {
    locLen = (i < mod) ? avl+1 : avl;
    this->oos[i] = new GCROp<Scalar>(this, idx, locLen);
    idx += locLen;
  }
}

template<class Scalar>
GenGCROrthoSet<Scalar>::~GenGCROrthoSet()
{
  if(allFPFiP) { delete [] allFPFiP; allFPFiP=0; }
  if(this->oos) {
    for(int i=0; i<this->numTasks; ++i)
      delete (GCROp<Scalar>*) this->oos[i];
    delete [] this->oos;
    this->oos=0;
  }
}

template<class Scalar>
void
GenGCROrthoSet<Scalar>::reset()
{
 this->numP        = 0;

 operation = &GCROp<Scalar>::reset;
 threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
void
GenGCROrthoSet<Scalar>::orthoAdd(Scalar *p, Scalar *Fp, Scalar FpFp)
{
  if(this->numP == this->maxP) this->numP--;
  allFPFiP[this->numP] = FpFp;
  this->numP ++;

  this->op1 = p;
  this->op2 = Fp;
  operation = &GCROp<Scalar>::addVec;
  threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
void
GenGCROrthoSet<Scalar>::orthoAddTimed(DistTimer &timer, Scalar *p, Scalar *Fp, 
                                      Scalar FpFp)
{
  double initTime = getTime();
  long initMem  = memoryUsed();

  if(this->numP == this->maxP) this->numP--;
  allFPFiP[this->numP] = FpFp;
  this->numP++;

  this->op1 = p;
  this->op2 = Fp;
  operation = &GCROp<Scalar>::addVec;
  threadManager->execTimedParal(timer,this->numTasks, this->oos);

  timer.addOverAll(memoryUsed()-initMem, getTime()-initTime);
}

template<class Scalar>
void
GenGCROrthoSet<Scalar>::orthogonalize(Scalar *r, Scalar *Fr, Scalar *p, Scalar *Fp)
{
  int i;
  Scalar *y = (Scalar *)dbg_alloca(this->numP*sizeof(Scalar));
   
  this->op2 = y;
  this->op1 = Fr;
  this->op3 = p;
  for(i = 0; i < this->numP; ++i) y[i] = 0.0;
  operation = &GCROp<Scalar>::Fdot;
  threadManager->execParal(this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom) this->fetiCom->globalSum(this->numP, y);
#endif

  for(i = 0; i < this->numP; ++i)
    y[i] /= allFPFiP[i];

  this->op1 = r;
  this->op3 = p;
  operation = &GCROp<Scalar>::multAdd;
  threadManager->execParal(this->numTasks, this->oos);

  this->op1 = Fr;
  this->op3 = Fp;
  operation = &GCROp<Scalar>::multFAdd;
  threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
void
GenGCROrthoSet<Scalar>::orthogonalizeTimed(DistTimer &timer, Scalar *r, 
                                           Scalar *Fr, Scalar *p, Scalar *Fp)
{
  double initTime = getTime();
  long initMem  = memoryUsed();

  int i;
  Scalar *y = (Scalar *) dbg_alloca(this->numP*sizeof(Scalar));

  this->op2 = y;
  this->op1 = Fr;
  this->op3 = p;
  for(i = 0; i < this->numP; ++i) y[i] = 0.0;
  operation = &GCROp<Scalar>::Fdot;
  threadManager->execTimedParal(timer,this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom) this->fetiCom->globalSum(this->numP, y);
#endif

  for(i = 0; i < this->numP; ++i)
    y[i] /= allFPFiP[i];

  this->op1 = r;
  this->op3 = p;
  operation = &GCROp<Scalar>::multAdd;
  threadManager->execTimedParal(timer,this->numTasks, this->oos);

  this->op1 = Fr;
  this->op3 = Fp;
  operation = &GCROp<Scalar>::multFAdd;
  threadManager->execTimedParal(timer,this->numTasks, this->oos);

  timer.addOverAll(memoryUsed()-initMem, getTime()-initTime);
}

template<class Scalar>
void
GenGCROrthoSet<Scalar>::predict(Scalar *r, Scalar *lambda)
{
  int i;
  Scalar *y = (Scalar *) dbg_alloca(this->numP*sizeof(Scalar));

  this->op2 = y;
  this->op1 = r; 
  this->op3 = lambda; 
  for(i = 0; i < this->numP; ++i) ScalarTypes::initScalar(y[i],0.0);
  operation = &GCROp<Scalar>::Fdot;
  threadManager->execParal(this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom) this->fetiCom->globalSum(this->numP, y);
#endif

  for(i = 0; i < this->numP; ++i)
    y[i] /= allFPFiP[i];

  operation = &GCROp<Scalar>::mult;
  threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
GCROp<Scalar>::GCROp(GenGCROrthoSet<Scalar> *_os, int _index, int _loclen)
{
  os = _os; this->idx = _index; this->loclen = _loclen;
  this->locAllP = locAllFiP = 0;
  this->numP = 0;
}

template<class Scalar>
GCROp<Scalar>::~GCROp()
{
 if(this->locAllP) { delete [] this->locAllP; this->locAllP = 0; }
 // don't delete locAllFiP (this points inside this->locAllP)
}

template<class Scalar>
void
GCROp<Scalar>::run()
{
  (this->*os->operation)();
}

template<class Scalar>
void
GCROp<Scalar>::addVec()
{
  if(this->numP == os->maxP)
    this->numP--;
  if(this->locAllP == 0) {
    //cerr << "allocating memory for orthoset, maxP = " << os->maxP << ", loclen = " << this->loclen << endl;
    this->locAllP   = new Scalar[2*os->maxP*this->loclen];
    locAllFiP = this->locAllP + os->maxP*this->loclen;
  }
  for(int i = 0; i < this->loclen; ++ i) {
    this->locAllP[i+ this->numP*this->loclen]   = os->op1[i+this->idx];
    locAllFiP[i+ this->numP*this->loclen] = os->op2[i+this->idx];
  }
  this->numP++;
}

template<class Scalar>
void
GCROp<Scalar>::Fdot()
{
  int i;
  char trans = 'C';
  Scalar *y = (Scalar *) dbg_alloca(os->maxP*sizeof(Scalar));

  Tgemv(trans, this->loclen, this->numP, Scalar(-1.0), locAllFiP,
        this->loclen, os->op1+this->idx, 1, Scalar(0.0), y, 1);
  os->lock.lock();
  for(i = 0; i < os->numP; ++i)
     os->op2[i] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
GCROp<Scalar>::Fdot_mod()
{
/*
  int i;
  char trans = 'C';
  Scalar *y = (Scalar *) dbg_alloca(os->maxP*sizeof(Scalar));

  Tgemv(trans, this->loclen, this->numP, Scalar(-1.0), locAllFiP,
        this->loclen, os->op1+this->idx, 1, Scalar(0.0), y, 1);
  os->lock.lock();
  for(i = 0; i < os->numP; ++i)
     os->op2[i] += ScalarTypes::conj(y[i]);
  os->lock.unlock();
*/
}


template<class Scalar>
void
GCROp<Scalar>::dot()
{
  int i;
  char trans = 'C';
  Scalar *y = (Scalar *) dbg_alloca(os->maxP*sizeof(Scalar));
  Tgemv(trans, this->loclen, this->numP, Scalar(-1.0), this->locAllP,
        this->loclen, os->op1+this->idx, 1, Scalar(0.0), y, 1);
  os->lock.lock();
  for(i = 0; i < os->numP; ++i)
     os->op2[i] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
GCROp<Scalar>::mult()
{
  char trans = 'N';
  Tgemv(trans, this->loclen, this->numP, Scalar(1.0), this->locAllP,
        this->loclen,  os->op2, 1, Scalar(1.0), os->op3+this->idx, 1);
}

template<class Scalar>
void
GCROp<Scalar>::multAdd()
{
  char trans = 'N';
  for(int i=0; i< this->loclen; ++i)
    os->op3[this->idx+i] = os->op1[this->idx+i];
  Tgemv(trans, this->loclen, this->numP, Scalar(1.0), this->locAllP,
        this->loclen,  os->op2, 1, Scalar(1.0), os->op3+this->idx, 1);
}

template<class Scalar>
void
GCROp<Scalar>::multFAdd()
{
  char trans = 'N';
  for(int i=0; i< this->loclen; ++i)
    os->op3[this->idx+i] = os->op1[this->idx+i];
  Tgemv(trans, this->loclen, this->numP, Scalar(1.0), locAllFiP,
        this->loclen,  os->op2, 1, Scalar(1.0), os->op3+this->idx, 1);
}

