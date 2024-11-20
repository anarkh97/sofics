#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>
#include <Feti.d/GMRESOrthoSet.h>
#include <Timers.d/DistTimer.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <iostream>

#define H(i,j) matrixH[(j)*(this->maxP+1) + (i)]
//#define v(i) (*allV[i])
//#define w(i) (*allW[i])

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
void GenGMRESOrthoSet<Scalar>::generateRotation(Scalar a, Scalar b, 
                                                Scalar &cs, Scalar &ss)
{
/*
  double temp = sqrt(ScalarTypes::sqNorm(a) + ScalarTypes::sqNorm(b));
  cs = ScalarTypes::conj(a) / temp;
  ss = b / temp;
*/
  double temp = sqrt(ScalarTypes::sqNorm(a) + ScalarTypes::sqNorm(b));
  cs = a / temp;
  ss = ScalarTypes::conj(b) / temp;
}

template<class Scalar>
void 
GenGMRESOrthoSet<Scalar>::applyRotation(Scalar &a, Scalar &b, 
                                        Scalar cs, Scalar ss)
{
/*
  Scalar temp1 = cs * a + ScalarTypes::conj(ss) *b ;
  Scalar temp2 = -ss* a + ScalarTypes::conj(cs) * b;
  a = temp1;
  b = temp2;
*/
  Scalar temp1 = ScalarTypes::conj(cs) * a + ss*b ;
  Scalar temp2 = -ScalarTypes::conj(ss) * a + cs * b;
  a = temp1;
  b = temp2;
}

template<class Scalar>
GenGMRESOrthoSet<Scalar>::GenGMRESOrthoSet(int _len, int maxSize, FSCommunicator *_fetiCom) 
{
  this->len = _len;
  this->numP = 0;
  this->fetiCom = _fetiCom;
  this->maxP = maxSize;

#ifdef DETERMINISTIC
 this->numTasks = 1;
#else
  this->numTasks = threadManager->numThr();
#endif

  this->oos = new TaskDescr*[this->numTasks];
  int locLen, idx=0;
  int avl = this->len/this->numTasks;
  int mod = this->len % this->numTasks;
  for(int i = 0; i < this->numTasks; ++i) {
    locLen = (i < mod) ? avl+1 : avl;
    this->oos[i] = new GMRESOp<Scalar>(this, idx, locLen);
    idx += locLen;
  }
  Scalar *p = new Scalar[4*(this->maxP+1)];
  givensC = p;
  givensS = p + this->maxP+1;
  g = p + 2 * (this->maxP+1);
  y = p + 3 * (this->maxP+1);
  matrixH = new Scalar[(this->maxP+1)*this->maxP];
  reInit();
}

template<class Scalar>
GenGMRESOrthoSet<Scalar>::~GenGMRESOrthoSet() 
{
  if(givensC) { delete [] givensC; givensC = 0; }
  if(matrixH) { delete [] matrixH; matrixH = 0; }
  // if(givensS) { delete [] givensS; givensS = 0; }
  // if(g) { delete [] g; g = 0; }
  // if(y) { delete [] y; y = 0; }
  if(this->oos) {
    int i;
    for(i=0; i<this->numTasks; ++i)
      delete (GMRESOp<Scalar>*) this->oos[i];
    delete [] this->oos;
    this->oos=0;
  }
}

template<class Scalar>
void 
GenGMRESOrthoSet<Scalar>::reInit() 
{
  this->numP = 0;
  for(int i=0; i < this->maxP+1 ; i++) { 
    g[i] = (0.0); 
    givensC[i] = Scalar(0.0); 
    givensS[i] = Scalar(0.0); 
    for(int j=0; j < this->maxP ; j++) H(i,j) = Scalar(0.0); 
  }

  operation = &GMRESOp<Scalar>::reInit;
  threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
void 
GenGMRESOrthoSet<Scalar>::init(Scalar *v0, double beta) 
{
  g[0] = beta;
  this->op1 = v0;
  operation = &GMRESOp<Scalar>::addVec;
  threadManager->execParal(this->numTasks, this->oos);
  this->numP++;
}

template<class Scalar>
double 
GenGMRESOrthoSet<Scalar>::orthoAdd(Scalar *Fv, Scalar *v) 
{
  int i;
  if(this->numP == this->maxP+1) { 
    fprintf(stderr, "Maximum number of directions exceeeded in GMRES.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  }

  this->op1 = Fv;
  this->op2 = y;
  for(i = 0; i < this->numP; ++i)
    y[i] = 0;
  operation = &GMRESOp<Scalar>::dot;
  threadManager->execParal(this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom)
    this->fetiCom->globalSum(this->numP, y);
#endif

  for(i = 0; i < this->numP; ++i) 
    H(i,this->numP-1) = y[i];

  this->op3 = v;
  operation = &GMRESOp<Scalar>::multAdd;
  threadManager->execParal(this->numTasks, this->oos);

  double nrm=0.0;
  for(i = 0; i < this->len; ++i) nrm += ScalarTypes::Real(ScalarTypes::conj(v[i])*v[i]);
#ifdef DISTRIBUTED
  nrm = this->fetiCom->globalSum(nrm);
#endif
  H(this->numP,this->numP-1) = Scalar(sqrt(nrm));

  this->op2 = & H(this->numP,this->numP-1);
  this->op1 = v;
  operation = &GMRESOp<Scalar>::addVecAndNorm;
  threadManager->execParal(this->numTasks, this->oos);

  for(i = 0; i < this->numP-1; i++)
    applyRotation(H(i,this->numP-1), H(i+1,this->numP-1), givensC[i], givensS[i]);
  generateRotation(H(this->numP-1,this->numP-1), H(this->numP,this->numP-1), givensC[this->numP-1], givensS[this->numP-1]);
  applyRotation(H(this->numP-1,this->numP-1), H(this->numP,this->numP-1), givensC[this->numP-1], givensS[this->numP-1]);
  applyRotation(g[this->numP-1], g[this->numP], givensC[this->numP-1], givensS[this->numP-1]);

  this->numP++;
  return ScalarTypes::norm(g[this->numP-1]);
}


template<class Scalar>
double 
GenGMRESOrthoSet<Scalar>::orthoAddTimed(DistTimer &timer, Scalar *Fv, Scalar *v) 
{
  double initTime = getTime();
  long initMem  = memoryUsed();

  if(this->len == 0) return 0.0;

  int i;
  if (this->numP == this->maxP+1) {
    fprintf(stderr, "Maximum number of directions exceeeded in GMRES.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  }

  if(this->len == 0) { this->numP++; return 0.0; }

  this->op1 = Fv;
  this->op2 = y;
  for(i = 0; i < this->numP; ++i)
    y[i] = 0;
  operation = &GMRESOp<Scalar>::dot;
  threadManager->execTimedParal(timer, this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if (this->fetiCom)
    this->fetiCom->globalSum(this->numP, y);
#endif

  for(i = 0; i < this->numP; ++i)
    H(i,this->numP-1) = y[i];

#ifdef DEBUGORTHO
  double dotNorm = 0.0;
  for(i = 0; i < this->numP; ++i)
    dotNorm += ScalarTypes::sqNorm(y[i]);
  fprintf(stderr, "DotNorm: %e\n", dotNorm);
#endif
  
  this->op3 = v;
  operation = &GMRESOp<Scalar>::multAdd;
  threadManager->execTimedParal(timer, this->numTasks, this->oos);

  double nrm=0.0;
  for(i = 0; i < this->len; ++i) nrm += ScalarTypes::Real(ScalarTypes::conj(v[i])*v[i]);
#ifdef DISTRIBUTED
  nrm = this->fetiCom->globalSum(nrm);
#endif
  H(this->numP,this->numP-1) = Scalar(sqrt(nrm));

  this->op2 = & H(this->numP,this->numP-1);
  this->op1 = v;
  operation = &GMRESOp<Scalar>::addVecAndNorm;
  threadManager->execTimedParal(timer, this->numTasks, this->oos);

  for(i = 0; i < this->numP-1; i++)
    applyRotation(H(i,this->numP-1), H(i+1,this->numP-1), givensC[i], givensS[i]);
  generateRotation(H(this->numP-1,this->numP-1), H(this->numP,this->numP-1), givensC[this->numP-1], givensS[this->numP-1]);
  applyRotation(H(this->numP-1,this->numP-1), H(this->numP,this->numP-1), givensC[this->numP-1], givensS[this->numP-1]);
  applyRotation(g[this->numP-1], g[this->numP], givensC[this->numP-1], givensS[this->numP-1]);

#ifdef DEBUGORTHO
  //XMLDEBUG
  this->op1 = v;
  this->op2 = y;
  for(i = 0; i < this->numP; ++i)
    y[i] = Scalar(0.0);
  operation = &GMRESOp<Scalar>::dot;
  threadManager->execTimedParal(timer, this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom)
    this->fetiCom->globalSum(this->numP, y);
#endif

  for(i = 0; i < this->numP; ++i) 
    H(i,this->numP-1) = y[i];
  
  dotNorm = 0.0;
  for(i = 0; i < this->numP; ++i)
    dotNorm += ScalarTypes::sqNorm(y[i]);
  fprintf(stderr, "After DotNorm: %e and %e\n", dotNorm, ScalarTypes::sqNorm(y[this->numP]));
  // END DEBUG
#endif  
  this->numP++;

  timer.addOverAll( memoryUsed()-initMem, getTime()-initTime );
  return ScalarTypes::norm(g[this->numP-1]);
}

template<class Scalar>
double
GenGMRESOrthoSet<Scalar>::ModorthoAddTimed(DistTimer &timer, Scalar *Fv, Scalar *v)
{
  double initTime = getTime();
  long initMem  = memoryUsed();
  int i;
  if (this->numP == this->maxP+1) {
    fprintf(stderr, "Maximum number of directions exceeeded in GMRES.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  }
  operation = &GMRESOp<Scalar>::reInitcurrent;
  threadManager->execParal(this->numTasks, this->oos);
  this->op1 = Fv;
  Scalar dotvalue[1];
  this->op2 = dotvalue;
  for(i = 0; i < this->numP; ++i) {
    dotvalue[0] = 0;
    operation = &GMRESOp<Scalar>::takecurrent;
    threadManager->execTimedParal(timer, this->numTasks, this->oos);
    operation = &GMRESOp<Scalar>::Moddot;
    threadManager->execTimedParal(timer, this->numTasks, this->oos);
#ifdef DISTRIBUTED
    if (this->fetiCom)
      this->fetiCom->globalSum(1, dotvalue);
#endif
    y[i] = dotvalue[0];
    H(i,this->numP-1) = y[i];
    this->op3 = v;
    operation = &GMRESOp<Scalar>::ModmultAdd;
    threadManager->execTimedParal(timer, this->numTasks, this->oos);
    this->op1 = v;
  }

  double nrm=0.0;
  for(i = 0; i < this->len; ++i) nrm += ScalarTypes::Real(ScalarTypes::conj(v[i])*v[i]);
#ifdef DISTRIBUTED
  nrm = this->fetiCom->globalSum(nrm);
#endif
 
  H(this->numP,this->numP-1) = Scalar(sqrt(nrm));
 
  this->op2 = & H(this->numP,this->numP-1);
  this->op1 = v;
  operation = &GMRESOp<Scalar>::addVecAndNorm;
  threadManager->execTimedParal(timer, this->numTasks, this->oos);
 
  for(i = 0; i < this->numP-1; i++)
    applyRotation(H(i,this->numP-1), H(i+1,this->numP-1), givensC[i], givensS[i]);
  generateRotation(H(this->numP-1,this->numP-1), H(this->numP,this->numP-1), givensC[this->numP-1], givensS[this->numP-1]);
  applyRotation(H(this->numP-1,this->numP-1), H(this->numP,this->numP-1), givensC[this->numP-1], givensS[this->numP-1]);
  applyRotation(g[this->numP-1], g[this->numP], givensC[this->numP-1], givensS[this->numP-1]);
 
  this->numP++;
 
  timer.addOverAll( memoryUsed()-initMem, getTime()-initTime );
  return ScalarTypes::norm(g[this->numP-1]);
}

template<class Scalar>
void GenGMRESOrthoSet<Scalar>::solution(Scalar *u) 
{
  int i,j;
  //Back substitute
  for(i=this->numP-2;i>=0;i--) {
    y[i] = g[i];
    for(j=this->numP-2; j>i;j--) y[i] -= H(i,j) * y[j]; 
    y[i] /= H(i,i);
  }
  this->op1 = u;
  this->op2 = y;

  operation = &GMRESOp<Scalar>::mult;  
  threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
GMRESOp<Scalar>::GMRESOp(GenGMRESOrthoSet<Scalar> *_os, int _index, int _loclen)
{
  os = _os; this->idx = _index; this->loclen = _loclen;
  this->locAllP = 0;
  this->numP = 0;
  currentindex = 0;
  currentvector = 0;
}

template<class Scalar>
void GMRESOp<Scalar>::reInit() 
{
 this->numP = 0;
}

template<class Scalar>
void GMRESOp<Scalar>::reInitcurrent() 
{
 currentindex = 0;
}

template<class Scalar>
void
GMRESOp<Scalar>::run()
{
 (this->*os->operation)();
}

template<class Scalar>
void
GMRESOp<Scalar>::addVec()
{
 if(this->locAllP == 0) {
   this->locAllP   = new Scalar[(os->maxP+1)*this->loclen];
   currentvector = new Scalar[this->loclen];
 }
 for(int i = 0; i < this->loclen; ++ i) {
   this->locAllP[i+ this->numP*this->loclen]   = os->op1[i+this->idx];
   currentvector[i] = os->op1[i+this->idx];
 }
 
 this->numP++;
}

template<class Scalar>
void
GMRESOp<Scalar>::addVecAndNorm()
{
 int i; 
 Scalar nrm = os->op2[0];
 for(i = 0; i < this->loclen; ++ i) {
   this->locAllP[i+ this->numP*this->loclen]   = os->op1[i+this->idx] / nrm;
   currentvector[i] = os->op1[i+this->idx] / nrm;
   os->op1[i+this->idx] = this->locAllP[i+ this->numP*this->loclen];
 }

 this->numP++;
}

template<class Scalar>
void
GMRESOp<Scalar>::takecurrent()
{
 int i;
 for(i = 0; i < this->loclen; ++ i)
   currentvector[i] = this->locAllP[i+ currentindex*this->loclen];
 currentindex++;
}

template<class Scalar>
void
GMRESOp<Scalar>::dot()
{
  if(this->loclen == 0) return;
  /* dot computes V^H w, so the dot of (w,v_i) */
  int i;
  char trans = 'C';
  Scalar *y = (Scalar *)dbg_alloca(os->maxP*sizeof(Scalar));

  Tgemv(trans, this->loclen, this->numP, Scalar(1.0), this->locAllP,
        this->loclen, os->op1+this->idx, 1, Scalar(0.0), y, 1);
  os->lock.lock();
  for(i = 0; i < this->numP; ++i)
    os->op2[i] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
GMRESOp<Scalar>::Moddot()
{
  if(this->loclen == 0) return;
  int i;
  char trans = 'C';
  Scalar y[1];
  Tgemv(trans, this->loclen, 1, Scalar(1.0), currentvector,
        this->loclen, os->op1+this->idx, 1, Scalar(0.0), y, 1);
  os->lock.lock();
  for(i = 0; i < 1; ++i)
     os->op2[i] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
GMRESOp<Scalar>::mult()
{
  if(this->loclen == 0) return;
  char trans = 'N';
 
  Tgemv(trans, this->loclen, this->numP-1, Scalar(-1.0), this->locAllP,
        this->loclen,  os->op2, 1, Scalar(0.0), os->op1+this->idx, 1);
}

template<class Scalar>
void
GMRESOp<Scalar>::multAdd()
{
  if(this->loclen == 0) return;
  char trans = 'N';

  for(int i=0; i< this->loclen; ++i)
    os->op3[this->idx+i] = os->op1[this->idx+i];
  Tgemv(trans, this->loclen, this->numP, Scalar(-1.0), this->locAllP,
        this->loclen,  os->op2, 1, Scalar(1.0), os->op3+this->idx, 1);
}

template<class Scalar>
void
GMRESOp<Scalar>::ModmultAdd()
{
  if(this->loclen == 0) return;
  char trans = 'N';

  for(int i=0; i< this->loclen; ++i)
    os->op3[this->idx+i] = os->op1[this->idx+i];
  Tgemv(trans, this->loclen, 1, Scalar(-1.0), currentvector,
        this->loclen,  os->op2, 1, Scalar(1.0), os->op3+this->idx, 1);
}

template<class Scalar>
GMRESOp<Scalar>::~GMRESOp()
{
  if(this->locAllP) { delete [] this->locAllP; this->locAllP = 0; }
  if(currentvector) { delete [] currentvector; currentvector = 0; }
}
