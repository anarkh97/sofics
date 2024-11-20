#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>
#include <Solvers.d/GmresOrthoSet.h>
#include <Timers.d/DistTimer.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <iostream>

#define _H(i,j) matrixH[(j)*(maxP+1) + (i)]
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
void GmresOrthoSet<Scalar>::generateRotation(Scalar a, Scalar b, 
                                             Scalar &cs, Scalar &ss)
{
  double temp = sqrt(ScalarTypes::sqNorm(a) + ScalarTypes::sqNorm(b));
  cs = a / temp;
  ss = ScalarTypes::conj(b) / temp;
}

template<class Scalar>
void 
GmresOrthoSet<Scalar>::applyRotation(Scalar &a, Scalar &b, 
                                     Scalar cs, Scalar ss)
{
  Scalar temp1 = ScalarTypes::conj(cs) * a + ss*b ;
  Scalar temp2 = -ScalarTypes::conj(ss) * a + cs * b;
  a = temp1;
  b = temp2;
}

template<class Scalar>
GmresOrthoSet<Scalar>::GmresOrthoSet(int _len, int maxSize, FSCommunicator *_com) 
{
  len = _len;
  numP = 0;
  com = _com;
  maxP = maxSize;

#ifdef DETERMINISTIC
 numTasks = 1;
#else
  numTasks = threadManager->numThr();
#endif

  oos = new TaskDescr*[numTasks];
  int locLen, idx=0;
  int avl = len/numTasks;
  int mod = len % numTasks;
  for(int i = 0; i < numTasks; ++i) {
    locLen = (i < mod) ? avl+1 : avl;
    oos[i] = new GmresOp<Scalar>(this, idx, locLen);
    idx += locLen;
  }
  Scalar *p = new Scalar[4*(maxP+1)];
  givensC = p;
  givensS = p + maxP+1;
  g = p + 2 * (maxP+1);
  y = p + 3 * (maxP+1);
  matrixH = new Scalar[(maxP+1)*maxP];
  reset();
}

template<class Scalar>
GmresOrthoSet<Scalar>::~GmresOrthoSet() 
{
  if(givensC) { delete [] givensC; givensC = 0; }
  if(matrixH) { delete [] matrixH; matrixH = 0; }
  // if(givensS) { delete [] givensS; givensS = 0; }
  // if(g) { delete [] g; g = 0; }
  // if(y) { delete [] y; y = 0; }
  if(oos) {
    int i;
    for(i=0; i<numTasks; ++i)
      delete (GmresOp<Scalar>*) oos[i];
    delete [] oos;
    oos=0;
  }
}

template<class Scalar>
void 
GmresOrthoSet<Scalar>::reset() 
{
  numP = 0;
  for(int i=0; i < maxP+1 ; i++) { 
    g[i] = (0.0); 
    givensC[i] = Scalar(0.0); 
    givensS[i] = Scalar(0.0); 
    for(int j=0; j < maxP ; j++) _H(i,j) = Scalar(0.0); 
  }

  operation = &GmresOp<Scalar>::reset;
  threadManager->execParal(numTasks, oos);
}

template<class Scalar>
void 
GmresOrthoSet<Scalar>::init(Scalar *v0, double beta) 
{
  g[0] = beta;
  op1 = v0;
  operation = &GmresOp<Scalar>::addVec;
  threadManager->execParal(numTasks, oos);
  numP++;
}

template<class Scalar>
double
GmresOrthoSet<Scalar>::init(Scalar *p)
{
  if(numP > 0) reset();

  double beta = 0.0;
  for(int i = 0; i < len; i++) beta += ScalarTypes::sqNorm(p[i]);
#ifdef DISTRIBUTED
  if(com)
    beta = com->globalSum(beta);
#endif
  beta = sqrt(beta);
  for(int i = 0; i < len; i++) p[i] /= beta;

  init(p, beta);

  return beta;
}

template<class Scalar>
double 
GmresOrthoSet<Scalar>::orthoAdd(Scalar *Fv, Scalar *v) 
{
  int i;
  if(numP == maxP+1) { 
    fprintf(stderr, "Maximum number of directions exceeeded in Gmres.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  }

  op1 = Fv;
  op2 = y;
  for(i = 0; i < numP; ++i)
    y[i] = 0;
  operation = &GmresOp<Scalar>::dot;
  threadManager->execParal(numTasks, oos);

#ifdef DISTRIBUTED
  if(com)
    com->globalSum(numP, y);
#endif

  for(i = 0; i < numP; ++i) 
    _H(i,numP-1) = y[i];

  op3 = v;
  operation = &GmresOp<Scalar>::multAdd;
  threadManager->execParal(numTasks, oos);

  double nrm=0.0;
  for(i = 0; i < len; ++i) nrm += ScalarTypes::Real(ScalarTypes::conj(v[i])*v[i]);
#ifdef DISTRIBUTED
  nrm = com->globalSum(nrm);
#endif
  _H(numP,numP-1) = Scalar(sqrt(nrm));

  op2 = & _H(numP,numP-1);
  op1 = v;
  operation = &GmresOp<Scalar>::addVecAndNorm;
  threadManager->execParal(numTasks, oos);

  for(i = 0; i < numP-1; i++)
    applyRotation(_H(i,numP-1), _H(i+1,numP-1), givensC[i], givensS[i]);
  generateRotation(_H(numP-1,numP-1), _H(numP,numP-1), givensC[numP-1], givensS[numP-1]);
  applyRotation(_H(numP-1,numP-1), _H(numP,numP-1), givensC[numP-1], givensS[numP-1]);
  applyRotation(g[numP-1], g[numP], givensC[numP-1], givensS[numP-1]);

  numP++;
  return ScalarTypes::norm(g[numP-1]);
}


/*
template<class Scalar>
double 
GmresOrthoSet<Scalar>::orthoAddTimed(DistTimer &timer, Scalar *Fv, Scalar *v) 
{
  double initTime = getTime();
  long initMem  = memoryUsed();

  if(len == 0) return 0.0;

  int i;
  if (numP == maxP+1) {
    fprintf(stderr, "Maximum number of directions exceeeded in Gmres.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  }

  if(len == 0) { numP++; return 0.0; }

  op1 = Fv;
  op2 = y;
  for(i = 0; i < numP; ++i)
    y[i] = 0;
  operation = &GmresOp<Scalar>::dot;
  threadManager->execTimedParal(timer, numTasks, oos);

#ifdef DISTRIBUTED
  if (com)
    com->globalSum(numP, y);
#endif

  for(i = 0; i < numP; ++i)
    _H(i,numP-1) = y[i];

#ifdef DEBUGORTHO
  double dotNorm = 0.0;
  for(i = 0; i < numP; ++i)
    dotNorm += ScalarTypes::sqNorm(y[i]);
  fprintf(stderr, "DotNorm: %e\n", dotNorm);
#endif
  
  op3 = v;
  operation = &GmresOp<Scalar>::multAdd;
  threadManager->execTimedParal(timer, numTasks, oos);

  double nrm=0.0;
  for(i = 0; i < len; ++i) nrm += ScalarTypes::Real(ScalarTypes::conj(v[i])*v[i]);
#ifdef DISTRIBUTED
  nrm = com->globalSum(nrm);
#endif
  _H(numP,numP-1) = Scalar(sqrt(nrm));

  op2 = & _H(numP,numP-1);
  op1 = v;
  operation = &GmresOp<Scalar>::addVecAndNorm;
  threadManager->execTimedParal(timer, numTasks, oos);

  for(i = 0; i < numP-1; i++)
    applyRotation(_H(i,numP-1), _H(i+1,numP-1), givensC[i], givensS[i]);
  generateRotation(_H(numP-1,numP-1), _H(numP,numP-1), givensC[numP-1], givensS[numP-1]);
  applyRotation(_H(numP-1,numP-1), _H(numP,numP-1), givensC[numP-1], givensS[numP-1]);
  applyRotation(g[numP-1], g[numP], givensC[numP-1], givensS[numP-1]);

#ifdef DEBUGORTHO
  //XMLDEBUG
  op1 = v;
  op2 = y;
  for(i = 0; i < numP; ++i)
    y[i] = Scalar(0.0);
  operation = &GmresOp<Scalar>::dot;
  threadManager->execTimedParal(timer, numTasks, oos);

#ifdef DISTRIBUTED
  if(com)
    com->globalSum(numP, y);
#endif

  for(i = 0; i < numP; ++i) 
    _H(i,numP-1) = y[i];
  
  dotNorm = 0.0;
  for(i = 0; i < numP; ++i)
    dotNorm += ScalarTypes::sqNorm(y[i]);
  fprintf(stderr, "After DotNorm: %e and %e\n", dotNorm, ScalarTypes::sqNorm(y[numP]));
  // END DEBUG
#endif  
  numP++;

  timer.addOverAll( memoryUsed()-initMem, getTime()-initTime );
  return ScalarTypes::norm(g[numP-1]);
}
*/

template<class Scalar>
double
GmresOrthoSet<Scalar>::ModorthoAdd(Scalar *Fv, Scalar *v)
{
  //double initTime = getTime();
  //long initMem  = memoryUsed();
  int i;
  if (numP == maxP+1) {
    fprintf(stderr, "Maximum number of directions exceeeded in Gmres.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  }
  operation = &GmresOp<Scalar>::reInitcurrent;
  threadManager->execParal(numTasks, oos);
  op1 = Fv;
  Scalar dotvalue[1];
  op2 = dotvalue;
  for(i = 0; i < numP; ++i) {
    dotvalue[0] = 0;
    operation = &GmresOp<Scalar>::takecurrent;
    threadManager->execParal(numTasks, oos);
    operation = &GmresOp<Scalar>::Moddot;
    threadManager->execParal(numTasks, oos);
#ifdef DISTRIBUTED
    if (com)
      com->globalSum(1, dotvalue);
#endif
    y[i] = dotvalue[0];
    _H(i,numP-1) = y[i];
    op3 = v;
    operation = &GmresOp<Scalar>::ModmultAdd;
    threadManager->execParal(numTasks, oos);
    op1 = v;
  }

  double nrm=0.0;
  for(i = 0; i < len; ++i) nrm += ScalarTypes::Real(ScalarTypes::conj(v[i])*v[i]);
#ifdef DISTRIBUTED
  if(com) nrm = com->globalSum(nrm);
#endif
 
  _H(numP,numP-1) = Scalar(sqrt(nrm));
 
  op2 = & _H(numP,numP-1);
  op1 = v;
  operation = &GmresOp<Scalar>::addVecAndNorm;
  threadManager->execParal(numTasks, oos);
 
  for(i = 0; i < numP-1; i++)
    applyRotation(_H(i,numP-1), _H(i+1,numP-1), givensC[i], givensS[i]);
  generateRotation(_H(numP-1,numP-1), _H(numP,numP-1), givensC[numP-1], givensS[numP-1]);
  applyRotation(_H(numP-1,numP-1), _H(numP,numP-1), givensC[numP-1], givensS[numP-1]);
  applyRotation(g[numP-1], g[numP], givensC[numP-1], givensS[numP-1]);
 
  numP++;
 
  //timer.addOverAll( memoryUsed()-initMem, getTime()-initTime );
  return ScalarTypes::norm(g[numP-1]);
}

template<class Scalar>
void GmresOrthoSet<Scalar>::solution(Scalar *u) 
{
  int i,j;
  //Back substitute
  for(i=numP-2;i>=0;i--) {
    y[i] = g[i];
    for(j=numP-2; j>i;j--) y[i] -= _H(i,j) * y[j]; 
    y[i] /= _H(i,i);
  }
  op1 = u;
  op2 = y;

  operation = &GmresOp<Scalar>::mult;  
  threadManager->execParal(numTasks, oos);
}

template<class Scalar>
GmresOp<Scalar>::GmresOp(GmresOrthoSet<Scalar> *_os, int _index, int _loclen)
{
  os = _os; idx = _index; loclen = _loclen;
  locAllP = 0;
  numP = 0;
  currentindex = 0;
  currentvector = 0;
}

template<class Scalar>
void GmresOp<Scalar>::reset() 
{
 numP = 0;
}

template<class Scalar>
void GmresOp<Scalar>::reInitcurrent() 
{
 currentindex = 0;
}

template<class Scalar>
void
GmresOp<Scalar>::run()
{
 (this->*os->operation)();
}

template<class Scalar>
void
GmresOp<Scalar>::addVec()
{
 if(locAllP == 0) {
   locAllP   = new Scalar[(os->maxP+1)*loclen];
   currentvector = new Scalar[loclen];
 }
 for(int i = 0; i < loclen; ++ i) {
   locAllP[i+ numP*loclen]   = os->op1[i+idx];
   currentvector[i] = os->op1[i+idx];
 }
 
 numP++;
}

template<class Scalar>
void
GmresOp<Scalar>::addVecAndNorm()
{
 int i; 
 Scalar nrm = os->op2[0];
 for(i = 0; i < loclen; ++ i) {
   locAllP[i+ numP*loclen]   = os->op1[i+idx] / nrm;
   currentvector[i] = os->op1[i+idx] / nrm;
   os->op1[i+idx] = locAllP[i+ numP*loclen];
 }

 numP++;
}

template<class Scalar>
void
GmresOp<Scalar>::takecurrent()
{
 int i;
 for(i = 0; i < loclen; ++ i)
   currentvector[i] = locAllP[i+ currentindex*loclen];
 currentindex++;
}

template<class Scalar>
void
GmresOp<Scalar>::dot()
{
  if(loclen == 0) return;
  /* dot computes V^H w, so the dot of (w,v_i) */
  int i;
  char trans = 'C';
  Scalar *y = (Scalar *)dbg_alloca(os->maxP*sizeof(Scalar));

  Tgemv(trans, loclen, numP, Scalar(1.0), locAllP,
        loclen, os->op1+idx, 1, Scalar(0.0), y, 1);
  os->lock.lock();
  for(i = 0; i < numP; ++i)
    os->op2[i] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
GmresOp<Scalar>::Moddot()
{
  if(loclen == 0) return;
  int i;
  char trans = 'C';
  Scalar y[1];
  Tgemv(trans, loclen, 1, Scalar(1.0), currentvector,
        loclen, os->op1+idx, 1, Scalar(0.0), y, 1);
  os->lock.lock();
  for(i = 0; i < 1; ++i)
     os->op2[i] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
GmresOp<Scalar>::mult()
{
  if(loclen == 0) return;
  char trans = 'N';
 
  Tgemv(trans, loclen, numP-1, Scalar(-1.0), locAllP,
        loclen,  os->op2, 1, Scalar(0.0), os->op1+idx, 1);
}

template<class Scalar>
void
GmresOp<Scalar>::multAdd()
{
  if(loclen == 0) return;
  char trans = 'N';

  for(int i=0; i< loclen; ++i)
    os->op3[idx+i] = os->op1[idx+i];
  Tgemv(trans, loclen, numP, Scalar(-1.0), locAllP,
        loclen,  os->op2, 1, Scalar(1.0), os->op3+idx, 1);
}

template<class Scalar>
void
GmresOp<Scalar>::ModmultAdd()
{
  if(loclen == 0) return;
  char trans = 'N';

  for(int i=0; i< loclen; ++i)
    os->op3[idx+i] = os->op1[idx+i];
  Tgemv(trans, loclen, 1, Scalar(-1.0), currentvector,
        loclen,  os->op2, 1, Scalar(1.0), os->op3+idx, 1);
}

template<class Scalar>
GmresOp<Scalar>::~GmresOp()
{
  if(locAllP) { delete [] locAllP; locAllP = 0; }
  if(currentvector) { delete [] currentvector; currentvector = 0; }
}
