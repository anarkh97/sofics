#include <Utils.d/Memory.h>
#include <iostream>
#include <Driver.d/Communicator.h>
#include <Feti.d/CGOrthoSet.h>
#include <Math.d/BLAS.h>
#include <Timers.d/GetTime.h>
#include <Timers.d/DistTimer.h>
#include <Utils.d/dbg_alloca.h>


#ifndef _TGEMV__
#define _TGEMV__
inline void Tgemv(const char &a, const int &b, const int &c,
                  const double &d, double *e, const int &f,
                  double *g, const int &h, const double &i, double *j, const int &k)
{
 _FORTRAN(dgemv)(a,b,c,d,e,f,g,h,i,j,k);
}

inline void Tgemv(const char &a, const int &b, const int &c,
                  const std::complex<double> &d, std::complex<double> *e, const int &f,
                  std::complex<double> *g, const int &h, const std::complex<double> &i, std::complex<double> *j, const int &k)
{
 _FORTRAN(zgemv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif

template<class Scalar>
GenCGOrthoSet<Scalar>::GenCGOrthoSet(int _len, int maxSize, FSCommunicator *_fetiCom)
{
 this->len  = _len;           // length of vectors
 this->maxP = maxSize;        // maximum number of vectors to be stored
 this->fetiCom = _fetiCom;    // PJSA
 this->numP = 0;              // initially zero vectors stored

 numOrthoSet = 1;       // number of ortho sets

 allPFiP = new Scalar[this->maxP];

 // this->maxP is large for the number of system indices we need...
 kindex  = new int[this->maxP+1];
 kindex[0] = 0;
 kindex[1] = 0;

#ifdef DETERMINISTIC
 this->numTasks = 1;
#else
 this->numTasks = threadManager->numThr();
#endif

 this->oos = new TaskDescr*[this->numTasks];
 int locLen, idx=0;
 int avl = this->len/this->numTasks;
 int mod = this->len % this->numTasks;
 int i;
 for(i = 0; i < this->numTasks; ++i) {
    locLen = (i < mod) ? avl+1 : avl;
    this->oos[i] = new CGOrthoOp<Scalar>(this, idx, locLen);
    idx += locLen;
 }
}

template<class Scalar>
GenCGOrthoSet<Scalar>::~GenCGOrthoSet()
{
  if(kindex) { delete [] kindex; kindex=0; } 
  if(allPFiP) { delete [] allPFiP; allPFiP=0; }
  if(this->oos) {
    int i;
    for(i=0; i<this->numTasks; ++i)
    delete (CGOrthoOp<Scalar>*) this->oos[i];
    delete [] this->oos;
    this->oos=0;
  }
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::clean_up()
{
 if(allPFiP) {
   delete [] allPFiP;
   allPFiP=nullptr;
 }
 if(kindex) {
   delete [] kindex;
   kindex= nullptr;
 }
 operation = &CGOrthoOp<Scalar>::clean_up;
 threadManager->execParal(this->numTasks, this->oos);

 if(this->oos) {
   delete [] this->oos;
   this->oos=nullptr;
 }
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::orthoAdd(Scalar *p, Scalar *Fp, Scalar pFp)
{
 if(this->numP == this->maxP) {
   this->numP--;
   kindex[numOrthoSet]--; // KHP: DEBUG
 }

 // Store dot product
 allPFiP[this->numP] = pFp;

 // Increase the number of stored directions
 this->numP += 1;
 kindex[numOrthoSet] += 1;

 this->op1 = p;
 this->op2 = Fp;

 operation = &CGOrthoOp<Scalar>::addVec;
 threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::orthoAddTimed(DistTimer &timer, Scalar *p, Scalar *Fp, Scalar pFp)
{
 double initTime = getTime();
 long initMem  = memoryUsed();

 if(this->numP == this->maxP) {
   this->numP--;
   kindex[numOrthoSet]--; // KHP: DEBUG
 }

 // Store dot product
 allPFiP[this->numP] = pFp;

 // Increase the number of stored directions
 this->numP += 1;
 kindex[numOrthoSet] += 1;

 this->op1 = p;
 this->op2 = Fp;

 operation = &CGOrthoOp<Scalar>::addVec;
 threadManager->execTimedParal(timer, this->numTasks, this->oos);

 timer.addOverAll( memoryUsed()-initMem, getTime()-initTime );
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::shiftOrthoSets()
{
 kSize = kindex[1];

 int i;
 for(i=kSize; i<this->numP; ++i)
   allPFiP[i-kSize] = allPFiP[i];

 operation = &CGOrthoOp<Scalar>::shift;
 threadManager->execParal(this->numTasks, this->oos);

 for(i=1; i<numOrthoSet; ++i)
   kindex[i] = kindex[i+1] - kSize;
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::newOrthoSet()
{
  if(numOrthoSet < this->maxP) { 
    // maximum number of ortho sets not exceeded, add one
    numOrthoSet += 1;
    if(this->numP == this->maxP) { 
      // make sure there is storage room for at least 
      // one more search direction
      kindex[numOrthoSet-1]--;
    }
  }
  // initialise the starting point of the next ortho set
  kindex[numOrthoSet] = kindex[numOrthoSet-1];
  this->numP = kindex[numOrthoSet];
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::reset()
{
 this->numP        = 0;
 kindex[0]   = 0;
 kindex[1]   = 0;
 numOrthoSet = 1;
 
 operation = &CGOrthoOp<Scalar>::reset;
 threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::orthogonalizeTimed(DistTimer &timer, Scalar *r, Scalar *p, bool hermitian)
{
  double initTime = getTime();
  long initMem  = memoryUsed();

  Scalar *y = (Scalar *)dbg_alloca(numDir()*sizeof(Scalar));
  
  this->op2 = y;
  this->op1 = r;
  this->op3 = p;

  int i;
  for(i = 0; i < numDir(); ++i)
    y[i] = 0.0;

  if(hermitian) operation = &CGOrthoOp<Scalar>::FdotH;
  else operation = &CGOrthoOp<Scalar>::Fdot;
  threadManager->execTimedParal(timer, this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom) this->fetiCom->globalSum(numDir(), y);
#endif

  for(i = 0; i < numDir(); ++i) 
    y[i] /= allPFiP[i+lastIndex()];
  
  operation = &CGOrthoOp<Scalar>::multAdd;
  threadManager->execTimedParal(timer, this->numTasks, this->oos);

  timer.addOverAll( memoryUsed()-initMem, getTime()-initTime );
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::orthogonalize(Scalar *r, Scalar *p)
{
  Scalar *y = (Scalar *)dbg_alloca(numDir()*sizeof(Scalar));
   
  this->op2 = y;
  this->op1 = r;
  this->op3 = p;

  int i;
  for(i = 0; i < numDir(); ++i)
    y[i] = 0.0;

  operation = &CGOrthoOp<Scalar>::Fdot;
  threadManager->execParal(this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom)
     this->fetiCom->globalSum(numDir(), y);
#endif

  for(i = 0; i < numDir(); ++i) 
    y[i] /= allPFiP[i+lastIndex()];
  
  operation = &CGOrthoOp<Scalar>::multAdd;
  threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::wtr()
{
  int i;
  Scalar *y = (Scalar *) dbg_alloca(sizeof(Scalar)*os->kSize);

  if(length() > 0) {  // PJSA
    // Compute y = W^t r
    Tgemv('T', length(), os->kSize,  1.0, 
          getAllP() + os->lastindex*length(), 
          length(),  os->op1+offset(), 1, 0.0, 
          y, 1);
  }
  else {
    // initialize y
    for(i=0; i < os->kSize; ++i) y[i] = 0.0;
  }

  os->lock.lock();
  for(i = 0; i < os->kSize; ++i)
     os->op2[i+os->lastindex] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::wtKpr()
{
  int i;
  Scalar *y = (Scalar *) dbg_alloca(sizeof(Scalar)*os->kSize);

  // Pr = preconditioned residual
  // 
  // Compute y = y - W^t K Pr
  // 
  // or now y = W^t r - W^t K Pr

  Scalar *Kw = getAllFiP() + os->lastindex*length();

  if(length() > 0) { 
    Tgemv('T', length(), os->kSize, -1.0, Kw,
          length(), os->op3+offset(), 1, 0.0, 
          y, 1);
  } 
  else {
    for(i=0; i < os->kSize; ++i) y[i] = 0.0;
  }

  os->lock.lock();
  for(i = 0; i < os->kSize; ++i)
    os->op2[i+os->lastindex] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::addWy()
{
  Scalar *w = getAllP() + os->lastindex*length();

  if(length() > 0) { 
    // Compute pr = pr + W y
    Tgemv('N', length(), os->kSize, 1.0, w,
          length(),  os->op2+os->lastindex, 1, 1.0,
          os->op3+offset(), 1);
  }
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::precondition(Scalar *r, Scalar *pr)
{
  // r  = residual
  // pr = preconditioned residual
  // numOrthoSet = number of stored Krylov spaces
  // this->numP        = number of total stored vectors
  // this->numTasks    = number of processors

  Scalar *y = (Scalar *)dbg_alloca(this->numP*sizeof(Scalar));

  int i, isys;

  for(i = 0; i < this->numP; ++i)
    y[i] = 0.0;

  this->op2 = y;
  this->op1 = r;
  this->op3 = pr;

  // KHP: testing an idea 
  int numO, firstSys;
  if( numOrthoSet <= 1 ) {
    firstSys = 0;
    numO     = 0;
  } else {
    firstSys = 0;
    numO = numOrthoSet - 1;
  }

  for(isys=firstSys; isys<numO; ++isys) {
    lastindex = kindex[isys];
    kSize     = kindex[isys+1] - lastindex;
    if(kSize == 0) continue;

    operation = &CGOrthoOp<Scalar>::wtr;
    threadManager->execParal(this->numTasks, this->oos);

    operation = &CGOrthoOp<Scalar>::wtKpr;
    threadManager->execParal(this->numTasks, this->oos);
    
#ifdef DISTRIBUTED
  if(this->fetiCom)
     this->fetiCom->globalSum(this->numP, y);
#endif
    // Compute y = y / W^tKW
    // so y = (W^t r - W^t K Pr)/(W^tKW)
    for(i=0; i<kSize; ++i)
      y[i+lastindex] /= allPFiP[i+lastindex];

    operation = &CGOrthoOp<Scalar>::addWy;
    threadManager->execParal(this->numTasks, this->oos);
  }
}

template<class Scalar>
void
GenCGOrthoSet<Scalar>::predict(Scalar *r, Scalar *lambda)
{
  int i;
  Scalar *y = (Scalar *)dbg_alloca(numDir()*sizeof(Scalar));

  this->op2 = y;
  this->op1 = r;
  this->op3 = lambda;

  for(i = 0; i < numDir(); ++i) y[i] = 0;

  operation = &CGOrthoOp<Scalar>::dot;
  threadManager->execParal(this->numTasks, this->oos);

#ifdef DISTRIBUTED
  if(this->fetiCom)
    this->fetiCom->globalSum(numDir(), y);
#endif
  for(i = 0; i < numDir(); ++i)
    y[i] /= allPFiP[i+lastIndex()];

  operation = &CGOrthoOp<Scalar>::mult;
  threadManager->execParal(this->numTasks, this->oos);
}

template<class Scalar>
CGOrthoOp<Scalar>::CGOrthoOp(GenCGOrthoSet<Scalar> *_os, int _index, int _loclen)
{
 os = _os; this->idx = _index; this->loclen = _loclen;
 this->locAllP = locAllFiP = 0;
 locAllD = 0;
 this->numP = 0;
}

template<class Scalar>
void
CGOrthoOp<Scalar>::clean_up()
{
 if(this->locAllP) {
   delete [] this->locAllP;
   this->locAllP=0;
 }
}

template<class Scalar>
void
CGOrthoOp<Scalar>::run()
{
 (this->*os->operation)();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::addVec()
{
 if(this->numP == os->maxP)
   this->numP--;

 if(this->locAllP == 0) {
   this->locAllP   = new Scalar[2*os->maxP*this->loclen];
   locAllFiP = this->locAllP + os->maxP*this->loclen;
 }
 int i;
 for(i = 0; i < this->loclen; ++ i) {
   this->locAllP[i+ this->numP*this->loclen]   = os->op1[i+this->idx];
   locAllFiP[i+ this->numP*this->loclen] = os->op2[i+this->idx];
 }
 this->numP++;
}

template<class Scalar>
void
CGOrthoOp<Scalar>::addVecFric()
{
 if(this->numP == os->maxP)
   this->numP--;

 if(this->locAllP == 0) {
   this->locAllP = new Scalar[ 3 * os->maxP * this->loclen];
   locAllD = this->locAllP + 2 * os->maxP*this->loclen;
   locAllFiP = this->locAllP + os->maxP*this->loclen;
 }
 int i;
 for(i = 0; i < this->loclen; ++ i) {
   this->locAllP[i+ this->numP*this->loclen]   = os->op1[i+this->idx];
   locAllFiP[i+ this->numP*this->loclen] = os->op2[i+this->idx];
   locAllD[ i + this->numP * this->loclen] = os->op3[ i + this->idx];
 }
 this->numP++;
}

template<class Scalar>
void
CGOrthoOp<Scalar>::reset()
{
 this->numP = 0;
}

template<class Scalar>
void
CGOrthoOp<Scalar>::Fdot()
{
  if(this->loclen == 0) return;
  Scalar *y = (Scalar *)dbg_alloca(os->numDir()*sizeof(Scalar));

  // Compute y = y - w^t K P r
  Tgemv('T', this->loclen, os->numDir(), Scalar(-1.0), 
        locAllFiP+os->lastIndex()*this->loclen, this->loclen, os->op1+this->idx, 1, Scalar(0.0), y, 1);

  os->lock.lock();

  int i;
  for(i = 0; i < os->numDir(); ++i)
     os->op2[i] += y[i];

  os->lock.unlock();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::FdotH()
{
  if(this->loclen == 0) return;
  Scalar *y = (Scalar *)dbg_alloca(os->numDir()*sizeof(Scalar));

  // Compute y = y - w^H K P r
  Tgemv('C', this->loclen, os->numDir(), Scalar(-1.0),
        locAllFiP+os->lastIndex()*this->loclen, this->loclen, os->op1+this->idx, 1, Scalar(0.0), y, 1);

  os->lock.lock();

  int i;
  for(i = 0; i < os->numDir(); ++i)
     os->op2[i] += y[i];

  os->lock.unlock();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::FdotSingleDir()
{
  // if the local vector length is zero, then no computation is required
  if(this->loclen == 0) return;

  // local declarations
  Scalar alpha;

  // Compute: - FiP^t * y for search direction i
  Tgemv('T', this->loclen, 1, -1.0, 
        locAllFiP + ( os->lastIndex() + os->op4) * this->loclen, 
        this->loclen, os->op1 + this->idx, 1, 0.0, &alpha, 1);

  // sum contributions to inproduct over threads
  os->lock.lock();
  *os->op3 += alpha;
  os->lock.unlock();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::multAddSingleD()
{
  // if the local vector length is zero, then no computation is required
  if(this->loclen == 0) return;

  // Compute: - FiP^t * y for search direction i
  Tgemv('N', this->loclen, 1, 1.0, 
        locAllD + ( os->lastIndex() + os->op4) * this->loclen, 
        this->loclen, os->op3, 1, 1.0, os->op1 + this->idx, 1);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::multAddSingleP()
{
  // if the local vector length is zero, then no computation is required
  if(this->loclen == 0) return;

  // Compute: - FiP^t * y for search direction i
  Tgemv('N', this->loclen, 1, 1.0, 
        this->locAllP + ( os->lastIndex() + os->op4) * this->loclen, 
        this->loclen, os->op3, 1, 1.0, os->op2 + this->idx, 1);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::RdotSingleD()
{
  // if the local vector length is zero, then no computation is required
  if(this->loclen == 0) return;

  // local declarations
  Scalar alpha;

  // Compute: di^t * r for search direction i
  Tgemv('T', this->loclen, 1, 1.0, 
        locAllD + ( os->lastIndex() + os->op4) * this->loclen, 
        this->loclen, os->op1 + this->idx, 1, 0.0, &alpha, 1);

  // sum contributions to inproduct over threads
  os->lock.lock();
  *os->op3 += alpha;
  os->lock.unlock();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::multSubtractSingleP()
{
  // if the local vector length is zero, then no computation is required
  if(this->loclen == 0) return;

  // Compute: lambda -= EtaJK * P for search direction i
  Tgemv('N', this->loclen, 1, -1.0, 
        this->locAllP + ( os->lastIndex() + os->op4) * this->loclen, 
        this->loclen, os->op3, 1, 1.0, os->op2 + this->idx, 1);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::multSubtractSingleFp()
{
  // if the local vector length is zero, then no computation is required
  if(this->loclen == 0) return;

  // Compute: r -= EtaJK * FiP for search direction i
  Tgemv('N', this->loclen, 1, -1.0, 
        locAllFiP + ( os->lastIndex() + os->op4) * this->loclen, 
        this->loclen, os->op3, 1, 1.0, os->op1 + this->idx, 1);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::dot()
{
  if(this->loclen == 0) return;
  Scalar *y = (Scalar *)dbg_alloca(os->numDir()*sizeof(Scalar));

  // Compute y = y - w^t r
  Tgemv('T', this->loclen, os->numDir(), -1.0, 
        this->locAllP+os->lastIndex()*this->loclen,
        this->loclen, os->op1+this->idx, 1, 0.0, y, 1);

  os->lock.lock();
  int i;
  for(i = 0; i < os->numDir(); ++i)
     os->op2[i] += y[i];
  os->lock.unlock();
}

template<class Scalar>
void
CGOrthoOp<Scalar>::mult()
{
  if(this->loclen == 0) return;
  Tgemv('N', this->loclen, os->numDir(), 1.0, 
        this->locAllP+os->lastIndex()*this->loclen,
        this->loclen,  os->op2, 1, 1.0, os->op3+this->idx, 1);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::multAdd()
{
  if(this->loclen == 0) return;

  // Compute precResid = precResid + W y
  int i;
  for(i=0; i< this->loclen; ++i)
    os->op3[this->idx+i] = os->op1[this->idx+i];

  Tgemv('N', this->loclen, os->numDir(), Scalar(1.0), 
        this->locAllP + os->lastIndex()*this->loclen,
        this->loclen,  os->op2, 1, Scalar(1.0), os->op3+this->idx, 1);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::multSubFp()
{
  if(this->loclen == 0) return;
  // Compute r = r - Fp * alpha
  int i;
  for(i=0; i< this->loclen; ++i)
    os->op3[this->idx+i] = os->op1[this->idx+i];

  Tgemv('N', this->loclen, os->numDir(), -1.0, 
        locAllFiP + os->lastIndex()*this->loclen,
        this->loclen,  os->op2, 1, 1.0, os->op3+this->idx, 1);
}

template<class Scalar>
void
CGOrthoOp<Scalar>::shift()
{
 int i,j;
 int shiftSize = os->kSize;

 for(i=shiftSize; i<os->numP; ++i)
   for(j=0; j<this->loclen; ++j) {
     this->locAllP[j+(i-shiftSize)*this->loclen] =   this->locAllP[j+i*this->loclen];
     locAllFiP[j+(i-shiftSize)*this->loclen] = locAllFiP[j+i*this->loclen];
     locAllD[j+(i-shiftSize)*this->loclen] =   locAllD[j+i*this->loclen];
   }
}

template<class Scalar>
CGOrthoOp<Scalar>::~CGOrthoOp()
{
 if(this->locAllP) { delete [] this->locAllP; this->locAllP = 0; }
 // don't delete locAllFiP or locAllD (these point inside this->locAllP)
}

template
class GenCGOrthoSet<double>;

template
class CGOrthoOp<double>;

template
class GenCGOrthoSet<std::complex<double>>;

template
class CGOrthoOp<std::complex<double>>;