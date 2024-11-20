#include <Utils.d/dbg_alloca.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

#include <Threads.d/Paral.h>
#include <Utils.d/MyComplex.h>

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif

template<class Scalar>
class GenVecOp : public TaskDescr 
{
	GenDistrVector<Scalar> *v1;
	const GenDistrVector<Scalar> *v2;
	const GenDistrVector<Scalar> *v3;
	const GenDistrVector<Scalar> *v4;
	Scalar *res;
	Scalar c;
	Scalar c1;
	Scalar c2;
	Scalar c3;
	void (GenVecOp<Scalar>::*f)(int);

   public:
     GenVecOp(void (GenVecOp<Scalar>::*_f)(int),
          GenDistrVector<Scalar> *_v1, const GenDistrVector<Scalar> *_v2=0);

     GenVecOp(void (GenVecOp<Scalar>::*_f)(int),
          GenDistrVector<Scalar> *_v1, Scalar *rc);

     GenVecOp(void (GenVecOp<Scalar>::*_f)(int), 
          GenDistrVector<Scalar> *_v1, const GenDistrVector<Scalar> *_v2, Scalar c);

     GenVecOp(void (GenVecOp<Scalar>::*_f)(int), 
          GenDistrVector<Scalar> *_v1, Scalar c);

     GenVecOp(void (GenVecOp<Scalar>::*_f)(int), 
          GenDistrVector<Scalar> *_v1, const GenDistrVector<Scalar> *_v2, Scalar *rc);


     GenVecOp(void (GenVecOp<Scalar>::*_f)(int),
          GenDistrVector<Scalar> *_v1, const GenDistrVector<Scalar> *_v2, Scalar c,
          const GenDistrVector<Scalar> *_v3);

     GenVecOp(void (GenVecOp<Scalar>::*_f)(int),
          GenDistrVector<Scalar> *_v1, Scalar c1, const GenDistrVector<Scalar> *_v2,
                            Scalar c2, const GenDistrVector<Scalar> *_v3);

     GenVecOp(void (GenVecOp<Scalar>::*_f)(int),
          GenDistrVector<Scalar> *_v1, Scalar c1, const GenDistrVector<Scalar> *_v2,
                            Scalar c2, const GenDistrVector<Scalar> *_v3, Scalar c3,
                            const GenDistrVector<Scalar> *_v4);

     void alloc(int);
     void dot(int);
     void dot_ignore_master_flag(int);
     void zero(int);
     void linAdd(int);
     void linAdd_inv(int);
     void linAdd1(int);
     void linAdd2(int);
     void linC(int);
     void linC1(int);
     void linC2(int);
     void linC3(int);
     void assign(int);
     void assign_times(int);
     void assign_div(int);
     void assign_plus(int);
     void assign_minus(int);
     void assign_divide(int);
     void assign_special(int);
     void swap(int);
     void run() override;
     void runFor(int) override;
     void negate(int);
     void compute_partial(int subNum);
     void partial_dot(int subNum);
     void linAdd_partial(int subNum);
     void partial_max(int subNum);
     void tdot(int);
};

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int), 
                           GenDistrVector<Scalar> *_v1, const GenDistrVector<Scalar> *_v2)
{
 f  = _f;
 v1 = _v1;
 v2 = _v2;
}

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int),
                           GenDistrVector<Scalar> *_v1, Scalar *_r)
{
 f  = _f;
 v1 = _v1;
 res = _r;
}

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int), 
                           GenDistrVector<Scalar> *_v1, const GenDistrVector<Scalar> *_v2, Scalar _c)
{
 f  = _f;
 v1 = _v1;
 v2 = _v2;
 c  = _c;
}

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int), GenDistrVector<Scalar> *_v1, Scalar _c)
             
{
 f  = _f;
 v1 = _v1;
 c  = _c;
}

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int), GenDistrVector<Scalar> *_v1, 
                           const GenDistrVector<Scalar> *_v2, Scalar _c, const GenDistrVector<Scalar> *_v3)
{
 f  = _f;
 v1 = _v1;
 v2 = _v2;
 c  = _c;
 v3 = _v3;
}

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int), GenDistrVector<Scalar> *_v1, 
                           Scalar _c1, const GenDistrVector<Scalar> *_v2,
                           Scalar _c2, const GenDistrVector<Scalar> *_v3)
{
 f  = _f;
 v1 = _v1;
 c1 = _c1;
 v2 = _v2;
 c2 = _c2;
 v3 = _v3;
}

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int), GenDistrVector<Scalar> *_v1,
                           Scalar _c1, const GenDistrVector<Scalar> *_v2,
                           Scalar _c2, const GenDistrVector<Scalar> *_v3,
                           Scalar _c3, const GenDistrVector<Scalar> *_v4)
{
 f  = _f;
 v1 = _v1;
 c1 = _c1;
 v2 = _v2;
 c2 = _c2;
 v3 = _v3;
 c3 = _c3;
 v4 = _v4;
}

template<class Scalar>
GenVecOp<Scalar>::GenVecOp(void (GenVecOp<Scalar>::*_f)(int), GenDistrVector<Scalar> *_v1, 
                           const GenDistrVector<Scalar> *_v2, Scalar *_r)
{
 f   = _f;
 v1  = _v1;
 v2  = _v2;
 res = _r;
}

template<class Scalar>
void
GenVecOp<Scalar>::compute_partial(int subNum)
{
 Scalar *d1 = v1->subData(subNum);
 res[subNum] = 0.0;
 for(int i = 0; i < v1->subLen(subNum); ++i)
   res[subNum] += d1[i]*ScalarTypes::conj(d1[i]);
}

template<class Scalar>
void
GenVecOp<Scalar>::partial_dot(int subNum)
{
 if(v1->getPartial(subNum) == 0.0) { res[subNum] = 0.0; return; }
 else {
   Scalar *d1 = v1->subData(subNum);
   Scalar *d2 = v2->subData(subNum);
   bool *masterFlag = v1->subMasterFlag(subNum);
   res[subNum] = 0.0;
   for(int i = 0; i < v1->subLen(subNum); ++i)
     if(masterFlag[i]) res[subNum] += d1[i]*ScalarTypes::conj(d2[i]);
 }
}

template<class Scalar>
void
GenVecOp<Scalar>::partial_max(int subNum)
{
  Scalar *d1 = v1->subData(subNum);
  res[subNum] = 0.0;
  for(int i = 0; i < v1->subLen(subNum); ++i) {
    double d_i = ScalarTypes::norm(d1[i]);
    if(d_i > ScalarTypes::Real(res[subNum])) res[subNum] = d_i;
  }
}

#ifdef DETERMINISTIC
template<class Scalar>
void
GenVecOp<Scalar>::dot(int subNum)
{
 Scalar *d1 = v1->subData(subNum);
 Scalar *d2 = v2->subData(subNum);
 bool *masterFlag = v1->subMasterFlag(subNum);
 int len = v1->subLen(subNum);
 Scalar r = 0;
 int i;
 for(i = 0; i < len; ++i)
   if(masterFlag[i]) 
     r += ScalarTypes::conj(d2[i])*d1[i];
 res[subNum] = r;
}
#else
template<class Scalar>
void
GenVecOp<Scalar>::dot(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 bool *masterFlag = v1->threadMasterFlag(threadNum);
 Scalar r = 0;
 int i;
 for(i = 0; i < len; ++i)
   if(masterFlag[i])
     r += ScalarTypes::conj(d2[i])*d1[i];
 res[threadNum] = r;
}
#endif

#ifdef DETERMINISTIC
template<class Scalar>
void
GenVecOp<Scalar>::dot_ignore_master_flag(int subNum)
{
 const Scalar *d1 = v1->subData(subNum);
 const Scalar *d2 = v2->subData(subNum);
 const int len = v1->subLen(subNum);

 Scalar r = Scalar();
 for(int i = 0; i < len; ++i) {
   r += ScalarTypes::conj(d2[i]) * d1[i];
 }
 res[subNum] = r;
}
#else
template<class Scalar>
void
GenVecOp<Scalar>::dot_ignore_master_flag(int threadNum)
{
 const Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 const int len = v1->threadLen(threadNum);

 Scalar r = Scalar();
 for(int i = 0; i < len; ++i) {
   r += ScalarTypes::conj(d2[i]) * d1[i];
 }
 res[threadNum] = r;
}
#endif

#ifdef DETERMINISTIC
template<class Scalar>
void
GenVecOp<Scalar>::tdot(int subNum)
{
 Scalar *d1 = v1->subData(subNum);
 Scalar *d2 = v2->subData(subNum);
 bool *masterFlag = v1->subMasterFlag(subNum);
 int len = v1->subLen(subNum);
 Scalar r = 0;
 int i;
 for(i = 0; i < len; ++i)
   if(masterFlag[i])
     r += d1[i]*d2[i];
 res[subNum] = r;
}
#else
template<class Scalar>
void
GenVecOp<Scalar>::tdot(int threadNum)
{
 const Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 bool *masterFlag = v1->threadMasterFlag(threadNum);
 Scalar r = 0;
 int i;
 for(i = 0; i < len; ++i)
   if(masterFlag[i])
     r += d1[i]*d2[i];
 res[threadNum] = r;
}
#endif

template<class Scalar>
void
GenVecOp<Scalar>::linAdd(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] += c*d2[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::linAdd_inv(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] += d2[i]/c;
}

template<class Scalar>
void
GenVecOp<Scalar>::linAdd_partial(int subNum)
{
 if(v2->getPartial(subNum) == 0.0) return;
 Scalar *d1 = v1->subData(subNum);
 const Scalar *d2 = v2->subData(subNum);
 for(int i = 0; i < v1->subLen(subNum); ++i) d1[i] += c*d2[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::linAdd1(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] += d2[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::linAdd2(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 const Scalar *d3 = v3->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] += c1*d2[i] + c2*d3[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::linC(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 const Scalar *d3 = v3->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] = d2[i] + c*d3[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::linC1(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);

 int len = v1->threadLen(threadNum);

 int i;
 for(i = 0; i < len; ++i)
    d1[i] = c*d2[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::linC2(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 const Scalar *d3 = v3->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] = c1*d2[i] + c2*d3[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::linC3(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 const Scalar *d3 = v3->threadData(threadNum);
 const Scalar *d4 = v4->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] = c1*d2[i] + c2*d3[i] + c3*d4[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::assign(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] = d2[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::assign_special(int subNum)
{
 Scalar *d1 = v1->subData(subNum);
 const Scalar *d2 = v2->subData(subNum);
 int len1 = v1->subLen(subNum);
 int len2 = v2->subLen(subNum);
 int minlen = std::min(len1,len2);
 int i;
 for(i = 0; i < minlen; ++i)
    d1[i] = d2[i];
 for(i = minlen; i < len1; ++i)
    d1[i] = 0;
}

template<class Scalar>
void
GenVecOp<Scalar>::swap(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 Scalar *d2 = const_cast<Scalar *>(v2->threadData(threadNum));
 int len = v1->threadLen(threadNum);
 Scalar tmp;
 int i;
 for(i = 0; i < len; ++i) {
    tmp   = d1[i];
    d1[i] = d2[i];
    d2[i] = tmp;
 }
}

template<class Scalar>
void 
GenDistrVector<Scalar>::addBlockSqr(int ii, Scalar c, GenDistrVector<Scalar> &x) 
{
 for(int i=0; i<len; ++i) v[i] += c*x[i]*x[i];
}

template<class Scalar>
void 
GenDistrVector<Scalar>::computeSqrt() 
{
 for(int i=0; i<len; ++i) v[i] = sqrt(v[i]);
}


template<class Scalar>
void
GenVecOp<Scalar>::assign_times(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] *= c;
}

template<class Scalar>
void
GenVecOp<Scalar>::assign_div(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] /= c;
}

template<class Scalar>
void
GenVecOp<Scalar>::assign_plus(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] += d2[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::assign_minus(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] -= d2[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::assign_divide(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 const Scalar *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i) {
    if(d2[i] != 0.0) d1[i] /= d2[i];
    else d1[i] *= -1.0e256;
 }
}

template<class Scalar>
void
GenVecOp<Scalar>::zero(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i = 0; i < len; ++i)
    d1[i] = 0.0;
}

template<class Scalar>
void
GenVecOp<Scalar>::negate(int threadNum)
{
 Scalar *d1 = v1->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i=0; i<len; ++i)
   d1[i] = -d1[i];
}

template<class Scalar>
void
GenVecOp<Scalar>::alloc(int threadNum)
{
 v1->thV[threadNum] = new Scalar[v1->threadLen(threadNum)];
}

template<class Scalar>
void
GenVecOp<Scalar>::run()
{
    throw "Unhandled case";
}

template<class Scalar>
void
GenVecOp<Scalar>::runFor(int threadNum)
{
 (this->*f)(threadNum);
}

template <class T, class Scalar>
class ExprVecAssign : public TaskDescr {
   GenDistrVector<Scalar> &lhs;
   const Expr<T,Scalar> &rhs;
 public:
   ExprVecAssign(GenDistrVector<Scalar> &v, const Expr<T,Scalar> &expr) :
     lhs(v), rhs(expr) { }
   void runFor(int);
};

template <class T, class Scalar>
void ExprVecAssign<T,Scalar>::runFor(int threadNum) {
  int len = lhs.threadLen(threadNum);
  int offset = lhs.threadOffset(threadNum);
  for(int i = 0; i < len; i++)
    lhs[offset+i] = rhs[offset+i];
}


template <class T, class Scalar>
class ExprVecIncr : public TaskDescr {
   GenDistrVector<Scalar> &lhs;
   const Expr<T,Scalar> &rhs;
 public:
   ExprVecIncr(GenDistrVector<Scalar> &v, const Expr<T,Scalar> &expr) :
     lhs(v), rhs(expr) { }
   void runFor(int);
};

template <class T, class Scalar>
void ExprVecIncr<T,Scalar>::runFor(int threadNum) {
  int len = lhs.threadLen(threadNum);
  int offset = lhs.threadOffset(threadNum);
  for(int i = 0; i < len; i++)
    lhs[offset+i] += rhs[offset+i];
}


template <class T, class Scalar>
class ExprVecDecr : public TaskDescr {
   GenDistrVector<Scalar> &lhs;
   const Expr<T,Scalar> &rhs;
 public:
   ExprVecDecr(GenDistrVector<Scalar> &v, const Expr<T,Scalar> &expr) :
     lhs(v), rhs(expr) { }
   void runFor(int);
};

template <class T, class Scalar>
void ExprVecDecr<T,Scalar>::runFor(int threadNum) {
  int len = lhs.threadLen(threadNum);
  int offset = lhs.threadOffset(threadNum);
  for(int i = 0; i < len; i++)
    lhs[offset+i] -= rhs[offset+i];
}

template<class Scalar>
Scalar
GenDistrVector<Scalar>::ident()
{
 Scalar res = 0;
 int i;
 for(i=0; i < len; ++i)
   res += ((i % 7)+1)*v[i];
 return res;
}

template<class Scalar>
void
GenDistrVector<Scalar>::initialize()
{
 // Set the overall length
 len = inf->len;

 // Set number of Domains
 numDom = inf->numDom;

 // Allocate memory for the values
 if(!v) v = new Scalar[len];

 // Array of Scalar pointers to hold sub vector values
 subV.resize(inf->numDom);

 // Array to hold sub vector lengths
 subVLen.resize(inf->numDom);
 subVOffset.resize(inf->numDom);

 if(inf->numDom>0)
 	subV[0] = v;
 Scalar *v2 = v;
 nT = std::min(threadManager->numThr(),inf->numDom); // PJSA: currently the number of threads being larger than
                                                     // the number of subdomains doesn't work
 thV.resize(nT);
 thLen.resize(nT);
 thOffset.resize(nT);

	// Distrinfo::computeOffsets should be explicitly called before this
	// constructor for interface vectors
	infoFlag = inf->subOffset != nullptr;
	if(infoFlag)
		masterFlag = new bool[len];
	else
		masterFlag = inf->masterFlag;

 int tLen;
 int md = 0;
 for(int iThread = 0; iThread < nT; ++iThread) {
   tLen = 0;
   thV[iThread] = v2;
   thOffset[iThread] = v2-v;
   for(int pseudomd = iThread; pseudomd < inf->numDom; pseudomd += nT) {
     subV[md] = v2;
     subVOffset[md] = v2-v;
     if(infoFlag) { // re-arrange masterFlag
       for(int i=0; i<inf->domLen[md]; ++i)
         masterFlag[i+subVOffset[md]] = inf->masterFlag[i+inf->subOffset[md]];
     }
     subVLen[md] = inf->domLen[md];
     v2 += inf->domLen[md];
     tLen += inf->domLen[md];
     ++md;
   }
   thLen[iThread] = tLen;
 }
 partial = 0;
}

template<class Scalar>
GenDistrVector<Scalar>::GenDistrVector(const DistrInfo &i) 
 : inf(&i)
{
 v = 0;
 myMemory = true;
 initialize();
}

template<class Scalar>
GenDistrVector<Scalar>::GenDistrVector(const GenDistrVector<Scalar> &x) 
 : inf(x.inf)
{
 v = 0;
 myMemory = true; 
 initialize();
 for(int i=0; i<len; ++i) v[i] = x.data()[i];
}

template<class Scalar>
void
GenDistrVector<Scalar>::resize(const DistrInfo &other)
{
  if(*inf == other) { inf = &other; return; }
  clean_up();
  inf = &other;
  v = 0;
  myMemory = true;
  initialize();
}

template<class Scalar>
void
GenDistrVector<Scalar>::conservativeResize(const DistrInfo &other)
{
  if(*inf == other) { inf = &other; return; }
  GenDistrVector<Scalar> copy(*this);
  clean_up();
  inf = &other;
  v = 0;
  myMemory = true;
  initialize();

  GenVecOp<Scalar> assignSpecialOp(&GenVecOp<Scalar>::assign_special, this, &copy);

  threadManager->execParal(numDom, &assignSpecialOp);
}

template<class Scalar>
GenDistrVector<Scalar>::~GenDistrVector()
{
  clean_up();
}

template<class Scalar>
void
GenDistrVector<Scalar>::clean_up()
{
 if(v) {
   if(myMemory) delete [] v;
   v=0;
 }

 if(infoFlag && masterFlag) {
   delete [] masterFlag;
   masterFlag = 0;
 }

 if(partial) { delete [] partial; partial = 0; }
}

template<class Scalar>
GenDistrVector<Scalar>::GenDistrVector(const DistrInfo &i, Scalar *_v, bool _myMemory)
 : inf(&i)
{
 myMemory = _myMemory;
 v = _v;
 initialize();
}

template<class Scalar>
void
GenDistrVector<Scalar>::zero()
{
 GenVecOp<Scalar> zeroAll(&GenVecOp<Scalar>::zero,this);
 threadManager->execParal(nT, &zeroAll);
}

template<class Scalar>
void
GenDistrVector<Scalar>::negate()
{
 GenVecOp<Scalar> negateAll(&GenVecOp<Scalar>::negate,this);
 threadManager->execParal(nT, &negateAll);
}

template<class Scalar>
double GenDistrVector<Scalar>::norm()
{
 Scalar product = (*this)*(*this); 
 double res = ScalarTypes::norm(product);
 return sqrt(res); 
}

template<class Scalar>
void
GenPartialDistrVector<Scalar>::computePartial()
{
 if(!this->partial) this->partial = new Scalar[this->numDom];
 GenVecOp<Scalar> op1(&GenVecOp<Scalar>::compute_partial, this, this->partial);
 threadManager->execParal(this->numDom, &op1);
}

template<class Scalar>
Scalar
GenPartialDistrVector<Scalar>::operator*(const GenDistrVector<Scalar> &x) const
{
 // partial dot product
 Scalar *p = (Scalar *) alloca(sizeof(Scalar)*this->numDom);
 GenVecOp<Scalar> op1(&GenVecOp<Scalar>::partial_dot, this, &x, p);
 threadManager->execParal(this->numDom, &op1);
 Scalar res = 0;
 for(int i=0; i < this->numDom; ++i) res += p[i];
#ifdef DISTRIBUTED
 if(structCom) res = structCom->globalSum(res);
#endif
 return res;
}

template<class Scalar>
double GenDistrVector<Scalar>::infNorm()
{
 Scalar *p = (Scalar *) alloca(sizeof(Scalar)*this->numDom);
 GenVecOp<Scalar> op1(&GenVecOp<Scalar>::partial_max, this, p);
 threadManager->execParal(this->numDom, &op1);
 double res = 0;
 for(int i=0; i < this->numDom; ++i) {
   double p_i = ScalarTypes::Real(p[i]);
   if(p_i > res) res = p_i;
 }
#ifdef DISTRIBUTED
 if(structCom) res = structCom->globalMax(res);
#endif
 return res;
}

namespace {
template <typename T>
auto lose_const(const T *t) { return const_cast<T *>(t); }
template <typename T>
auto lose_const(T *t) { return t; }
}

template<class Scalar>
Scalar
GenDistrVector<Scalar>::operator*(const GenDistrVector<Scalar> &x) const
{
 Scalar *partial = (Scalar *) dbg_alloca(sizeof(Scalar)*numDom);
 GenVecOp<Scalar> dotAll(&GenVecOp<Scalar>::dot,
                         lose_const(this), lose_const(&x), partial);
#ifdef DETERMINISTIC
 threadManager->execParal(numDom, &dotAll);
#else
 threadManager->execParal(nT, &dotAll);
#endif
 Scalar res = 0;
 int i;
#ifdef DETERMINISTIC
 for(i=0; i < numDom; ++i)
#else
 for(i=0; i < nT; ++i)
#endif
    res += partial[i];
#ifdef DISTRIBUTED
 if(structCom)
   res = structCom->globalSum(res);
#endif
 return res;
}

template<class Scalar>
Scalar
dot_ignore_master_flag(const GenDistrVector<Scalar> &v1, const GenDistrVector<Scalar> &v2)
{
 const int numDom = v1.num();
 Scalar *partial = (Scalar *) dbg_alloca(sizeof(Scalar) * numDom);
 GenVecOp<Scalar> dotAll(&GenVecOp<Scalar>::dot_ignore_master_flag,
                         const_cast<GenDistrVector<Scalar> *>(&v1),
                         const_cast<GenDistrVector<Scalar> *>(&v2),
                         partial);
#ifdef DETERMINISTIC
 threadManager->execParal(numDom, &dotAll);
#else
 const int nT = v1.numThreads();
 threadManager->execParal(nT, &dotAll);
#endif
 Scalar res = Scalar();
#ifdef DETERMINISTIC
 for(int i = 0; i < numDom; ++i)
#else
 for(int i = 0; i < nT; ++i)
#endif
    res += partial[i];
#ifdef DISTRIBUTED
 if(structCom)
   res = structCom->globalSum(res);
#endif
 return res;
}

template<class Scalar>
Scalar
GenDistrVector<Scalar>::operator^(const GenDistrVector<Scalar> &x) const
{
 Scalar *partial = (Scalar *) dbg_alloca(sizeof(Scalar)*numDom);
 GenVecOp<Scalar> dotAll(&GenVecOp<Scalar>::tdot,
                         const_cast<GenDistrVector<Scalar>*>(this),
                         const_cast<GenDistrVector<Scalar>*>(&x), partial);
#ifdef DETERMINISTIC
 threadManager->execParal(numDom, &dotAll);
#else
 threadManager->execParal(nT, &dotAll);
#endif
 Scalar res = 0;
 int i;
#ifdef DETERMINISTIC
 for(i=0; i < numDom; ++i)
#else
 for(i=0; i < nT; ++i)
#endif
    res += partial[i];
#ifdef DISTRIBUTED
 if(structCom)
   res = structCom->globalSum(res);
#endif
 return res;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator=(const GenDistrVector<Scalar> &x)
{
 if(*inf != x.info()) {
   resize(x.info());
 }
 GenVecOp<Scalar> assignOp(&GenVecOp<Scalar>::assign, this, const_cast<GenDistrVector<Scalar> *>(&x));
 threadManager->execParal(nT, &assignOp);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator=(Scalar c)
{
 int i;
 for(i=0; i < len; ++i)
   v[i] = c;
 return *this;
}

template<class Scalar>
void
GenDistrVector<Scalar>::initRand()
{
 for(int i=0; i < len; ++i)
   v[i] = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator*=(Scalar c)
{
 GenVecOp<Scalar> assign_times(&GenVecOp<Scalar>::assign_times, this, c);
 threadManager->execParal(nT, &assign_times);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator/=(Scalar c)
{
 GenVecOp<Scalar> assign_div(&GenVecOp<Scalar>::assign_div, this, c);
 threadManager->execParal(nT, &assign_div);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator+=(GenDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in +=\n");
 }
 GenVecOp<Scalar> assign_plus(&GenVecOp<Scalar>::assign_plus,this, &x);
 threadManager->execParal(nT, &assign_plus);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator+=(const GenDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in +=\n");
 }
 GenVecOp<Scalar> assign_plus(&GenVecOp<Scalar>::assign_plus,this, const_cast<GenDistrVector<Scalar>* >(&x));
 threadManager->execParal(nT, &assign_plus);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator-=(GenDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in -=\n");
 }
 GenVecOp<Scalar> assign_minus(&GenVecOp<Scalar>::assign_minus,this, &x);
 threadManager->execParal(nT, &assign_minus);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator/=(GenDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in /=\n");
 }
 GenVecOp<Scalar> assign_divide(&GenVecOp<Scalar>::assign_divide,this, &x);
 threadManager->execParal(nT, &assign_divide);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linAdd(GenDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in linAdd\n");
 }
 GenVecOp<Scalar> add(&GenVecOp<Scalar>::linAdd1, this, &x);
 threadManager->execParal(nT, &add);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linAdd(Scalar c, GenDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in linAdd\n");
 }
 GenVecOp<Scalar> add(&GenVecOp<Scalar>::linAdd, this, &x, c);
 threadManager->execParal(nT, &add);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linAdd_inv(Scalar c, GenDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in linAdd_inv\n");
 }
 GenVecOp<Scalar> add(&GenVecOp<Scalar>::linAdd_inv, this, &x, c);
 threadManager->execParal(nT, &add);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linAdd(Scalar c, GenPartialDistrVector<Scalar> &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in linAdd\n");
 }
 GenVecOp<Scalar> add(&GenVecOp<Scalar>::linAdd_partial, this, &x, c);
 threadManager->execParal(numDom, &add);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linAdd(Scalar c1, GenDistrVector<Scalar> &x, Scalar c2, 
                               GenDistrVector<Scalar> &y)
{
 if(x.len != y.len) {
  fprintf(stderr, "Length error in linAdd\n");
 }
 GenVecOp<Scalar> add2(&GenVecOp<Scalar>::linAdd2, this, c1, &x, c2, &y);
 threadManager->execParal(nT, &add2);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linC(const GenDistrVector<Scalar> &x, Scalar c, const GenDistrVector<Scalar> &y)
{
 if(x.len != y.len) {
  fprintf(stderr, "Length error in linC\n");
 }
 GenVecOp<Scalar> linC(&GenVecOp<Scalar>::linC, this, &x, c, &y);
 threadManager->execParal(nT, &linC);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linC(const GenDistrVector<Scalar> &v, Scalar c)
{
 if(*inf != v.info()) {
   resize(v.info());
 }
 GenVecOp<Scalar> linC1(&GenVecOp<Scalar>::linC1, this, &v, c);
 threadManager->execParal(nT, &linC1);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linC( Scalar c1, const GenDistrVector<Scalar> &x,
                              Scalar c2, const GenDistrVector<Scalar> &y)
{
 GenVecOp<Scalar> linC2(&GenVecOp<Scalar>::linC2, this, c1, &x, c2, &y);
 threadManager->execParal(nT, &linC2);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::linC( Scalar c1, const GenDistrVector<Scalar> &x,
                              Scalar c2, const GenDistrVector<Scalar> &y,
                              Scalar c3, const GenDistrVector<Scalar> &z)
{
 GenVecOp<Scalar> linC3(&GenVecOp<Scalar>::linC3, this, c1, &x, c2, &y, c3, &z);
 threadManager->execParal(nT, &linC3);
 return *this;
}

template<class Scalar>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::swap( GenDistrVector<Scalar> &x)
{
 GenVecOp<Scalar> linC2(&GenVecOp<Scalar>::swap,this,&x);
 threadManager->execParal(nT, &linC2);
 return *this;
}

template<class Scalar>
void
GenDistrVector<Scalar>::doubleUp(int threadNum, Scalar coef, GenDistrVector<Scalar> *x, 
                                 GenDistrVector<Scalar> *y, GenDistrVector<Scalar> *z)
{
 Scalar *d  = threadData(threadNum);
 Scalar *d1 = x->threadData(threadNum);
 Scalar *d2 = y->threadData(threadNum);
 Scalar *d3 = z->threadData(threadNum);
 int len1 = threadLen(threadNum);
 int len2 = y->threadLen(threadNum);
 int i;
 for(i = 0; i < len1; ++i) {
    d[i] += coef*d1[i];
 }
 for(i = 0; i < len2; ++i) {
    d2[i] += coef*d3[i];
 }
}

template<class Scalar>
void
GenDistrVector<Scalar>::tripleUp(int threadNum, Scalar coef, GenDistrVector<Scalar> *x,
                                 GenDistrVector<Scalar> *y, GenDistrVector<Scalar> *z,
                                 GenDistrVector<Scalar> *a, GenDistrVector<Scalar> *b)
{
 Scalar *d  = threadData(threadNum);
 Scalar *d1 = x->threadData(threadNum);
 Scalar *d2 = y->threadData(threadNum);
 Scalar *d3 = z->threadData(threadNum);
 Scalar *d4 = a->threadData(threadNum);
 Scalar *d5 = b->threadData(threadNum);
 int len1 = threadLen(threadNum);
 int len2 = y->threadLen(threadNum);
 int len3 = a->threadLen(threadNum);
 int i;
 for(i = 0; i < len1; ++i) {
    d[i] += coef*d1[i];
 }
 for(i = 0; i < len2; ++i) {
    d2[i] += coef*d3[i];
 }
 for(i = 0; i < len3; ++i) {
    d4[i] += coef*d5[i];
 }
}

#include <Threads.d/PHelper.h>

template<class Scalar>
void doubleUpdate(Scalar nu, GenDistrVector<Scalar> &x1, GenDistrVector<Scalar> &dx1,
                  GenDistrVector<Scalar> &x2, GenDistrVector<Scalar> &dx2)
{
 execParal(x1.numThreads(), &x1, &GenDistrVector<Scalar>::doubleUp, nu, &dx1, &x2, &dx2);
}

template<class Scalar>
void tripleUpdate(Scalar nu, GenDistrVector<Scalar> &x1, GenDistrVector<Scalar> &dx1,
                  GenDistrVector<Scalar> &x2, GenDistrVector<Scalar> &dx2,
                  GenDistrVector<Scalar> &x3, GenDistrVector<Scalar> &dx3)
{
 execParal(x1.numThreads(), &x1, &GenDistrVector<Scalar>::tripleUp, nu, &dx1, &x2, &dx2, &x3, &dx3);
}

template<class Scalar>
void
GenDistrVector<Scalar>::print()
{
 int i;
 for(i=0; i < len; ++i) std::cerr << "v[" << i << "] = " << v[i] << std::endl;
 std::cerr << std::endl;
}

template<class Scalar>
void
GenDistrVector<Scalar>::printNonZeroTerms()
{
 int i;
 for(i=0; i < len; ++i) 
   if(fabs(v[i]) > 0.00000001) std::cerr << "v[" << i << "] = " << v[i] << std::endl;
  std::cerr << std::endl;
}


template<class Scalar>
void
GenDistrVector<Scalar>::printAll()
{
 std::cerr << "Length of GenDistrVector<Scalar> = " << len << std::endl;
 std::cerr << "Number of Domains = " << numDom << std::endl;
 int i,j;
 for(i=0; i<numDom; ++i) {
   std::cerr << "--- Sub Vector Length = " << subVLen[i] << std::endl;
   for(j=0; j<subVLen[i]; ++j)
     std::cerr << "v(" << j+1 << ") = " << subV[i][j] << std::endl;
 }
 std::cerr << std::endl;
}

template <class Scalar>
template <class T>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator=(const Expr<T,Scalar> &x) {
 ExprVecAssign<T, Scalar> assign(*this,x);
 threadManager->execParal(nT, &assign);
 return *this;
}

template <class Scalar>
template <class T>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator+=(const Expr<T,Scalar> &x) {
 ExprVecIncr<T, Scalar> assign(*this,x);
 threadManager->execParal(nT, &assign);
 return *this;
}

template <class Scalar>
template <class T>
GenDistrVector<Scalar> &
GenDistrVector<Scalar>::operator-=(const Expr<T,Scalar> &x) {
 ExprVecDecr<T, Scalar> assign(*this,x);
 threadManager->execParal(nT, &assign);
 return *this;
}

template<class Scalar>
double norm(const GenDistrVector<Scalar> &v) {
  return v.norm();
}

template<class T1, class Scalar>
class NormCompute : public TaskDescr {
   const DistrInfo &info;
   const Expr<T1,Scalar,const DistrInfo &> &expr;
  public:
   NormCompute(const Expr<T1,Scalar,const DistrInfo &> &e)
    : info(e.info()), expr(e) { }
   double *partialDot;
   void runFor(int threadNum) {
     partialDot[threadNum] = 0;
     for(int i = info.threadOffset[threadNum];
                i < info.threadOffset[threadNum+1]; ++i)
       partialDot[threadNum] += ScalarTypes::sqNorm(expr[i]);
   }
   double total() {
     double res = 0;
     for(int i = 0; i < info.numLocThreads; ++i)
       res += partialDot[i];
     #ifdef DISTRIBUTED
     if(structCom)
       res = structCom->globalSum(res);
     #endif
     return res;
   }
};

template<class T1, class Scalar>
double norm(const Expr<T1,Scalar,const DistrInfo &> &e) {
    int nT = e.info().numLocThreads;
    NormCompute<T1,Scalar> nc(e);
    nc.partialDot = (double *)dbg_alloca(sizeof(double)*nT);
    threadManager->execParal(nT, &nc);
    return sqrt(nc.total());
}
