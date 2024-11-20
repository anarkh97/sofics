#ifndef _FDIST_VECTOR_H_
#define _FDIST_VECTOR_H_

#include <Utils.d/dbg_alloca.h>
#include <Math.d/FVector.h>
#include <Feti.d/DistrVector.h>

#define MPI_OMP_REDUCTION
//------------------------------------------------------------------------------
namespace NewVec {

template<class Scalar>
class DistVec : public Vec<Scalar> {

public:

  typedef DistrInfo InfoType;

private:

  const DistrInfo &distInfo;

  Vec<Scalar> **subVec;

public:
  
  DistVec(const DistrInfo &);
  DistVec(const DistVec<Scalar> &);
  DistVec(const DistrInfo &, Scalar *, Vec<Scalar> **);
  ~DistVec();

  template<class T>
  DistVec<Scalar> &operator=(const T);

  template<class T>
  DistVec<Scalar> &operator*=(const T);

  DistVec<Scalar> &operator=(const DistVec<Scalar> &);
  DistVec<Scalar> &operator+=(const DistVec<Scalar> &);
  DistVec<Scalar> &operator-=(const DistVec<Scalar> &);

  Scalar operator*(const DistVec<Scalar> &);

  template<class T>
  DistVec<Scalar> &operator=(const Expr<T, Scalar> &);

  template<class T>
  DistVec<Scalar> &operator+=(const Expr<T, Scalar> &);

  template<class T>
  DistVec<Scalar> &operator-=(const Expr<T, Scalar> &);
  
  template<class T>
  Scalar operator*(const Expr<T, Scalar> &);

  Vec<Scalar> &operator() (int i) const { return *subVec[i]; }

  void createSubVec();

  DistVec<Scalar> *alias() const;

  Scalar sum() const;

  Scalar min() const;

  Scalar max() const;

  const InfoType &info() const { return distInfo; }

  int getNumLocSub() const { return distInfo.numDom; }

  int subSize(int isub) const { return distInfo.domLen[isub]; }

  Scalar *subData(int isub) const { return this->v+distInfo.subOffset[isub]; }

  int *getMasterFlag(int i) const { return distInfo.getMasterFlag(i); }

  void reduce(DistVec<Scalar> &, int **, int *);

  void setNewData(Scalar *d); 
};

//------------------------------------------------------------------------------

template <class Scalar>
DistVec<Scalar>::DistVec(const DistrInfo &dI) : 
  Vec<Scalar>(dI.len), distInfo(dI), subVec(0)
{

  createSubVec();

}

//------------------------------------------------------------------------------

template <class Scalar>
DistVec<Scalar>::DistVec(const DistVec<Scalar> &v2) :
  Vec<Scalar>(v2.size()), distInfo(v2.distInfo), subVec(0)
{

// #pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

     int locOffset = distInfo.threadOffset[iThread];
     int locLen = distInfo.threadLen[iThread];

     for (int i = 0; i < locLen; ++i) 
       this->v[locOffset+i] = v2[locOffset+i];

  }

  createSubVec();

}

//------------------------------------------------------------------------------

template<class Scalar>
DistVec<Scalar>::DistVec(const DistrInfo &dI, Scalar *vv, Vec<Scalar> **sv) : 
  Vec<Scalar>(dI.len, vv), distInfo(dI), subVec(sv)
{

}

//------------------------------------------------------------------------------

template <class Scalar>
DistVec<Scalar>::~DistVec() 
{ 

  if (/*this->locAlloc &&*/ subVec) {

    //#pragma omp parallel for BUG alloc
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread)
      if (subVec[iThread]) { delete subVec[iThread]; subVec[iThread] = 0; }

    delete [] subVec; subVec = 0; 

  }

}

//------------------------------------------------------------------------------

template <class Scalar>
inline
void 
DistVec<Scalar>::createSubVec()
{
  subVec = new Vec<Scalar>*[distInfo.numLocSub];

  //#pragma omp parallel for BUG alloc
  for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread)
    subVec[iThread] = new Vec<Scalar>(distInfo.subLen[iThread], this->v+distInfo.subOffset[iThread]);

}

template <class Scalar>
inline
void
DistVec<Scalar>::setNewData(Scalar *d)
{
  if(this->locAlloc && this->v) delete [] this->v; this->v=d;
  if(subVec==0) createSubVec();
  else {
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread)
      subVec[iThread]->setNewData(this->v+distInfo.subOffset[iThread]); 
  }
}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> *
DistVec<Scalar>::alias() const
{

  return new DistVec<Scalar>(distInfo, this->v, subVec); 

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator=(const T y)
{

// #pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] = y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator*=(const T y)
{

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] *= y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator=(const DistVec<Scalar> &v2)
{

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] = v2.v[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator+=(const DistVec<Scalar> &v2)
{

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] += v2.v[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator-=(const DistVec<Scalar> &v2)
{

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread)  {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] -= v2.v[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
inline
Scalar 
DistVec<Scalar>::operator*(const DistVec<Scalar> &x)
{

  int iThread;

  Scalar res = 0;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(dbg_alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) allres[iThread] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  locres += this->v[locOffset+i] * x.v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	locres += this->v[locOffset+i] * x.v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif
      
    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) res += allres[iThread];
#endif

  return res;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;

// #pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] = x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
DistVec<Scalar> &
DistVec<Scalar>::operator+=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] += x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator-=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] -= x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
template<class T>
inline
Scalar 
DistVec<Scalar>::operator*(const Expr<T, Scalar> &expr)
{

  int iThread;

  Scalar res = 0;

  const T &x = expr.x;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(dbg_alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) allres[iThread] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  locres += this->v[locOffset+i] * x[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	locres += this->v[locOffset+i] * x[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif
      
    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) res += allres[iThread];
#endif

  return res;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
inline
Scalar 
DistVec<Scalar>::sum() const
{
  
  int iThread;

  Scalar res = 0;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(dbg_alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) allres[iThread] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  locres += this->v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	locres += this->v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) res += allres[iThread];
#endif

  return res;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
inline
Scalar 
DistVec<Scalar>::min() const
{
  
  int iThread;

  Scalar *allmin = reinterpret_cast<Scalar *>(dbg_alloca(sizeof(Scalar) * distInfo.numLocSub));

//#pragma omp parallel for
  for (iThread = 0; iThread < distInfo.numLocSub; ++iThread)
    allmin[iThread] = subVec[iThread]->min();

  Scalar vmin = allmin[0];
  for (iThread = 1; iThread < distInfo.numLocSub; ++iThread)
    vmin = vmin < allmin[iThread] ? vmin : allmin[iThread];

  distInfo.com->globalMin(1, &vmin);

  return vmin;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
inline
Scalar 
DistVec<Scalar>::max() const
{
  
  int iThread;

  Scalar *allmax = reinterpret_cast<Scalar *>(dbg_alloca(sizeof(Scalar) * distInfo.numLocSub));

//#pragma omp parallel for
  for (iThread = 0; iThread < distInfo.numLocSub; ++iThread)
    allmax[iThread] = subVec[iThread]->max();

  Scalar vmax = allmax[0];
  for (iThread = 1; iThread < distInfo.numLocSub; ++iThread)
    vmax = vmax > allmax[iThread] ? vmax : allmax[iThread];

  distInfo.com->globalMax(1, &vmax);

  return vmax;

}

//--------------------------------------------------------------------------

template<class Scalar, int dim>
class DistSVec : public SVec<Scalar,dim> {

public:

  typedef DistrInfo InfoType;

private:

  const DistrInfo &distInfo;

  SVec<Scalar,dim> **subVec;

public:

  DistSVec(const DistrInfo &);
  DistSVec(const DistSVec<Scalar,dim> &);
  DistSVec(const DistrInfo &, Scalar (*)[dim], SVec<Scalar,dim> **);
  ~DistSVec();

  DistSVec<Scalar,dim> &operator=(const Scalar);
  DistSVec<Scalar,dim> &operator*=(const Scalar);

  DistSVec<Scalar,dim> &operator=(const DistSVec<Scalar,dim> &);
  DistSVec<Scalar,dim> &operator+=(const DistSVec<Scalar,dim> &);
  DistSVec<Scalar,dim> &operator-=(const DistSVec<Scalar,dim> &);

  DistSVec<Scalar,dim> &operator=(const DistVec<Scalar> &);

  Scalar operator*(const DistSVec<Scalar,dim> &);

  template<class T>
  DistSVec<Scalar,dim> &operator=(const Expr<T, Scalar> &);

  template<class T>
  DistSVec<Scalar,dim> &operator+=(const Expr<T, Scalar> &);

  template<class T>
  DistSVec<Scalar,dim> &operator-=(const Expr<T, Scalar> &);

  template<class T>
  Scalar operator*(const Expr<T, Scalar> &);

  SVec<Scalar,dim> &operator() (int i) const { return *subVec[i]; }

  void createSubVec();

  void set(const Scalar *);
  
  void restrict();

  void average();

  int nonOverlapSize() const;

  DistSVec<Scalar,dim> *alias() const;

  const InfoType &info() const { return distInfo; }

  int getNumLocSub() const { return distInfo.numLocSub; }

  int subSize(int i) const { return distInfo.subLen[i]; }

  Scalar (*subData(int i) const)[dim] { return this->v+distInfo.subOffset[i]; }

  int *getMasterFlag(int i) const { return distInfo.getMasterFlag(i); }

  void reduce(DistSVec<Scalar, dim> &, int **, int *);

};

//---------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::DistSVec(const DistrInfo &dI) : 
  SVec<Scalar,dim>(dI.len), distInfo(dI), subVec(0)
{

  createSubVec();

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::DistSVec(const DistSVec<Scalar,dim> &y) :
  SVec<Scalar,dim>(y.size()), distInfo(y.distInfo), subVec(0)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = yy[locOffset+i];

  }

  createSubVec();

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::DistSVec(const DistrInfo &dI, Scalar (*vv)[dim], 
				SVec<Scalar,dim> **sv) : 
  SVec<Scalar,dim>(dI.len, vv), distInfo(dI), subVec(sv)
{

}

//----------------------------------------------------------------------

template<class Scalar>
void DistVec<Scalar>::reduce(DistVec<Scalar> &reducedDV, int **mFlag, 
				int *numFlags)  {

  for (int iThread = 0; iThread < distInfo.numLocSub; iThread++)  {

    Scalar *v1 = subData(iThread);
    Scalar *v2 = reducedDV.subData(iThread);

    int count = 0;
    for (int iData = 0; count < numFlags[iThread] ; iData++)
      if (mFlag[iThread][iData] >= 0)
	v2[count++] = v1[iData];
  }
}

//----------------------------------------------------------------------

template<class Scalar, int dim>
void DistSVec<Scalar,dim>::reduce(DistSVec<Scalar, dim> &reducedDSV, 
				int **mFlag, int *numFlags)  {

  for (int iThread = 0; iThread < distInfo.numLocSub; iThread++)  {

    Scalar (*v1)[dim] = subData(iThread);
    Scalar (*v2)[dim] = reducedDSV.subData(iThread);

    int count = 0;
    for (int iData = 0; count < numFlags[iThread]; iData++)  {

      if (mFlag[iThread][iData] >= 0)  {
        for (int iComp = 0; iComp < dim; iComp++)
          v2[count][iComp] = v1[iData][iComp];
        count++;
      }
    }
  }

}

//----------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::~DistSVec() 
{ 

  if (/*this->locAlloc &&*/ subVec) {

    //#pragma omp parallel for BUG alloc
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread)
      if (subVec[iThread]) { delete subVec[iThread]; subVec[iThread] = 0; }

    delete [] subVec; subVec = 0; 

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void 
DistSVec<Scalar,dim>::createSubVec()
{

  subVec = new SVec<Scalar,dim>*[distInfo.numLocSub];

  //#pragma omp parallel for BUG alloc
  for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread)
    subVec[iThread] = new SVec<Scalar,dim>(distInfo.subLen[iThread], this->v+distInfo.subOffset[iThread]);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void
DistSVec<Scalar,dim>::set(const Scalar *y)
{

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      for (int j = 0; j < dim; ++j)
	this->v[locOffset+i][j] = y[j];
    
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> *
DistSVec<Scalar,dim>::alias() const
{

  return new DistSVec<Scalar,dim>(distInfo, this->v, subVec); 

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const Scalar y)
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator*=(const Scalar y)
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] *= y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const DistSVec<Scalar,dim> &y)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator+=(const DistSVec<Scalar,dim> &y)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] += yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator-=(const DistSVec<Scalar,dim> &y)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] -= yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const DistVec<Scalar> &y)
{

  const Scalar *yy = y.data();

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = distInfo.threadOffset[iThread];
    int locLen = distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) 
      for (int j = 0; j < dim; ++j)
	this->v[locOffset+i][j] = yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------
    
template<class Scalar, int dim>
inline
Scalar 
DistSVec<Scalar,dim>::operator*(const DistSVec<Scalar,dim> &x)
{

  int iThread;

  Scalar res = 0;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(dbg_alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) allres[iThread] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  for (int j = 0; j < dim; ++j)
	    locres += this->v[locOffset+i][j] * x.v[locOffset+i][j];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	for (int j = 0; j < dim; ++j)
	  locres += this->v[locOffset+i][j] * x.v[locOffset+i][j];
      
#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) res += allres[iThread];
#endif

  return res;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator+=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] += x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator-=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocThreads; ++iThread) {

    int locOffset = dim * distInfo.threadOffset[iThread];
    int locLen = dim * distInfo.threadLen[iThread];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] -= x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------
    
template<class Scalar, int dim>
template<class T>
inline
Scalar 
DistSVec<Scalar,dim>::operator*(const Expr<T, Scalar> &expr)
{

  int iThread;

  Scalar res = 0;

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(dbg_alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) allres[iThread] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  for (int j = 0; j < dim; ++j)
	    locres += vv[(locOffset+i)*dim + j] * x[(locOffset+i)*dim + j];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
//#pragma omp parallel for reduction(+: res)
#else
//#pragma omp parallel for
#endif
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      int locOffset = distInfo.subOffset[iThread];
      int locLen = distInfo.subLen[iThread];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	for (int j = 0; j < dim; ++j)
	  locres += vv[(locOffset+i)*dim + j] * x[(locOffset+i)*dim + j];
      
#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iThread]] = locres;
#endif

    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iThread=0; iThread<distInfo.numGlobSub; ++iThread) res += allres[iThread];
#endif

  return res;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
int 
DistSVec<Scalar,dim>::nonOverlapSize() const
{

  int tsize = 0;

  if (distInfo.masterFlag) {

//#pragma omp parallel for reduction(+: tsize)
    for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

      for (int i=0; i<distInfo.subLen[iThread]; ++i)
	if (distInfo.masterFlag[ distInfo.subOffset[iThread]+i ]) ++tsize;

    }

  } 
  else
    tsize = distInfo.len;

  return tsize;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void
DistSVec<Scalar,dim>::restrict()
{

  if (!distInfo.masterFlag) return;

//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread) {
 
    int locOffset = distInfo.subOffset[iThread];
    
    for (int i=0; i<distInfo.subLen[iThread]; ++i)
      if (distInfo.masterFlag[locOffset+i] == 0)
	for (int j=0; j<dim; ++j) this->v[locOffset+i][j] = 0;
      
  }

}

//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
void
DistSVec<Scalar,dim>::average()
{
   if (!distInfo.invNdWeight) return;
//#pragma omp parallel for
  for (int iThread = 0; iThread < distInfo.numLocSub; ++iThread) {

    int locOffset = distInfo.subOffset[iThread];
   
    for (int i=0; i<distInfo.subLen[iThread]; ++i) {
       for (int j=0; j<dim; ++j) 
          this->v[locOffset+i][j] *= distInfo.invNdWeight[locOffset+i];
    }
   }

}
*/
//------------------------------------------------------------------------------

} // end of namespace
#endif
