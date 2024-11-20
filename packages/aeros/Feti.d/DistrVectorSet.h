#ifndef _DISTRVECTORSET_H_
#define _DISTRVECTORSET_H_

#include <Feti.d/DistrVector.h>

template<class Scalar>
class GenDistrVectorSet 
{
 protected:
   int    numVectors;
   int    length;
   GenDistrVector<Scalar> **vecSet;
   bool   myMemory;

   void copy(const GenDistrVector<Scalar>& v);
   void copy(const GenDistrVectorSet<Scalar>& vset);

  public:

   // Constructors
   GenDistrVectorSet() { numVectors = 0; length = 0; vecSet = 0; }
   GenDistrVectorSet(int numVec, DistrInfo &dinfo, Scalar initialValue = 0.0);
   GenDistrVectorSet(int numVec, DistrInfo &dinfo, Scalar *vectors, bool myMemory = true);

   virtual ~GenDistrVectorSet();

   // Copy constructor
   GenDistrVectorSet(const GenDistrVectorSet<Scalar> &vecset);

   GenDistrVector<Scalar>& operator [] (int i) const;
   GenDistrVector<Scalar>* operator +  (int i) const;

   int numVec() const { return numVectors; }
   int size()   const { return length;     }

   void print(const char *msg = "");
   void zero();
};

template<class Scalar>
inline GenDistrVector<Scalar>&
GenDistrVectorSet<Scalar>::operator [] (int i) const
{
  return *(vecSet[i]);
}

template<class Scalar>
inline GenDistrVector<Scalar>*
GenDistrVectorSet<Scalar>::operator +  (int i) const
{
 return *(vecSet+i);
}

template<class Scalar>
inline
GenDistrVectorSet<Scalar>::~GenDistrVectorSet()
{
  if(vecSet) { 
    for(int i=0; i<numVectors; ++i) delete vecSet[i];
    delete [] vecSet; vecSet = 0; 
  }
}

template<class Scalar>
class GenStackDistrVectorSet: public GenDistrVectorSet<Scalar>
{
  public:
    GenStackDistrVectorSet() { this->numVectors = 0; this->length = 0; this->vecSet = 0; }
    GenStackDistrVectorSet(int numVec, DistrInfo &dinfo, Scalar *ptr);
    GenStackDistrVectorSet(int numVec, DistrInfo &dinfo, Scalar *ptr, Scalar initval);
    ~GenStackDistrVectorSet() { }
};

template<class Scalar>
GenStackDistrVectorSet<Scalar>::GenStackDistrVectorSet(int _numVec, DistrInfo &dinfo, Scalar *ptr)
{
  this->numVectors = _numVec;
  this->vecSet     = new GenStackDistVector<Scalar> * [this->numVectors];
  this->length     = this->numVectors*dinfo.totLen();

  for(int i=0; i<this->numVectors; i++)
    this->vecSet[i] = new GenStackDistVector<Scalar>(dinfo,ptr+i*dinfo.totLen());
}

template<class Scalar>
GenStackDistrVectorSet<Scalar>::GenStackDistrVectorSet(int _numVec, DistrInfo &dinfo, Scalar *ptr, Scalar initval)
{
  this->numVectors = _numVec;
  this->vecSet     = new GenStackDistVector<Scalar> * [this->numVectors];
  this->length     = this->numVectors*dinfo.totLen();

  for(int i=0; i<this->numVectors; i++) {
    this->vecSet[i] = new GenStackDistVector<Scalar>(dinfo,ptr+i*dinfo.totLen());
    (*this->vecSet[i]) = initval;
  }
}

#ifdef _TEMPLATE_FIX_
#include <Feti.d/DistrVectorSet.C>
#endif

#endif
