#ifndef _VECTORSET_H_
#define _VECTORSET_H_

#include <iostream>
#include <Math.d/Vector.h>

template <class Scalar> class GenVector;
template <class Scalar> class GenStackVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;

template<class Scalar>
class GenVectorSet {
 protected: 
   int    numVectors;
   int    length;
   GenVector<Scalar> *vecSet;
   bool myMemory;
   
   void copy(const GenVector<Scalar>&    v);
   void copy(const GenVectorSet<Scalar>& vset);
  
  public:
  
   // Constructors
   GenVectorSet() { numVectors = 0; length = 0; vecSet = 0; }
   GenVectorSet(int numVec, int length);
   GenVectorSet(int numVec, int length, Scalar initialValue);
   GenVectorSet(int numVec, int length, Scalar *vectors, bool myMemory = true);
   
   virtual ~GenVectorSet(); //HB: made the destructor virtual to deal with 
                            //    the derived class GenStackVectorSet
   
   // Copy constructor
   GenVectorSet(const GenVectorSet<Scalar> &vecset);

   GenVector<Scalar>& operator [] (int i) const;
   GenVector<Scalar>* operator +  (int i) const;

   int numVec() const { return numVectors; }
   int size()   const { return length;     }
   
   void print(const char *msg = "");
   void zero();
};


template<class Scalar>
inline GenVector<Scalar>&
GenVectorSet<Scalar>::operator [] (int i) const
{
  return vecSet[i];
}

template<class Scalar>
inline GenVector<Scalar>*
GenVectorSet<Scalar>::operator +  (int i) const
{
 return vecSet+i;
}

template<class Scalar>
inline 
GenVectorSet<Scalar>::~GenVectorSet()
{
  if(vecSet) { delete [] vecSet; vecSet = 0; }
}

//HB
template<class Scalar>
class GenStackVectorSet: public GenVectorSet<Scalar> 
{
  public:
    GenStackVectorSet() { this->numVectors = 0; this->length = 0; this->vecSet = 0; }
    GenStackVectorSet(int numVec, int length, Scalar *ptr);
    GenStackVectorSet(int numVec, int length, Scalar *ptr, Scalar initval);
    ~GenStackVectorSet() { }
};
                                                                                                                    
template<class Scalar>
GenStackVectorSet<Scalar>::GenStackVectorSet(int _numVec, int _length, Scalar *ptr)
{
  this->numVectors = _numVec;
  this->length     = _length;
  this->vecSet     = new GenStackVector<Scalar>[this->numVectors];

  for(int i=0; i<this->numVectors; i++)
    this->vecSet[i].setData(&ptr[i*this->length],this->length);
}

template<class Scalar>
GenStackVectorSet<Scalar>::GenStackVectorSet(int _numVec, int _length, Scalar *ptr, Scalar initval)
{
  this->numVectors = _numVec;
  this->length     = _length;
  this->vecSet     = new GenStackVector<Scalar>[this->numVectors];
                                                                                                                    
  for(int i=0; i<this->numVectors; i++){
    this->vecSet[i].setData(&ptr[i*this->length],this->length);
    this->vecSet[i] = initval;
  }
}

#ifdef _TEMPLATE_FIX_
#include <Math.d/VectorSet.C>
#endif

#endif
