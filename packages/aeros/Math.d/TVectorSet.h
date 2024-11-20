#ifndef _VECTOR_SET_H_
#define _VECTOR_SET_H_

#include <memory>
#include <iostream>

//------------------------------------------------------------------------------

template<class VecType>
class VecSet {

  int numVec;
  const typename VecType::InfoType &len;

  std::allocator<VecType> alloc;
  VecType *vecSet;

public:

  VecSet(int, const typename VecType::InfoType &);
  VecSet(const VecSet<VecType> &);
  ~VecSet();

  VecType &operator[] (int i) const { return vecSet[i]; }
  
  int numVectors() const { return numVec; }
  typename VecType::InfoType &size() const { return len; }

  void resize(int);

  void print(const char *msg = "") {
    if (msg) std::cerr << msg << std::endl;

    for (int i=0; i<numVec; ++i) {
      std::cerr << "vector " << i << ":";
      vecSet[i].print();
    }
  }

};

//------------------------------------------------------------------------------

template<class VecType>
VecSet<VecType>::VecSet(int _numVec, const typename VecType::InfoType &_len) : len(_len)
{

  numVec = _numVec;

  vecSet = static_cast<VecType*>(alloc.allocate(numVec));

  for (int i = 0; i < numVec; ++i)
    new (static_cast<void *>(vecSet+i)) VecType(len);

}

//------------------------------------------------------------------------------

template<class VecType>
void VecSet<VecType>::resize(int n)
{
 int i;
 for (i = 0; i < numVec; ++i)
   vecSet[i].~VecType();

 alloc.deallocate(vecSet, numVec);

 numVec = n;

 vecSet = static_cast<VecType*>(alloc.allocate(numVec));

 for (i = 0; i < numVec; ++i)
    new (static_cast<void *>(vecSet+i)) VecType(len);

}

//------------------------------------------------------------------------------

template<class VecType>
VecSet<VecType>::VecSet(const VecSet<VecType> &vectorSet)
{
 
  numVec = vectorSet.numVectors();
  len = vectorSet.size();

  vecSet = static_cast<VecType*>(alloc.allocate(numVec));

  for (int i = 0; i < numVec; ++i) {
    new (static_cast<void *>(vecSet+i)) VecType(len);
    vecSet[i] = vectorSet[i];
  }

}

//------------------------------------------------------------------------------

template<class VecType>
VecSet<VecType>::~VecSet() 
{

  if (vecSet) {

    for (int i = 0; i < numVec; ++i) vecSet[i].~VecType();

    alloc.deallocate(vecSet, numVec);

    vecSet = 0;

  }

}

//------------------------------------------------------------------------------

#endif
