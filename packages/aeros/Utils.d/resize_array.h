#ifndef _RESIZE_ARRAY_H_
#define _RESIZE_ARRAY_H_

template <class Type> class ResizeArray {
   Type *d;
   Type init_v;
   int  csize;  // current size
   bool myData;
   ResizeArray(const ResizeArray &); // should never be used
   ResizeArray(ResizeArray &); // should never be used
 public:
   ResizeArray(Type i, int ini_size=16);
   ~ResizeArray();
    Type &operator[] (int i);
    const Type &operator[](int i) const { return d[i]; }
   void resize(int ns);
   int max_size() { return csize; }
   Type *operator + (int i); //dangerous operator to use with caution.
    // It will not check that the returned pointer is not used beyond
    // the existing values
   Type *yield(); // When we want to keep the array and remove the resizing
   void deleteArray()  { delete [] d; d = 0; csize = 0; }
   void restartArray() { d = new Type[csize]; };
   Type *data(bool _myData = true) { myData = _myData; return d; }
};

template <class Type>
inline
Type &ResizeArray<Type>::operator[](int i) {
 if(i >= csize)
    resize(i);
 return d[i];
}

template <class Type>
inline
Type *ResizeArray<Type>::operator+(int i) {
  if(i >= csize)
    resize(i);
 return d+i;
}

#ifdef _TEMPLATE_FIX_
#include <Utils.d/resize_array.C>
#endif

#endif
