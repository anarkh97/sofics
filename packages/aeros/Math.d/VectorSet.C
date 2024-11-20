#include <cstdio>

#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>

template<class Scalar>
GenVectorSet<Scalar>::GenVectorSet(int _numVectors, int _length)
{
 // This constructor creates a vector set of numVectors
 // each intialized to zero

 numVectors = _numVectors;
 length     = _length;

 vecSet = new GenVector<Scalar>[numVectors];

 GenVector<Scalar> v(length, 0.0);

 copy(v);
}

template<class Scalar>
GenVectorSet<Scalar>::GenVectorSet(int _numVectors, int _length, Scalar initialValue)
{
 // This constructor creates a vector set of numVectors
 // each intialized to initialValue 

 numVectors = _numVectors;
 length     = _length;

 vecSet = new GenVector<Scalar>[numVectors];

 GenVector<Scalar> v(length, initialValue);

 copy(v);
}

template<class Scalar>
GenVectorSet<Scalar>::GenVectorSet(const GenVectorSet<Scalar> &vectorSet)
{
 // copy constructor constructs an idential copy of the incoming
 // Vectorset object
 
 numVectors = vectorSet.numVec();
 length     = vectorSet.size();

 vecSet = new GenVector<Scalar>[numVectors];

 copy(vectorSet);

}

template<class Scalar>
GenVectorSet<Scalar>::GenVectorSet(int _numVectors, int _length, Scalar* vectors, bool _myMemory)
{
 // Set number of vectors and length of those vectors
 numVectors = _numVectors;
 length     = _length;
 myMemory   = _myMemory;

 // Allocate
 vecSet = new GenVector<Scalar>[numVectors];

 for(int i=0; i<numVectors; ++i) {
//   if(myMemory) {
//     vecSet[i].init(length);
//     for(int j=0; j<length; ++j)
//        vecSet[i][j] = vectors[j + i*length];
//   } else {
     vecSet[i].setData(&vectors[i*length],length);
     vecSet[i].setmyMemory((myMemory && i==0) ? true : false);
//   }
 }
 zero();

/*
 if(myMemory) {
    vecSet = new GenVector<Scalar>[numVectors];
    int i,j;
    for(i=0; i<numVectors; ++i)
      for(j=0; j<length; ++j) 
         vecSet[i][j] = vectors[j + i*length];
 } else {
    int i;
    vecSet = new GenStackVector<Scalar>[numVectors];
    for(i=0; i<numVectors; ++i)
      vecSet[i].setData(&vectors[i*length],length);
 }
 zero();
*/
}

template<class Scalar>
void
GenVectorSet<Scalar>::copy(const GenVector<Scalar> &v)
{
 int i;
 for(i=0; i<numVectors; ++i)
   vecSet[i] = v;
}

template<class Scalar>
void
GenVectorSet<Scalar>::copy(const GenVectorSet<Scalar> &vectorSet)
{
 int i;
 for(i=0; i<numVectors; ++i)
   vecSet[i] = vectorSet[i];
}

template<class Scalar>
void 
GenVectorSet<Scalar>::print(const char *message)
{
   if(message) fprintf(stderr,"%s\n",message);
   int i,j;
   for(i=0; i<numVectors; ++i) {
     for(j=0; j<length; ++j)
       fprintf(stderr,"%e\n",vecSet[i][j]);
     fprintf(stderr,"vector %d\n",i+1);
   }

}

template<class Scalar>
void 
GenVectorSet<Scalar>::zero()
{
   for(int i=0;i<numVectors; ++i)
      vecSet[i].zero();
}     
