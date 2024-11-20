#include <cstdio>

#include <Feti.d/DistrVector.h>
#include <Feti.d/DistrVectorSet.h>

template<class Scalar>
GenDistrVectorSet<Scalar>::GenDistrVectorSet(int _numVectors, DistrInfo &dinfo, Scalar initialValue)
{
 // This constructor creates a distributed vector set of numVectors
 // each initialized to initialValue 
 numVectors = _numVectors;
 vecSet = new GenDistrVector<Scalar> * [numVectors];
 length = numVectors*dinfo.totLen();
 for(int i=0; i<numVectors; ++i) {
   vecSet[i] = new GenDistrVector<Scalar>(dinfo);
   (*vecSet[i]) = initialValue;
 }
}

template<class Scalar>
GenDistrVectorSet<Scalar>::GenDistrVectorSet(const GenDistrVectorSet<Scalar> &vectorSet)
{
 // copy constructor constructs an idential copy of the incoming
 // distributed vectorset object
 numVectors = vectorSet.numVec();
 vecSet = new GenDistrVector<Scalar> * [numVectors];
 length = 0;
 for(int i=0; i<numVectors; ++i) {
   vecSet[i] = new GenDistrVector<Scalar>(vectorSet[i]);
   length += vecSet[i]->size();
 }
}

template<class Scalar>
GenDistrVectorSet<Scalar>::GenDistrVectorSet(int _numVectors, DistrInfo &dinfo, Scalar* vectors, bool _myMemory)
{
 numVectors = _numVectors;
 myMemory = _myMemory;
 vecSet = new GenDistrVector<Scalar> * [numVectors];
 length = numVectors*dinfo.totLen();
 int i,j;
 for(i=0; i<numVectors; ++i) {
   if(myMemory) {
     vecSet[i] = new GenDistrVector<Scalar>(dinfo);
     for(j=0; j<dinfo.totLen(); ++j) (*vecSet[i])[j] = vectors[j + i*dinfo.totLen()];
   }
   else {
     vecSet[i] = new GenStackDistVector<Scalar>(dinfo,vectors+i*dinfo.totLen());
   }
 }
 zero();
}

template<class Scalar>
void
GenDistrVectorSet<Scalar>::copy(const GenDistrVector<Scalar> &v)
{
 int i;
 for(i=0; i<numVectors; ++i)
   (*vecSet[i]) = v;
}

template<class Scalar>
void
GenDistrVectorSet<Scalar>::copy(const GenDistrVectorSet<Scalar> &vectorSet)
{
 int i;
 for(i=0; i<numVectors; ++i)
   (*vecSet[i]) = vectorSet[i];
}

template<class Scalar>
void 
GenDistrVectorSet<Scalar>::print(const char *message)
{
   if(message) fprintf(stderr,"%s\n",message);
   int i,j;
   for(i=0; i<numVectors; ++i) {
     for(j=0; j<length; ++j)
       fprintf(stderr,"%e\n",(*vecSet[i])[j]);
     fprintf(stderr,"vector %d\n",i+1);
   }

}

template<class Scalar>
void 
GenDistrVectorSet<Scalar>::zero()
{
   for(int i=0;i<numVectors; ++i)
      (*vecSet[i]) = 0.0;
}     
