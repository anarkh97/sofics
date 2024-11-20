#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Math.d/SparseMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/SparseSet.h>

template<class Scalar>
GenSparseSet<Scalar>::~GenSparseSet()
{
 // if(sm) { delete [] sm; sm=0; }
}

template<class Scalar> 
GenSparseSet<Scalar>::GenSparseSet(int num)
 : sm(nullptr, num)
{
 numSM = 0;
}

template<class Scalar>
int
GenSparseSet<Scalar>::addSparseMatrix(std::shared_ptr<GenSparseMatrix<Scalar>> _sm)
{ 
 sm[numSM++] = _sm;
 return numSM-1;
}

template<class Scalar> 
void
GenSparseSet<Scalar>::setSparseMatrix(int i, std::shared_ptr<GenSparseMatrix<Scalar>>_sm)
{
 sm[i] = _sm;
 if(i > numSM) numSM = i+1; // this shouldn't happen
}

template<class Scalar> 
int
GenSparseSet<Scalar>::numRow()
{
 return (sm[0]) ? sm[0]->numRow() : 0;
}

template<class Scalar> 
int
GenSparseSet<Scalar>::numCol()
{
 int i,nc=0;
 for(i=0; i<numSM; ++i) 
   nc += (sm[i]) ? sm[i]->numCol(): 0;
 return nc;
}

template<class Scalar> 
void
GenSparseSet<Scalar>::mult(Scalar *rhs, Scalar *result)
{
 int i;
 int colOffset=0;
 for(i=0; i<numSM; ++i) {
   if(sm[i]) sm[i]->mult(rhs, result+colOffset);
   colOffset += (sm[i]) ? sm[i]->numCol() : 0;
 }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::multAdd(Scalar *rhs, Scalar *result)
{
 int i;
 int colOffset=0;
 for(i=0; i<numSM; ++i) {
   if(sm[i]) sm[i]->multAdd(rhs, result+colOffset);
   colOffset += (sm[i]) ? sm[i]->numCol() : 0;
 }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::multIdentity(Scalar **result)
{
  int i;
  int colOffset=0;
  for(i=0; i<numSM; ++i) {
    if(sm[i]) sm[i]->multIdentity(result+colOffset);
    colOffset += (sm[i]) ? sm[i]->numCol() : 0;
  }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::multIdentity(Scalar *result)
{
  int i;
  int colOffset=0;
  for(i=0; i<numSM; ++i) {
    if(sm[i]) sm[i]->multIdentity(result+colOffset);
    colOffset += (sm[i]) ? sm[i]->numCol() : 0;
  }
}

template<class Scalar>
void
GenSparseSet<Scalar>::transposeMult(Scalar *rhs, Scalar *result)
{
  int i;
  int colOffset=0;
  for(i=0; i<numSM; ++i) {
    if(sm[i]) sm[i]->transposeMult(rhs+colOffset,result);
    colOffset += (sm[i]) ? sm[i]->numCol() : 0;
  }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::transposeMultAdd(Scalar *rhs, Scalar *result)
{
  int i;
  int colOffset=0;
  for(i=0; i<numSM; ++i) {
    if(sm[i]) sm[i]->transposeMultAdd(rhs+colOffset,result);
    colOffset += (sm[i]) ? sm[i]->numCol() : 0;
  }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::transposeMultSubtract(Scalar *rhs, Scalar *result)
{
  int i;
  int colOffset=0;
  for(i=0; i<numSM; ++i) {
    if(sm[i]) sm[i]->transposeMultSubtract(rhs+colOffset,result);
    colOffset += (sm[i]) ? sm[i]->numCol() : 0;
  }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::multSub(int nRHS, Scalar **rhs, Scalar **result)
{
  int i,n;
  int colOffset=0;
  for(i=0; i<numSM; ++i) {
    for(n=0; n<nRHS; ++n)
      if(sm[i]) sm[i]->multSub(rhs[n],result[n]+colOffset);
    colOffset += (sm[i]) ? sm[i]->numCol() : 0;
  }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::multIdentity(DComplex **result, int col1, int col2) {

 int i;
 int colOffset = 0;
 int colDone = 0;

 for (i=0; i<numSM; ++i) {

   if (sm[i]==0)
     continue;

   int buff = sm[i]->numCol(); 

   if ((col1 < colOffset+buff) && (colOffset < col2)) {

     int start = col1-colOffset;
     start = (start<0) ? 0 : start;
     int stop = col2-colOffset;
     stop = (stop>buff) ? buff : stop;

     sm[i]->multIdentity(result+colDone, start, stop);
     colDone += stop-start;
   }
   colOffset += buff;
 }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::transposeMultSubtract(int nRHS, DComplex **rhs, DComplex **result) 
{
  // Transpose - multiplication of a set of vectors & subtraction
  for(int n=0; n<nRHS; ++n)
     transposeMultSubtract(rhs[n], result[n]);
}

template<class Scalar> 
void
GenSparseSet<Scalar>::transposeMultSubtract(int nRHS, DComplex **rhs, FullMC *result, int col1) 
{
  // Transpose - multiplication of a set of vectors & subtraction
  int i, n;
  int nrow = result->numRow();
  DComplex *temp = (DComplex *) dbg_alloca(sizeof(DComplex)*nrow);
  for (n=col1; n<nRHS; ++n) {
     for (i=0; i<nrow; ++i)
        temp[i] = (*result)[i][n];
     transposeMultSubtract(rhs[n-col1], temp);
     for (i=0; i<nrow; ++i)
        (*result)[i][n] = temp[i]; 
  }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::transposeMultSubtract(int nRHS, DComplex **rhs, DComplex *result, int col1) 
{
  // Transpose - multiplication of a set of vectors & subtraction
  int n;
  int nrow = numCol();
  for (n=col1; n<nRHS; ++n) {
    transposeMultSubtract(rhs[n-col1], result+(n-col1)*nrow);
  }
}

template<class Scalar> 
void
GenSparseSet<Scalar>::multSub(DComplex *rhs, DComplex *result) 
{
  // Multiplication of one vector & subtraction
  int i;
  int colOffset=0;
  for(i=0; i<numSM; ++i) {
    if(sm[i]) sm[i]->multSubtract(rhs+colOffset,result);
    colOffset += (sm[i]) ? sm[i]->numCol() : 0;
  }
}

