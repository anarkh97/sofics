// Diagonal Preconditioner
// used with the dof based sparse matrix (DBSparseMatrix)
template<class AnyVector, class AnyMatrix>
DiagPrec<AnyVector, AnyMatrix>::DiagPrec(AnyMatrix *A) : diag(A->dim())
{
// int size = A->dim();
 int size = A->neqs();
 int i;
 for(i=0; i<size; ++i)  
 diag[i] = (1.0/A->diag(i));
}

template<class AnyVector, class AnyMatrix>
void
DiagPrec<AnyVector, AnyMatrix>::apply(AnyVector &r, AnyVector &pr)
{
 this->time -= getTime();
 int i;
 int size = r.size();
 for(i=0; i<size; ++i)  
   pr[i] = r[i]*diag[i]; 
 this->time += getTime();
}

/*
// Block Diagonal Preconditioner
// used with the node based sparse matrix (NBSparseMatrix)
template<class Scalar, class AnyVector, class AnyMatrix>
BlockDiagPrec<Scalar, AnyVector, AnyMatrix>::BlockDiagPrec(AnyMatrix *A)
{
 firstDof  = A->getFirstDof();
 numnod    = A->numNodes();
 diagBlock = new GenFullM<Scalar>*[numnod];

 int i;
 for(i=0; i<numnod; ++i) {
   GenFullM<Scalar> *diagMat = A->getDiagMatrix(i);
   diagBlock[i]   = (diagMat == 0) ? 0 : new GenFullM<Scalar>(*diagMat);
   if(diagBlock[i])
     diagBlock[i]->factor();
 } 
}

template<class Scalar, class AnyVector, class AnyMatrix>
BlockDiagPrec<Scalar, AnyVector, AnyMatrix>::~BlockDiagPrec()
{
 if(diagBlock && numnod) {
   for(int i=0; i<numnod; ++i)
     if(diagBlock[i]) { delete diagBlock[i]; diagBlock[i] = 0; }
   delete [] diagBlock; diagBlock = 0;
 }
}

template<class Scalar, class AnyVector, class AnyMatrix>
void
BlockDiagPrec<Scalar, AnyVector, AnyMatrix>::apply(AnyVector &r, AnyVector &pr)
{
 this->time -= getTime();
 pr = r;
 Scalar *prd = pr.data();
 int i;
 for(i=0; i<numnod; ++i) {
   diagBlock[i]->reSolve(prd+firstDof[i]);
 }
 this->time += getTime();
}
*/

// Sfem Block Diagonal Preconditioner
template<class Scalar, class AnyVector, class AnyMatrix>
ScalarBlockDiagPrec<Scalar, AnyVector, AnyMatrix>::ScalarBlockDiagPrec(AnyMatrix *_A)
{
 A = _A;
 firstDof  = A->getFirstDof();
 numnod    = A->numNodes();
 scalars = A->getBlockScalarMultipliers();
}

template<class Scalar, class AnyVector, class AnyMatrix>
ScalarBlockDiagPrec<Scalar, AnyVector, AnyMatrix>::~ScalarBlockDiagPrec()
{
  delete firstDof;
  delete scalars;
}

template<class Scalar, class AnyVector, class AnyMatrix>
void
ScalarBlockDiagPrec<Scalar, AnyVector, AnyMatrix>::apply(AnyVector &r, AnyVector &pr)
{
 this->time -= getTime(); 
 for(int i=0; i<numnod; ++i) {
   if(r.isnnz(i)==1) {   
     if(r.getBlock(i).sqNorm()==0)  pr.getBlock(i).zero(); 
     else {
       A->getMeanSolver()->solve(r.getBlock(i),pr.getBlock(i));
       pr.scaleBlock(i,scalars[i]); 
     }
   }
 }
 this->time += getTime(); 
}

