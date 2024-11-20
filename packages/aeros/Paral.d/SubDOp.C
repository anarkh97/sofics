#include <Feti.d/DistrVector.h>
#include <Math.d/SparseMatrix.h>
#include <Paral.d/Assembler.h>
#include <iostream>

template<class Scalar>
GenSubDOp<Scalar>::~GenSubDOp() 
{ 
  if(sops) { 
    for(int i=0; i<numSub; ++i) if(sops[i]) delete sops[i]; 
    delete [] sops; 
  } 
}

template<class Scalar>
void
GenSubDOp<Scalar>::mult(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &opV)
{
  SubDOpMult<Scalar> sapply(sops, &v, &opV);
  threadManager->execParal(numSub, &sapply);
  if(assembler) assembler->assemble(opV);
}

template<class Scalar>
void
GenSubDOp<Scalar>::transposeMult(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &opV)
{
  SubDOpTransposeMult<Scalar> sapply(sops, &v, &opV);
  threadManager->execParal(numSub, &sapply);
  if(assembler) assembler->assemble(opV);
}

template<class Scalar>
void
GenSubDOp<Scalar>::squareRootMult(GenDistrVector<Scalar> &v)
{
  SubDOpSquareRootMult<Scalar> sapply(sops, &v);
  threadManager->execParal(numSub, &sapply);
  if(assembler) assembler->assemble(v);
}

template<class Scalar>
void
GenSubDOp<Scalar>::inverseSquareRootMult(GenDistrVector<Scalar> &v)
{
  SubDOpInverseSquareRootMult<Scalar> sapply(sops, &v);
  threadManager->execParal(numSub, &sapply);
  if(assembler) assembler->assemble(v);
}


template<class Scalar>
void
GenSubDOp<Scalar>::multAdd(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &opV)
{
  SubDOpMultAdd<Scalar> sapply(sops, &v, &opV);
  threadManager->execParal(numSub, &sapply);
  if(assembler) assembler->assemble(opV);
}

template<class Scalar>
void
GenSubDOp<Scalar>::multNoAssemble(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &opV)
{
  SubDOpMult<Scalar> sapply(sops, &v, &opV);
  threadManager->execParal(numSub, &sapply);
}

template<class Scalar>
void
GenSubDOp<Scalar>::assemble(GenDistrVector<Scalar> &v)
{
  if(assembler) assembler->assemble(v);
}

template<class Scalar>
void 
GenSubDOp<Scalar>::invertDiag()
{
  if(assembler) {
    //use multInvertDiag with a DistrVector that has only 1 as entries
    // or do something else
    std::cerr << "pb in SubDOp::invertDiag() " << std::endl;
  }
  else {
    SubDOpInvertDiag<Scalar> sapply(sops);
    threadManager->execParal(numSub, &sapply);
  }
}

template<class Scalar>
void
GenSubDOp<Scalar>::multInvertDiag(GenDistrVector<Scalar> &z)
{
  SubDOpMultInvertDiag<Scalar> sapply(sops, &z);
  threadManager->execParal(numSub, &sapply);
  if(assembler) assembler->assemble(z, 1);
}

template<class Scalar>
void
GenSubDOp<Scalar>::zeroAll()
{
  SubDOpZeroAll<Scalar> sapply(sops);
  threadManager->execParal(numSub, &sapply);
}

template<class Scalar>
void
SubDOpMult<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  // Get the pointer to the part of the vector s corresponding to subdomain sNum
  Scalar *src = s->subData(sNum);
  // Get the pointer to the part of the vector tg corresponding to subdomain sNum
  Scalar *tg  = d->subData(sNum);
  // multiply for subdomain sNum
  sm[sNum]->mult(src,tg);
}

template<class Scalar>
void
SubDOpTransposeMult<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  // Get the pointer to the part of the vector s corresponding to subdomain sNum
  Scalar *src = s->subData(sNum);
  // Get the pointer to the part of the vector tg corresponding to subdomain sNum
  Scalar *tg  = d->subData(sNum);
  // multiply for subdomain sNum
  sm[sNum]->transposeMult(src,tg);
}

template<class Scalar>
void
SubDOpSquareRootMult<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  // Get the pointer to the part of the vector s corresponding to subdomain sNum
  Scalar *src = s->subData(sNum);
  // multiply for subdomain sNum
  sm[sNum]->squareRootMult(src);
}

template<class Scalar>
void
SubDOpInverseSquareRootMult<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  // Get the pointer to the part of the vector s corresponding to subdomain sNum
  Scalar *src = s->subData(sNum);
  // multiply for subdomain sNum
  sm[sNum]->inverseSquareRootMult(src);
}


template<class Scalar>
void
SubDOpMultAdd<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  // Get the pointer to the part of the vector s corresponding to subdomain sNum
  Scalar *src = s->subData(sNum);
  // Get the pointer to the part of the vector tg corresponding to subdomain sNum
  Scalar *tg  = d->subData(sNum);
  // multiply for subdomain sNum
  sm[sNum]->multAdd(src,tg);
}

template<class Scalar>
void
SubDOpMultInvertDiag<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  // Get the pointer to the part of the vector s corresponding to subdomain sNum
  Scalar *src = z->subData(sNum);
  // multiply for subdomain sNum
  for(int j=0; j< sm[sNum]->dim(); ++j)
    src[j] /= sm[sNum]->diag(j);
}

template<class Scalar>
void
SubDOpInvertDiag<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  sm[sNum]->invertDiag();
}

template<class Scalar>
void
SubDOpZeroAll<Scalar>::runFor(int sNum)
{
  if(!sm[sNum]) return;
  sm[sNum]->zeroAll();
}

