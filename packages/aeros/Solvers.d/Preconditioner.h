#ifndef _PRECONDITIONER_H_
#define _PRECONDITIONER_H_

template<class AnyVector>
class Preconditioner {
 public:
   double time;
   Preconditioner() { time = 0.0; }
   virtual ~Preconditioner() {}
   virtual void apply(AnyVector &r, AnyVector &pr) = 0; 
};

template<class AnyVector>
class NullPreconditioner : public Preconditioner<AnyVector> 
{
 public:
  void apply(AnyVector &r, AnyVector &pr) { pr = r; }
};

template<class AnyVector, class AnyMatrix>
class DiagPrec : public Preconditioner<AnyVector> 
{
  AnyVector diag; 
 public:
  DiagPrec(AnyMatrix *A);
  virtual ~DiagPrec() { }
  void apply(AnyVector &r, AnyVector &pr);
};

/*
template<class Scalar, class AnyVector, class AnyMatrix>
class BlockDiagPrec : public Preconditioner<AnyVector> 
{
  int     numnod;	// number of nodes
  int    *firstDof;	// int array of the nodes first dof
  GenFullM<Scalar> **diagBlock;    // array of pointers to FullM
 public:
  BlockDiagPrec(AnyMatrix *); // currently NBSparseMatrix/SfemBlockMatrix can be used here
  virtual ~BlockDiagPrec();
  void apply(AnyVector &r, AnyVector &pr);
};
*/

template<class Scalar, class AnyVector, class AnyMatrix>
class ScalarBlockDiagPrec : public Preconditioner<AnyVector> // Each block diagonal is scalar multiple of one matrix
{
  int     numnod;       // number of nodes
  int    *firstDof;     // int array of the nodes first dof
  int blocksize;
  AnyMatrix *A;
  //GenSolver<Scalar> *diagBlock; 
  Scalar *scalars;
 public:
  ScalarBlockDiagPrec(AnyMatrix *); 
  virtual ~ScalarBlockDiagPrec();
  void apply(AnyVector &r, AnyVector &pr);
};


#ifdef _TEMPLATE_FIX_
#include <Solvers.d/Preconditioner.C>
#endif

#endif
