#include <cstdio>
#include <iostream>
#include <Solvers.d/KProject.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Preconditioner.h>
#include <Math.d/Vector.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Math.d/NBSparseMatrix.h>

#include <Solvers.d/PCGSolver.h>

template<class Scalar,
         class AnyVector,
         class AnyOperator>
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::GenPCGSolver(AnyOperator *_A, int _precno, int _maxitr, double _tolpcg, int _maxVecStorage, int _verbose) 
: BasePCG<Scalar, AnyVector, AnyOperator, KrylovProjector<Scalar,AnyVector>, Preconditioner<AnyVector> >  
  (_maxitr, _tolpcg, _verbose, _A)
{
  precno     = _precno;
  kryflg     = 0;
  initflg    = 0;
  reorthoflg = 0;

  this->solveTime = 0.0;

  if(_maxVecStorage)
    this->proj = new KrylovProjector<Scalar,AnyVector>(_A->neqs(),_maxVecStorage); 
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
void
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::initPrec()
{
  if(this->prec == 0) {
    if(precno == 0)
      this->prec = new NullPreconditioner<AnyVector>();
    else if(precno == 1)
      this->prec = new DiagPrec<AnyVector,AnyOperator>(this->A);
    else if(precno == 2)
      this->prec = new ScalarBlockDiagPrec<Scalar,AnyVector,AnyOperator>(this->A);
  }
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
void
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::reSolve(AnyVector &rhs)
{
    const_cast<GenPCGSolver<Scalar, AnyVector, AnyOperator> *>(this)->initPrec();
    this->solveTime -= getTime(); this->memUsed -= memoryUsed();
    AnyVector sol(rhs.size());
    this->doSolve(rhs,sol);
    rhs = sol;
    this->solveTime += getTime(); this->memUsed += memoryUsed();
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
void
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::solve(const Scalar *rhs, Scalar *solution)
{
 std::cerr << "GenPCGSolver::solve(Scalar *rhs, Scalar *solution) is not implemented\n";
/*
 this->solveTime -= getTime(); this->memUsed -= memoryUsed();
 AnyVector sol(solution, this->A->dim());
 AnyVector f(rhs, this->A->dim());
 this->doSolve(f,sol);
 for(int i=0; i<this->A->dim(); ++i) solution[i] = sol[i];
 this->solveTime += getTime(); this->memUsed += memoryUsed();
*/
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
void
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::solve(const AnyVector &rhs, AnyVector &solution)
{
    const_cast<GenPCGSolver<Scalar, AnyVector, AnyOperator> *>(this)->initPrec();
    this->solveTime -= getTime(); this->memUsed -= memoryUsed();
    this->doSolve(rhs,solution);
    this->solveTime += getTime(); this->memUsed += memoryUsed();
    this->times.precond = this->prec->time;
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
void
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::reSolve(Scalar *rhs)
{
  std::cerr << "GenPCGSolver::reSolve(Scalar *rhs) const is not implemented\n";
/*
  this->solveTime -= getTime(); this->memUsed -= memoryUsed();
  AnyVector sol(rhs, this->A->dim());
  reSolve(sol);
  for(int i=0; i<this->A->dim(); ++i) rhs[i] = sol[i];
  this->solveTime += getTime(); this->memUsed += memoryUsed();
*/
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
void
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::reSolve(int nRHS, Scalar **RHS)
{
 std::cerr << "GenPCGSolver::reSolve(int nRHS, Scalar **RHS) const is not implemented\n";
/*
 this->solveTime -= getTime(); this->memUsed -= memoryUsed();
 int i,n;
 for(n=0; n<nRHS; ++n) {
   AnyVector sol(RHS[n], this->A->dim());
   AnyVector rhs(RHS[n], this->A->dim());
   this->doSolve(rhs, sol);
   for(i=0; i<this->A->dim(); ++i)
     RHS[n][i] = sol[i];
 }
 this->solveTime += getTime(); this->memUsed += memoryUsed();
*/
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
void
GenPCGSolver<Scalar, AnyVector, AnyOperator>
::reSolve(int nRHS, AnyVector *RHS)
{
    const_cast<GenPCGSolver<Scalar, AnyVector, AnyOperator> *>(this)->initPrec();
    this->solveTime -= getTime(); this->memUsed -= memoryUsed();
    AnyVector sol(RHS[0].size());
    for(int n=0; n<nRHS; ++n) {
        sol = RHS[n];
        this->doSolve(RHS[n],sol);
        RHS[n] = sol;
    }
    this->solveTime += getTime(); this->memUsed += memoryUsed();
}

template<class Scalar,
         class AnyVector,
         class AnyOperator>
int
GenPCGSolver<Scalar, AnyVector, AnyOperator>::neqs() const
{
  return BasePCG<Scalar,AnyVector,AnyOperator,KrylovProjector<Scalar,AnyVector>,Preconditioner<AnyVector> >::neqs();
}

