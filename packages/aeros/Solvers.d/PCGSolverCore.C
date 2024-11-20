#include <cstdio>
#include <Solvers.d/KProject.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Solvers.d/Preconditioner.h>
#include <Math.d/Vector.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Solvers.d/PCGSolver.h>

template<>
void
GenPCGSolver<double>::reSolve(GenVector<double>& rhs) const
{
  solveTime -= getTime();
  Vector sol(rhs.size());
  doSolve(rhs,sol);
  rhs = sol;
  solveTime += getTime();
  memUsed += memoryUsed();
}

template<>
GenPCGSolver<double>::GenPCGSolver(SparseMatrix *_dynamK, int _precno, int _maxitr, 
                                   double _tolpcg, Rbm *_rbm)
: BasePCG<Vector, SparseMatrix, KrylovProjector, Preconditioner>  
  (_maxitr, _tolpcg, _dynamK, 0,0)
{
  if(_precno == 0)
    prec = new NullPreconditioner;
  else if(_precno == 1) {
    prec = new DiagPrec(_dynamK);
  }

  kryflg     = 0;
  initflg    = 0;
  reorthoflg = 0;
  numrbm     = 0;

  // initialize rbm information
  rbm    = _rbm;
  if(rbm)
    numrbm = rbm->numRBM();
  else
    numrbm = 0;

  solveTime = 0.0;

  memUsed = -memoryUsed();
}

template<>
GenPCGSolver<DComplex>::GenPCGSolver(ComplexSparseMatrix *_dynamK, int _precno, int _maxitr,
                                     double _tolpcg, Rbm *_rbm)
: BasePCG<ComplexVector, ComplexSparseMatrix, KrylovProjector, Preconditioner>
  (_maxitr, _tolpcg, _dynamK, 0,0)
{ 
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

/*
template<>
GenPCGSolver<double>::GenPCGSolver(NBSparseMatrix *AA, int _precno, int _maxitr, 
                                  double _tolpcg, Rbm *_rbm)
: BasePCG<Vector, SparseMatrix, KrylovProjector,Preconditioner> 
  (_maxitr, _tolpcg, AA, 0,0)
{
  if(_precno == 0)
    prec = new NullPreconditioner;
  else if(_precno == 1)
//    prec = new BlockDiagPrec(AA); // BlockDiagPrec no longer in use
    prec = new DiagPrec(AA);

  kryflg     = 0;
  initflg    = 0;
  reorthoflg = 0;
  numrbm     = 0;

  // initialize rbm information
  rbm    = _rbm;
  if(rbm)
    numrbm = rbm->numRBM();
  else
    numrbm = 0;

  solveTime = 0.0;
  memUsed = -memoryUsed();
}

template<>
GenPCGSolver<DComplex>::GenPCGSolver(NBComplexSparseMatrix *AA, int _precno, int _maxitr,
                                     double _tolpcg, Rbm *_rbm)
: BasePCG<ComplexVector, ComplexSparseMatrix, KrylovProjector, Preconditioner>
  (_maxitr, _tolpcg, AA, 0,0)
{ 
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}
*/

template<>
void
GenPCGSolver<double>::solve(double *rhs, double *solution)
{
 StackVector sol(solution, A->dim());
 StackVector f(rhs, A->dim());

 solveTime -= getTime();
 doSolve(f,sol);
 solveTime += getTime();
 memUsed += memoryUsed();
}

template<>
void
GenPCGSolver<DComplex>::solve(DComplex *rhs, DComplex *solution)
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<double>::solve(Vector &rhs, Vector &solution)
{
 solveTime -= getTime();
 doSolve(rhs,solution);
 solveTime += getTime();
 memUsed += memoryUsed();
}

template<>
void
GenPCGSolver<DComplex>::solve(ComplexVector &rhs, ComplexVector &solution)
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<double>::reSolve(double *rhs) const
{
  solveTime -= getTime();
  StackVector sol( rhs, A->dim() );
  reSolve(sol);
  memUsed += memoryUsed();
  solveTime += getTime();
  memUsed += memoryUsed();
}

template<>
void
GenPCGSolver<DComplex>::reSolve(DComplex *rhs) const
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<DComplex>::reSolve(ComplexVector& rhs) const
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<double>::reSolve(int nRHS, double **RHS) const
{
 solveTime -= getTime();

 int i,n;
 for(n=0; n<nRHS; ++n) {
   Vector sol(RHS[n], A->dim());
   Vector rhs(RHS[n], A->dim());
   doSolve(rhs, sol);
   for(i=0; i<A->dim(); ++i)
     RHS[n][i] = sol[i];
 }

 solveTime += getTime();
 memUsed += memoryUsed();
}

template<>
void
GenPCGSolver<DComplex>::reSolve(int nRHS, DComplex **RHS) const
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<double>::reSolve(int nRHS, Vector *RHS) const
{
 fprintf(stderr,"Begin pcg algorithm\n");
 solveTime -= getTime();
 Vector sol(RHS[0].size());
 int n;
 for(n=0; n<nRHS; ++n) {
   sol = RHS[n];
   doSolve(RHS[n],sol);
   RHS[n] = sol;
 }
 solveTime += getTime();
 memUsed += memoryUsed();
}

template<>
void
GenPCGSolver<DComplex>::reSolve(int nRHS, ComplexVector *RHS) const
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<double>::getRBMs(double *rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

template<>
void
GenPCGSolver<DComplex>::getRBMs(double *rigidBodyModes)
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<double>::getRBMs(Vector *rigidBodyModes)
{
  rbm->getRBMs(rigidBodyModes);
}

template<>
void
GenPCGSolver<DComplex>::getRBMs(Vector *rigidBodyModes)
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
void
GenPCGSolver<double>::getRBMs(VectorSet& rigidBodyModes)
{
 rbm->getRBMs(rigidBodyModes);
}

template<>
void
GenPCGSolver<DComplex>::getRBMs(VectorSet& rigidBodyModes)
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n"); 
}

template<>
int
GenPCGSolver<double>::neqs() const
{
  return BasePCG<Vector,SparseMatrix,KrylovProjector,Preconditioner>::neqs();
}

template<>
int
GenPCGSolver<DComplex>::neqs() const
{
  fprintf(stderr, "GenPCGSolver<DComplex> is not implemented \n");
  return 0;
}
