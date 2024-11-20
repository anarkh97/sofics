#include <cstdio>

#include <Solvers.d/ComplexSolver.h>
#include <Math.d/ComplexVector.h>

int
ComplexSolver::numRBM()
{
 return 0;
}

void
ComplexSolver::reSolve(FullMC *) const {

 fprintf(stderr,"Selected reSolve does not support FullMC matrix input.\n");

}


void
ComplexSolver::reSolve(int nRHS, DComplex **RHS) const
{
 // nRHS = number of rhs
 int i;
 for (i = 0; i < nRHS; ++i)
    reSolve(RHS[i]);
}

void
ComplexSolver::reSolve(int nRHS, DComplex *RHS) const {

 fprintf(stderr,"Selected reSolve does not support DComplex* input.\n");

}

void
ComplexSolver::reSolve(int nRHS, ComplexVector *RHS) const
{
 int i;
 for (i = 0; i < nRHS; ++i)
   reSolve(RHS[i]);
}

void
ComplexSolver::reSolve(ComplexVector &v) const
{
 reSolve(v.data());
}

void
ComplexSolver::reSolve(DComplex*) const
{
 fprintf(stderr,"Selected reSolve does not support complex.\n");
}

void
ComplexSolver::solve(DComplex *, DComplex *)
{
 fprintf(stderr,"Selected Solver does not support complex.\n");
}

void
ComplexSolver::solve(ComplexVector &, ComplexVector &)
{
 fprintf(stderr,"Selected Solver does not support a Complex Vector.\n");
}

void
ComplexSolver::factor()
{
 
}

double
ComplexSolver::getSolutionTime()
{
 return 0.0; 
}

long
ComplexSolver::size()
{
 return 0;
}
