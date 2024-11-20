#include <cstdio>
#include <Solvers.d/Solver.h>
template<>
void
GenSolver<DComplex>::reSolve(ComplexVector &v)
{
 reSolve(v.data());
}

template<>
void
GenSolver<double>::reSolve(ComplexVector &v)
{
  fprintf(stderr, "WARNING: GenSolver<double>::reSolve(ComplexVector &v) const not implemented \n");
}

template<>
void
GenSolver<double>::reSolve(Vector &v)
{
 reSolve(v.data());
}

template<>
void
GenSolver<DComplex>::reSolve(Vector &v)
{
  fprintf(stderr, "WARNING: GenSolver<DComplex>::reSolve(Vector &v) const not implemented \n");
}

