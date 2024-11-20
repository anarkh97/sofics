#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/Vector.h>

template<class Scalar> 
void
GenSolver<Scalar>::reSolve(int nRHS, Scalar **RHS)
{
 int i;
 for (i = 0; i < nRHS; ++i)
    reSolve(RHS[i]);
}

template<class Scalar> 
void
GenSolver<Scalar>::reSolve(int nRHS, Scalar *rhs)
{

  int numUncon=neqs();
  Scalar **rhsP = (Scalar**) dbg_alloca(sizeof(Scalar*)*nRHS);
  int i;
  for(i=0; i<nRHS; ++i)
    rhsP[i] = rhs+i*numUncon;
  reSolve(nRHS,rhsP);
}

template<class Scalar> 
void
GenSolver<Scalar>::reSolve(int nRHS, GenVector<Scalar> *RHS)
{ 
 int i;
 for (i = 0; i < nRHS; ++i)
   reSolve(RHS[i]);
}

template<>
void
GenSolver<DComplex>::reSolve(ComplexVector &v);

template<>
void
GenSolver<double>::reSolve(ComplexVector &v);

template<>
void
GenSolver<double>::reSolve(Vector &v);

template<>
void
GenSolver<DComplex>::reSolve(Vector &v);

template<class Scalar> 
void
GenSolver<Scalar>::reSolve(Scalar*)
{
 fprintf(stderr,"Selected solver does not support reSolve(Scalar*) function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::reSolve(GenFullM<Scalar> *)
{
 fprintf(stderr,"Selected solver does not support reSolve(GenFullM<Scalar> *) function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::forward(GenVector<Scalar> &rhs)
{
 fprintf(stderr,"Selected solver does not support forward(GenVector<Scalar> &) function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::forward(Scalar*)
{
 fprintf(stderr,"Selected solver does not support forward(Scalar*) function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::backward(GenVector<Scalar> &rhs)
{
 fprintf(stderr,"Selected solver does not support backward(GenVector<Scalar> &) function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::backward(Scalar*)
{
 fprintf(stderr,"Selected solver does not support backward(Scalar*) function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::upperMult(Scalar*)
{
 fprintf(stderr,"Selected solver does not support upperMult(Scalar*) function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::lowerMult(Scalar*)
{
 fprintf(stderr,"Selected solver does not support lowerMult(Scalar*) function\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::reBuild(FullSquareMatrix *, int, int)
{
 fprintf(stderr,"Selected Solver does not support non-linear analysis.\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::reBuild(FullSquareMatrix *, FullSquareMatrix *, Scalar)
{
 fprintf(stderr,"Selected Solver does not support non-linear analysis.\n");
}

template<class Scalar>
int
GenSolver<Scalar>::dim() const
{
 fprintf(stderr,"Selected Solver does not support dim() function\n");
 return 0;
}
#include <typeinfo>
template<class Scalar> 
int
GenSolver<Scalar>::numRBM() const
{
    std::cerr << "Solver " << typeid(*this).name() << " does not support numRBM() function" << std::endl;
// fprintf(stderr,"Selected Solver does not support numRBM() function\n");
 return 0;
}

template<class Scalar> 
void
GenSolver<Scalar>::getRBMs(double *)
{
 fprintf(stderr,"Selected Solver does not support getRBMs(double *) function\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::getRBMs(Vector *)
{
 fprintf(stderr,"Selected Solver does not support getRBMs(Vector *) function\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::getRBMs(VectorSet &)
{
 fprintf(stderr,"Selected Solver does not support getRBMs(VectorSet &) function\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::solve(const GenVector<Scalar>  &rhs, GenVector<Scalar>  &sol)
{
 solve(rhs.data(), sol.data());
}

template<class Scalar> 
void
GenSolver<Scalar>::solve(const Scalar *, Scalar *)
{
 fprintf(stderr,"Selected Solver does not support solve(Scalar*, Scalar*) function\n");
}

// This function does nothing, it is supposed to do nothing.
// When a solver does not have a factoring step, i.e. pcg, bcg
// or frontal, factoring is not used.

template<class Scalar> 
void
GenSolver<Scalar>::factor()
{
 fprintf(stderr,"Selected Solver does not implement factor()\n"); 
}

template<class Scalar>
void
GenSolver<Scalar>::parallelFactor()
{
 factor();
}


template<class Scalar> 
void
GenSolver<Scalar>::reBuildGeometricRbms(GeomState *)
{
 fprintf(stderr,"Selected Solver does not support reBuilding Geomtric Rbms\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::clean_up()
{
}

template<class Scalar> 
void 
GenSolver<Scalar>::addBoeing(int, const int *, const int *, const double *, const int *, Scalar multiplier)
{
 fprintf(stderr,"Selected Solver does not support addBoeing function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::addone(Scalar d, int dofi, int dofj)
{
  fprintf(stderr,"Selected Solver does not support addone function\n");
}

template<class Scalar>
Scalar
GenSolver<Scalar>::getone(int dofi, int dofj)
{
  fprintf(stderr,"Selected Solver does not support getone function\n");
  return 0;
}

template<class Scalar>
void
GenSolver<Scalar>::unify(FSCommunicator *)
{
  fprintf(stderr,"Selected Solver does not support unify function\n");
}

template<class Scalar>
void
GenSolver<Scalar>::add(Scalar *d)
{
  fprintf(stderr,"Selected Solver does not support add(double *) function\n");
}

template<class Scalar>
Scalar*
GenSolver<Scalar>::getData()
{
  fprintf(stderr,"Selected Solver does not support getData() function\n");
  return 0;
}
