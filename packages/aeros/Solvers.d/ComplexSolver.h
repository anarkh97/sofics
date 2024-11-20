#ifndef _COMPLEXSOLVER_H_
#define _COMPLEXSOLVER_H_

#include <Utils.d/MyComplex.h>

class ComplexVector;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class FullMC;
class Rbm;

class ComplexSolver {

  public:

    virtual int neqs() const = 0;
    virtual int dim() const = 0;

    virtual int numRBM();

    // complex solve functions
    virtual void solve(DComplex *rhs, DComplex *solution);
    virtual void solve(ComplexVector &rhs, ComplexVector &solution);

    // reSolve functions overwriting the rhs vector with the solution
    virtual void reSolve(DComplex *rhs);
    virtual void reSolve(ComplexVector &v);

    // Multiple rhs reSolve functions
    virtual void reSolve(int nRHS, DComplex **);
    virtual void reSolve(int nRHS, DComplex *);
    virtual void reSolve(int nRHS, ComplexVector *);
    virtual void reSolve(FullMC *);

    virtual void factor();

    // function to return solution time
    virtual double getSolutionTime();

    // function to return memory used
    virtual long size();

};

#endif
