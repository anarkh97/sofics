#ifndef __ADDEDMASSMATRIX_H__
#define __ADDEDMASSMATRIX_H__

#include <Math.d/DBSparseMatrix.h>
#include <Math.d/BLKSparseMatrix.h>

template<class Scalar, class ConstraintOperator>
class AddedMassMatrix : public GenDBSparseMatrix<Scalar>
{
  private: 
    GenBLKSparseMatrix<Scalar> *fluidSolver;
    ConstraintOperator *op;
    void (ConstraintOperator::*multC)(const GenVector<Scalar> &, GenVector<Scalar> &);
    void (ConstraintOperator::*trMultC)(const GenVector<Scalar> &, GenVector<Scalar> &);

  public:
    // Constructor
    AddedMassMatrix(const Connectivity* con, const DofSetArray* dsa, const ConstrainedDSA* c_dsa,
                    ConstraintOperator *_op,
                    void (ConstraintOperator::*_multC)(const GenVector<Scalar> &, GenVector<Scalar> &),
                    void (ConstraintOperator::*_trMultC)(const GenVector<Scalar> &, GenVector<Scalar> &)) 
     : op(_op), multC(_multC), trMultC(_trMultC), GenDBSparseMatrix<Scalar>(con, dsa, c_dsa) {}

    // Destructor
    ~AddedMassMatrix() { }

    void setFluidSolver(GenBLKSparseMatrix<Scalar> *_fluidSolver) { fluidSolver = _fluidSolver; }

    void mult(const GenVector<Scalar> &x, GenVector<Scalar> &y) { 
      GenVector<Scalar> tmp(fluidSolver->neqs());
      GenVector<Scalar> y_added(y);

      GenDBSparseMatrix<Scalar>::mult(x, y); // y = M*x

      (op->*trMultC)(x, tmp);               // tmp = C^T*x
      fluidSolver->reSolve(tmp);            // tmp = F^{-1}*C^T*x
      (op->*multC)(tmp, y_added);           // y_added = C*F^{-1}*C^T*x            

      y += y_added;                         // y = (M + C*F^{-1}*C^T)*x
    }
};

#endif
