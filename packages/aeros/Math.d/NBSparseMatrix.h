#ifndef NB_SPARSEMATRIX_H_
#define NB_SPARSEMATRIX_H_

#include <Utils.d/MyComplex.h>
#include <Math.d/SparseMatrix.h>

template <class Scalar> class GenSparseMatrix;
class Connectivity;
template <class Scalar> class GenFullM;
class ConstrainedDSA;
template <class Scalar> class GenVector;
template <class Scalar> class GenFullSquareMatrix;

// Node Based Sparse Matrix

template<class Scalar>
class GenNBSparseMatrix : public GenSparseMatrix<Scalar>
{
   const Connectivity *con;     //!< \brief Node to node connectivity.
   GenFullM<Scalar> *allM;
   int *firstDof;               //!< \brief Array of the first DOF number for a node.
   const ConstrainedDSA *dsa;
   int *dofToNode;
   int numnodes;
   int *bc;

 public:
   GenNBSparseMatrix(const Connectivity *con, const ConstrainedDSA *c_dsa);
   virtual ~GenNBSparseMatrix();   
 
   void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result ) const override; // Matrix-Vector multiply
   void mult(const Scalar *rhs, Scalar *result ) const override;   // Matrix-Vector multiply
   void multAdd(const Scalar *rhs, Scalar *result) const override;
   Scalar diag(int n) const override;                 // returns diagonal values of matrix
   Scalar &diag(int n) override;                 // returns diagonal values of matrix
   void add(const FullSquareMatrix &m, const int *dofs) override;
   void add(const FullSquareMatrixC &m, const int *dofs) override;
   void addImaginary(const FullSquareMatrix &m, const int *dofs) override;
   void zeroAll() override;

   // returns number of unconstrained dof
   int  dim() const override { return dsa->size(); }
   int  neqs() const override { return dsa->size(); }

   // returns a pointer to a diagonal matrix of node i
   GenFullM<Scalar>* getDiagMatrix(int i) override;
   int*   getFirstDof() override { return firstDof; }
   int    numNodes() const override {return numnodes; }
};

typedef GenNBSparseMatrix<double> NBSparseMatrix;
typedef GenNBSparseMatrix<DComplex> NBComplexSparseMatrix;

#endif
