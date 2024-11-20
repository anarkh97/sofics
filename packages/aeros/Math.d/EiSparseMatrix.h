#ifndef EI_SPARSEMATRIX_H_
#define EI_SPARSEMATRIX_H_

#ifdef USE_EIGEN3
#include <complex>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#ifdef EIGEN_CHOLMOD_SUPPORT
#include <Eigen/CholmodSupport>
#endif
#ifdef EIGEN_UMFPACK_SUPPORT
#include <Eigen/UmfPackSupport>
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif
#ifdef EIGEN_SPQR_SUPPORT
#include <Eigen/SPQRSupport>
#endif
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>

template<class Scalar, class SolverClass = Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >
class GenEiSparseMatrix : public SparseData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
 protected:
   // this is a symmetric sparse matrix using CSR storage (upper triangluar part only in self-adjoint case) 
   // and Eigen 3 implementation via MappedSparseMatrix
   bool selfadjoint;
   int nnz;
   Scalar *unonz;
   Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor, int> M;
   SolverClass solver;
   Eigen::SparseMatrix<Scalar> *M_copy;

 public:
   GenEiSparseMatrix(const Connectivity *, const DofSetArray *, const DofSetArray *, bool= true);
   GenEiSparseMatrix(const Connectivity *, const DofSetArray *, const int *, bool=true);
   GenEiSparseMatrix(const Connectivity *, const EqNumberer *, bool=true);
   virtual ~GenEiSparseMatrix();

   Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor, int>& getEigenSparse() { return M; }
   SolverClass& getEigenSolver() { return solver; }

   // GenSparseMatrix assembly
   void add(const FullSquareMatrix &, const int *dofs) override;
   void addCoef(int, int, Scalar);
   void add(const GenAssembledFullM<Scalar> &, const int *) override;
   void addImaginary(const FullSquareMatrix &, const int *dofs) override;
   void add(const FullSquareMatrixC &, const int *dofs) override;

   // GenSparseMatrix matrix-vector multiplications
   void mult(const GenVector<Scalar> &, GenVector<Scalar> &) const override;
   void mult(const GenVector<Scalar> &rhs, Scalar *result) const override;
   void mult(const Scalar *, Scalar *) const override;
   void multAdd(const Scalar *, Scalar *) const override;
   void transposeMult(const GenVector<Scalar> & rhs, GenVector<Scalar> & result) const override;
   void transposeMult(const Scalar *, Scalar *) const override;
   void upperMult(Scalar* result) override;
   void backward(Scalar* result) override;

   // GenSparseMatrix miscellaneous
   void zeroAll() override;
   int dim() const override { return numUncon; }
   double getMemoryUsed() const override;
   int numRow() const override { return numUncon; }
   int numCol() const override { return numUncon; }
   Scalar diag(int dof) const override;
   Scalar &diag(int dof) override;

   // GenSolver factor
   void factor() override;

    int numRBM() const override { return 0; }

   // GenSolver solve
   void solve(const Scalar *rhs, Scalar *solution) override;
   void solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution ) override;
   void reSolve(Scalar *rhs) override ;
   void reSolve(GenVector<Scalar> &rhs) override;
/* TODO
   void reSolve(int nRHS, Scalar **RHS);
   void reSolve(int nRHS, GenVector<Scalar> * RHS);
   void reSolve(GenFullM<Scalar> *);
*/
   // GenSolver miscellaneous
   int  neqs() const override { return numUncon; }
   long size() const override;
   void unify(FSCommunicator *communicator) override;
   void print() override;
   void printSparse(const std::string& filename) override;
};


template<class Scalar, class SolverClass = Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >
class WrapEiSparseMat : public GenEiSparseMatrix<Scalar, SolverClass>
{
  public:
    struct CtorData {
      const Connectivity *cn;
      const DofSetArray *dsa;
      const DofSetArray *cdsa;
      bool flg;
      CtorData(const Connectivity *c, const DofSetArray *d, const DofSetArray *dc, bool f=true) {
        cn = c;
        dsa = d;
        cdsa = dc;
        flg = f;
      }
    };

    WrapEiSparseMat(CtorData &ctd) : GenEiSparseMatrix<Scalar,SolverClass>(ctd.cn, ctd.dsa, ctd.cdsa, ctd.flg) {}
};

typedef GenEiSparseMatrix<double> EiSparseMatrix;
typedef GenEiSparseMatrix<DComplex> EiComplexSparseMatrix;


#ifdef _TEMPLATE_FIX_
#include <Math.d/EiSparseMatrix.C>
#endif

#endif
#endif
