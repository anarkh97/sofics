#ifndef BLOCK_SPARSE_MATRIX_H_
#define BLOCK_SPARSE_MATRIX_H_

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class Rbm;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
class DofSetArray;
class DofSetArray;
class Connectivity;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class SolverCntl;

#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h> 
#include <Solvers.d/Solver.h>

template<class Scalar>
class GenBLKSparseMatrix :
        public SparseData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
protected:
   const SolverCntl& scntl;

   Scalar *unonz;

   int nsuper;
   int nsub; 
   int nnzl; 
   int tmpsiz;

   int *colcnt;
   int *snode;
   int *xsuper;
   int *invsuper;

   int *adj; 
   int *xadj;

   int *iwork;
   int iwsiz;

   int *xlindx;
   int *lindx;

   int  *xlnz;
   Scalar *lnz;

   int *perm;
   int *invp;

   double tol;
   
   Rbm *rbm;   // pointer to rigid body modes
   int numrbm; // number of rigid body modes
   int ngrbm;  // number of geometric rbms
   int defblk; // size of last deficient block
   int lbdef; 
   int *iprow;
   int *ipcol;
   int *def;

   bool myRbm;

 public:

   // Stiffness matrix
   GenBLKSparseMatrix(const Connectivity *, const DofSetArray *, const DofSetArray *,
                      double tolerance, const SolverCntl &_scntl, Rbm *rbm = 0);
   // Kii
   GenBLKSparseMatrix(const Connectivity *, const DofSetArray *, int *dofmap,
                      double tolerance, const SolverCntl &_scntl);
   // GtG, Kcc
   GenBLKSparseMatrix(const Connectivity *, const EqNumberer *, double tolerance, const SolverCntl &_scntl,
                      int ngrbm = 0);

   virtual ~GenBLKSparseMatrix();

   void allocateMemory();

   Scalar  diag(int dof) const override;
   Scalar  &diag(int dof) override;
   void    unify(FSCommunicator *communicator) override;
   void    addBoeing(int nlines, const int *Kai, const int *Kaj,
                     const double *nz, const int *map, Scalar multiplier) override;
   void    add(const FullSquareMatrix &, const int *dofs) override;
   void    add(const FullSquareMatrixC &, const int *dofs) override;
   void    add(const FullM &knd, int fRow, int fCol);
   void    add(const GenAssembledFullM<Scalar> &kel, const int *dofs) override;
   void    add(int row_dof, int col_dof, Scalar s) override { addone(s, row_dof, col_dof); }
   void    addone(Scalar d, int dofi, int dofj) override;
   void    add(Scalar *_lnz) override;
   Scalar  getone(int row, int col) override;
   void    zeroAll() override;
   void    clean_up() override;
   int     dim() const override { return numUncon; }
   int     neqs() const override { return numUncon; }
   void    printAll();
   void    print() override;
   Scalar* getData() override { return lnz; }

   void    factor() override;

   void    solve(const Scalar *rhs, Scalar *solution) override;
   void    solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution) override;

   void    reSolve(Scalar *rhs) override;
   void    reSolve(GenVector<Scalar> &rhs) override;

   void    reSolve(int nRHS, Scalar **RHS) override;
   void    reSolve(int nRHS, GenVector<Scalar> * RHS) override;

   void    reSolve(GenFullM<Scalar> *) override;

   double getMemoryUsed() const override;
   long size() const override { return (numUncon) ? xlnz[numUncon] : 0; }

   void    getRBMs(double *) override;
   void    getRBMs(Vector *) override;
   void    getRBMs(VectorSet &) override;
   int     numRBM() const override;

   void    addDiscreteMass(int dof, Scalar mass) override;
   void    addImaginary(const FullSquareMatrix &kel, const int *dofs) override;

   void mult(const Scalar *rhs, Scalar *result) const override;
   void getNullSpace(Scalar *ns) override;

 private:
   void init();
   void computeRBMs();
};

template<class Scalar>
class WrapSparseMat : public GenBLKSparseMatrix<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa, *cdsa;
      double trbm;
      SolverCntl& scntl;
      Rbm *rbm;
      CtorData(Connectivity *c, DofSetArray *d, DofSetArray *dc, double t, SolverCntl& _scntl, Rbm *r)
       : cn(c), dsa(d), cdsa(dc), trbm(t), scntl(_scntl), rbm(r) {} 
    };

    WrapSparseMat(CtorData &ctd) : GenBLKSparseMatrix<Scalar>(ctd.cn, ctd.dsa, ctd.cdsa,
        ctd.trbm, ctd.scntl, ctd.rbm) {}
};

typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
typedef GenBLKSparseMatrix<DComplex> BLKSparseMatrixC;

#endif
