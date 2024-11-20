#ifndef _DSC_SOLVER_H_
#define _DSC_SOLVER_H_

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class Connectivity;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#ifdef DISTRIBUTED
#include <mpi.h>
#endif

class DSCsolver :
      public SparseData, public SparseMatrix, public Solver {

   double *unonz;   // storage of A

   int  *adj; 
   int *xadj;

   int scheme_number;
   int numNodes;

   double tol;

#ifdef DISTRIBUTED
   MPI_Comm dscComm;
#endif
   int color;
   int maxNum;

 public:

   DSCsolver(Connectivity *cn, EqNumberer *eqNums, int sch_number);
   virtual ~DSCsolver();

   void    add(const FullSquareMatrix &knd, const int *dofs) override {};
   void    add(const FullM &knd, int fRow, int fCol) override;
   void    zeroAll() override;

   void    factor() override;

   void    reSolve(double *rhs) override;
   void    reSolve(Vector &rhs) override;

   void    reSolve(int nRHS, double **RHS) override;
   void    reSolve(int nRHS, Vector * RHS) override;
   void    unify(FSCommunicator *communicator) override;
   void    print() override;

   int     dim()  const override { return numUncon;  }
   int     neqs() const override { return numUncon;  }

   // Functions that need to be written
   int numRBM() const override { return -1; }

   // Functions not needed
   double    diag(int) const override { return 1.0; }
   double &diag(int) override { throw "Crazy programmers\n"; }
   long size()    { return 0; }
};

#endif
