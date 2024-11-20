#ifndef _SKYMATRIX_H_
#define _SKYMATRIX_H_

#include <cstdio>

#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/MyComplex.h>
#if defined(sgi) && ! defined(_OPENMP)
#include <ulocks.h>
#endif

template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;
class CoordSet;

class SkyData {
 protected:
   int *lacol;          //<! \brief last active column (lacol)
   int *pivot;          //<! \brief stores pivot information for the solver
   int *seqid;          //<! \brief stores pivot information for geometric sky solver
   int *dlp;            //<! \brief diagonal location pointer

                        //<! \brief in some instances, rowColNum is allocated
                        //<! \brief from a SkyData constructor and other times,
                        //<! \brief it is allocated in the ConstrainedDSA class.
   int myRCN;           //<! \brief 0=not allocated, 1=allocated memory
   const int *rowColNum;//<! \brief Mapping vector from full DOF indices to unconstrained.

   int neq;             //<! \brief # of equations
   int numUncon;        //<! \brief # unconstrained dofs
   double TOLERANCE;    //<! \brief Tolerance value for factoring

   int nzem;            //<! \brief # zero energy modes
   Rbm *rbm;            //<! \brief stores rigid body mode information
   int isTRBM;          //<! \brief whether we are TRBM or GRBM

   int myRbm;

 public:
   SkyData();
   SkyData(const Connectivity *cn, const EqNumberer *dsa, double trbm, const int *bc);
   SkyData(const Connectivity *cn, const EqNumberer *dsa, double trbm);
   SkyData(const Connectivity *cn, const DofSetArray *c_dsa, double trbm, Rbm *rigid=0);
   SkyData(const EqNumberer *_dsa, const Connectivity *cn, double trbm, const int *rCN);
   SkyData(int n, double tolerance = 1.0E-4);
   virtual ~SkyData();
 private:
   void initialize();
};

// Sky Matrix:

template<class Scalar>
class GenSkyMatrix : public SkyData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
protected:
   Scalar *skyA = nullptr;  // Sky array
   int isScaled = false;    // whether to scale the matrix or not.
   Scalar *scale = nullptr; // vector to store the matrix scaling
   bool wasScaled = false; // HB: indicate if the scaling has already been applied before
                            //     getting in Factor/parallelFactor method
   // Timing data members
   mutable double solveTime = 0.0;
   double constructTime = 0.0;
 public:
   // Constructors
   GenSkyMatrix() = default;
   GenSkyMatrix(const Connectivity *cn, const EqNumberer *dsa, double trbm, const int *bc);
   GenSkyMatrix(const Connectivity *cn, const EqNumberer *dsa, double trbm, int isScaled = 0);
   GenSkyMatrix(const Connectivity *cn, const DofSetArray *, double trbm, Rbm *rigid = 0);
   GenSkyMatrix(const Connectivity *cn, const EqNumberer *dsa, double trbm, const int *rCN, int dummy);
   GenSkyMatrix(GenFullM<Scalar> *mat, double tolerance = 1.0E-4);
   GenSkyMatrix(const Connectivity *, const EqNumberer *, const ConstrainedDSA *, double trbm);
   GenSkyMatrix(int, double);
   // Destructor
   virtual ~GenSkyMatrix();

   // returns memory used in megabytes
   double getMemoryUsed() const override { return 8*dlp[numUncon-1]/(1024.0*1024.0); }
   long size() const override { return (numUncon) ? dlp[numUncon - 1] : 0; }
   Scalar* getData() override { return skyA; }

   void mult(const GenVector<Scalar> &, GenVector<Scalar> &) const override;
   void mult(const Scalar *rhs, Scalar *result) const override;
   Scalar diag(int) const override;                           // returns diagonal entry
   Scalar &diag(int) override;

   // assembly
   void add(const FullSquareMatrix &, const int *dofs) override;
   void add(const FullSquareMatrixC &, const int *dofs) override;
   void addImaginary(const FullSquareMatrix &, const int *dofs) override;
   void add(const GenAssembledFullM<Scalar> &kel, const int *dofs) override;
   void add(const GenFullM<Scalar> &, int rowStart, int colStart) override;
   void addone(Scalar d, int dofi, int dofj) override;
   void addBoeing(int, const int *, const int *, const double *, const int *, Scalar multiplier) override;
   void add(Scalar *_skyA) override;
   void addPoint(Scalar, int, int);

   void Factor();
   void Factor(Rbm *rigid);
   void factor() override;
   // Parallel factorization. It's used for Coarse problems only
   virtual void parallelFactor() override;

   void solve(const Scalar *rhs, Scalar *solution) override;
   void solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution) override;

   void reSolve(GenVector<Scalar> &rhs) override;
   void reSolve(Scalar *rhs) override;
   void reSolve(int numRHS, Scalar **RHS) override;
   void reSolve(int numRHS, GenVector<Scalar> *RHS) override;
   void reSolve(int numRHS, Scalar *RHS) override;
   void reSolve(GenFullM<Scalar> *mat) override;

   void forward(GenVector<Scalar> &rhs) override;
   void backward(GenVector<Scalar> &rhs) override;

   void unify(FSCommunicator *communicator) override;

   GenVector<Scalar>* getNullSpace();   // retrieve ZERO ENERGY MODES
   void getNullSpace(Scalar *rbm) override;      // retrieve ZERO ENERGY MODES

   int  numRBM() const override { return nzem; } // retrieve the number of rigid body modes
   void getRBMs(double *) override;         // retrieve the rigid body modes
   void getRBMs(Vector* vs) override;       // retrieve the rigid body modes
   void getRBMs(VectorSet& vs) override;    // retrieve the rigid body modes

   void addDiscreteMass(int cdof, Scalar diMass) override;
   void add(int row_dof, int col_dof, Scalar s) override { addone(s, row_dof, col_dof); }
   void zeroAll() override;
   int  dim() const override { return numUncon; }
   int  neqs()const override { return numUncon; }
   void allMult(Scalar x); // multiply the whole skyline matrix by x

   void printMemory();
   void printConstructTime();
   void print1(int dof);
   void print(FILE * =stderr);
   void printMatlab(int subNumber);
   void printMatlab(char *fileName);
   void printDiagonals();

   void clean_up() override;

   Scalar getone(int row, int col) override;

   void applyScaling(Scalar *vector) const;
   void symmetricScaling();
   bool IsScaled() { return (isScaled)? true : false; }
   void setIsScaled(int _isScaled) { isScaled = _isScaled; }

   //HB
   double rmsBandwidth() const;

#if defined(sgi) && ! defined(_OPENMP)
   void pfact(int, int, barrier_t *, Scalar *);
#else
   void pfact(int, int, Scalar *);
#endif
};

template<class Scalar>
class WrapSkyMat : public GenSkyMatrix<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa;
      double trbm;
      Rbm *rbm;
      CtorData(Connectivity *c, DofSetArray *d, double t, Rbm *r) {
        cn = c;
        dsa = d;
        trbm = t;
        rbm = r;
        if(r && r->numComponents() > 1) {
          std::cerr << " *** WARNING: multi-component GRBM with skyline solver and direct constraint method is not supported\n";
        }
      }
    };

    WrapSkyMat(CtorData &ctd) : GenSkyMatrix<Scalar>(ctd.cn, ctd.dsa, ctd.trbm, ctd.rbm) {}
};


typedef GenSkyMatrix<double> SkyMatrix;
typedef GenSkyMatrix<DComplex> SkyMatrixC;

#endif
