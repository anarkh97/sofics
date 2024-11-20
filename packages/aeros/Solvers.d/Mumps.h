#ifndef MUMPS_H_
#define MUMPS_H_

// Mumps include files
#ifdef USE_MUMPS
// include MUMPS librairies here or in header file
#include "dmumps_c.h" // double precision mumps header
#include "zmumps_c.h" // complex double precision
#endif

#include <gsl/span>
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/MultiDomainSolver.h>
#include <Utils.d/MyComplex.h>
#include <Comm.d/Communicator.h>

class EqNumberer;
class ConstrainedDSA;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class Connectivity;
class SolverCntl;
class FSCommunicator;

#ifdef USE_MUMPS
template<class Scalar>
class MumpsId {
};

template <>
class MumpsId<double> {
 public:
   DMUMPS_STRUC_C id;
   typedef double MumpsType;
};

template <>
class MumpsId<complex<double> > {
 public:
   ZMUMPS_STRUC_C id;
   typedef mumps_double_complex MumpsType;
};
#endif

template<class Scalar>
class GenMumpsSolver : public GenSolver<Scalar>, public GenSparseMatrix<Scalar>, public SparseData, public MultiDomainSolver<Scalar> 
{
   const SolverCntl& scntl;
   int 	neq;        // number of equations = id.n for Mumps
   Scalar *unonz;   // matrix of elements = id.a for mumps
   int nNonZero;    // number of non zero entries = id.nz for Mumps	
   int nrbm;        // number of zero pivots detected

#ifdef USE_MUMPS
   MumpsId<Scalar> mumpsId;
#endif

   FSCommunicator *mpicomm;
   bool host;
   Timings times;

   bool mumpsCPU; // JAT 052214
   Communicator *groupcomm; // JAT 072616

public:
	GenMumpsSolver(const Connectivity *nToN, const EqNumberer *dsa, const SolverCntl& _scntl, int *map=0, FSCommunicator *_mpicomm = 0);
	GenMumpsSolver(const Connectivity *nToN, const DofSetArray *dsa, const ConstrainedDSA *c_dsa, const SolverCntl& _scntl,
	               FSCommunicator *_mpicomm = 0);
	GenMumpsSolver(const Connectivity *nToN, const DofSetArray *dsa, const ConstrainedDSA *c_dsa, int nsub,
	               gsl::span<GenSubDomain<Scalar> *> sd,
	               const SolverCntl& _scntl, FSCommunicator *_mpicomm = 0);

   virtual ~GenMumpsSolver();

   void add(const FullSquareMatrix &, const int *dofs);
   void addImaginary(const FullSquareMatrix &, const int *dofs);
   void add(const FullSquareMatrixC&, const int *dofs);
   void add(const GenFullM<Scalar> &, int, int);
   void add(const GenAssembledFullM<Scalar> &, const int *);
   void addDiscreteMass(int dof, Scalar);
   void addone(Scalar d, int dofi, int dofj) { GenSparseMatrix<Scalar>::add(dofi, dofj, d); }

   void unify(FSCommunicator *);
   void factor();

   void solve(const Scalar *rhs, Scalar *solution);
   void reSolve(Scalar *rhs);
   void reSolve(int nRHS, Scalar *rhs);
   void reSolve(int nRHS, Scalar **rhs);
   void reSolve(int nRHS, GenVector<Scalar> *rhs);
   void getNullSpace(Scalar *rbm);

   void mult(const Scalar *rhs, Scalar *result);

   int dim() const { return neq; }
   int neqs() const { return neq; }
   long size() const;
   int numRBM() const { return nrbm; }
   void getRBMs(Vector *rbms);
   void getRBMs(VectorSet& rbms);
   int* getPivnull_list();

   void print();
   Scalar  diag(int dof) const;
   Scalar &diag(int dof);

   void zeroAll();

   // for parallel solver
   double getSolutionTime() const override { return ((GenSolver<Scalar> *)this)->getSolutionTime(); }
   Timings& getTimers() { return times; }
   void refactor() { factor(); }
   
 private:
   void init();
   void printStatistics();
#ifdef USE_MUMPS
   void copyToMumpsLHS(mumps_double_complex *&m, DComplex *&d, int len);
   void copyToMumpsLHS(double *&m, double *&d, int len);
   void copyToMumpsRHS(mumps_double_complex *&m, const DComplex *d, int len);
   void copyToMumpsRHS(double *&m, const double *d, int len);
   void copyFromMumpsRHS(DComplex *d, mumps_double_complex *m, int len);
   void copyFromMumpsRHS(double *d, double *m, int len);
   void copyToMumpsRHS(mumps_double_complex *&m, DComplex **d, int len, int nRHS);
   void copyToMumpsRHS(double *&m, double **d, int len, int nRHS);
   void copyFromMumpsRHS(DComplex **d, mumps_double_complex *m, int len, int nRHS);
   void copyFromMumpsRHS(double **d, double *m, int len, int nRHS);
   void copyToMumpsRHS(mumps_double_complex *&m, const GenVector<DComplex> *d, int len, int nRHS);
   void copyToMumpsRHS(double *&m, const GenVector<double> *d, int len, int nRHS);
   void copyFromMumpsRHS(GenVector<DComplex> *d, mumps_double_complex *m, int len, int nRHS);
   void copyFromMumpsRHS(GenVector<double> *d, double *m, int len, int nRHS);
#endif
};

template<class Scalar>
class WrapMumps : public GenMumpsSolver<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa;
      ConstrainedDSA *cdsa;
      SolverCntl& scntl;
      FSCommunicator *com;
      CtorData(Connectivity *c, DofSetArray *d, ConstrainedDSA *dc, SolverCntl& _scntl, FSCommunicator *_com)
        : cn(c), dsa(d), cdsa(dc), scntl(_scntl), com(_com) {}
    };

    WrapMumps(CtorData &ctd) : GenMumpsSolver<Scalar>(ctd.cn, ctd.dsa, ctd.cdsa, ctd.scntl, ctd.com) {}
};


typedef GenMumpsSolver<double> MumpsSolver;
typedef GenMumpsSolver<DComplex> ComplexMumpsSolver;

#endif
