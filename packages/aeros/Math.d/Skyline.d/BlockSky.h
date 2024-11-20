#ifndef _BLOCK_SKY_H_
#define _BLOCK_SKY_H_

/*****************************************************************************
 *                   Copyright (C) 1999 CMSoft                               *
 *                                                                           *
 *  These lines of code and declarations contain unpublished proprietary     *
 *  information of CMSoft. They may not be copied or duplicated in whole     *
 *  or part without prior authorization from CMSoft.                         *
 *                                                                           *
 *****************************************************************************/

#ifdef sgi
#include <ulocks.h>
#endif

class Connectivity;
class EqNumberer;
class DofSetArray;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
template <class Scalar> class GenAssembledFullM;
typedef GenAssembledFullM<double> AssembledFullM;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

#include <cstdio>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Utils.d/MyComplex.h>

template<class Scalar>
class GenBlockSky : public GenSolver<Scalar>, public GenSparseMatrix<Scalar> 
{
     int nBlocks;   // number of blocks or supernodes
     int *blWeight; // multiplicity of each supernode
     int *blHeight; // height of each supernodal rectangle
     int *blTop; // begining of each supernodal rectangle
     int *firstDof;
     int *lastCol;
     int *dlp; // diagonal location pointer (non-FORTRAN style)
     int myRCN;
     const int *rowColNum;
     int maxBlockSize; // maximum area of a block
     Scalar *invDiag;
     int *perm;

     double tol;

     int neq;
     int nzem;
     int *sing; // PJSA: location of singularities (used for GtGstar in FetiDPSolver)
   protected:
     Scalar *skyA;
     #ifdef sgi
       void pFactor(int iThread, int numThread, barrier_t *b, 
                    Scalar *ltmp, Scalar *invDiag, Scalar *origDiag);
     #endif
   public:
     GenBlockSky(const Connectivity *nodeToNode, const EqNumberer *eqnums, double tol);
     GenBlockSky(const Connectivity *nodeToNode, const DofSetArray *eqnums, double tol);
     GenBlockSky(const Connectivity *nodeToNode, const DofSetArray *dsa, double tol,
              int *glInternalMap);
     virtual ~GenBlockSky();
     void factor() override;
     void unify(FSCommunicator *communicator) override;
     void parallelFactor() override;
     void solve(const Scalar *rhs, Scalar *solution) override;
     void solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution ) override;
     void reSolve(Scalar *f) override;
     void reSolve(GenVector<Scalar> &f) override;
     void reSolve(int numRHS, Scalar **RHS) override;
     // assembly
     void add(const FullSquareMatrix &, const int *dofs) override;
     void add(const FullSquareMatrixC &, const int *dofs) override;
     void addImaginary(const FullSquareMatrix &kel, const int *dofs) override;
     void add(const AssembledFullM &, const int *dofs);
     void add(const GenAssembledFullM<complex<double> > &, const int *dofs);
     void add(const FullM &, int rowStart, int colStart);
     void add(int row_dof, int col_dof, Scalar s) override;
     void addBoeing(int, const int *, const int *, const double *, const int *, Scalar multiplier) override;
     void addDiscreteMass(int dof, Scalar dmass) override;

     void print(FILE * = stderr);
     int neqs() const override { return neq; }
     long size() const override { return dlp[neq-1]; }
     Scalar diag(int i) const override { return skyA[dlp[i]]; }  // Jing Li's problem
     Scalar &diag(int i) override { return skyA[dlp[i]]; }  // Jing Li's problem
     void zeroAll() override;
     void clean_up() override;
     int dim() const override { return neq;  }
     int numRBM() const  override { return nzem; }
     int* getSingularities() { return sing; }
   private:
     void initialize();
};

typedef GenBlockSky<double> BlockSky;

#endif
