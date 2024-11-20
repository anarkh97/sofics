#ifndef SPOOLES_H_
#define SPOOLES_H_

// Spooles include files
#ifdef USE_SPOOLES
extern "C" {
  #include <misc.h>
  #include <FrontMtx.h>
  #include <SymbFac.h>
}
#endif

#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h>


class ConstrainedDSA;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class Connectivity;
class SolverCntl;

class FSCommunicator;

template<class Scalar>
class GenSpoolesSolver : public GenSolver<Scalar>, public GenSparseMatrix<Scalar>,
                         public SparseData 
{
   const SolverCntl& scntl;
   int neq;               // number of equations
   int *constrNum;        // constrained equation numbers
   Scalar *unonz;
   double cpus[22];
   int stats[7];
   int pivotingflag;
   int numThreads;
   int nNonZero;
   int isScaled;        // whether to scale the matrix or not
   Scalar *scale;       // vector to store the matrix scaling
   int msglvl;
   FILE *msgfile;
   long _size;

#ifdef USE_SPOOLES
   InpMtx   *inpMtx;	 // input Matrix
   FrontMtx *frontMtx;    // front Matrix
   IV       *newToOldIV;  // new to old permutation table
   IV       *oldToNewIV;  // old to new permutation table
   SubMtxManager *mtxManager;
   IV *ownersIV;
   DenseMtx *mtxB, *mtxX;
   IVL *symbfacIVL;
   ETree *frontETree;
   DV *cumopsDV;
   Graph *graph;
#endif

 public:
   GenSpoolesSolver(const Connectivity *nToN, const EqNumberer *dsa, const SolverCntl& _scntl, int *map=0);
   GenSpoolesSolver(const Connectivity *nToN, const DofSetArray *_dsa, const ConstrainedDSA *c_dsa, const SolverCntl& _scntl);

   virtual void clean_up() override {
     cleanUp();
     if(unonz) { delete [] unonz; unonz = 0; }
     if(scale) { delete [] scale; scale = 0; }
   }

   virtual ~GenSpoolesSolver();

   void add(const FullSquareMatrix &, const int *dofs) override;
   void addImaginary(const FullSquareMatrix &, const int *dofs) override;
   void add(const FullSquareMatrixC&, const int *dofs) override; // RT addded to support PML, DGM
   void add(const GenFullM<Scalar> &, const int *dofs);
   void add(const GenFullM<Scalar> &, int, int) override;
   void add(const GenAssembledFullM<Scalar> &, const int *) override;
   void addDiscreteMass(int dof, Scalar) override;
   void add(int dofi, int dofj, Scalar d) override; //HB: add upper part only
   void addone(Scalar d, int dofi, int dofj) override { add(dofi, dofj, d); }
   //void addBoeing(int nlines, const int *Kai, const int *Kaj,
   //               const double *nz, int *map, Scalar multiplier);

   void unify(FSCommunicator *) override;
   void parallelFactor() override;
   void factor() override;
   void allFactor(bool fctIsParal);

   void reSolve(Scalar *rhs) override;
   void solve(const Scalar *rhs, Scalar *solution) override;
   
   void print() override;
   int dim() const override { return neq;      }
   int neqs() const override { return neq; }
   double getMemoryUsed() const override;
   long size() const override;
   Scalar  diag(int dof) const override;
   Scalar &diag(int dof) override;

   void    zeroAll() override;
   void    cleanUp();

   int numRBM() const override { return 0; } // note: spooles should not be used for singular matrices.
 
 private:
   void init();
   void applyScaling(Scalar *vector);
   void symmetricScaling();
};

template<class Scalar>
class WrapSpooles : public GenSpoolesSolver<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa;
      ConstrainedDSA *cdsa;
      SolverCntl& scntl;
      CtorData(Connectivity *c, DofSetArray *d, ConstrainedDSA *dc, SolverCntl& _scntl)
       : cn(c), dsa(d), cdsa(dc), scntl(_scntl) {}
    };

    WrapSpooles(CtorData &ctd) : GenSpoolesSolver<Scalar>(ctd.cn, ctd.dsa, ctd.cdsa, ctd.scntl) {}
};

typedef GenSpoolesSolver<double> SpoolesSolver;
typedef GenSpoolesSolver<DComplex> ComplexSpoolesSolver;

#endif
