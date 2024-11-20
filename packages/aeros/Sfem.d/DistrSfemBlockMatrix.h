#ifndef DISTR_SFEM_BLOCK_MATRIX_H_
#define DISTR_SFEM_BLOCK_MATRIX_H_

#include <Utils.d/MyComplex.h>
#include <Math.d/matrix.h>
#include <Feti.d/Feti.h>
#include <Paral.d/SubDOp.h>
#include <iostream>

struct DistrBlockInfo {
  DistrInfo *blockinfo;
  int nblocks;
  int* nnzblkindex;
};

template<class Scalar> 
class DistrBlockVector {
    GenDistrVector<Scalar> **v;
//    const DistrBlockInfo &inf;
    DistrBlockInfo inf;
//    Scalar *v_scalar; // Scalar version of v
    double* blocknorms;
  public:
    DistrBlockVector(DistrBlockInfo _inf) : inf(_inf){ 
      v = new GenDistrVector<Scalar> * [inf.nblocks]; 
      for(int i=0; i<inf.nblocks; ++i) v[i] = new GenDistrVector<Scalar>(inf.blockinfo[i]);
      blocknorms = new double[inf.nblocks];
    }
    DistrBlockVector(DistrBlockVector &v2) : inf(v2.info()) {
      v = new GenDistrVector<Scalar> * [inf.nblocks];
      for(int i=0; i<inf.nblocks; ++i) v[i] = new GenDistrVector<Scalar>(*v2.getBlock(i));
      blocknorms = new double[inf.nblocks];
    }

	GenDistrVector<Scalar>** getv() {return v;}
	GenDistrVector<Scalar>** const getv() const {return v;}

    GenDistrVector<Scalar>&  getBlock(int iblock) { return *(v[iblock]); }
    const GenDistrVector<Scalar>&  getBlock(int iblock) const { return *(v[iblock]); }
//    GenDistrVector<Scalar>&  getBlock(int iblock, int junk) { return *(v[iblock]); } // YYY remove later
    DistrBlockInfo info() { return inf; }

    void scaleBlock(int iblock, Scalar s) { *(v[iblock]) *= (1/s); } // YYY change to /= operator
    void setBlock(int iblock, GenDistrVector<Scalar> &vi) { *(v[iblock]) = vi; }

    ~DistrBlockVector(); 

    int size() const;
    void zero();
    Scalar operator[](int i) const { int n = inf.blockinfo[0].len; return (*(v[i/n]))[i%n]; } 
    Scalar & operator[](int i) { int n = inf.blockinfo[0].len; return (*(v[i/n]))[i%n]; }
    Scalar operator * (DistrBlockVector &);
    double norm(); 
    double sqNorm(); 
    DistrBlockVector& operator=(const DistrBlockVector<Scalar> &);
    DistrBlockVector& operator=(Scalar c);
    DistrBlockVector& operator*=(Scalar c);
    DistrBlockVector& operator+=(DistrBlockVector<Scalar> &);
    DistrBlockVector& operator-=(DistrBlockVector<Scalar> &);
    DistrBlockVector& linAdd(DistrBlockVector<Scalar> &);
    DistrBlockVector& linAdd(Scalar, DistrBlockVector<Scalar> &);
    DistrBlockVector& linAdd(Scalar, DistrBlockVector<Scalar> &, Scalar, DistrBlockVector<Scalar> &);
    DistrBlockVector& linC(const DistrBlockVector<Scalar> &, Scalar);
    DistrBlockVector& linC(const DistrBlockVector<Scalar> &, Scalar, const DistrBlockVector<Scalar> &);
    DistrBlockVector& linC(Scalar, const DistrBlockVector<Scalar> &, Scalar, const DistrBlockVector<Scalar> &);
    DistrBlockVector& swap(DistrBlockVector<Scalar> &);
    Scalar sum();
   
    void updateBlock(int ii, Scalar c, DistrBlockVector<Scalar> &);
    void copyBlock(DistrBlockVector<Scalar> &, int ii);
    void addBlockSqr(int ii, Scalar c, DistrBlockVector<Scalar> &);
    void computeSqrt();
    void computeRealz(int ii, Scalar c, DistrBlockVector<Scalar> &);

    Scalar *data() const  { std::cerr << "Warning :: DistrBlockVector::data() not implemented"; return 0;  } 
    void print();
    void printNonZeroTerms();
    void printAll();
    void setn(int _n) {};
    void computeBlockNorms(); 
    double* getBlockNorms() { return blocknorms;}
    void printBlockNorms(); 
    void setNnzBlocks(int* bl) {for (int i=0; i<inf.nblocks; ++i) inf.nnzblkindex[i]=bl[i];} // YYY DG
    int isnnz(int i) {return inf.nnzblkindex[i];}
    void printBlockDetails();  
};

template<class Scalar>
class DistrSfemBlockMatrix {
    int n, L, P, ndim;   // comes from class Sfem
    GenDistrVector<Scalar>*** kiuj;
    SubDOp** allK;
    double ***inpsi;
    Connectivity *blockToBlock;
    Scalar *diags;
    int *firstdof;
    GenFetiSolver<Scalar>* meansolver;
    Scalar* scalarfactors;
    DistrBlockInfo &inf; 
  public:
    DistrSfemBlockMatrix(int L, int P, int _ndim, int output_order, DistrBlockInfo & info); 
    ~DistrSfemBlockMatrix() {}; // YYY DG delete Ki, allK  all of the above // also call delete sfbm
    void setKi(SubDOp* Ki, int i);
    void matvec_sfem_block(DistrBlockVector<Scalar> &u);
    void mult(DistrBlockVector<Scalar> &u, DistrBlockVector<Scalar> &ku);
    DistrBlockInfo &dim(); 
    int neqs() const {return n*P;} //  in some cases like dim()
    Scalar diag(int i) {Scalar r =0.0; std::cerr <<"DistrSfemBlockMatrix::diag not implemented" << std::endl; return r;}
    int* getFirstDof();
    int numNodes() const;
    GenFullM<Scalar>* getDiagMatrix(int i);
    void setMeanSolver(GenFetiSolver<Scalar> *prc);
    GenFetiSolver<Scalar>* getMeanSolver();
    Scalar* getBlockScalarMultipliers();
    int getBlockSize() {std::cerr << "DistrSfemBlockMatrix::getBlockSize() not implemented" << std::endl; return 0;}
};

#ifdef _TEMPLATE_FIX_
#include <Sfem.d/DistrSfemBlockMatrix.C>
#endif

#endif
