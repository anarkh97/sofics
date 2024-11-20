#ifndef _SUBDOP_H_
#define _SUBDOP_H_
#include <iostream>
#include <Threads.d/Paral.h>
#include <iostream>
#include <Paral.d/Assembler.h>

template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
class DistrInfo;

template<class Scalar>
class GenSubDOp {
       int numSub;                      // Number of subdomains
       GenSparseMatrix<Scalar> **sops;  // List of subdomain matrices
       GenAssembler<Scalar> *assembler; // Assembler object. 
                                        // If NULL, no assembly is performed
    public:
       GenSubDOp(int _numSub, GenSparseMatrix<Scalar> **_sops, GenAssembler<Scalar> *_assembler = 0) {
          numSub = _numSub; sops =_sops; assembler =_assembler;
       }
       GenSubDOp(int _numSub) { numSub = _numSub; sops = new GenSparseMatrix<Scalar> * [numSub]; assembler = 0; }
// RT
       GenSubDOp(int _numSub, GenSparseMatrix<Scalar> ***_sops, int i, GenAssembler<Scalar> *_assembler = 0) {
          numSub = _numSub;
          sops = new GenSparseMatrix<Scalar> *[numSub];
          for(int j=0;j<numSub;j++) {
            if (_sops[j])  {
                *(sops+j) = _sops[j][i];
            }
            else *(sops+j) = 0;
          }
          assembler =_assembler;
       }
// RT end
       ~GenSubDOp();
       void setAssembler(GenAssembler<Scalar> *_a) { assembler = _a; }
       void mult(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Opv); // multiplication routine
       void squareRootMult(GenDistrVector<Scalar> &v); //multiplication by square root
       void inverseSquareRootMult(GenDistrVector<Scalar> &v); //multiplication by inverse square root
       void transposeMult(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Opv); // multiplication routine
       void multAdd(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Opv);
       void multNoAssemble(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &opV); // multiplication routine without assembling
       void assemble(GenDistrVector<Scalar> &v);
       void invertDiag();
       void multInvertDiag(GenDistrVector<Scalar> &);
       GenAssembler<Scalar>* getAssembler() { return assembler; }
       void partialClean() { if(sops) delete [] sops; }
       void zeroAll();
       GenSparseMatrix<Scalar> * operator[](int i) const { return sops[i]; }
       GenSparseMatrix<Scalar> *& operator[](int i) { return sops[i]; }
       int getNumSub() { return numSub; }
       int neqs() const { return 0; } // TODO
};

template<class Scalar>
class SubDOpMult : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
     GenDistrVector<Scalar> *s, *d;
  public:
     SubDOpMult(GenSparseMatrix<Scalar> **_sm, GenDistrVector<Scalar> *_s, GenDistrVector<Scalar> *_d) {
       sm = _sm; s = _s; d = _d;
     }
     void runFor(int);
};

template<class Scalar>
class SubDOpTransposeMult : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
     GenDistrVector<Scalar> *s, *d;
  public:
     SubDOpTransposeMult(GenSparseMatrix<Scalar> **_sm, GenDistrVector<Scalar> *_s, GenDistrVector<Scalar> *_d) {
       sm = _sm; s = _s; d = _d;
     }
     void runFor(int);
};

template<class Scalar>
class SubDOpSquareRootMult : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
     GenDistrVector<Scalar> *s;
  public:
     SubDOpSquareRootMult(GenSparseMatrix<Scalar> **_sm, GenDistrVector<Scalar> *_s) {
       sm = _sm; s = _s;
     }
     void runFor(int);
};

template<class Scalar>
class SubDOpInverseSquareRootMult : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
     GenDistrVector<Scalar> *s;
  public:
     SubDOpInverseSquareRootMult(GenSparseMatrix<Scalar> **_sm, GenDistrVector<Scalar> *_s) {
       sm = _sm; s = _s;
     }
     void runFor(int);
};

template<class Scalar>
class SubDOpMultAdd : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
     GenDistrVector<Scalar> *s, *d;
  public:
     SubDOpMultAdd(GenSparseMatrix<Scalar> **_sm, GenDistrVector<Scalar> *_s, GenDistrVector<Scalar> *_d) {
       sm = _sm; s = _s; d = _d;
     }
     void runFor(int);
};

template<class Scalar>
class SubDOpMultInvertDiag : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
     GenDistrVector<Scalar> *z;
  public:
     SubDOpMultInvertDiag(GenSparseMatrix<Scalar> **_sm, GenDistrVector<Scalar> *_z) {
       sm = _sm; z = _z;
     }
     void runFor(int);
};

template<class Scalar>
class SubDOpInvertDiag : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
  public:
     SubDOpInvertDiag(GenSparseMatrix<Scalar> **_sm) {
       sm = _sm;
     }
     void runFor(int);
};

template<class Scalar>
class SubDOpZeroAll : public TaskDescr {
     GenSparseMatrix<Scalar> **sm;
  public:
     SubDOpZeroAll(GenSparseMatrix<Scalar> **_sm) {
       sm = _sm;
     }
     void runFor(int);
};


typedef GenSubDOp<double> SubDOp;

#ifdef _TEMPLATE_FIX_
  #include <Paral.d/SubDOp.C>
#endif

#endif
