#ifndef _FETI_MPC_SPARSE_H_
#define _FETI_MPC_SPARSE_H_

#include <cstdio>
#include <Math.d/SparseMatrix.h>
#include <Driver.d/Mpc.h>
#include <Utils.d/MyComplex.h>

// Container class that holds the MPC equation information

template<class Scalar>
class GenMpcSparse : public GenSparseMatrix<Scalar> 
{
	const std::vector<std::unique_ptr<SubLMPCons<Scalar> > > &mpc;
    int NumRow;
    int NumCol;
    const DofSetArray *dsa;
    const DofSetArray *DSA;
    int *wetInterfaceMap;
    int mpcOffset;
   
  public:
    GenMpcSparse(int numMPC, const std::vector<std::unique_ptr<SubLMPCons<Scalar> > > &_mpc, const DofSetArray *_dsa);
    GenMpcSparse(int numMPC, const std::vector<std::unique_ptr<SubLMPCons<Scalar> > > &_mpc, const DofSetArray *_dsa, const DofSetArray *_DSA,
                 int *_wetInterfaceMap, int mpcOffset);
    virtual ~GenMpcSparse();

    double getMemoryUsed() const override;

    void add(const double *const *kel, const int *dofs, int kndof);
    void add(const FullSquareMatrix &mel, const int *dofs);
    void add(const FullM&, int, int);
    void add(const GenAssembledFullM<Scalar> &, const int*) override;

    void addBoeing(int, const int *, const int *, const double *, const int *, Scalar multiplier);

    void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const override ;
    void mult(const Scalar *rhs, Scalar *result) const override;
    
    void multIdentity(int dof, Scalar *result) const;
    void multIdentity(Scalar **result) const override;
    void multIdentity(Scalar *result) const override;
    void multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const override;
    void multSubtract(const Scalar *rhs, Scalar *result) const override;

    void zeroAll() override;
    int  dim() const override { return 0; }
    int  neqs() const override { return 0; }
    int  numRow() const override { return NumRow; }
    int  numCol() const override { return NumCol; }

    Scalar diag(int) const override { throw "GenMpcSparse::diag - 1 - should never be called"; }
    Scalar &diag(int) override { throw "GenMpcSparse::diag - 2 - should never be called"; }
    long size() const { return 0; }
    void print(FILE *file=stderr, const char* msg="A");
    void negate();

    void multSub(const Scalar *rhs, Scalar *result) const override;
    void multSub(int numRHS, Scalar **rhs, Scalar **result) const;
    void multSubWI(const Scalar *rhs, Scalar *result) const;
    void transposeMultAdd(const Scalar *rhs, Scalar *result) const override;
    void multAdd(const Scalar *rhs, Scalar *result) const override;
    void transposeMultSubtract(const Scalar *rhs, Scalar *result) const override;
    void transposeMultSubtractWI(const Scalar *rhs, Scalar *wi_result) const;
    void transposeMult(const Scalar *rhs, Scalar *result) const override;

  private:
    GenMpcSparse<Scalar>& operator = (const GenMpcSparse<Scalar> &);
    GenMpcSparse<Scalar>(const GenMpcSparse<Scalar> &);
};

typedef GenMpcSparse<double> MpcSparse;
typedef GenMpcSparse<DComplex> MpcSparseC;

#ifdef _TEMPLATE_FIX_
#include <Math.d/MpcSparse.C>
#endif


#endif
