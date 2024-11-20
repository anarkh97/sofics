#ifndef CUC_SPARSE_H_
#define CUC_SPARSE_H_

#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h>

class ControlLawInfo;

class LMPCons;

/**\brief Constrained to Unconstrained Sparse Matrix.
 * \details This is a rectangular matrix that needs to be stored in sparse format) */
template<class Scalar>
class GenCuCSparse : public SparseData, public GenSparseMatrix<Scalar> {
    std::vector<Scalar> Kuc;
  public:
    GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const int *bc);
    GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const DofSetArray *c_dsa);
    GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const int *glbmap, const int *glimap);

	GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const std::vector<int> &glbmap, const int *glimap)
		: GenCuCSparse(con, dsa, glbmap.data(), glimap)
	{
	}
	GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const std::vector<int> &glbmap, const std::vector<int> &glimap)
		: GenCuCSparse(con, dsa, glbmap.data(), glimap.data())
	{
	}
	GenCuCSparse(const Connectivity *con, const DofSetArray *dsa, const int *glbmap, const std::vector<int> &glimap)
		: GenCuCSparse(con, dsa, glbmap, glimap.data())
	{
	}

    GenCuCSparse(LMPCons **mpc, int numMPC, const DofSetArray *c_dsa);
    GenCuCSparse(int numInterface, const int *glbmap, int numRBM, Scalar *rbms, int lda = -1);
    /**
     *
     * @param num Number of constraint DOFs.
     * @param neq Number of unconstrained DOFs.
     * @param count Count of terms per constrained DOF.
     * @param list Unconstrained indices of each term.
     * @param val Values of the terms.
     */
	GenCuCSparse(int num, int neq, gsl::span<int> count, gsl::span<int> list, std::vector<Scalar> val);

	~GenCuCSparse() override = default;

    double getMemoryUsed() const override;

    void add(const FullSquareMatrix &mel, const int *dofs) override;
    void addImaginary(const FullSquareMatrix &kel, const int *dofs) override;
    void add(const FullSquareMatrixC &mel, const int *dofs) override;
    void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
    void addDiscreteMass(int dof, Scalar dimass) override;
    void add(int row_dof, int col_dof, Scalar s) override;

    void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const override;
    void mult(const Scalar *rhs, Scalar *result) const override;

    void multIdentity(Scalar **result) const override;
    /// \copydoc
    void multIdentity(Scalar *result) const override;
    void multIdentity(Scalar **v, int start, int stop) const override;

    void multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const override;
    void multSubtract(const Scalar *rhs, Scalar *result) const override;

    void multSubAdd(const Scalar *rhs, Scalar *result) const;

    Scalar diag(int) const override { throw "GenCuCSparse::diag - 1 - should never be called"; }
    Scalar &diag(int) override { throw "GenCuCSparse::diag - 2 - should never be called"; }
    void zeroAll() override;
    int  dim() const override { return numConstrained; }
    int  neqs() const override { return numConstrained; }
    int  numRow() const override { return neq-numConstrained; }
    int  numCol() const override { return numConstrained; }
    void print() override;
    long size() const;

    void negate();
    void clean_up() override;

    void multSub(const Scalar *rhs, Scalar *result) const override;
    void multSub(int numRHS, const Scalar **rhs, Scalar **result) const override;
    void transposeMultAdd(const Scalar *rhs,Scalar *result) const override;
    void multAdd(const Scalar *rhs,Scalar *result) const override;
    void transposeMultSubtract(const Scalar *rhs,Scalar *result) const override;
    void transposeMultSubtractClaw(const Scalar *rhs,Scalar *result, int ncl,
                                   int *clawfDofs) const override;
    void transposeMult(const Scalar *rhs,Scalar *result) const override;

    void doWeighting(Scalar *weight);  // PJSA 4-3-03
    void doWeighting(int *weight);

    void multAddNew(const Scalar *rhs, Scalar *result) const;
    void transposeMultNew(const Scalar *rhs, Scalar *result) const;
    void transposeMultAddNew(const Scalar *rhs, Scalar *result) const;
    void transposeMultSubNew(const Scalar *rhs, Scalar *result) const;

    void  mult(  const Scalar *rhs, Scalar *result, Scalar alpha, Scalar beta) const; //HB 12/15/04
    void  trMult(const Scalar *rhs, Scalar *result, Scalar alpha, Scalar beta) const; //HB 12/15/04
};

typedef GenCuCSparse<double> CuCSparse;
typedef GenCuCSparse<DComplex> CuCComplexSparse;

#endif
