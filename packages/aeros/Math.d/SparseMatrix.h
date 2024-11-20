#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <gsl/span>
#include <Utils.d/dofset.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/MyComplex.h>

#include <string>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
class Connectivity;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenAssembledFullM;
typedef GenAssembledFullM<DComplex> AssembledFullMC;
template <class Scalar> class GenFullM;
typedef GenFullM<DComplex> FullMC;
template <class Scalar> class GenSolver;
template <class Scalar> class DistrBlockVector;


template<class Scalar>
class GenSMV {
 GenSparseMatrix<Scalar> &M;
 GenVector<Scalar>       &v;
 public:
   GenSMV(GenSparseMatrix<Scalar> &_M, GenVector<Scalar> &_v) : M(_M), v(_v) {}
};

template<class Scalar>
class GenSparseMatrix {
		// TODO Get this out of here! Single use object that should be built
		// by the class using it (ScalarBlockDiagPrec)
        GenSolver<Scalar>* meansolver;
        int *firstdof;
        Scalar* scalarfactors;
public:
        virtual void zeroAll() = 0;
        virtual int dim() const = 0;
        virtual double norm() const;
        virtual void clean_up();
        virtual double getMemoryUsed() const;
        virtual int  numRow() const;
        virtual int  numCol() const;
        virtual Scalar diag(int dof) const = 0;
        virtual Scalar &diag(int dof) = 0;
        virtual void invertDiag();
        virtual void printSparse(const std::string& filename);

	/** \brief Assemble an element's contribution into the matrix.
	 *
	 * @param kel Elemental matrix.
	 * @param dofs DOFs of the elemental matrix.
	 */
	// TODO Make this the real add virtual method
	void add(const FullSquareMatrix &kel, gsl::span<const int> dofs) {
		add(kel, dofs.data());
	}
	void add(const FullSquareMatrix &kel, const gsl::span<int> &dofs) {
		add(kel, dofs.data());
	}
	void add(const FullSquareMatrixC &kel, gsl::span<const int> dofs) {
		add(kel, dofs.data());
	}
	void add(const GenAssembledFullM<Scalar> &kel, gsl::span<const int> dofs) {
		add(kel, dofs.data());
	}
    void addImaginary(const FullSquareMatrix &kel, gsl::span<const int> dofs) {
        add(kel, dofs.data());
    }
        virtual void add(const FullSquareMatrix &, const int *dofs) = 0;
        virtual void addImaginary(const FullSquareMatrix &, const int *dofs);
        virtual void add(const FullSquareMatrixC &, const int *dofs);
        virtual void add(const GenFullM<Scalar> &knd, int fRow, int fCol);
        /// \brief Assemble skipping the unconstrNum mapping.
        virtual void add(const GenAssembledFullM<Scalar> &kel, const int *dofs) ;
        virtual void addDiscreteMass(int dof, Scalar mass);
        virtual void add(int row_dof, int col_dof, Scalar s);

        virtual void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const;
        virtual void mult(const GenVector<Scalar> &rhs, Scalar *result) const;
        virtual void mult(const Scalar *rhs, Scalar *result) const;
        virtual void multAdd(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const;
        virtual void multAdd(const Scalar *rhs, Scalar *result) const;
        virtual void multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const;
        virtual void multSubtract(const Scalar *rhs, Scalar *result) const;
        virtual void transposeMult(const GenVector<Scalar> & rhs, GenVector<Scalar> & result) const;
        virtual void transposeMult(const Scalar *, Scalar *) const;
        virtual void transposeMultAdd(const Scalar *, Scalar *) const;
        virtual void transposeMultSubtract(const Scalar *, Scalar *) const;
        virtual void transposeMultSubtractClaw(const Scalar *, Scalar *, int, int *) const;
        virtual void multSub(const Scalar *, Scalar *) const;
        virtual void multSub(int, const Scalar **, Scalar **) const;
        virtual void multDiag(const Scalar *, Scalar *) const;
        /** \brief Compute \f$ R = R + M\f$.
         *
         * @param R Pointer to result matrix to modify in dense row-wise matrix storage.
         */
        virtual void multIdentity(Scalar *R) const;
        virtual void multIdentity(Scalar **) const;
        virtual void multIdentity(Scalar **v, int start, int stop) const;
        virtual void squareRootMult(Scalar * result);
        virtual void inverseSquareRootMult(Scalar * result);
        virtual GenFullM<Scalar> * getFullMatrix();
        virtual int* getFirstDof();
        virtual int numNodes() const;
        virtual GenFullM<Scalar>* getDiagMatrix(int i);
        virtual ~GenSparseMatrix() = 0;
        virtual Scalar* getBlockScalarMultipliers();
        virtual void setMeanSolver(GenSolver<Scalar> *prc);
        virtual GenSolver<Scalar>* getMeanSolver();
        virtual int getBlockSize();
        virtual int  neqs() const = 0;
        virtual void print();

        void mult(DistrBlockVector<double>&, DistrBlockVector<double>&) { }; // hack to get code to compile

        virtual void matvec(GenVector<Scalar> &rhs, GenVector<Scalar> &result) { mult(rhs, result); }
};


class LMPCons;

/** \brief Structure data of a sparse matrix arranged in a column by column storage.
 *
 * \details The column by column storage is only a matter of vocabulary.
 * Effectively xunonz, rowu have the same structure as a Connectivity.
 *
 * Some users of the structure use a zero based indexing. Some others a one based indexing.
 * The following convention can be assumed and must be respected by all users:
 * xunonz[0] is either 1 for 1 based indexing (Fortran) or 0 for 0 based indexing (C/C++).
 */
class SparseData {
 protected:
  std::vector<int> unconstrNum;
  std::vector<int> constrndNum;
  std::vector<int> xunonz; //!< \brief Pointer into rowu for the start of each column.
  std::vector<int> rowu; //!< Row index of each term.
  std::vector<int> colu; //!< Column index for  an entry
  int numConstrained = 0;
  int numUncon = 0;
  int neq = 0;
 public:
    SparseData() = default;
    // Constructors for data structures of type CuCSparse and CuCComplexSparse

    /// \brief Form the sparse data for constrained to unconstrained DOFs. Generates 0-based indexing.
    SparseData(const Connectivity *con, const DofSetArray *dsa, const int *bc);
	/// \brief Form the sparse data for constrained to unconstrained DOFs. Generates 0-based indexing.
    SparseData(const Connectivity *con, const DofSetArray *dsa, const DofSetArray *c_dsa);
    SparseData(const Connectivity *con, const DofSetArray *dsa, const int *glbmap, const int *glimap);

    // Constructors for data structures of type DBSparseMatrix 
    // and ComplexDBSparseMatrix and Spooles/Mumps
    SparseData(const DofSetArray *dsa, const DofSetArray *c_dsa, const Connectivity *con,
               int expand = 0, int make_colu = 0, bool unsym = false);
    SparseData(const EqNumberer *dsa, const Connectivity *con, const int *rcn, int expand = 0, int make_colu = 0);

    SparseData(const DofSetArray *_dsa, const int *glInternalMap,
               const Connectivity *cn, int expand);

	// This constructor is for the Esmond sparse solver (BLKSparseMatrix)
	/** \brief Form the sparse data for a square matrix in 1-based indexing.
	 *
	 * @param cn Node to node connectivity.
	 * @param eqn Equation numbering for the nodes.
	 * @param expand If false only the upper triangular part is formed.
	 */
	SparseData(const Connectivity *cn, const EqNumberer *eqn, double trbm, bool expand = true);

	/** \brief Form the sparse data for a square matrix in 0-based indexing.
	 *
	 * @param eqn Equation numbering for the nodes.
	 * @param cn Node to node connectivity.
	 */
    SparseData(const EqNumberer *eqn, const Connectivity *cn);

    // KHP: for storing mpcs.
    // ML: Can we make the LMPConst const? Can't implicitely cast LMPCons ** to const LMPConst **
    SparseData(LMPCons **mpc, int numColumns, const DofSetArray *cdsa);
    SparseData(int num, const int *xyzCount, const int *xyzList);

    // for storing G for use in FETI-DP
    SparseData(int numInterface,
               const int *glbmap, int numModes, int ldm);

    virtual ~SparseData();

	size_t nnz() const { return rowu.size(); }

	size_t numCol() const { return xunonz.size()-1; }

	bool usesOneBasedIndexing() const { return xunonz[0] == 1; }

	const std::vector<int> &colPointers() const { return xunonz; }

	const std::vector<int> &rowIndices() const { return rowu; }

    void clean_up();

};

typedef GenSMV<double> SMV; 
typedef GenSparseMatrix<double> SparseMatrix; 
typedef GenSparseMatrix<DComplex> ComplexSparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/SparseMatrix.C>
#endif

#endif
