//
// Created by Michel Lesoinne on 3/28/18.
//

#ifndef FEM_SPARSEMATRIX_H
#define FEM_SPARSEMATRIX_H

#include "Types.h"
#include "VectorReference.h"

namespace FetiLib {


template <typename Scalar>
class SparseMatrix {
public:
	/** \brief Create a Sparse Matrix from Compressed Sparse Row format data.
	 *
	 * @param nrow The number of rows in the matrix.
	 * @param ncol The number of columns in the matrix.
	 * @param rowStart Vector of pointers indicating where each row starts.
	 * This vector is of size nrow+1 with the last entry being the total number of non-zero elements.
	 * @param columnIndices Vector of column indices. The size of this vector is rowStart[nrow].
	 * @param nzEntries Vector of non-zero entries in the matrix. The size of this vector is rowStart[nrow].
	 * @return A matrix useable by FetiLib.
	 * The user can free the storage used to pass the data immediately after calling the routine.
	 */
	static SparseMatrix fromCSR(local_dof_index nrow, local_dof_index ncol,
	                            VectorReference<const size_t> rowStart,
	                            VectorReference<const local_dof_index> columnIndices, VectorReference<Scalar> nzEntries
	);

private:

};

} // namespace FetiLib


template<typename Scalar>
FetiLib::SparseMatrix<Scalar>
FetiLib::SparseMatrix<Scalar>::fromCSR(FetiLib::local_dof_index nrow, FetiLib::local_dof_index ncol,
                                       FetiLib::VectorReference<const size_t> rowStart,
                                       FetiLib::VectorReference<const FetiLib::local_dof_index> columnIndices,
                                       FetiLib::VectorReference<Scalar> nzEntries) {
	return SparseMatrix();
}

#endif //FEM_SPARSEMATRIX_H
