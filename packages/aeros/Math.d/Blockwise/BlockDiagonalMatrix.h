//
// Created by Michel Lesoinne on 1/9/18.
//

#ifndef FEM_BLOCKDIAGONALMATRIX_H
#define FEM_BLOCKDIAGONALMATRIX_H


#include <vector>

/** \brief Block diagonal matrix with each block being of the same matrix type. */
template <typename BlockMatrix>
class BlockDiagonalMatrix {
public:
	BlockDiagonalMatrix(std::vector<BlockMatrix> blocks) : blocks(std::move(blocks)) {}
private:
	/** \brief Vector of matrices that form a diagonal matrix. */
	std::vector<BlockMatrix> blocks;
};


#endif //FEM_BLOCKDIAGONALMATRIX_H
