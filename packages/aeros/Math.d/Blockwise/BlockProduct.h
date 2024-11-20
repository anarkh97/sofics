//
// Created by Michel Lesoinne on 1/11/18.
//

#ifndef FEM_BLOCKPRODUCT_H
#define FEM_BLOCKPRODUCT_H

#include <tuple>
#include "BlockDiagonalMatrix.h"

template <typename ... MT>
struct BlockProduct {
	std::tuple<const MT &...> m;
};

template <typename ... ML, typename TR>
auto operator*(const BlockProduct<ML...> &PL, const BlockDiagonalMatrix<TR> &MR) {
	return BlockProduct<ML..., BlockDiagonalMatrix<TR>>
			{ std::tuple_cat(PL.m, std::tuple<const BlockDiagonalMatrix<TR> &>{MR}) };
};

template <typename ML, typename TR>
auto operator*(const BlockDiagonalMatrix<ML> &PL, const BlockDiagonalMatrix<TR> &MR) {
	return BlockProduct<BlockDiagonalMatrix<ML>,BlockDiagonalMatrix<TR>>
			{ std::tuple<const BlockDiagonalMatrix<ML> &, const BlockDiagonalMatrix<TR> &>{PL, MR} };
};

#endif //FEM_BLOCKPRODUCT_H
