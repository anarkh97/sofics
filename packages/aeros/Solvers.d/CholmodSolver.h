//
// Created by Michel Lesoinne on 2018-12-14.
//

#ifndef FEM_CHOLDMOD_H
#define FEM_CHOLDMOD_H

#include <Utils.d/Connectivity.h>
#include <Math.d/SparseMatrix.h>
#include "Solver.h"

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
    getCholmod(const Connectivity *cn, const EqNumberer *_dsa);

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getCholmod(const Connectivity *cn, const DofSetArray &dsa, const ConstrainedDSA &cdsa);

#endif //FEM_CHOLDMOD_H
