//
// Created by Michel Lesoinne on 2019-01-14.
//

#ifndef FEM_PARDISOSOLVER_H
#define FEM_PARDISOSOLVER_H

#include <Utils.d/Connectivity.h>
#include <Math.d/SparseMatrix.h>
#include "Solver.h"

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getPardiso(const Connectivity *cn, const EqNumberer *_dsa);


#endif //FEM_PARDISOSOLVER_H
