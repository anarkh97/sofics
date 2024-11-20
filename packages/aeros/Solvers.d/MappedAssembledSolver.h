/*
 * MappedAssembledSolver.h
 *
 *  Created on: Apr 15, 2009
 *      Author: michel
 */

#ifndef MAPPEDASSEMBLEDSOLVER_H_
#define MAPPEDASSEMBLEDSOLVER_H_
#include <iostream>

#include <Utils.d/dbg_alloca.h>
#include <Utils.d/resize_array.h>
#include <Math.d/FullMatrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/dofset.h>

typedef GenVectorSet<double> VectorSet;

struct DOFMap {
  int ndofs;
  int *dofs;
  double *coefs;
  double rhs; // PJSA extension for non-zero lmpc rhs
  DOFMap() : ndofs(0), dofs(NULL), coefs(NULL), rhs(0.0) {}
  ~DOFMap() { if(dofs) delete [] dofs; if(coefs) delete [] coefs; }
};

DOFMap *getDofMaps(int size);

class CoordinateMap
{
    long step;
    ResizeArray<long> flag;
    ResizeArray<int> subScript;

  protected:
    int numEqs;
    DOFMap *dofMaps;
    int numMappedEqs;
    DOFMap *eqMaps;

    template <class Scalar, class MatrixType >
    MatrixType map(const MatrixType &mat, const int *dofs, ResizeArray<int> &mappedDofs,
                   ResizeArray<Scalar> &m) {
      // TODO this function is not thread safe
      // Create the set of DOFs that are used.
      step += 1;
      int nMappedDOFs = 0;
      int ndofs = mat.dim();
      // Whether the mapping is a one to one renumbering only
      bool isRenumbering = true;
      for(int i = 0; i < mat.dim(); ++i) {
        if(dofMaps[dofs[i]].ndofs != 1 || flag[dofMaps[dofs[i]].dofs[0]] == step ||
            dofMaps[dofs[i]].coefs[0] != 1.0)
          isRenumbering = false;
        //std::cerr << "mapping " << dofs[i] << " to" ;
        for(int j = 0; j < dofMaps[dofs[i]].ndofs; ++j) {
          //std::cerr << " " << dofMaps[dofs[i]].dofs[j]
          //          << " with flag " << flag[dofMaps[dofs[i]].dofs[j]];
          if(flag[dofMaps[dofs[i]].dofs[j]] != step) {
            flag[dofMaps[dofs[i]].dofs[j]] = step;
            mappedDofs[nMappedDOFs] = dofMaps[dofs[i]].dofs[j];
            subScript[dofMaps[dofs[i]].dofs[j]] = nMappedDOFs;
            nMappedDOFs++;
          }
        }
        //std::cerr << std::endl;
      }

/* TODO this is ok if there is no rhs
      // If we only have a renumbering, we are done and can return;
      if(isRenumbering)
        return mat;
*/
      // If we have a more complex mapping, we do the transformation of the matrix
      // The intermediate matrix needs to be nMappedDOFs x ndofs. Allocate memory for it.
      Scalar *dataSpace = (Scalar *) dbg_alloca(nMappedDOFs*ndofs*sizeof(Scalar));
      GenStackFSFullMatrix<Scalar> mTmp(nMappedDOFs, ndofs, dataSpace);
      // Multiply by the mapping matrix on the left
      mTmp = (Scalar) 0;
      // mTmp = M^T * mat
      // mTmp[i][j] = sum_k ( M[k][i] * mat[k][j] )
      for(int k = 0; k < ndofs; ++k) {
        DOFMap &dMapK = dofMaps[dofs[k]];
        for(int iptr = 0; iptr < dMapK.ndofs; ++iptr) {
          int i = subScript[dMapK.dofs[iptr]];
          for(int j = 0; j < ndofs; ++j)
            mTmp[i][j] += dMapK.coefs[iptr]*mat[k][j];
        }
      }

      // Multiply the mapping on the right.
      m[ndofs+nMappedDOFs*nMappedDOFs] = 0; // resize the workspace, include space for vector containing contribution of non-zero lmpc rhs
      GenStackVector<Scalar> vec(ndofs, m.data());
      vec.zero();
      MatrixType res(nMappedDOFs, m.data()+ndofs);
      res.zero();
      // res = mTmp*M
      // res[i][j] =  sum_k ( mTmp[i][k] * M[k][j] )
      for(int k = 0; k < ndofs; ++k) {
        DOFMap &dMapK = dofMaps[dofs[k]];
        for(int jptr = 0; jptr < dMapK.ndofs; ++jptr) {
          int j = subScript[dMapK.dofs[jptr]];
          for(int i = 0; i < nMappedDOFs; ++i)
            res[i][j] += mTmp[i][k]*dMapK.coefs[jptr];
        }
        for(int i = 0; i < ndofs; ++i)
          vec[i] += mat[i][k]*dMapK.rhs;
      }

      return res;
    }

  public:
    CoordinateMap(int nEq, DOFMap *baseMap, int nMappedEq, DOFMap *eqMap):
      numEqs(nEq), dofMaps(baseMap), numMappedEqs(nMappedEq), eqMaps(eqMap),
      step(0), flag(0), subScript(0) {}
    ~CoordinateMap() { delete [] dofMaps; delete [] eqMaps; }
};


template<class BaseSolver, class Scalar, class Map = CoordinateMap>
class MappedAssembledSolver : public BaseSolver, public Map
{
    ResizeArray<int> mappedDofs;
    ResizeArray<double> dd;
    ResizeArray<complex<double> > dz;
    ResizeArray<Scalar> ds;
    GenVector<Scalar> f;
    ConstrainedDSA *cdsa;
    Rbm *rbm;

  public:
    template <class BaseArgs>
    MappedAssembledSolver(BaseArgs &ba, int nEq, DOFMap *baseMap, int nMappedEq, DOFMap *eqMap, ConstrainedDSA *_cdsa, Rbm *_rbm=0) :
      BaseSolver(ba), Map(nEq, baseMap, nMappedEq, eqMap), mappedDofs(0), dd(0), dz(0), ds(0), f(nMappedEq), cdsa(_cdsa), rbm(_rbm) {
        f.zero();
    }
    ~MappedAssembledSolver() {}

    // Note: None of the add routines are thread safe
    void add(const FullSquareMatrix &mat, const int *dofs) {
      FullSquareMatrix m = Map::map(mat, dofs, mappedDofs, dd);
      m.setMyval(0);
      BaseSolver::add(m, mappedDofs.data());
      for(int i = 0; i < mat.dim(); ++i) {  // PJSA extension for non-zero lmpc rhs
        int j = cdsa->getRCN(dofs[i]);      // the first m.dim() entries in the workspace dd now store the
        if(j > -1) f[j] += dd[i];           // contribution to f 
      }                                     // TODO f should be zero'd when solver is rebuilt
    }                                       // TODO consider just building Kcc/Kuc and then making f inside solve
                                            //      this would be better for the linear case if just the rhs is changing
                                            //      with time, but the element matrices are the same
    void addImaginary(const FullSquareMatrix &mat, const int *dofs) {
      FullSquareMatrix m = Map::map(mat, dofs, mappedDofs, dd);
      BaseSolver::addImaginary(m, mappedDofs.data());
    }

    void add(const FullSquareMatrixC &mat, const int *dofs) {
      FullSquareMatrixC m = Map::map(mat, dofs, mappedDofs, dz);
      BaseSolver::add(m, mappedDofs.data());
    }

    void add(const GenFullM<Scalar> &mat, const int *dofs) {
      GenFullM<Scalar> m = Map::map(mat, dofs, mappedDofs, dz);
      BaseSolver::add(m, mappedDofs.data());
    }

    void add(const GenFullM<Scalar> &mat, int fi, int fj) {
      std::cerr << "MappedAssembledSolver::add(const GenFullM<Scalar>&, int, int) is not implemented\n";
    }

    void add(const GenAssembledFullM<Scalar> &mat, const int *dofs) {
      std::cerr << "MappedAssembledSolver::add(GenAssembledFullM<Scalar>&, int*) is not implemented\n";
    }

    void addDiscreteMass(int dof, Scalar s) {
      GenSparseMatrix<Scalar>::addDiscreteMass(dof, s);
    }

    void add(int dofi, int dofj, Scalar s) {
      GenSparseMatrix<Scalar>::add(dofi, dofj, s);
    }

    void factor() {
      BaseSolver::factor();
    } 

    void reSolve(Scalar *rhs);
    void reSolve(GenVector<Scalar> &rhs) { reSolve(rhs.data()); }
    void reSolve(int nRHS, Scalar **rhs) { for(int i=0; i<nRHS; ++i) reSolve(rhs[i]); }
    void reSolve(int nRHS, Scalar  *rhs) { for(int i=0; i<nRHS; ++i) reSolve(rhs + i*Map::numMappedEqs); }
    void reSolve(int nRHS, GenVector<Scalar> *rhs) { for(int i=0; i<nRHS; ++i) reSolve(rhs[i].data()); }

    void solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution) {
      solution = rhs;
      reSolve(solution.data());
    }
    void solve(const Scalar *rhs, Scalar *solution) {
      std::cerr << "MappedAssembledSolver::solve(Scalar *, Scalar *) is not implemented\n";
    }

    double getResidualNorm(const GenVector<Scalar> &v) {
      GenVector<Scalar> v1(Map::numMappedEqs);
      GenVector<Scalar> v2(BaseSolver::neqs());
      v2.zero();
      for(int i = 0; i < Map::numMappedEqs; ++i) {
        v1[i] = Map::eqMaps[i].rhs;
        for(int j = 0; j < Map::eqMaps[i].ndofs; ++j) {
          v2[Map::eqMaps[i].dofs[j]] += Map::eqMaps[i].coefs[j]*(v[i]-f[i]);
        }
      }
      return sqrt(v1.sqNorm()+v2.sqNorm());
    }

    void getRBMs(double *rbms) {
      std::cerr << "MappedAssembledSolver::getRBMs(double *) is not implemented\n";
    }

    void getRBMs(Vector *rbms);

    void getRBMs(VectorSet &rbms);
};

template<class BaseSolver, class Scalar, class Map>
void MappedAssembledSolver<BaseSolver, Scalar, Map>::reSolve(Scalar *s)
{
  GenVector<Scalar> s2(BaseSolver::neqs());
  s2.zero();
  for(int i = 0; i < Map::numMappedEqs; ++i) {
    Scalar delta = s[i] - f[i];
    for(int j = 0; j < Map::eqMaps[i].ndofs; ++j) {
      s2[Map::eqMaps[i].dofs[j]] += Map::eqMaps[i].coefs[j]*delta; // PJSA extension for non-zero lmpc rhs
    }
  }
  if(BaseSolver::neqs() > 0)
    BaseSolver::reSolve(s2.data());
  for(int i = 0; i < Map::numMappedEqs; ++i) {
    s[i] = Map::eqMaps[i].rhs; // PJSA extension for non-zero lmpc rhs
    for(int j = 0; j < Map::eqMaps[i].ndofs; ++j)
      s[i] += s2[Map::eqMaps[i].dofs[j]]*Map::eqMaps[i].coefs[j];
  }
}

template<class BaseSolver, class Scalar, class Map>
void MappedAssembledSolver<BaseSolver, Scalar, Map>::getRBMs(Vector *rbms)
{
  if(rbm && rbm->numRBM() == BaseSolver::numRBM()) {
    rbm->getRBMs(rbms);
    return;
  }
  Vector *rbms2 = new Vector[BaseSolver::numRBM()];
  for(int iRBM = 0; iRBM < BaseSolver::numRBM(); ++iRBM) {
    rbms2[iRBM].initialize(BaseSolver::neqs());
    rbms2[iRBM].zero();
  }
  BaseSolver::getRBMs(rbms2);
  for(int iRBM = 0; iRBM < BaseSolver::numRBM(); ++iRBM) {
    const Vector &s = rbms[iRBM];
    const Vector &s2 = rbms2[iRBM];
    for(int i = 0; i < Map::numMappedEqs; ++i) {
      s[i] = 0;
      for(int j = 0; j < Map::eqMaps[i].ndofs; ++j)
        s[i] += s2[Map::eqMaps[i].dofs[j]]*Map::eqMaps[i].coefs[j];
    }
  }
  delete [] rbms2;
}

template<class BaseSolver, class Scalar, class Map>
void MappedAssembledSolver<BaseSolver, Scalar, Map>::getRBMs(VectorSet &rbms)
{
  if(rbm && rbm->numRBM() == BaseSolver::numRBM()) {
    rbm->getRBMs(rbms);
    return;
  }
  VectorSet rbms2(BaseSolver::numRBM(), BaseSolver::neqs(), 0.0);
  BaseSolver::getRBMs(rbms2);
  for(int iRBM = 0; iRBM < BaseSolver::numRBM(); ++iRBM) {
    const Vector &s = rbms[iRBM];
    const Vector &s2 = rbms2[iRBM];
    for(int i = 0; i < Map::numMappedEqs; ++i) {
      s[i] = 0;
      for(int j = 0; j < Map::eqMaps[i].ndofs; ++j)
        s[i] += s2[Map::eqMaps[i].dofs[j]]*Map::eqMaps[i].coefs[j];
    }
  }
}

#endif /* MAPPEDASSEMBLEDSOLVER_H_ */
