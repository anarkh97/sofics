#include <cstdio>
#include <iostream>

template<class Scalar>
GenDiagMatrix<Scalar>::GenDiagMatrix(const DofSetArray *_dsa)
{
  dsa = _dsa;
  neq = dsa->size();
  v = new Scalar[neq];
  int i;
  for(i = 0; i < neq; ++i)
    v[i] = 0.0;
}

template<class Scalar>
GenDiagMatrix<Scalar>::~GenDiagMatrix()
{
  if(v) { delete [] v; v = 0; }
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::addBoeing(int nl, const int *Kai, const int *Kaj, 
                                 const double *nz, const int *map, Scalar multiplier)
{
  if(neq == 0) return;
  int i, j;
  for(i = 0; i < nl; ++i) {
    if(map[i] == -1) continue;
    int dofI = dsa->getRCN(map[i]);
    if(dofI < 0) continue;
    for(j = Kai[i]; j < Kai[i+1]; ++j) {
      if(map[Kaj[j-1]-1] == -1) continue;
      int dofJ = dsa->getRCN(map[Kaj[j-1]-1]);
      if(dofJ == dofI)  {
        if(dofI > neq) { fprintf(stderr, "Inconsistent data diagmatrix\n");
                         return; }
        v[dofI] += (nz[j-1]*multiplier);
        break;
      }
    }
  }
}

template<class Scalar>
Scalar
GenDiagMatrix<Scalar>::diag(int i) const
{
 return v[i];
}

template<class Scalar>
Scalar &
GenDiagMatrix<Scalar>::diag(int i)
{
 return v[i];
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::add(const FullSquareMatrix &Mat, const int *map)
{
  // local objects
  int MatSize, i, dofI;

  // return for ZERO diagonal matrix size
  if(neq == 0) return;

  // Get the size of the square input matrix 
  MatSize = Mat.dim();

  // add the diagonal elements of the input matrix to the diagonal matrix
  for(i = 0; i < MatSize; i++) {
    // get the DOF in Krr
    dofI = dsa->getRCN( map[ i]);

    // skip constrained DOF's
    if ( dofI < 0) continue;

    // add the diagonal element
    v[dofI] += Mat[i][i];
  }
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::add(const GenFullM<Scalar> &, int, int)
{
  fprintf(stderr,"WARNING: GenDiagMatrix<Scalar>::add(const GenFullM<Scalar> &, int , int) is not implemented\n");
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::add(const GenAssembledFullM<Scalar> &, const int *dofs)
{
  fprintf(stderr,"WARNING: GenDiagMatrix<Scalar>::add(const GenAssembledFullM<Scalar> &, const int *dofs) is not implemented\n");
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::solve(const Scalar *rhs, Scalar *sol)
{
  this->solveTime -= getTime();
  for(int i=0; i<neq; ++i)
    sol[i] = rhs[i]/v[i];
  this->solveTime += getTime();
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::reSolve(Scalar *rhsSol)
{
  this->solveTime -= getTime();
  for(int i = 0 ; i < neq ; ++i)
    rhsSol[i] = rhsSol[i]/v[i];
  this->solveTime += getTime();
}

/*
template<class Scalar>
void
GenDiagMatrix<Scalar>::print()
{
  fprintf(stderr," ***** PRINTING GenDiagMatrix:\n");
  for (int i = 0 ; i < neq ; ++i)
    fprintf(stderr,"  %f",v[i]);
  fprintf(stderr,"\n");
}
*/

template<class Scalar>
void
GenDiagMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const
{
  mult(rhs.data(), result.data());
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result) const
{
 for (int i = 0 ; i < neq ; i++)
   result[i] = v[i]*rhs[i];
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::squareRootMult(Scalar *result)
{
  for (int i = 0 ; i < neq ; i++){
    result[i] *= std::sqrt(v[i]);
  }
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::inverseSquareRootMult(Scalar *result)
{
  for (int i = 0 ; i < neq ; i++){
    result[i] /= std::sqrt(v[i]);
  }
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::addDiscreteMass(int dof, Scalar mass)
{
  int dofI = dsa->getRCN(dof);
  if(dofI > -1) // skip constrained dofs
    v[dofI] += mass; 
}

template<class Scalar>
void
GenDiagMatrix<Scalar>::factor()
{
  // check for zero
  int count = 0;
  double small = 1e-5;
  for (int i = 0 ; i < neq ; i++) {
    if(v[i] == 0.0) { 
      v[i] = small;
      count++;
    }
  }
  if(count > 0) std::cerr << " *** WARNING: " << count << " zero diagonal/s detected in mass matrix set to " << small << std::endl;
}
