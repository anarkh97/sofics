#include <complex>
#include <cstdlib>
#include <Math.d/DBSparseMatrix.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/Vector.h>
#include <Driver.d/Communicator.h>

extern "C"      {
void    _FORTRAN(micspsmvp)(const int&, const double*, const int*, const int*, const double*, double*);
void    _FORTRAN(sptmv)(const double*, const int*, const int*, const int&, const double*, double*);
void    _FORTRAN(cspsmvp)(const int&, const complex<double> *, const int*, const int*, const complex<double> *,
                          complex<double> *);
void    _FORTRAN(cdspsmvp)(const int&, const double*, const int*, const int*, const complex<double>*, complex<double>*);
}

inline void Tcspsmvp(const int& a, const double* b, const int* c, const int* d, const double* e, double* f)
{
_FORTRAN(micspsmvp)(a,b,c,d,e,f);
}

inline void Tcspsmvp(const int& a, const complex<double> * b, const int* c, const int* d, const complex<double> * e,
                     complex<double> * f)
{
_FORTRAN(cspsmvp)(a,b,c,d,e,f);
}


template<class Scalar>
GenDBSparseMatrix<Scalar>::~GenDBSparseMatrix()
{
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::clean_up()
{
 unonz.clear();
 unonz.shrink_to_fit();
 scale.clear();
 scale.shrink_to_fit();
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::zeroAll()
{
  int i;
  for(i=0; i < xunonz[numUncon]; ++i)
    ScalarTypes::initScalar(unonz[i], 0.0, 0.0);
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::print()
{
  for (int dof = 0; dof < numUncon; ++dof)
    for (int i = xunonz[dof]; i<xunonz[dof+1]; ++i)
      std::cerr << " " << rowu[i-1] - 1 << " " << dof << " " << unonz[i-1] << "\n";
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::print(char *fileName)
{
 FILE *file;
 if(fileName)
   file=fopen(fileName,"w");
 else
   file=stderr;
 int i, m;
 for(i=0; i<numUncon; ++i)
   print1(i,file);
 if(fileName)
   fclose(file);
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::print1(int dof, FILE *fid)
{
 FILE *file = (fid) ? fid : stderr;

 int i;
 for(i=xunonz[dof]; i<xunonz[dof+1]; ++i)
   fprintf(file," %d %d %24.16e\n", rowu[i-1]-1, dof, ScalarTypes::Real(unonz[i-1]));
}

template<class Scalar>
GenFullM<Scalar> *
GenDBSparseMatrix<Scalar>::getFullMatrix()
{
  GenFullM<Scalar> *ret = new GenFullM<Scalar>(dim());
  ret->zero();
  for(int i=0; i<numUncon; ++i) {
    for(int j=xunonz[i]; j<xunonz[i+1]; ++j) {
      (*ret)[rowu[j-1]-1][i] = unonz[j-1];
    }
  }
  return ret;
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::add(int idof, int jdof, Scalar s)
{
 if((idof < 0) || (jdof < 0)) return;
 int dofsi = (jdof > idof) ? idof : jdof;
 int dofsj = (jdof > idof) ? jdof : idof;

 if(unconstrNum[dofsi] == -1 || unconstrNum[dofsj] == -1) return;
 int mstart = xunonz[unconstrNum[dofsj]];
 int mstop  = xunonz[unconstrNum[dofsj]+1];
 for(int m=mstart; m<mstop; ++m) {
   if(rowu[m-1] == (unconstrNum[dofsi] + 1)) {
     unonz[m-1] += s;
     break;
   }
 }
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::add(const FullSquareMatrix &kel, const int *dofs)
{
 int i, j, m;
 int kndof = kel.dim();                           // Dimension of element stiffness.
 for(i = 0; i < kndof; ++i) {                     // Loop over rows.
    if(dofs[i] == -1) continue;                   // Skip constrained dofs
    if(unconstrNum[dofs[i]] == -1) continue;      // Skip irrelevant dofs, EC, 20070820
    for(j = 0; j < kndof; ++j) {                       // Loop over columns.
       if(dofs[i] > dofs[j]) continue;                 // Work with upper symmetric half.
       if(unconstrNum[dofs[j]] == -1) continue;        // Skip constrained dofs
       if(dofs[j] == -1) continue;                     // Skip irrelevant dofs, EC, 20070820
       int mstart = xunonz[unconstrNum[dofs[j]]];
       int mstop  = xunonz[unconstrNum[dofs[j]]+1];
       for(m=mstart; m<mstop; ++m) {
          if( rowu[m-1] == (unconstrNum[dofs[i]] + 1) ) {
            unonz[m-1] += kel[i][j];
            break;
          }
       }
    }
 }
}

template<class Scalar>
void GenDBSparseMatrix<Scalar>::add(const GenFullM<Scalar> &kel, gsl::span<const int> dofs) {
    int kndof = kel.dim();
    if(dofs.size() != kndof)
        throw std::logic_error("Not the right number of dofs in GenDBSparseMatrix add method.");
    for (int i = 0; i < kndof; ++i) {                     // Loop over rows.
        if (dofs[i] == -1) continue;                   // Skip constrained dofs
        if (unconstrNum[dofs[i]] == -1) continue;      // Skip irrelevant dofs, EC, 20070820
        for (int j = 0; j < kndof; ++j) {                       // Loop over columns.
            if(dofs[i] > dofs[j]) continue;                 // Work with upper symmetric half.
            if(unconstrNum[dofs[j]] == -1) continue;        // Skip constrained dofs
            int mstart = xunonz[unconstrNum[dofs[j]]];
            int mstop  = xunonz[unconstrNum[dofs[j]]+1];
            for (int m=mstart; m<mstop; ++m) {
                if( rowu[m-1] == (unconstrNum[dofs[i]] + 1) ) {
                    unonz[m-1] += kel[i][j];
                    break;
                }
            }
        }
    }
}

template<class Scalar>
void GenDBSparseMatrix<Scalar>::add(const GenAssembledFullM<Scalar> &kel, const int *dofs) {
    int kndof = kel.dim();
    for (int i = 0; i < kndof; ++i) {                     // Loop over rows.
        if (dofs[i] == -1) continue;                   // Skip constrained dofs
        for (int j = 0; j < kndof; ++j) {                       // Loop over columns.
            if(dofs[i] > dofs[j]) continue;                 // Work with upper symmetric half.
            if(dofs[j] == -1) continue;                     // Skip irrelevant dofs
            int mstart = xunonz[dofs[j]];
            int mstop  = xunonz[dofs[j]+1];
            for (int m=mstart; m<mstop; ++m) {
                if( rowu[m-1] == (dofs[i] + 1) ) {
                    unonz[m-1] += kel[i][j];
                    break;
                }
            }
        }
    }
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::add(const GenFullM<Scalar> &knd, int fRow, int fCol)
{
 int i, j, m;

 int nrow  = knd.numRow();
 int nCol  = knd.numCol();

 for(i = 0; i < nrow; ++i) {
    int rowi = fRow + i;
    for(j = 0; j < nCol; ++j) {
       int colj = fCol + j;
//     if( rowi > colj ) continue; // Work with upper symmetric half.
       int mstart = xunonz[colj] - 1;
       int mstop  = xunonz[colj+1] - 1;
       for(m=mstart; m<mstop; ++m) {
          if(rowu[m] == (rowi+1)) {
            unonz[m] += knd[i][j];
            break;
          }
       }
    }
 }
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::addBoeing(int nl, const int *Kai, const int *Kaj,
                                     const double *nz, int *map, Scalar multiplier)
{
 int i, j;
 for(i = 0; i < nl; ++i) {
   if(map[i] >= neq) {
      fprintf(stderr, "Out of bounds, %d %d\n",map[i],numUncon);
      return ;
   }
   if(map[i] == -1) continue;
   int firstDof = unconstrNum[map[i]];
   if(firstDof < 0) continue;
   for(j = Kai[i]; j < Kai[i+1]; ++j) {
     if(map[Kaj[j-1]-1] == -1) continue;
     int secondDof = unconstrNum[map[Kaj[j-1]-1]];
     if(secondDof < 0) continue;
     int dofI, dofJ;
     if(firstDof >= secondDof) {
        dofI = firstDof;
        dofJ = secondDof;
     } else {
        dofJ = firstDof;
        dofI = secondDof;
     }
     int mstart = xunonz[dofI];
     int mstop  = xunonz[dofI+1];
     int m;
     for(m=mstart; m<mstop; ++m) {
       if(m > xunonz[numUncon]) { fprintf(stderr, "Exceeded\n"); return; }
       if(m <= 0) { fprintf(stderr, "Underbound\n");  return; }
       if(rowu[m-1] == dofJ+1) {
          unonz[m-1] += (nz[j-1]*multiplier);
          break;
       }
     }
     if(m==mstop) {
       std::cerr << "Warning: m = mstop in DBSparseMatrix::addBoeing(), i = " << i << ", map[i] = " << map[i] << std::endl;
       //exit(-1);
     }
   }
 }
}

template<class Scalar>
GenDBSparseMatrix<Scalar>::GenDBSparseMatrix(const Connectivity *cn, const DofSetArray *dsa, const int *rCN)
: SparseData(dsa,cn,rCN)
{
  // ... Allocate memory for unonz & initialize to zero
  unonz.resize(xunonz[numUncon]);

  // ... INITIALIZE THE VECTOR CONTAINING THE SPARSE MATRIX TO ZERO
  zeroAll();

  isScaled=0;
}

template<class Scalar>
GenDBSparseMatrix<Scalar>::GenDBSparseMatrix(const Connectivity *cn, const DofSetArray *_dsa, const ConstrainedDSA *c_dsa) :
  SparseData(_dsa,c_dsa,cn)
{
  // ... Allocate memory for matrix value vector unonz
  unonz.resize(xunonz[numUncon]);
  // ... INITIALIZE THE VECTOR CONTAINING THE SPARSE MATRIX TO ZERO
  zeroAll();

  isScaled=0;
}

template<class Scalar>
GenDBSparseMatrix<Scalar>::GenDBSparseMatrix(const Connectivity *cn, const EqNumberer *eqNums)
 : SparseData(cn, eqNums,1.0E-6,false)
{
  // ... Allocate memory for matrix value vector unonz
  unonz.resize(xunonz[numUncon]-xunonz[0]);

  // ... INITIALIZE THE VECTOR CONTAINING THE SPARSE MATRIX TO ZERO
  zeroAll();

  isScaled=0;
}

template<class Scalar>
double
GenDBSparseMatrix<Scalar>::getMemoryUsed() const
{
 return 8.0*xunonz[numUncon]/(1024.0*1024.0);
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result) const
{
 int nn = numUncon;
 Tcspsmvp(nn, unonz.data(), xunonz.data(), rowu.data(), rhs, result);
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, Scalar *result) const
{
  Tcspsmvp(numUncon, unonz.data(), xunonz.data(), rowu.data(), rhs.data(), result);
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::transposeMult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const
{
  // PJSA 11/2/09 (same as mult since this class is for symmetric matrices)
  mult(rhs, result);
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::transposeMult(const Scalar * rhs, Scalar * result) const
{
   mult(rhs,result);
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::multAdd(const Scalar *rhs, Scalar *result) const
{
 int nn = numUncon;
 Scalar *tmp = (Scalar *) dbg_alloca(sizeof(Scalar)*nn);
 Tcspsmvp(nn, unonz.data(), xunonz.data(), rowu.data(), rhs, tmp);
 for(int i=0; i<nn; ++i) result[i] += tmp[i];
}


template<class Scalar>
void
GenDBSparseMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const
{
 Tcspsmvp(numUncon, unonz.data(), xunonz.data(), rowu.data(), rhs.data(), result.data());
}

template<class Scalar>
Scalar
GenDBSparseMatrix<Scalar>::diag(int dof) const
{
  int m, mstart, mstop;

  mstart = xunonz[dof]-1;
  mstop  = xunonz[dof+1]-1;

  for(m=mstart; m<mstop; ++m) {
    if(rowu[m]-1 == dof) {
      if(unonz[m] == 0.0)
        return (1.0);
      else
        return unonz[m];
    }
  }

  throw "GenDBSparseMatrix<Scalar>::diag - 1 - this should never be reached";
}

template<class Scalar>
Scalar &
GenDBSparseMatrix<Scalar>::diag(int dof)
{
  int m, mstart, mstop;

  mstart = xunonz[dof]-1;
  mstop  = xunonz[dof+1]-1;

  for(m=mstart; m<mstop; ++m) {
    if(rowu[m]-1 == dof) {
        return unonz[m];
    }
  }

  throw "GenDBSparseMatrix<Scalar>::diag - 2 - this should never be reached";
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::makeIdentity()
{
 for(int dof = 0; dof < numUncon; ++dof) {
   int m, mstart, mstop;
   mstart = xunonz[dof]-1;
   mstop  = xunonz[dof+1]-1;

   for(m=mstart; m<mstop; ++m) {
     if(rowu[m]-1 == dof)
       unonz[m] = 1.0;
     else unonz[m] = 0.0;
   }
 }
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::invertDiag()
{
  int m, mstart, mstop;

  int dof;
  for(dof=0; dof<numUncon; ++dof) {
    mstart = xunonz[dof]-1;
    mstop  = xunonz[dof+1]-1;

    for(m=mstart; m<mstop; ++m) {
      if(rowu[m]-1 == dof)
        unonz[m] = 1.0/unonz[m];
    }
  }
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  // dof is in unconstrained dof numbering
  int cdof = unconstrNum[dof];
  if(cdof < 0) return;
  int diagLocator = xunonz[cdof+1]-2; // This should be the diagonal
  unonz[diagLocator] += dmass;
}

template<class Scalar>
long
GenDBSparseMatrix<Scalar>::size() const
{
 return (numUncon) ? xunonz[numUncon] : 0;
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::multDiag(const Scalar *x, Scalar *b) const
{
  int i;
  for(i=0; i<numUncon; ++i) {
    Scalar dv = diag(i);
    b[i] = dv*x[i];
  }
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::multDiag(int numRHS, const Scalar **x, Scalar **b) const
{
  int i;
  int n = 0;

  int multiple = 2;
  switch(multiple) {
    default:
    case 2:
      {
       for( ; n<numRHS-1; n += 2)
        for(i=0; i<numUncon; ++i) {
          Scalar dv = diag(i);
          b[n][i]   = dv*x[n][i];
          b[n+1][i] = dv*x[n+1][i];
        }
      }
    case 1:
      {
       for( ; n<numRHS; ++n)
        for(i=0; i<numUncon; ++i) {
          Scalar dv = diag(i);
          b[n][i] = dv*x[n][i];
        }
      }
      break;
  }

}

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
 communicator->globalSum(xunonz[numUncon], unonz.data());
#endif
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::symmetricScaling()
{
 if(dim() == 0) return;
 isScaled=1;

 scale.resize(numUncon);
 int i,m;
 for(i=0; i<numUncon; ++i)
   scale[i] = Scalar(1.0)/ScalarTypes::sqrt(diag(i));

 unonz[0] *= scale[0]*scale[0];

 for(i=1; i<numUncon; ++i) {
   for(m=begin(i); m<end(i); ++m)
     unonz[m-1] *= scale[rowu[m-1]-1]*scale[i];
 }
}

template<class Scalar>
void
GenDBSparseMatrix<Scalar>::applyScaling(Scalar *vector)
{
 if(isScaled) {
   int i;
   for(i=0; i<numUncon; ++i)
     vector[i] *= (scale[i]*scale[i]);
 }
}
