#include <cstdio>
#include <Math.d/NBSparseMatrix.C>

template<>
void
GenNBSparseMatrix<DComplex>::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
 int i, j, dof;

 // For each dof figure out to which node it belongs and
 // its unconstrained number 
 
 int kndof = kel.dim();

 // memory for dof to node tables
 int *dton = (int *) dbg_alloca(sizeof(int)*kndof);
 int *dn   = (int *) dbg_alloca(sizeof(int)*kndof);

 for(i = 0; i < kndof; ++i) {
   dof = dsa->getRCN(dofs[i]);  
   if(dof == -1)
     dn[i] = -1;
   else {
     dton[i] = dofToNode[ dof ];
     dn[i]   = dof - dsa->firstdof( dton[i] );
   }
 }
  
// Now add the contribution to the matrix
  for(i = 0; i < kndof; ++i) {
    if(dn[i] == -1) continue;
    for(j=0; j < kndof; ++j) {
      if(dn[j] == -1) continue;
      int matid = con->offset(dton[i],dton[j]);
      allM[matid][dn[i]][dn[j]] += complex<double>(0.0,kel[i][j]);
    }
  }
}

template<>
void
GenNBSparseMatrix<double>
   ::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  fprintf(stderr, "GenNBSparseMatrix<double> cannot addImaginary\n");
}

template<>
void
GenNBSparseMatrix<DComplex>::add(const FullSquareMatrixC &kel, const int *dofs)
{
 int i, j, dof;

 // For each dof figure out to which node it belongs and
 // its unconstrained number 
 
 int kndof = kel.dim();

 // memory for dof to node tables
 int *dton = (int *) alloca(sizeof(int)*kndof);
 int *dn   = (int *) alloca(sizeof(int)*kndof);

 for(i = 0; i < kndof; ++i) {
   dof = dsa->getRCN(dofs[i]);  
   if(dof == -1)
     dn[i] = -1;
   else {
     dton[i] = dofToNode[ dof ];
     dn[i]   = dof - dsa->firstdof( dton[i] );
   }
 }
  
// Now add the contribution to the matrix
  for(i = 0; i < kndof; ++i) {
    if(dn[i] == -1) continue;
    for(j=0; j < kndof; ++j) {
      if(dn[j] == -1) continue;
      int matid = con->offset(dton[i],dton[j]);
      allM[matid][dn[i]][dn[j]] += kel[i][j];
    }
  }
}

template<>
void
GenNBSparseMatrix<double>
   ::add(const FullSquareMatrixC &kel, const int *dofs)
{
  fprintf(stderr, "GenNBSparseMatrix<double> cannot add FullSquareMatrixC\n");
}

#define NBSPARSE_INSTANTIATION_HELPER(Scalar) \
template \
GenFullM<Scalar> * \
GenNBSparseMatrix<Scalar>::getDiagMatrix(int); \
template \
void \
GenNBSparseMatrix<Scalar>::zeroAll(); \
template \
void \
GenNBSparseMatrix<Scalar>::add(const FullSquareMatrix&, const int*); \
template \
GenNBSparseMatrix<Scalar>::GenNBSparseMatrix(const Connectivity*, const ConstrainedDSA*); \
template \
GenNBSparseMatrix<Scalar>::~GenNBSparseMatrix(); \
template \
void \
GenNBSparseMatrix<Scalar>::mult(const Scalar*, Scalar*) const; \
template \
void \
GenNBSparseMatrix<Scalar>::multAdd(const Scalar*, Scalar*) const; \
template \
void \
GenNBSparseMatrix<Scalar>::mult(const GenVector<Scalar>&, GenVector<Scalar>&) const; \
template \
Scalar \
GenNBSparseMatrix<Scalar>::diag(int dof) const; \
template \
Scalar & \
GenNBSparseMatrix<Scalar>::diag(int dof);

NBSPARSE_INSTANTIATION_HELPER(double);
NBSPARSE_INSTANTIATION_HELPER(complex<double>);

