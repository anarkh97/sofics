#include <cstdio>
#include <Math.d/DBSparseMatrix.C>
#include "DBSparseMatrix.h"


template<>
void
GenDBSparseMatrix<complex<double> >::add(const FullSquareMatrixC &kel, const int *dofs)
{
 int i, j, m;
 int kndof = kel.dim();                 // Dimension of element stiffness.
 for( i = 0; i < kndof; ++i ) {          // Loop over rows.
    if( unconstrNum[dofs[i]] == -1 ) continue;      // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {              // Loop over columns.
       if( dofs[i] > dofs[j] ) continue; // Work with upper symmetric half.
       if( unconstrNum[dofs[j]] == -1 ) continue;   // Skip constrained dofs
       int mstart = xunonz[unconstrNum[dofs[j]]];
       int mstop  = xunonz[unconstrNum[dofs[j]]+1];
       for(m=mstart; m<mstop; ++m)
       {
          if( rowu[m-1] == (unconstrNum[dofs[i]] + 1) )
          {
            unonz[m-1] += kel[i][j];
            break;
          }
       }
       // if(m == mstop) fprintf(stderr," *** ERROR: GenDBSparseMatrix<Scalar>::add\n");
    }
 }
}

template<>
void
GenDBSparseMatrix<double>::add(const FullSquareMatrixC &kel, const int *dofs)
{
 fprintf(stderr, "Error: adding complex to real matrix\n");
}

template<>
void
GenDBSparseMatrix<complex<double> >::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
 int i, j, m;
 int kndof = kel.dim();                 // Dimension of element stiffness.
 for( i = 0; i < kndof; ++i ) {          // Loop over rows.
    if( unconstrNum[dofs[i]] == -1 ) continue;      // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {              // Loop over columns.
       if( dofs[i] > dofs[j] ) continue; // Work with upper symmetric half.
       if( unconstrNum[dofs[j]] == -1 ) continue;   // Skip constrained dofs
       int mstart = xunonz[unconstrNum[dofs[j]]];
       int mstop  = xunonz[unconstrNum[dofs[j]]+1];
       for(m=mstart; m<mstop; ++m)
       {
          if( rowu[m-1] == (unconstrNum[dofs[i]] + 1) )
          {
            unonz[m-1] += complex<double>(0.0, kel[i][j]);
            break;
          }
       }
       // if(m == mstop) fprintf(stderr," *** ERROR: GenDBSparseMatrix<Scalar>::add\n");
    }
 }
}

template<>
void
GenDBSparseMatrix<double>
   ::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  fprintf(stderr, "GenDBSparseMatrix<double> cannot addImaginary\n");
}

/*
template<>
void
GenDBSparseMatrix<double>::transposeMult(const Vector &rhs, Vector &result)
{
  _FORTRAN(sptmv)(unonz, xunonz, rowu, numUncon, rhs.data(), result.data() );
}

template<>
void
GenDBSparseMatrix<DComplex>::transposeMult(const ComplexVector &rhs, ComplexVector &result)
{
  fprintf(stderr, "GenDBSparseMatrix<DComplex>::transposeMult(const ComplexVector &rhs, ComplexVector &result) not implemented \n");
}
*/

template<>
void
GenDBSparseMatrix<double>::multcomplex(const DComplex *rhs, DComplex *result) const
{
 int nn = numUncon;
  _FORTRAN(cdspsmvp)(nn, unonz.data(), xunonz.data(), rowu.data(), rhs, result );
}

template<>
void
GenDBSparseMatrix<DComplex>::multcomplex(const DComplex *rhs, DComplex *result) const
{
  mult(rhs, result);
}

template class GenDBSparseMatrix<double>;
template class GenDBSparseMatrix<complex<double>>;
