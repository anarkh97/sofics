#include <cstdio>
#include <Math.d/Skyline.d/BlockSky.C>

template<>
void
GenBlockSky<complex<double> >::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
 // Construct stiffness matrix K (skyA)

 int i, j, ri, rj;

 int kndof = kel.dim();                 // Element stiffness dimension

 for( i = 0; i < kndof; ++i ) {          // Loop over rows.

   if( (ri = rowColNum[dofs[i]]) == -1 ) continue;// Skip constrained dofs

   for( j = 0; j < kndof; ++j ) {          // Loop over columns.

     if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.

     if( (rj = rowColNum[dofs[j]]) == -1 ) continue; // Skip constrained dofs

     skyA[dlp[rj] - rj + ri ] += complex<double>(0.0, kel[i][j]);
   }
 }

}

template<>
void
GenBlockSky<double>
   ::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  fprintf(stderr, "GenBlockSky<double> cannot addImaginary\n");
}

template<>
void
GenBlockSky<complex<double> >::add(const FullSquareMatrixC &kel, const int *dofs)
{
 // Construct stiffness matrix K (skyA)

 int i, j, ri, rj;

 int kndof = kel.dim();                 // Element stiffness dimension

 for( i = 0; i < kndof; ++i ) {          // Loop over rows.

   if( (ri = rowColNum[dofs[i]]) == -1 ) continue;// Skip constrained dofs

   for( j = 0; j < kndof; ++j ) {          // Loop over columns.

     if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.

     if( (rj = rowColNum[dofs[j]]) == -1 ) continue; // Skip constrained dofs

     skyA[dlp[rj] - rj + ri ] += kel[i][j];
   }
 }

}

template<>
void
GenBlockSky<double>
   ::add(const FullSquareMatrixC &kel, const int *dofs)
{
  fprintf(stderr, "GenBlockSky<double> cannot add FullSquareMatrixC\n");
}

template<>
void
GenBlockSky<double>::add(const GenAssembledFullM<complex<double> > &kel, const int *dofs)
{
 fprintf(stderr, "ERROR: add(const GenAssembledFullM<complex<double> > &, int*) is being called"
                 "on an inappropriate matrix\n");
}



template class GenBlockSky<double>;
template class GenBlockSky<complex<double>>;

