#include <Math.d/BLKSparseMatrix.C>

template<> 
void
GenBLKSparseMatrix<complex<double> >
   ::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
 int i, j, k, rowi, colj, offset, p1, position, csuper, fstcol, lxbeg, lxend;

 int kndof = kel.dim();

 for(i = 0; i < kndof; ++i ) {
   if((rowi = unconstrNum[dofs[i]]) == -1) continue;
   p1 = invp[rowi] - 1;
   position = xlnz[p1+1];
   csuper = invsuper[p1] - 1;
   fstcol = xsuper[csuper] - 1;
   lxbeg = xlindx[csuper]-1;
   lxend = xlindx[csuper+1]-1;
   for(j = 0; j < kndof; ++j ) {
     if((colj = unconstrNum[dofs[j]]) == -1) continue;
     int irow   = invp[colj] - 1;
     if(irow >= fstcol) {
       offset = lxend - lxbeg;
       for(k=lxbeg; k<lxend; ++k) {
         offset -= 1;
         if(lindx[k]-1 == irow) {
           lnz[position - 2 - offset] += complex<double>(0.0, kel[i][j]);
           break;
         }
       }
     }
   }
 }
}

template<> 
void
GenBLKSparseMatrix<complex<double> >
   ::add(const FullSquareMatrixC &kel, const int *dofs)
{

 int i, j, k, rowi, colj, offset, p1, position, csuper, fstcol, lxbeg, lxend;

 int kndof = kel.dim();

 for( i = 0; i < kndof; ++i ) {
   if( (rowi = unconstrNum[dofs[i]]) == -1 ) continue;
     p1     = invp[rowi] - 1;
     position = xlnz[p1+1];
     csuper = invsuper[p1] - 1;
     fstcol = xsuper[csuper] - 1;
     lxbeg  = xlindx[csuper]-1;
     lxend  = xlindx[csuper+1]-1;
     for( j = 0; j < kndof; ++j ) {
       if( (colj = unconstrNum[dofs[j]]) == -1 ) continue;
       int irow   = invp[colj] - 1;
       if ( irow >= fstcol ) {
         offset = lxend - lxbeg;
         for(k=lxbeg; k<lxend; ++k) {
           offset -= 1;
           if(lindx[k]-1 == irow) {
             lnz[ position - 2 - offset] += kel[i][j];
             break;
           }
         }
       }

     }
 }

}

template<> 
void
GenBLKSparseMatrix<double>
   ::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  fprintf(stderr, "GenBLKSparseMatrix<double> cannot addImaginary\n");
}

template<> 
void
GenBLKSparseMatrix<double>
   ::add(const FullSquareMatrixC &kel, const int *dofs)
{
  fprintf(stderr, "GenBLKSparseMatrix<double> cannot add FullSquareMatrixC\n");
}

template class GenBLKSparseMatrix<double>;
template class GenBLKSparseMatrix<complex<double>>;

