#include <Math.d/CuCSparse.C>

template<>
void
GenCuCSparse<complex<double> >::add(const FullSquareMatrixC &kel, const int *dofs)
{
 int i, j, m, ri;
 int kndof = kel.dim();
 for(i=0; i<kndof; ++i) {
   if((ri = unconstrNum[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(constrndNum[dofs[j]] == -1) continue;
     for(m=xunonz[constrndNum[dofs[j]]]; m<xunonz[constrndNum[dofs[j]]+1]; ++m) {
       if(rowu[m] == ri) {
         Kuc[m] += kel[i][j];
         break;
       }
     }
   }
 }
}

template<>
void
GenCuCSparse<double>::add(const FullSquareMatrixC &kel, const int *dofs)
{
  fprintf(stderr, "GenCuCSparse<double> cannot add FullSquareMatrixC \n");
}

template<>
void
GenCuCSparse<complex<double> >::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
 int i, j, m, ri;
 int kndof = kel.dim();

 for(i=0; i<kndof; ++i) {
   if((ri = unconstrNum[dofs[i]]) == -1) continue;
   for(j=0; j<kndof; ++j) {
     if(constrndNum[dofs[j]] == -1) continue;
     for(m=xunonz[constrndNum[dofs[j]]]; m<xunonz[constrndNum[dofs[j]]+1]; ++m) {
       if(rowu[m] == ri) {
         Kuc[m] += complex<double>(0.0, kel[i][j]);
         break;
       }
     }
   }
 }
}

template<>
void
GenCuCSparse<double>::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  fprintf(stderr, "GenCuCSparse<double> cannot addImaginary\n");
}

template class GenCuCSparse<double>;
template class GenCuCSparse<complex<double>>;
