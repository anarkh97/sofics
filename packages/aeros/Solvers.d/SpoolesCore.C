#include <Solvers.d/Spooles.C>

template<>
void
GenSpoolesSolver<double>::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  fprintf(stderr, "GenSpoolesSolver<double> cannot addImaginary\n");
}

template<>
void
GenSpoolesSolver<complex<double> >
   ::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  int i, j, m, mstart, mstop;
  int kndof = kel.dim();                       // Dimension of element stiff.

  for(i = 0; i < kndof; ++i) {                 // Loop over rows.
    if(unconstrNum[dofs[i]] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {               // Loop over columns.
      if(unconstrNum[dofs[j]] == -1) continue; // Skip constrained dofs
      if(unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue;
      mstart = xunonz[unconstrNum[dofs[j]]];
      mstop  = xunonz[unconstrNum[dofs[j]]+1];
      for(m=mstart; m<mstop; ++m) {
        // if(rowu[m-1] > unconstrNum[dofs[j]]+1)
        //   fprintf(stderr, "Bigger: %d %d\n", rowu[m-1]-1, unconstrNum[dofs[j]]);
        if(rowu[m-1] == (unconstrNum[dofs[i]] + 1)) {
          unonz[m-1] += complex<double>(0.0, kel[i][j]);
          break;
        }
      }
    }
  }
}

template<>
void
GenSpoolesSolver<double>::add(const FullSquareMatrixC &kel, const int *dofs)
{
  fprintf(stderr,"GenSpoolesSolver<double>::add(const FullSquareMatrixC &kel, const int *dofs) is not implemented.\n");
}

template<>
void
GenSpoolesSolver<complex<double> >::add(const FullSquareMatrixC &kel, const int *dofs)
{
  int i, j, m, mstart, mstop;
  int kndof = kel.dim();                       // Dimension of element stiff.

  for(i = 0; i < kndof; ++i) {                 // Loop over rows.
    if(unconstrNum[dofs[i]] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {               // Loop over columns.
      if(unconstrNum[dofs[j]] == -1) continue; // Skip constrained dofs
      if(unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue; 
      mstart = xunonz[unconstrNum[dofs[j]]];
      mstop  = xunonz[unconstrNum[dofs[j]]+1];
      for(m=mstart; m<mstop; ++m) {
        if(rowu[m-1] == (unconstrNum[dofs[i]] + 1)) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

template class GenSpoolesSolver<double>;
template class GenSpoolesSolver<complex<double>>;

