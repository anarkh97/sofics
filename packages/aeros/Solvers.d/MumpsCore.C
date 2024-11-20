#include <Solvers.d/Mumps.C>

template<>
void
GenMumpsSolver<double>::addImaginary(const FullSquareMatrix &kel, const int *dofs)
{
  fprintf(stderr, "GenMumpsSolver<double> cannot addImaginary\n");
}

template<>
void
GenMumpsSolver<complex<double> >::addImaginary(const FullSquareMatrix &kel, const int *dofs)
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
GenMumpsSolver<double>::add(const FullSquareMatrixC &kel, const int *dofs)
{
  fprintf(stderr, "GenMumpsSolver<double>::add(const FullSquareMatrixC &kel, int *dofs is not implemented\n");
}

template<>
void
GenMumpsSolver<complex<double> >::add(const FullSquareMatrixC &kel, const int *dofs)
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
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

#ifdef USE_MUMPS
template<>
void
GenMumpsSolver<DComplex>::copyToMumpsLHS(mumps_double_complex *&m, DComplex *&d, int len)
{
/*
  m = new mumps_double_complex[len]; // this is inefficient ... doubles memory use
  for(int i=0; i<len; ++i) {
    m[i].r = d[i].real();
    m[i].i = d[i].imag();
  }
  delete [] d; d = 0;
*/
  m = reinterpret_cast<mumps_double_complex *>(d);
}

template<>
void
GenMumpsSolver<double>::copyToMumpsLHS(double *&m, double *&d, int len)
{
  m = d;
}

template<>
void
GenMumpsSolver<DComplex>::copyToMumpsRHS(mumps_double_complex *&m, const DComplex *d, int len)
{
  m = new mumps_double_complex[len];
  for(int i=0; i<len; ++i) {
    m[i].r = d[i].real();
    m[i].i = d[i].imag();
  }
}

template<>
void
GenMumpsSolver<DComplex>::copyFromMumpsRHS(DComplex *d, mumps_double_complex *m, int len)
{
  for(int i=0; i<len; ++i) {
    DComplex di(m[i].r,m[i].i);
    d[i] = di;
  }
  delete [] m;
}

template<>
void
GenMumpsSolver<double>::copyToMumpsRHS(double *&m, const double *d, int len)
{
  m = new double[len];
  for(int i=0; i<len; ++i) m[i] = d[i];
}

template<>
void
GenMumpsSolver<double>::copyFromMumpsRHS(double *d, double *m, int len)
{
  for(int i=0; i<len; ++i) d[i] = m[i];
  delete [] m;
}

template<>
void
GenMumpsSolver<DComplex>::copyToMumpsRHS(mumps_double_complex *&m, DComplex **d, int len, int nRHS)
{
  m = new mumps_double_complex[len*nRHS];
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      m[I].r = d[i][j].real();
      m[I].i = d[i][j].imag();
      I++;
    }
  }
}

template<>
void
GenMumpsSolver<DComplex>::copyFromMumpsRHS(DComplex **d, mumps_double_complex *m, int len, int nRHS)
{
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      DComplex dij(m[I].r,m[I].i);
      d[i][j] = dij;
      I++;
    }
  }
  delete [] m;
}

template<>
void
GenMumpsSolver<double>::copyToMumpsRHS(double *&m, double **d, int len, int nRHS)
{
  m = new double[len*nRHS];
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      m[I] = d[i][j];
      I++;
    }
  }
}

template<>
void
GenMumpsSolver<double>::copyFromMumpsRHS(double **d, double *m, int len, int nRHS)
{
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      d[i][j] = m[I];
      I++;
    }
  }
  delete [] m;
}

template<>
void
GenMumpsSolver<DComplex>::copyToMumpsRHS(mumps_double_complex *&m, const GenVector<DComplex> *d, int len, int nRHS)
{
  m = new mumps_double_complex[len*nRHS];
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      m[I].r = d[i][j].real();
      m[I].i = d[i][j].imag();
      I++;
    }
  }
}

template<>
void
GenMumpsSolver<DComplex>::copyFromMumpsRHS(GenVector<DComplex> *d, mumps_double_complex *m, int len, int nRHS)
{
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      DComplex dij(m[I].r,m[I].i);
      d[i][j] = dij;
      I++;
    }
  }
  delete [] m;
}

template<>
void
GenMumpsSolver<double>::copyToMumpsRHS(double *&m, const GenVector<double> *d, int len, int nRHS)
{
  m = new double[len*nRHS];
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      m[I] = d[i][j];
      I++;
    }
  }
}

template<>
void
GenMumpsSolver<double>::copyFromMumpsRHS(GenVector<double> *d, double *m, int len, int nRHS)
{
  int I = 0;
  for(int i=0; i<nRHS; ++i) {
    for(int j=0; j<len; ++j) {
      d[i][j] = m[I];
      I++;
    }
  }
  delete [] m;
}
#endif

template class GenMumpsSolver<double>;
template class GenMumpsSolver<complex<double>>;

