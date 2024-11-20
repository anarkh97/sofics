extern "C" {
  // Lapack 3.2: Symmetric diagonal rescaling
  void _FORTRAN(dsyequb)(const char * uplo, const int* n, const double* a, const int* lda,
                         double* s, double* scond, double* amax, double* work, int* info);

}

namespace Pita {

void equilibrateSym(int matrixSize, const double * data, double relTol, double * scaling, double * scond, double * amax);

} // end namespace Pita
