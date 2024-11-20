extern "C" {
  // Lapack 3.2: Cholesky factorization with complete diagonal pivoting
  void _FORTRAN(dpstrf)(const char* uplo, const int* n, double* a, const int* lda, int* piv,
                        int* rank, const double* tol, double* work, int* info);

}
