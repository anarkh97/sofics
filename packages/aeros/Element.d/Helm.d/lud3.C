
#define _LUD3_C

void lud3(double **A, int n)
{

  int i,j,k;
  double m;

  for (k=0;k<n;k++) {
    for (i=k+1;i<n;i++) {
      m = A[i][k]/A[k][k];
      A[i][k] = m;
      for (j=k+1;j<n;j++) 
        A[i][j] = A[i][j] - m*A[k][j];
    }
  }

}

