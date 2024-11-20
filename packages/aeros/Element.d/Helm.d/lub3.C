
#define _LUB3_C

void lub3(double **A, int n, double *b)
{

  int i,j, k;
  double sum;

  // Ly=b
  for (i=0;i<n;i++) {
      sum = b[i];
      for (j=0;j<i;j++) 
        sum -= A[i][j]*b[j];
      b[i] = sum;
  }

  // Ux=y
  for (i=n-1;i>=0;i--) {
    sum = b[i];
    for (j=i+1;j<n;j++) 
      sum -= A[i][j]*b[j];
    b[i] = sum/A[i][i];
  }

}
