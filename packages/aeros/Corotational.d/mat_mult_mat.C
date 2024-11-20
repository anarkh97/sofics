
void
mat_mult_mat(const double A[3][3], const double B[3][3], double C[3][3], int transflag)
{
 int i,j;
 if(transflag == 0) {
   // C = A*B
   for(i=0; i<3; ++i)
     for(j=0; j<3; ++j)
       C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];
 }
 else if(transflag == 1) {
   // C = A^t*B
   for(i=0; i<3; ++i)
     for(j=0; j<3; ++j)
       C[i][j] = A[0][i]*B[0][j] + A[1][i]*B[1][j] + A[2][i]*B[2][j];

 } 
 else if(transflag == 2) {
   // C = A*B^t
   for(i=0; i<3; ++i)
     for(j=0; j<3; ++j)
       C[i][j] = A[i][0]*B[j][0] + A[i][1]*B[j][1] + A[i][2]*B[j][2];

 }
}

void
mat_mult_vec(double A[3][3], double b[3], double c[3], int transflag)
{
  if(transflag == 0) {
    // c = A*b
    for(int i=0; i<3; ++i)
      c[i] = A[i][0]*b[0]+A[i][1]*b[1]+A[i][2]*b[2];
  }
  else if(transflag == 1) {
    // c = A^t*b
    for(int i=0; i<3; ++i)
      c[i] = A[0][i]*b[0]+A[1][i]*b[1]+A[2][i]*b[2];
  }
}
