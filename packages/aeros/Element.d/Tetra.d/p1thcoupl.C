#include <cstdio>
 
typedef double Coord[3];
typedef double Matrix33[3][3];
 
void compute_transf(Coord& p1, Coord& p2, Coord& p3, Coord& p4, Matrix33& a);
void compute_cofactors(Matrix33& a, Matrix33& b);
void compute_coupling(Coord *p, double elC[12][4]);
 
extern double p1_dNdxiN[3][4][4];
 
void compute_coupling(Coord *p, double elC[12][4]) {
 
  int i,j,k,m;
  Matrix33 transf_matrix, cofactor_matrix;
  int perm[30] = {
     1,    4,   10,    7,    2,    5,   11,    8,    3,    6,   12,    9};
  int perm4[4] = { 1,    2,    4,    3 };
 
  compute_transf(p[0], p[1], p[3], p[2], transf_matrix);
  compute_cofactors(transf_matrix, cofactor_matrix);
 
  // Compute \int dN_i/dx_k N_j
  double dNdxN[3][4][4];
  for(k=0;k<3;k++) for(i=0;i<4;i++) for(j=0;j<4;j++) {
    dNdxN[k][i][j] = 0.0;
    for(m=0;m<3;m++) dNdxN[k][i][j] += cofactor_matrix[k][m] * p1_dNdxiN[m][i][j];
  }
 
  // Assign C(k*4+i,j) = \int N_i dN_j/dx_k
  for(i=0;i<4;i++) for(j=0;j<4;j++)
    for(k=0;k<3;k++) {
     elC[perm[k*4+i]-1][perm4[j]-1] = dNdxN[k][i][j];
    }
}


void compute_transf(Coord& p1, Coord& p2, Coord& p3, Coord& p4, Matrix33& a) {

  a[0][0] = p1[0] - p4[0]; 
  a[1][0] = p1[1] - p4[1]; 
  a[2][0] = p1[2] - p4[2]; 

  a[0][1] = p2[0] - p4[0]; 
  a[1][1] = p2[1] - p4[1]; 
  a[2][1] = p2[2] - p4[2]; 

  a[0][2] = p3[0] - p4[0]; 
  a[1][2] = p3[1] - p4[1]; 
  a[2][2] = p3[2] - p4[2]; 
}


void compute_cofactors(Matrix33& a, Matrix33& b) {

  b[0][0] = (1.0) * (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
  b[0][1] = (-1.0) * (a[1][0] * a[2][2] - a[1][2] * a[2][0]);
  b[0][2] = (1.0) * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

  b[1][0] = (-1.0) * (a[0][1] * a[2][2] - a[0][2] * a[2][1]);
  b[1][1] = (1.0) * (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
  b[1][2] = (-1.0) * (a[0][0] * a[2][1] - a[0][1] * a[2][0]);

  b[2][0] = (1.0) * (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
  b[2][1] = (-1.0) * (a[0][0] * a[1][2] - a[0][2] * a[1][0]);
  b[2][2] = (1.0) * (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
}
