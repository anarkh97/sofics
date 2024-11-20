#include <cstdio>

typedef double Coord[3];
typedef double Matrix33[3][3];

void compute_transf(Coord& p1, Coord& p2, Coord& p3, Coord& p4, Matrix33& a);
void compute_cofactors(Matrix33& a, Matrix33& b);
void compute_coupling(Coord *p, double elC[30][10]);
int check_if_flat(Coord *p);

extern double p2_dNdxiN[3][10][10];

void compute_coupling(Coord *p, double elC[30][10]) { 

  int i,j,k,m;
  Matrix33 transf_matrix, cofactor_matrix;
  int perm[30] = { 
     1,    4,   10,    7,   13,   25,   28,   19,   22,   16,    2,    5,
    11,    8,   14,   26,   29,   20,   23,   17,    3,    6,   12,    9,
    15,   27,   30,   21,   24,   18 };
  int perm10[10] = { 1,    2,    4,    3,    5,    9,   10,    7,    8,    6 };

/*  if (!check_if_flat(p)) { 
    fprintf(stderr, "Warning: The element has curved faces. The result may be inaccurate.\n");
  }
*/
  compute_transf(p[0], p[1], p[3], p[2], transf_matrix);
  compute_cofactors(transf_matrix, cofactor_matrix);

  // Compute \int dN_i/dx_k N_j
  double dNdxN[3][10][10];
  for(k=0;k<3;k++) for(i=0;i<10;i++) for(j=0;j<10;j++) {
    dNdxN[k][i][j] = 0.0;
    for(m=0;m<3;m++) dNdxN[k][i][j] += cofactor_matrix[k][m] * p2_dNdxiN[m][i][j];
  }

  // Assign C(k*10+i,j) = \int N_i dN_j/dx_k
  for(i=0;i<10;i++) for(j=0;j<10;j++)
    for(k=0;k<3;k++) { 
     elC[perm[k*10+i]-1][perm10[j]-1] = dNdxN[k][i][j];
    }
}


double rel_middle_point_distance_sq(Coord& p1, Coord& p2, Coord &p3) {

  double dist = 0.0, len = 0.0;
  int i;
  for(i=0; i<3 ; i++) {
    double x = (p3[i]-(p1[i]+p2[i])/2.0);
    dist += x * x;
    double y = p2[i] - p1[i];
    len += y * y;
  }
  return dist/len;
}


int check_if_flat(Coord *p) {

  if (rel_middle_point_distance_sq(p[0],p[1],p[4]) > 1e-6) return 0; 
  if (rel_middle_point_distance_sq(p[0],p[2],p[6]) > 1e-6) return 0; 
  if (rel_middle_point_distance_sq(p[0],p[3],p[7]) > 1e-6) return 0; 
  if (rel_middle_point_distance_sq(p[1],p[2],p[5]) > 1e-6) return 0; 
  if (rel_middle_point_distance_sq(p[1],p[3],p[8]) > 1e-6) return 0; 
  if (rel_middle_point_distance_sq(p[2],p[3],p[9]) > 1e-6) return 0; 
  return 1;
}


/*
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
*/
