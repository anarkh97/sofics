#ifndef _B2UTIL_H_
#define _B2UTIL_H_

/*
 * include file for declaration of tensor utility functions
 * 
 */
 
/* compiler includes */
#include <cmath>
 
/* function declarations */
inline void   b2MatrixAddition      ( int k, int l, double cA, double cB,  
                                      double *A, double *B, double *C);
inline void   b2MatrixProduct       ( int k, int l, int m, 
                                      double *A, double *B, double *C);
inline void   b2SkewMatrix          ( double w[3], double wskew[9]);
inline void   b2CrossProduct        ( double a[3], double b[3], 
                                      double c[3]);
inline void   b2DyadProduct         ( double a[3], double b[3], 
                                      double c[9]);
inline double b2VectorNorm          ( int n, double *w);

inline void b2MatrixAddition( int k, int l, double cA, double cB,  
                              double *A, double *B, double *C)
/*
 * programmer : G.Rebel, 03-10-1995,  TU-Delft, The Netherlands
 *
 * function   : Compute matrix C from : C = cA * A + cB * B, 
 *              in which cA, cB are scalars and A, B are 
 *              matrices
 *
 * precision  : double
 *
 * dependency : none
 *
 * input      : k  = number of rows of matrices A, B and C
 *              l  = number of columns of matrices A, B and C
 *              cA = constant with which matrix A is multiplied
 *              cB = constant with which matrix B is multiplied
 *              A  = matrix A (dimension of A : k x l)
 *              B  = matrix B (dimension of B : k x l) 
 *
 * output     : C  = cA * A + cB * B, (dimension of C : k x l)
 *
 */
{
    int i, j;
    
    for ( i = 0; i < k; i++)
    {
	for ( j = 0; j < l; j++)
	{
	    *(C+(i*l+j)) = cA * *(A+(i*l+j)) + cB * *(B+(i*l+j));
	}
    }
}

inline void b2MatrixProduct( int k, int l, int m, 
                             double *A, double *B, double *C) 
/*
 * programmer : G.Rebel, 02-10-1995,  TU-Delft, The Netherlands
 *
 * function   : Compute the matrix product C of the matrices A
 *              and B as follows : C = A B
 *
 * precision  : double
 *
 * dependency : none
 *
 * input      : k = number of rows of matrix A
 *              l = number of columns of matrix A
 *                   (= number of rows of matrix B)
 *              m = number of columns of matrix B
 *              A = matrix A (dimension of A : k x l)
 *              B = matrix B (dimension of B : l x m) 
 *
 * output     : C = matrix product : A B
 *                  (dimension of C : k x m)
 *
 */
{
    int i, j, n;
    
    for ( i = 0; i < k; i++)
    {
	for ( j = 0; j < m; j++)
	{
	    C[i*m+j] = 0.0;
	    for ( n = 0; n < l; n++)
	    {
		*(C+(i*m+j)) += *(A+(i*l+n)) * *(B+(n*m+j));
	    }
	}
    }
}

inline void b2SkewMatrix( double w[3], double wskew[9])
/*
 * programmer : G.Rebel, 06-10-1995,  TU-Delft, The Netherlands
 *
 * function   : Compute skew symmetric matrix W from its axial
 *              vector w as follows :
 * 
 *                                    |  0  -w3  w2 |
 *                        W = w x I = |  w3  0  -w1 |
 *                                    | -w2  w1  0  |
 *
 * precision  : double
 *
 * dependency : none
 *
 * input      : w     = axial vector of skew symmetric matrix W
 *
 * output     : wskew = elements of matrix W stored as follows :
 *
 *                               | wskew[0] wskew[1] wskew[2] |                     
 *                          W =  | wskew[3] wskew[4] wskew[5] | 
 *                               | wskew[6] wskew[7] wskew[8] |
 *
 */
{
    wskew[0] = 0.0;
    wskew[1] = -w[2];
    wskew[2] = w[1];
    
    wskew[3] = w[2];
    wskew[4] = 0.0;
    wskew[5] = -w[0];
    
    wskew[6] = -w[1];
    wskew[7] = w[0];
    wskew[8] = 0.0;
}
   
inline void b2CrossProduct( double a[3], double b[3], double c[3])
/*
 * programmer : G.Rebel, 31-10-1995, TU-Delft, The Netherlands
 *
 * function   : Compute the cross-product c of the vectors a
 *              and b as follows : c = a x b
 *
 * precision  : double
 *
 * dependency : none
 *
 * input      : a = vector a
 *              b = vector b 
 *
 * output     : c = cross product : c = a x b
 *
 */
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

inline void b2DyadProduct( double a[3], double b[3], double c[9])
/*
 * programmer : G.Rebel, 31-10-1995, TU-Delft, The Netherlands
 *
 * function   : Compute the dyad-product c of the vectors a
 *              and b ; c = a (x) b
 *
 * precision  : double
 *
 * dependency : none
 *
 * input      : a = vector a
 *              b = vector b 
 *
 * output     :                                  | a b   a b   a b  |
 *                                               |  1 1   1 2   1 3 |
 *                                               |                  |
 *              c = dyad product : c = a (x) b = | a b   a b   a b  |
 *                                               |  2 1   2 2   2 3 |
 *                                               |                  |
 *                                               | a b   a b   a b  |
 *                                               |  3 1   3 2   3 3 |
 *
 *                  stored as : c[(i-1)*3+(j-1)] = c
 *                                                  ij
 */
{ 
    /* variable definitions */
    int i, j;
 
    /* compute dyad-product c = a (x) b */
    for ( i = 1; i < 4; i++)
    {
	for ( j = 1; j < 4; j++)
	{
	    *(c+((i-1)*3+j-1)) = *(a+(i-1)) * *(b+(j-1));
	}
    }
}

inline double b2VectorNorm( int n, double *w) 
/*
 * programmer : G.Rebel, 02-10-1995,  TU-Delft, The Netherlands
 *
 * function   : Compute the length of vector w. The length of w
 *              is defined as its Euclidean norm.
 *
 * precision  : double
 *
 * dependency : sqrt
 *
 * input      : n = number of elements of vector w
 *              w = vector w (dimension of w : n x 1) 
 *
 * output     : function returns the length of w
 *
 */
{
    int i;
    double length=0.0;
    
    for ( i = 0; i < n; i++)
    {
        length += w[i] * w[i];
    }
    return ( sqrt( length));
}

#endif /* EOF, do NOT add anything below this line ! */
