#ifndef _B2ROTATION_C_
#define _B2ROTATION_C_

/* user defined includes */
#include "b2Util.h"
#include "b2Rotation.h"

/* define M_PI if not allready defined */
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

int b2QEtoA( double E[3], double A[3], double Q[9])
/*
 * programmer : G.Rebel, 3-6-1997, TU-Delft, The Netherlands
 *
 * function   : compute the rotation tensor that maps vector E
 *              onto vector A without a drill, i.e., compute Q in
 *
 *                 A = Q E
 *
 *              in which Q is a rotation tensor that represents a
 *              rotation about a vector that is perpendicular to the
 *              plane spanned by E and A. This Q is given by
 *
 *                 Q = (E . A) I  + (E x A) x I  + 
 *                              3              3
 *
 *                                 1
 *                           + --------- (E x A) (x) (E x A)
 *                             1 + E . A
 *
 * dependency : b2MatrixAddition, b2MatrixProduct, b2SkewMatrix,
 *              b2CrossProduct, b2DyadProduct
 *
 * input      : E = 3-dimensional vector E
 *              A = 3-dimensional vector A
 *
 * output     : on success function computes rotation tensor Q
 *              as specified above
 *
 * return    : status flag is returned :
 *               0 -> no errors
 *              -2 -> division by zero, since : E . A  = -1.0
 */
{
    /* variable definitions */
    int    k;
    double fac, wrk3[3], wrk9[9];
	
    /* compute : E . A */
    b2MatrixProduct( 1, 3, 1, &E[0], &A[0], &fac);
    if ( fac == -1.0)
    {
        /* division by zero detected */
	return -2;
    }
	
    /*************
     *           *
     * compute Q *
     *           *
     *************/
	
    /*
     * compute : (E . A) I
     *                    3
     */
    for ( k = 0; k < 9; k++)
    {
        Q[k] = 0.0;
    }
    Q[0] = fac;
    Q[4] = fac;
    Q[8] = fac;
	
    /*
     * compute : (E . A) I  + (E x A) x I
     *                    3              3
     */
    b2CrossProduct( &E[0], &A[0], &wrk3[0]);
    b2SkewMatrix( &wrk3[0], &wrk9[0]);
    b2MatrixAddition( 3, 3, 1.0, 1.0, &wrk9[0], &Q[0], &Q[0]);
	
    /*
     * compute : Q  = (E . A) I  + (E x A) x I  +
     *                         3              3
     * 
     *                         1
     *                   + ---------- (E x A) (x) (E x A)
     *                     1 + E . A
     */
    fac = 1.0 / ( 1.0 + fac );
    b2DyadProduct( &wrk3[0], &wrk3[0], &wrk9[0]);
    b2MatrixAddition( 3, 3, fac, 1.0, &wrk9[0], &Q[0], &Q[0]);
	
    /******************************
     *                            *
     * computation of Q  finished *
     *                            *
     ******************************/
    
    /* no errors */
    return 0;
}

void b2RotationMatrix( double w[3], double work[9], double q[9])
/*
 * programmer : G.Rebel, 06-10-1995,  TU-Delft, The Netherlands
 *
 * function   : Compute the rotation matrix Q from its corresponding
 *              finite rotation vector w from :
 *
 *                      sin(|w|)           1 - cos(|w|)
 *              Q = I + -------- (w x I) + ------------ (w x (w x I)) 
 *                        |w|                |w|*|w| 
 *
 * precision  : double
 *
 * dependency : b2SkewMatrix, b2MatrixProduct, b2SinNwDivNw, 
 *              b21MCosNwDivNw2, b2MatrixAddition
 *
 * input      : w = finite rotation vector w (dimension of w : 3 x 1)
 *
 * work space : work
 *
 * output     : q = elements of matrix Q stored as follows :
 *
 *                          | q[0] q[1] q[2] |                     
 *                     Q =  | q[3] q[4] q[5] | 
 *                          | q[6] q[7] q[8] |
 */
{
/* define matrix from array */
#define Q(i, j) q[(i)*3+(j)]
    
    /* variable definitions */
    int i;
    double cA, cB;
    
    /* Compute the skew symmetric matrix : w x I */
    b2SkewMatrix( &w[0], &work[0]);
    
    /* Compute the square of skew symmetric matrix w x I */
    b2MatrixProduct( 3, 3, 3, &work[0], &work[0], &q[0]);
    
    /*
     *           sin(|w|)           1 - cos(|w|)
     * compute : -------- (w x I) + ------------ (w x (w x I)) 
     *             |w|                |w|*|w|
     */
    cA = b2SinNwDivNw( &w[0]);
    cB = b21MCosNwDivNw2( &w[0]);           
    b2MatrixAddition( 3, 3, cA, cB, &work[0], &q[0], &q[0]);    
    
    /* Add the identity matrix to Q */
    for ( i = 0; i < 3; i++)
    {
	Q(i, i) += 1.0; 
    }
}

double b2SinNwDivNw( double w[3]) 
/*
 * programmer : G.Rebel, 06-10-1995,  TU-Delft, The Netherlands
 *
 *                         sin(|w|)
 * function   : Compute :  --------
 *                           |w|
 *
 * precision  : double
 *
 * dependency : sin, b2VectorNorm, b2RestrictAngle
 *
 * input      : w = vector w (dimension of w : 3 x 1)
 *                  (|w|=0 is allowed!!!)
 *
 *                                  sin(|w|)
 * output     : function returns :  --------
 *                                    |w|
 */
{   
    double Normw, Normw2, value;
    
    Normw = b2VectorNorm( 3, &w[0]);
    
    if ( Normw > 0.5)
    {
        /* direct computation */
	return ( sin( b2RestrictAngle( Normw)) / Normw);
    }
    else
    {
        /* computation using 16 digit accurate 
	   Taylor expansion. This is done for
	   preventing division by small or zero
	   values of Normw */
	Normw2 = Normw * Normw;
	value = 1.0 - Normw2 / 6.0;
	Normw = Normw2 * Normw2;
	value += Normw / 120.0;
	Normw *= Normw2;
	value -= Normw / 5040.0;
	Normw *= Normw2;
	value += Normw / 362880.0;
	Normw *= Normw2;
	value -= Normw / 39916800.0;
	Normw *= Normw2;
	value += Normw / 6227020800.0;
	return ( value);	 
    }
}

double b21MCosNwDivNw2( double w[3]) 
/*
 * programmer : G.Rebel, 06-10-1995,  TU-Delft, The Netherlands
 *
 *                         1 - cos(|w|)
 * function   : Compute :  ------------
 *                           |w|*|w|
 *
 * precision  : double
 *
 * dependency : cos, b2VectorNorm, b2RestrictAngle
 *
 * input      : w = vector w (dimension of w : 3 x 1)
 *                  (|w|=0 is allowed!!!)
 *
 *                                  1 - cos(|w|)
 * output     : function returns :  ------------
 *                                    |w|*|w|
 */ 
{   
    double Normw, Normw2, value;
    
    Normw = b2VectorNorm( 3, &w[0]);
    
    if ( Normw > 0.5)
    {
        /* direct computation */
	return ( ( 1.0 - cos( b2RestrictAngle( Normw))) / 
	         ( Normw * Normw));
    }
    else
    {
        /* computation using 16 digit accurate 
	   Taylor expansion. This is done for
	   preventing division by small or zero
	   values of Normw */
	Normw2 = Normw * Normw;
	value = 0.5 - Normw2 / 24.0;
	Normw = Normw2 * Normw2;
	value += Normw / 720.0;
	Normw *= Normw2;
	value -= Normw / 40320.0;
	Normw *= Normw2;
	value += Normw / 3628800.0;
	Normw *= Normw2;
	value -= Normw / 479001600.0;
	Normw *= Normw2;
	value += Normw / 87178291200.0;
	return ( value);	 
    }
}

double b2RestrictAngle( double angle) 
/*
 * programmer : G.Rebel, 06-10-1995,  TU-Delft, The Netherlands
 *
 * function   : Add 2*pi*k (k=.., -2, -1, 0, 1, 2, ...) to the
 *              value of angle such that the result lies in
 *              between 0 and 2*pi
 *
 * precision  : double
 *
 * dependency : rint
 *
 * input      : angle = current angle
 *
 * output     : function returns the value of angle that lies
 *              in between 0 and 2*pi 
 *
 */
{
    /* variable definitions */   
    double lower, upper, k;
    
    /* compute the lower and upper bounds of k */
    lower = - angle / ( 2.0 * M_PI);
    upper = ( 2.0 * M_PI - angle) / ( 2.0 * M_PI);
    
    /* determine the value of k */
    k = (double)(rint ( lower));
    if ( k < lower)
    {
	k = (double)(rint ( upper));
    }
    
    /* return value in between 0 and 2*pi */
    return ( angle + 2.0 * M_PI * k);
}

#endif /* EOF, do NOT add anything below this line ! */
