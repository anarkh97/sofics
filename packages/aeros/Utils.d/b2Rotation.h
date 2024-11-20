#ifndef _B2ROTATION_H_
#define _B2ROTATION_H_

/*
 * include file for declaration of the functions used
 * to construct rotation tensors
 * 
 */
 
/* function declarations */
int    b2QEtoA         ( double E[3], double A[3], double Q[9]);
void   b2RotationMatrix( double w[3], double work[9], double q[9]);
double b2SinNwDivNw    ( double w[3]);
double b21MCosNwDivNw2 ( double w[3]);
double b2RestrictAngle ( double angle);

#endif /* EOF, do NOT add anything below this line ! */
