#include <Corotational.d/utilities.h>
#include <cmath>

void vec_to_mat( double rvec[3], double rten[3][3] )

/*****************************************************************
 *  Compute the rotation tensor from a rotation matrix
 *
 *  Input:
 *  rten:    rotation tensor
 *
 *  Output:
 *  rvec:    rotation vector
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{

   double q[4];

   vec_to_quat( rvec, q );

   quat_to_mat( q, rten );

}

void vec_to_mat_new( double rvec[3], double rten[3][3] )
/*****************************************************************
 *  Compute the rotation tensor from a rotation matrix
 *
 *  Input:
 *  rten:    rotation tensor
 *
 *  Output:
 *  rvec:    rotation vector
 *
 *  Coded by Teymour Manzouri
 *****************************************************************/
{
   double th, sth, cth,alpha,beta;

   th = sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]);

   sth = sin( th );
   cth = cos( th );

   if ( th < 1.0E-15 ) {
      alpha = 1.0;
      beta  = 0.5;
   } else {
      alpha = sth/th;
      beta  = (1.0 - cth)/th/th;
   }

   rten[0][0] = 1.0 - beta*(rvec[1]*rvec[1]+rvec[2]*rvec[2]);
   rten[1][1] = 1.0 - beta*(rvec[0]*rvec[0]+rvec[2]*rvec[2]);
   rten[2][2] = 1.0 - beta*(rvec[1]*rvec[1]+rvec[0]*rvec[0]);

   rten[0][1] = -alpha*rvec[2] + beta*rvec[0]*rvec[1] ;
   rten[0][2] =  alpha*rvec[1] + beta*rvec[0]*rvec[2] ;
   rten[1][2] = -alpha*rvec[0] + beta*rvec[1]*rvec[2] ;

   rten[1][0] =  alpha*rvec[2] + beta*rvec[0]*rvec[1] ;
   rten[2][0] = -alpha*rvec[1] + beta*rvec[0]*rvec[2] ;
   rten[2][1] =  alpha*rvec[0] + beta*rvec[1]*rvec[2] ;

}

void dvec_to_mat( double rvec[3], double drvec[3], double rten[3][3],
                  double drten[3][3] )

/*****************************************************************
 *  Compute the rotation tensor from a rotation matrix
 *
 *  Input:
 *  rten:    rotation tensor
 * drten:    derivative of rotation tensor
 *
 *  Output:
 *  rvec:    rotation vector
 *  drvec:  derivative of rotation vector
 *  Coded by Bjorn Haugen
 *****************************************************************/
{

   double q[4];
   double dq[4];
   
   dvec_to_quat( rvec, drvec, q, dq );

   dquat_to_mat( q, dq, rten, drten );

}
