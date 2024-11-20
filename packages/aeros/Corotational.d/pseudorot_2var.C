#include <cmath>
#include <Corotational.d/utilities.h>

void pseudorot_2var(const double r[3], const double f[3], double scndvar[3][3])
/******************************************************************
 *
 *  Purpose: Compute the variation of a contraction between
 *      Rmat'*f  where f is constant vector and Rmat is the
 *      rotation variation of a pseudovector r with respect
 *      to instantanious rotations
 *      
 *  Input:
 *     r : initial rotation vector ( vector to be varied );
 *     f : constant vector that is contracted with the varied
 *         rotation varitation matrix
 *
 *  Output:
 *     scndvar: matrix of the variation of the contraction with respect to
 *             w = [eps, 0, 0], w = [0, eps, 0] and w = [0, 0, eps]
 *
 *  Coded by: Bjorn Haugen; adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int i, j;
   double eta, nu, th, sth, sthh, cthh, fac;
   double mat[3][3], sr[3][3], fr[3][3], spsq[3][3], fstvar[3][3];

   th = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );

/* Check small value of theta : th = 0.01 gives small error between ex.
 * and approximate solution for eta
 */
 
   if ( th < 0.05 ) {
      eta = (1.0/12.0)  + (th*th)/720.0  + ((th*th)*(th*th))/30240.0;
      nu  = (1.0/360.0) + (th*th)/7560.0 + ((th*th)*(th*th))/201600.0;
   } else {
      sth  = sin( th );
      sthh = sin( 0.5*th );
      cthh = cos( 0.5*th );
      eta  = (sthh - 0.5*th*cthh)/( (th*th)*sthh );
      nu   = (th*(th +sth) - 8.0*sthh*sthh )/(4.0*((th*th)*(th*th))*sthh*sthh);
   }

/* Compute the spin of the rotation vectors */
   sr[0][0] = 0.0;
   sr[1][1] = 0.0;
   sr[2][2] = 0.0;

   sr[1][2] = -r[0];
   sr[0][2] =  r[1];
   sr[0][1] = -r[2];

   sr[2][1] =  r[0];
   sr[2][0] = -r[1];
   sr[1][0] =  r[2];


// Compute first part of mat =  -SpinF/2
   mat[0][0] = 0.0;
   mat[1][1] = 0.0;
   mat[2][2] = 0.0;

   mat[1][2] =  0.5*f[0];
   mat[0][2] = -0.5*f[1];
   mat[0][1] =  0.5*f[2];

   mat[2][1] = -0.5*f[0];
   mat[2][0] =  0.5*f[1];
   mat[1][0] = -0.5*f[2];

/* Add the terms  mat += eta*( (r'*f)I +r*f' -2f*r') */
   fac = r[0]*f[0] + r[1]*f[1] + r[2]*f[2];

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         spsq[i][j] = sr[i][0]*sr[0][j]
                     +sr[i][1]*sr[1][j]
                     +sr[i][2]*sr[2][j];
         fr[i][j] = f[i]*r[j];
         mat[i][j] += eta*( r[i]*f[j] - 2.0*f[i]*r[j] );
      }
      mat[i][i] += eta*fac;
   }

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         mat[i][j] += nu*( spsq[i][0]*fr[0][j]
                          +spsq[i][1]*fr[1][j]
                          +spsq[i][2]*fr[2][j] );
      }
   }

// Multiply with first variation
   pseudorot_var( r, fstvar );

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         scndvar[i][j] =  mat[i][0]*fstvar[0][j]
                         +mat[i][1]*fstvar[1][j]
                         +mat[i][2]*fstvar[2][j];
      }
   }

}
