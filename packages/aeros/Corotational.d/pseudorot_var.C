#include <cmath>
#include <Corotational.d/utilities.h>

void pseudorot_var(const double rvec[3], double varmat[3][3] )
/******************************************************************
 *
 *  Purpose: Compute the variation of a rotation pseudo vector r with
 *  respect to incremental rotation w where the rotations obey the rule
 *  
 *      R(r+w) = R(w)*R(r)
 *
 *  The variation is given as V = I -Spin/2 +eta*Spin*Spin
 *
 *  Input:
 *     rvec : initial rotation vector
 *
 *  Output:
 *     varmat: matrix of pseudovector variations with repspect to
 *             w = [eps, 0, 0], w = [0, eps, 0] and w = [0, 0, eps]
 *
 *  Coded by: Bjorn Haugen; adjusted for C++ by Teymour Manzouri
 *****************************************************************/
{
   int i, j;
   double eta, th, sthh, cthh;
   double spin[3][3];

   th = sqrt( rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]);

// Check small value of theta : th = 0.05 gives small error between ex.
// and approximate solution

   double th2 = th*th;
   if ( th < 0.05 ) {
      eta = (1.0/12.0) + (1.0/720.0)*th2 + (1.0/30240.0)*(th2*th2);
   } else {
      sthh = sin( 0.5*th );
      cthh = cos( 0.5*th );
      eta  = (sthh - 0.5*th*cthh)/( th2*sthh );
   }

// Compute the spin of the vector
   spin[0][0] = 0.0;
   spin[1][1] = 0.0;
   spin[2][2] = 0.0;

   spin[1][2] = -rvec[0];
   spin[0][2] =  rvec[1];
   spin[0][1] = -rvec[2];

   spin[2][1] =  rvec[0];
   spin[2][0] = -rvec[1];
   spin[1][0] =  rvec[2];

   // Compute first part of var = I - Spin/2

   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) 
        varmat[i][j] = -0.5*spin[i][j];
      varmat[i][i] = 1.0;
   }

// Add the term eta*Spin*Spin
   for( i=0; i<3; i++ ) {
      for( j=0; j<3; j++ ) {
         varmat[i][j] += eta*(spin[i][0]*spin[0][j]
                             +spin[i][1]*spin[1][j]
                             +spin[i][2]*spin[2][j]);
      }
   }

}
