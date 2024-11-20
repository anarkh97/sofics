#include <Corotational.d/utilities.h>
#include <iostream>

void orthonorm3( double rten[3][3] )
/*****************************************************************
 *  Orthonormalize the rotation tensor using the Euler
 *  parameter definitions
 *
 *  Input:
 *  rten:    3x3 rotation matrix
 *
 *  Output:
 *  rten:    3x3 rotation matrix orthonormalized
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{

   double q[4];

   mat_to_quat( rten, q );

   quat_to_mat( q, rten );
   
}
