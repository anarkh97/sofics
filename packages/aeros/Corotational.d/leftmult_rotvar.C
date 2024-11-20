#include <Corotational.d/utilities.h>
#include <Math.d/FullSquareMatrix.h>

void
leftmult_rotvar(int num_nodes, int itrans, double rotvar[][3][3], 
                FullSquareMatrix &stiff)

/*****************************************************************
 *
 *  Purpose:
 *     Multiply | I  0 |*stiff or | I  0 |*stiff 
 *              | 0  V'|          | 0  V |
 *     where V is rotvar or variation of pseudovector rotation 
 *     with respect to instantanious rotation
 *
 *  Input:
 *     num_nodes : number of nodes
 *     itrans    : = 1 use transpose V
 *     rotvar    : rotvar[i][j][k] is component j k of rotation
 *                 gradient matrix for node i
 *     stiff     : stiffness
 *  Output:
 *     stiff     : stiff updated with the rotation gradient product
 *
 *  Coded by: Bjorn Haugen
 *******************************************************************/
{

   int i, j, inod, ndofs, ir;
   double scr[3];

// Transform groups of 3 degrees of freedom

   ndofs = 6*num_nodes;

   if ( itrans == 1 ) {                       /* Use Transpose rotvar */
      for( inod=0; inod<num_nodes; inod++ ) {
         ir = 6*inod + 3; 
         for( j=0; j<ndofs; j++ ) {
            for( i=0; i<3; i++ ) {
                  scr[i] =  rotvar[inod][0][i]*stiff[ir +0][j]
                           +rotvar[inod][1][i]*stiff[ir +1][j]
                           +rotvar[inod][2][i]*stiff[ir +2][j];
            }
	    stiff[ir+0][j] = scr[0];
	    stiff[ir+1][j] = scr[1];
            stiff[ir+2][j] = scr[2];
         }
      }
   }

   else {                                     /*  Use rotvar */
      for( inod=0; inod<num_nodes; inod++ ) {
         ir = 6*inod + 3;
         for( j=0; j<ndofs; j++ ) {
            for( i=0; i<3; i++ ) {
                  scr[i] =  rotvar[inod][i][0]*stiff[ir +0][j]
                           +rotvar[inod][i][1]*stiff[ir +1][j]
                           +rotvar[inod][i][2]*stiff[ir +2][j];
            }
            stiff[ir+0][j] = scr[0];
            stiff[ir+1][j] = scr[1];
            stiff[ir+2][j] = scr[2];
         }
      }
   }

}
