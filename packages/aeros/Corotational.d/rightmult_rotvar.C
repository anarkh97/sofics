#include <Corotational.d/utilities.h>
#include <Math.d/FullSquareMatrix.h>

void rightmult_rotvar(int  num_nodes, int itran, double  rotvar[][3][3]
                                             ,FullSquareMatrix  &stiff )
/*****************************************************************
 *
 *  Purpose:
 *     Multiply stiff*| I  0 | or stiff*| I  0 |
 *                    | 0  V |          | 0  V'|
 *     where V is rotvar  or variation of pseudovector rotation 
 *     with respect to instantanious rotation
 *
 *  Input:
 *     num_nodes : number of nodes
 *     itran     : =1 user transpose rotvar
 *     rotvar    : rotvar[i][j][k] is compinent j k of rotation
 *                 gradient matrix for node i
 *     stiff     : stiffness
 *  Output:
 *     stiff     : stiff updated with the rotation gradient product
 *
 *  Coded by: Bjorn Haugen
 *******************************************************************/
{

   int i, j, ndofs, jnod, jr;
   double scr[3];

// Transform groups of 3 degrees of freedom

   ndofs = 6*num_nodes;

   if( itran == 1 ) {                       /*  Use transpose rotvar */
      
      for( jnod=0; jnod<num_nodes; jnod++ ) {
         jr = 6*jnod +3;
         for( i=0; i<ndofs; i++ ) {
            for( j=0; j<3; j++ ) {
                  scr[j] =  stiff[i][jr +0]*rotvar[jnod][j][0]
                           +stiff[i][jr +1]*rotvar[jnod][j][1]
                           +stiff[i][jr +2]*rotvar[jnod][j][2];
            }
            stiff[i][jr+0] = scr[0];
            stiff[i][jr+1] = scr[1];
            stiff[i][jr+2] = scr[2];
         }
      }

   }

   else {                                     /*  Use rotvar */
      for( jnod=0; jnod<num_nodes; jnod++ ) {
         jr = 6*jnod +3;
         for( i=0; i<ndofs; i++ ) {
            for( j=0; j<3; j++ ) {
                  scr[j] =  stiff[i][jr +0]*rotvar[jnod][0][j]
                           +stiff[i][jr +1]*rotvar[jnod][1][j]
                           +stiff[i][jr +2]*rotvar[jnod][2][j];
            }
            stiff[i][jr+0] = scr[0];
            stiff[i][jr+1] = scr[1];
            stiff[i][jr+2] = scr[2];
         }
      }
   }

}
