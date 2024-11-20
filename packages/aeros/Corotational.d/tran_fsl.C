#include <Corotational.d/utilities.h>
#include <Math.d/FullSquareMatrix.h>

void tran_fsl(  double force[], FullSquareMatrix &stiff, 
                double t0n[3][3], int num_nodes )
/*****************************************************************
 *
 *  Purpose:
 *     Transform the force and stiffness matrix to global system
 *     or local predefined system if this is defined
 *
 *  Method:
 *            force_global = T'*force
 *            stiff_global = T'*stiff*T
 *
 *  Input:
 *     force : local coordinate internal force vector
 *     stiff : local coordinate stiffness matrix
 *     t0n   : transformation matrix T to global system
 *
 *  Output:
 *     stiff & force : now with respect to global coordinate system
 *
 *  Coded by: Bjorn Haugen, adjusted for C++ by Teymour Manzouri
 *******************************************************************/
{

   int i, j, ii, jj, k, inod, jnod, imat, jmat;
   double lstiff[3][3], stt[3], tf0, tf1, tf2;
   double tmi[3][3], tmj[3][3];

   for( inod=0; inod<num_nodes; inod++ ) {

   // Form transformation matrix tmi

      for( i=0; i<3; i++ ) 
        for( j=0; j<3; j++ ) 
          tmi[i][j] = t0n[i][j];

   // Transform force vector f = T'*f

      for( imat = 0; imat<2; imat ++ ) {
         ii  = inod*6 +imat*3;
         tf0 = tmi[0][0]*force[ii  ]
              +tmi[1][0]*force[ii+1]
              +tmi[2][0]*force[ii+2];
         tf1 = tmi[0][1]*force[ii  ]
              +tmi[1][1]*force[ii+1]
              +tmi[2][1]*force[ii+2];
         tf2 = tmi[0][2]*force[ii  ]
              +tmi[1][2]*force[ii+1]
              +tmi[2][2]*force[ii+2];
         force[ii  ] = tf0;
         force[ii+1] = tf1;
         force[ii+2] = tf2;
      }


      for( jnod=0; jnod<num_nodes; jnod++ ) {

      // Form transformation matrix tmj
         for( i=0; i<3; i++ ) 
           for( j=0; j<3; j++ ) 
             tmj[i][j] = t0n[i][j];

         for( imat = 0; imat<2; imat ++ ) {
            ii = inod*6 +imat*3;
            for( jmat = 0; jmat<2; jmat ++ ) {
               jj = jnod*6 +jmat*3;
               for( j=0; j<3; j++ ) {
                  for( i=0; i<3; i++ ) {
                     for( k=0; k<3; k++ ) 
		     stt[k] = stiff[ii+k][jj  ]*tmj[0][j]
                             +stiff[ii+k][jj+1]*tmj[1][j]
                             +stiff[ii+k][jj+2]*tmj[2][j];
			     
                     lstiff[i][j] = tmi[0][i]*stt[0]
                                   +tmi[1][i]*stt[1]
                                   +tmi[2][i]*stt[2];
                  }
               }
               for( i=0; i<3; i++ ) {
                  for( j=0; j<3; j++ ) 
		    stiff[ii+i][jj+j] = lstiff[i][j];
               }
            }
         }
      }
   }

}
