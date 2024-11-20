#include <Corotational.d/utilities.h>
#include <iostream>

void inc_rottensor( double rvec[3], double rten[3][3] )

/*****************************************************************
 *  Update rotation tensor with an additional rotation 
 *  according to a (spatial) rotation vector
 *
 *  Input:
 *  rvec:   rotation vector
 *
 *  Output:
 *  rten:   updated rotation tensor 
 *****************************************************************/
{

      int    i, j;
      double drten[3][3], nrten[3][3];

   // Compute incremental rotation tensor

      form_rottensor( rvec, drten );

   // Compute updated rotation tensor
      for( i=0; i<3; i++ ) {
         for( j=0; j<3; j++ ) {
            nrten[i][j] = drten[i][0]*rten[0][j]
                         +drten[i][1]*rten[1][j]
                         +drten[i][2]*rten[2][j];
         }
      }

   // Orthonormalize the rotation in order to avoid round off

      orthonorm3( nrten );

   // Store updated rotation tensor in old

      for( i=0; i<3; i++ )
        for( j=0; j<3; j++ ) 
          rten[i][j] = nrten[i][j];

}

void inc_rottensor( double rten[3][3], double rvec[3] )

/*****************************************************************
 *  Update rotation tensor with an additional rotation 
 *  according to a (material) rotation vector
 *
 *  Input:
 *  rvec:   rotation vector
 *
 *  Output:
 *  rten:   updated rotation tensor 
 *****************************************************************/
{

      int    i, j;
      double drten[3][3], nrten[3][3];

   // Compute incremental rotation tensor

      form_rottensor( rvec, drten );

   // Compute updated rotation tensor
      for( i=0; i<3; i++ ) {
         for( j=0; j<3; j++ ) {
            nrten[i][j] = rten[i][0]*drten[0][j]
                         +rten[i][1]*drten[1][j]
                         +rten[i][2]*drten[2][j];
         }
      }

   // Orthonormalize the rotation in order to avoid round off

      orthonorm3( nrten );

   // Store updated rotation tensor in old

      for( i=0; i<3; i++ )
        for( j=0; j<3; j++ ) 
          rten[i][j] = nrten[i][j];

}
