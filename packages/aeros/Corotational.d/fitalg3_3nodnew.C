#include <cmath>
#include <Corotational.d/utilities.h>

double fitalg3_3nodnew( double x0[3][3], double xn[3][3])
/*****************************************************************
 *
 *  Purpose:
 *     Find continum mechanics rotation from x0 coordinates to
 *     xn coordinates.
 *
 *  Method:
 *     Local z-axis is computed normal to the element both for t0 and t0n.
 *         Local x-axis is long side 1-2 of the x0 coordinates for t0
 *     tranformation matrix and xl0 are the coordinates of the elemtent
 *     in this local coordinate system.
 *         Local x-axis is rotated relative to side 1-2 for the deformed
 *     element given by xn coordinates.
 *     xln are the coordinates of the elemtent in this local coordinate system.
 *         The x-axis is rotated so that xl0 are a best fit of the undeformed
 *     element on top of the deformed elemtent given by xln. fitalg
 *     gives which algorithm to use for this rotation.
 *
 *
 *  Input;
 *     fitalg : == 3 x axis rotated so that xl0 is rotated relative
 *                   to xln with the continuum mechanics definition
 *                   of rotations
 *     x0[i][j] : global coordinate component j of node i in undef. conf.
 *     xn[i][j] : global coordinate component j of node i in def. conf.
 *
 *  Output:
 *     t0       : tranformation matrix to local undef coord
 *     t0n      : tranformation matrix to local def coord
 *     xl0[i][j]: local coordinate component j of node i in undef. conf.
 *                also best fit of undef. element on def element.
 *     xln[i][j]: local coordinate component j of node i in def. conf.
 *
 *  Coded by: Bjorn Haugen; Adjusted for C++ by Teymour Manzouri
 *******************************************************************/
{
   int    i, nod;
   double xc0[3], xcn[3], u3, u3p, h3;
   double t0[3][3], t0n[3][3], xl0[3][3], xln[3][3];

// Compute coordinates of centroid for undef. and def. conf.

   for( i=0; i<3; i++ ) {
      xc0[i] = ( x0[0][i] + x0[1][i] + x0[2][i] )/3.0;
      xcn[i] = ( xn[0][i] + xn[1][i] + xn[2][i] )/3.0;
   }

   int nodi = 0;
   int nodj = 1;
   int nodk = 2;

// Compute t0 tranformation matrix with x axis along side j-k
   for( i=0; i<3; i++ ) t0[0][i] = x0[nodj][i] -x0[nodi][i];
   normalize( t0[0] );

   for( i=0; i<3; i++ ) t0[1][i] = x0[nodk][i] - x0[nodi][i];
   crossprod( t0[0], t0[1], t0[2] );
   normalize( t0[2] );

   crossprod( t0[2], t0[0], t0[1] );
   normalize( t0[1] );

// Compute local coordinates of undeformed element
   for( nod=0; nod<3; nod++) {
      for( i=0; i<3; i++ ) {
         xl0[nod][i] = t0[i][0]*(x0[nod][0] -xc0[0])
                      +t0[i][1]*(x0[nod][1] -xc0[1])
                      +t0[i][2]*(x0[nod][2] -xc0[2]);
     }
   }

// Compute t0n tranformation matrix with x axis along side 1-2
   for( i=0; i<3; i++ ) t0n[0][i] = xn[nodj][i] -xn[nodi][i];
   normalize( t0n[0] );

   for( i=0; i<3; i++ ) t0n[1][i] = xn[nodk][i] - xn[nodi][i];
   crossprod( t0n[0], t0n[1], t0n[2] );
   normalize( t0n[2] );

   crossprod( t0n[2], t0n[0], t0n[1] );
   normalize( t0n[1] );

// Compute local coordinates of deformed element
   for( nod=0; nod<3; nod++) {
      for( i=0; i<3; i++ ) {
         xln[nod][i] = t0n[i][0]*(xn[nod][0] -xcn[0])
                      +t0n[i][1]*(xn[nod][1] -xcn[1])
                      +t0n[i][2]*(xn[nod][2] -xcn[2]);
     }
   }

   u3  =  (xln[nodk][0] -xln[nodi][0]) -(xl0[nodk][0] -xl0[nodi][0]);

   u3p = ((xln[nodj][0] -xln[nodi][0]) -(xl0[nodj][0] -xl0[nodi][0]))
        *(xl0[nodk][0]  -xl0[nodi][0])/(xln[nodj][0] -xln[nodi][0]);

   h3  = (xln[nodk][1] - xln[nodi][1]);

   /* talpha = -(u3 - u3p)/h3; */

   // double alphacorr = 0.5*atan( -(u3 - u3p)/h3 );

   return (0.5*atan( -(u3 - u3p)/h3 ));

}
