#include <cmath>
#include <Corotational.d/utilities.h>
#include <iostream>
#include <Utils.d/linkfc.h>

// PJSA: 5/2/2011 I observe much better results with the eigen version
// which uses algorithm from "Quaternion Calculus and Fast Animation",
// Ken Shoemake, 1987 SIGGRAPH course notes
#ifdef USE_EIGEN3
#define USE_EIGEN_MAT_TO_QUAT
#ifdef USE_EIGEN_MAT_TO_QUAT
#include <Eigen/Core>
#include <Eigen/Geometry>
#endif
#endif

using std::cerr;

extern "C"      {
   void _FORTRAN(dsyev)(const char &, const char &, const int &,
                        double *, const int &, double *, double *,
                        const int &, int &);
}


void mat_to_quat( double rten[3][3], double q[4] )
/*****************************************************************
 *  Compute the quaternion representation of a rotation tensor
 *  by means of equations from Richard A. Spurrier comments on
 *  "Singularity free Extraction of a Quaternion from a 
 *   Direction-Cosine Matrix" in J. of Rockets and Spacecraft 1978.
 *
 *  Input:
 *  rten:    rotation tensor
 *
 *  Output:
 *  q:    rotation quaternion ( Euler parameters )
 *
 *  Return:
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
#ifdef USE_EIGEN_MAT_TO_QUAT
      // using this version causes reparameterization of rotations at -2/3*pi,4/3*pi
      Eigen::Matrix3d rm;
      for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) rm(i,j) = rten[i][j];
      Eigen::Quaterniond qq(rm);
      q[0] = qq.w();
      q[1] = qq.x();
      q[2] = qq.y();
      q[3] = qq.z();
#else
      // using this routine causes reparameterization of rotations at -1/2*pi,3/2*pi
      static int p[5] = {0,1,2,0,1};
      int    imax, i, j, k;
      double trace;

      trace = rten[0][0] +rten[1][1] +rten[2][2];

      imax = 0;
      if ( rten[1][1] > rten[imax][imax] ) imax = 1;
      if ( rten[2][2] > rten[imax][imax] ) imax = 2;

      if ( trace > rten[imax][imax] ) {
         q[0] = std::sqrt( 1.0 +trace )/2.0;
         q[1] = ( rten[2][1] -rten[1][2] )/(4.0*q[0]);
         q[2] = ( rten[0][2] -rten[2][0] )/(4.0*q[0]);
         q[3] = ( rten[1][0] -rten[0][1] )/(4.0*q[0]);
      }
      else {
         i = p[imax];
         j = p[imax+1];
         k = p[imax+2];
         q[i+1] = std::sqrt( rten[i][i]/2.0 +(1.0 -trace)/4.0 );
         q[  0] = ( rten[k][j] -rten[j][k] )/(4.0*q[i+1]);
         q[j+1] = ( rten[j][i] +rten[i][j] )/(4.0*q[i+1]);
         q[k+1] = ( rten[k][i] +rten[i][k] )/(4.0*q[i+1]);
      }

/*
      // PJSA 2-12-09 this is a variant due to Markley (Journal of Guidance and Control Vol 31 No 2 pp 440-442, 2008)
      // using this version causes a reparameterization of the rotation at -pi,+pi
      if ( trace > rten[imax][imax] ) {
         q[0] = ( 1.0 +trace );
         q[1] = ( rten[2][1] -rten[1][2] );
         q[2] = ( rten[0][2] -rten[2][0] );
         q[3] = ( rten[1][0] -rten[0][1] );
      }
      else {
         i = p[imax];
         j = p[imax+1];
         k = p[imax+2];
         q[i+1] = 1.0 + rten[i][i] - rten[j][j] - rten[k][k];
         q[  0] = ( rten[k][j] -rten[j][k] );
         q[j+1] = ( rten[j][i] +rten[i][j] );
         q[k+1] = ( rten[k][i] +rten[i][k] );
      }
      double norm = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
      for(int l=0; l<4; ++l) q[l] /= norm;
*/

/*
      // PJSA 3-12-09 method due to Bar-Itzhack (Journal of Guidance and Control Vol 23 No 6 pp 1085-1087, 2000)
      // which always gives the quaternion that corresponds to the orthogonal matrix closest to rten
      // this version causes a reparameterization of the rotation at -3/2*pi, 1/2*pi
      double A[16] = {
        (rten[0][0]-rten[1][1]-rten[2][2])/3.0, (rten[0][1]+rten[1][0])/3.0,            (rten[0][2]+rten[2][0])/3.0,            (rten[2][1]-rten[1][2])/3.0,
        0.0,                                    (rten[1][1]-rten[0][0]-rten[2][2])/3.0, (rten[1][2]+rten[2][1])/3.0,            (rten[0][2]-rten[2][0])/3.0,
        0.0,                                    0.0,                                    (rten[2][2]-rten[0][0]-rten[1][1])/3.0, (rten[1][0]-rten[0][1])/3.0,
        0.0,                                    0.0,                                    0.0,                                    (rten[0][0]+rten[1][1]+rten[2][2])/3.0 };

      int info;
      double w[4], work[11];
      _FORTRAN(dsyev)('V', 'L', 4, A, 4, w, work, 11, info);

      q[0] = A[15];
      q[1] = A[12];
      q[2] = A[13];
      q[3] = A[14];

      //double norm = std::sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
      //for(int l=0; l<4; ++l) q[l] /= norm;
*/
#endif

  // always return the positive quaternion.
  // see: "INTERPOLATION OF ROTATIONAL VARIABLES IN NONLINEAR DYNAMICS OF 3D BEAMS" 
  // by Jelenic and Crisfield, Int. J. Numer. Meth. Engng. 43, 1193–1222 (1998)
  // section 3.5.2. for discussion of uniqueness of the extracted rotational pseudovector.
  // also appendix B1 of "A THREE-DIMENSIONAL FINITE-STRAIN ROD MODEL
  // PART II: COMPUTATIONAL ASPECTS" by Simo and Vu-Quoc, Computer Methods in Applied 
  // Mechanics and Engineering, Volume 58, Issue 1, October 1986
  if(q[0] < 0) {
    //cerr << "negative quaternion in mat_to_quat\n";
    for(int l=0; l<4; ++l) q[l] = -q[l]; // this amounts to reparameterizing at -pi,pi
  }

}
