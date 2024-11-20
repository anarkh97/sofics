#include <cmath>
#include <Corotational.d/utilities.h>
#include <iostream>

#ifdef USE_EIGEN3
#define USE_EIGEN_QUAT_TO_MAT
#ifdef USE_EIGEN_QUAT_TO_MAT
#include <Eigen/Core>
#include <Eigen/Geometry>
#endif
#endif

void quat_to_mat( double q[4], double rten[3][3] )
/*****************************************************************
 *  Compute the rotation matrix from a quaternion representation
 *  of the rotation
 *
 *  Input:
 *  q:    rotation quaternion ( Euler parameters )
 *
 *  Output:
 *  rten:    rotation tensor
 *
 *  Return:
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
#ifdef USE_EIGEN_QUAT_TO_MAT
   Eigen::Quaterniond qq;
   qq.w() = q[0];
   qq.x() = q[1];
   qq.y() = q[2];
   qq.z() = q[3];
   qq.normalize();
   Eigen::Matrix3d rm = qq.toRotationMatrix();
   for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) rten[i][j] = rm(i,j);
#else
   // Compute norm of q
   double norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
   if(norm == 0) {
     std::cerr << "WARNING: Skipping quaternion normalize( null )\n";
   }
   else {
     // Normalize q
     q[0] /= norm;
     q[1] /= norm;
     q[2] /= norm;
     q[3] /= norm;
   }

   rten[0][0] = 2.0*(q[1]*q[1] + q[0]*q[0]) - 1.0;
   rten[1][1] = 2.0*(q[2]*q[2] + q[0]*q[0]) - 1.0;
   rten[2][2] = 2.0*(q[3]*q[3] + q[0]*q[0]) - 1.0;

   rten[0][1] = 2.0*(q[1]*q[2] - q[3]*q[0]);
   rten[0][2] = 2.0*(q[1]*q[3] + q[2]*q[0]);
   rten[1][2] = 2.0*(q[2]*q[3] - q[1]*q[0]);

   rten[1][0] = 2.0*(q[2]*q[1] + q[3]*q[0]);
   rten[2][0] = 2.0*(q[3]*q[1] - q[2]*q[0]);
   rten[2][1] = 2.0*(q[3]*q[2] + q[1]*q[0]);

/*
   rten[0][0] = q[1]*q[1] - q[2]*q[2] - q[3]*q[3] + q[0]*q[0];
   rten[1][0] = 2.0*(q[1]*q[2] + q[3]*q[0]);
   rten[2][0] = 2.0*(q[1]*q[3] - q[2]*q[0]);
   rten[0][1] = 2.0*(q[1]*q[2] - q[3]*q[0]);
   rten[1][1] = -q[1]*q[1] + q[2]*q[2] - q[3]*q[3] + q[0]*q[0];
   rten[2][1] = 2.0*(q[2]*q[3] + q[1]*q[0]);
   rten[0][2] = 2.0*(q[1]*q[3] + q[2]*q[0]);
   rten[1][2] = 2.0*(q[2]*q[3] - q[1]*q[0]);
   rten[2][2] = -q[1]*q[1] - q[2]*q[2] + q[3]*q[3] + q[0]*q[0];
*/
#endif
}

void dquat_to_mat( double q[4], double dq[4], double rten[3][3], 
                   double drten[3][3] )
/*****************************************************************
 *  Compute the rotation matrix from a quaternion representation
 *  of the rotation
 *
 *  Input:
 *  q:    rotation quaternion ( Euler parameters )
 *
 *  Output:
 *  rten:    rotation tensor
 *
 *  Return:
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
   // Compute norm of q
   double norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
   
   // if(norm == 0.0) fprintf(stderr,"\n\n ERROR IN QUAT_TO_MAT \n\n");
  
   rten[0][0] = 2.0*(q[1]*q[1] + q[0]*q[0]) - 1.0;
  drten[0][0] = 4.0*(q[1]*dq[1]+ q[0]*dq[0]);
   rten[1][1] = 2.0*(q[2]*q[2] + q[0]*q[0]) - 1.0;
  drten[1][1] = 4.0*(q[2]*dq[2] + q[0]*dq[0]);
   rten[2][2] = 2.0*(q[3]*q[3] + q[0]*q[0]) - 1.0;
  drten[2][2] = 4.0*(q[3]*dq[3] + q[0]*dq[0]);

   rten[0][1] = 2.0*(q[1]*q[2] - q[3]*q[0]);
  drten[0][1] = 2.0*(dq[1]*q[2]+ q[1]*dq[2] - dq[3]*q[0] - q[3]*dq[0]);
   rten[0][2] = 2.0*(q[1]*q[3] + q[2]*q[0]);
  drten[0][2] = 2.0*(dq[1]*q[3]+ q[1]*dq[3] + dq[2]*q[0] + q[2]*dq[0]);
   rten[1][2] = 2.0*(q[2]*q[3] - q[1]*q[0]);
  drten[1][2] = 2.0*(dq[2]*q[3]+ q[2]*dq[3] - dq[1]*q[0] - q[1]*dq[0]);

   rten[1][0] = 2.0*(q[2]*q[1] + q[3]*q[0]);
  drten[1][0] = 2.0*(dq[2]*q[1] + q[2]*dq[1] + dq[3]*q[0] + q[3]*dq[0]);
   rten[2][0] = 2.0*(q[3]*q[1] - q[2]*q[0]);
  drten[2][0] = 2.0*(dq[3]*q[1] + q[3]*dq[1] - dq[2]*q[0] - q[2]*dq[0]);
   rten[2][1] = 2.0*(q[3]*q[2] + q[1]*q[0]);
  drten[2][1] = 2.0*(dq[3]*q[2] + q[3]*dq[2] + dq[1]*q[0] + q[1]*dq[0]);
}
