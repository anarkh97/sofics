#include <cmath>
#include <Corotational.d/utilities.h>
#include <iostream>

double form_rottensor( double rvec[3], double rten[3][3] )

/*****************************************************************
 *  Compute the rotation tensor based on a rotation vector
 *
 *  Input:
 *  rvec:   rotation vector
 *
 *  Output:
 *  rten:   rotation tensor
 *
 *  Return:
 *
 *  th = sqrt(rvec[0]*rvec[0] +rvec[1]*rvec[1] + rvec[2]*rvec[2])
 *     = rotation angle
 *****************************************************************/
{
      double q[4];

      bool print = false;
      if(print) std::cerr << "rvec = " << rvec[0] << " " << rvec[1] << " " << rvec[2] << std::endl;

      vec_to_quat( rvec, q );
      if(print) std::cerr << "q = " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << std::endl;

      quat_to_mat( q, rten );
      if(print) std::cerr << "rten = " << rten[0][0] << " " << rten[0][1] << " " << rten[0][2] << std::endl
                          << "       " << rten[1][0] << " " << rten[1][1] << " " << rten[2][2] << std::endl
                          << "       " << rten[2][0] << " " << rten[2][1] << " " << rten[2][2] << std::endl;

      return sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]); 
}
