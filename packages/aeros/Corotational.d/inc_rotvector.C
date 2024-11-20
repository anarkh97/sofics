#include <Element.d/Function.d/utilities.hpp>

void inc_rotvector( double rvec[3], double rten[3][3] )

/*****************************************************************
 *  Update rotation vector without rescaling
 *
 *  Input/Output:
 *  rvec:   rotation vector
 *
 *  Input:
 *  rten:   updated rotation tensor 
 *****************************************************************/
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Vector3d> Psi_n(&rvec[0]);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&rten[0][0]);

  // first, compute the updated rescaled rotation vector
  Eigen::Vector3d Psi;
  mat_to_vec<double>(R, Psi);
 
  if(Psi_n.norm() > M_PI) {
    Psi_n = complement_rot_vec(unscale_rotvec(Psi, complement_rot_vec<double>(Psi_n)));
  }
  else {
    Psi_n = unscale_rotvec<double>(Psi, Psi_n);
  }
#endif
}

