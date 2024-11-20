#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/utilities.hpp>
#include <iostream>
#include <cmath>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

void tran_rvec(double rten[3][3], double rvec[3], bool rescale, int rotvecouttype, double th[3])
/*****************************************************************
 *
 *  Purpose:
 *     Transform the rotation vector
 *
 *  Input:
 *     rten           : rotation tensor
 *     rvec           : Euler rotation vector (not necessarily rescaled)
 *     rotvecouttype  : type of the output rotation vector
 *                      0 -> Euler, 1 -> Complement, 2 -> Linear, 
                        3 -> ReducedEulerRodrigues, 4 -> CayleyGibbsRodrigues,
                        5 -> WienerMilenkovic, 6 -> BauchauTrainelli
 *     rescale        : true -> apply transformation to the "rescaled"
 *                      rotation vector, i.e. mat_to_vec(rten, Psi), ||Psi|| <= pi
 *                      false -> apply transformation to the rotation
 *                      vector rvec which is not necessarily rescaled
 *
 *  Output:
 *     th : transformed rotation vector
 *
 *******************************************************************/
{
#ifdef USE_EIGEN3
  using std::sin;
  using std::tan;
  using std::pow;

  Eigen::Vector3d Psi;
  if(rescale) {
    mat_to_vec(rten, Psi.data());
  }
  else {
    Psi << rvec[0], rvec[1], rvec[2];
  }

  Eigen::Map<Eigen::Vector3d> PsiOut(&th[0]);
  switch(rotvecouttype) {
    case 0 : // Euler
      PsiOut = Psi;
      break;
    case 1 : // Complement
      PsiOut = complement_rot_vec<double>(Psi);
      break;
    case 2 : { // Linear
      double psi = Psi.norm();
      if(psi != 0) PsiOut = sin(psi)*Psi.normalized();
      else PsiOut.setZero();
    } break;
    case 3 : { // ReducedEulerRodrigues
      double psi = Psi.norm();
      if(psi != 0) PsiOut = 2*sin(psi/2)*Psi.normalized();
      else PsiOut.setZero();
    } break;
    case 4 : { // CayleyGibbsRodrigues
      double psi = Psi.norm();
      if(psi != 0) PsiOut = 2*tan(psi/2)*Psi.normalized();
      else PsiOut.setZero();
    } break;
    case 5 : { // WienerMilenkovic
      double psi = Psi.norm();
      if(psi != 0) PsiOut = 4*tan(psi/4)*Psi.normalized();
      else PsiOut.setZero();
    } break;
    case 6 : { // BauchauTrainelli
      double psi = Psi.norm();
      if(psi != 0) PsiOut = pow(6*(psi-sin(psi)),1/3.)*Psi.normalized();
      else PsiOut.setZero();
    } break;
    default :
      std::cerr << " *** Warning: requested rotation vector transformation to "
                << " rotvecouttype = " << rotvecouttype << " is not supported\n";
      for(int i = 0; i < 3; ++i) th[i] = 0;
  }
#else
  if(rotvecouttype == 0) {
    if(rescale) {
      mat_to_vec(rten, th);
    }
    else {
      for(int i = 0; i < 3; ++i) th[i] = rvec[i];
    }
  }
  else {
    std::cerr << " *** Warning: requested rotation vector transformation to "
              << " rotvecouttype = " << rotvecouttype << " is not supported\n";
    for(int i = 0; i < 3; ++i) th[i] = 0;
  }
#endif
}
