#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/utilities.hpp>
#include <iostream>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

void tran_accel(double rten[3][3], double rvec[3], double V[3], double A[3],
                int angularintype, int angularouttype, bool rescalein,
                bool rescaleout, double a[3])
/*****************************************************************
 *
 *  Purpose:
 *     Transform the angular acceleration
 *
 *  Input:
 *     rten           : rotation tensor
 *     rvec           : rotation vector
 *     V              : input angular velocity
 *     A              : input angular acceleration
 *     angularintype  : type of the input angular velocity/accel.
 *                      0 -> spatial, 1 -> convected, 2 -> total
 *     angularouttype : type of the transformed angular accel. 
 *                      0 -> spatial, 1 -> convected, 2 -> total
 *     rescalein      : true -> if angularintype == 2, V is the 
 *                      first time derivative of the "rescaled"
 *                      rotation vector, i.e. mat_to_vec(rten, Psi), ||Psi|| <= pi
 *                      false -> if angularintype == 2, V is the
 *                      first time derivative of the rotation
 *                      vector rvec which is not necessarily rescaled
 *     rescaleout     : true -> if angularouttype == 2, v is the
 *                      first time derivative of the "rescaled"
 *                      rotation vector, i.e. mat_to_vec(rten, Psi), ||Psi|| <= pi
 *                      false -> if angularouttype == 2, v is the
 *                      first time derivative of the rotation
 *                      vector rvec which is not necessarily rescaled
 *
 *  Output:
 *     v : transformed angular acceleration
 *
 *******************************************************************/
{
  switch (angularintype) {

    case 0 : {
      switch(angularouttype) {
        case 0 : // spatial to spatial: a = A
          for(int j = 0; j < 3; ++j) a[j] = A[j];
          break;
        case 1 : // spatial to convected: a = R^T*A
          mat_mult_vec(rten, A, a, 1); 
          break;
#ifdef USE_EIGEN3
        case 2 : { // spatial to total: a = T^{-T}*(A - Tdot^T*Psidot);
          Eigen::Vector3d Psi;
          if(rescaleout) {
            mat_to_vec(rten, Psi.data());
          }
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          }
          Eigen::Matrix3d T, Tdot;
          tangential_transf(Psi, T);
          Eigen::Map<Eigen::Vector3d> omega(V), alpha(A), Psiddot(a);
          Eigen::Vector3d Psidot = T.transpose().inverse()*omega;
          tangential_transf_dot(Psi, Psidot, Tdot);
          Psiddot = T.transpose().inverse()*(alpha - Tdot.transpose()*Psidot);
        } break;
#endif
        default :
          std::cerr << " *** Warning: requested angular acceleration transformation from " 
                    << " angularintype = " << angularintype << " to "
                    << " angularouttype = " << angularouttype << " is not supported\n";
          for(int j = 0; j < 3; ++j) a[j] = 0;
      }
    } break;

    case 1 : {
      switch(angularouttype) {
        case 0 : // convected to spatial: a = R*A
          mat_mult_vec(rten, A, a, 0);
          break;
        case 1 : // convected to convected: a = A
          for(int j = 0; j < 3; ++j) a[j] = A[j];
          break;
#ifdef USE_EIGEN3
        case 2 : { // convected to total: a = T.inverse()*(A - Tdot*Psidot)
          Eigen::Vector3d Psi;
          if(rescaleout) {
            mat_to_vec(rten, Psi.data());
          } 
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          } 
          Eigen::Matrix3d T, Tdot;
          tangential_transf(Psi, T);
          Eigen::Map<Eigen::Vector3d> Omega(V), Alpha(A), Psiddot(a);
          Eigen::Vector3d Psidot = T.inverse()*Omega;
          tangential_transf_dot(Psi, Psidot, Tdot);
          Psiddot = T.inverse()*(Alpha - Tdot*Psidot);
        } break;
#endif
        default :
          std::cerr << " *** Warning: requested angular acceleration transformation from " 
                    << " angularintype = " << angularintype << " to "
                    << " angularouttype = " << angularouttype << " is not supported\n";
          for(int j = 0; j < 3; ++j) a[j] = 0;
      }
    } break;

    case 2 : {
      switch(angularouttype) {
#ifdef USE_EIGEN3
        case 0 : { // total to spatial: a = T^T*Psiddot + Tdot^T*Psidot;
          Eigen::Vector3d Psi;
          if(rescalein) {
            mat_to_vec(rten, Psi.data());
          }
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          }
          Eigen::Vector3d Psidot = Eigen::Map<Eigen::Vector3d>(V),
                          Psiddot = Eigen::Map<Eigen::Vector3d>(A);
          Eigen::Map<Eigen::Vector3d> alpha(a);
          Eigen::Matrix3d T, Tdot;
          tangential_transf(Psi, T);
          tangential_transf_dot(Psi, Psidot, Tdot);
          alpha = T.transpose()*Psiddot + Tdot.transpose()*Psidot;
        } break;
        case 1 : { // total to convected: a = T*Psiddot + Tdot*Psidot;
          Eigen::Vector3d Psi;
          if(rescalein) {
            mat_to_vec(rten, Psi.data());
          }
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          }
          Eigen::Vector3d Psidot = Eigen::Map<Eigen::Vector3d>(V),
                          Psiddot = Eigen::Map<Eigen::Vector3d>(A);
          Eigen::Map<Eigen::Vector3d> Alpha(a);
          Eigen::Matrix3d T, Tdot;
          tangential_transf(Psi, T);
          tangential_transf_dot(Psi, Psidot, Tdot);
          Alpha = T*Psiddot + Tdot*Psidot;
        } break;
#endif
        case 2 : // total to total:
          if(rescalein == rescaleout) {
            for(int j = 0; j < 3; ++j) a[j] = A[j];
          }
          else {
#ifdef USE_EIGEN3
            Eigen::Vector3d Psi, PsiOut;
            if(!rescalein && rescaleout) {
              Psi << rvec[0], rvec[1], rvec[2];
              mat_to_vec(rten, PsiOut.data());
            }
            else {
              PsiOut << rvec[0], rvec[1], rvec[2];
              mat_to_vec(rten, Psi.data());
            }
            if(Psi.isApprox(PsiOut)) {
              for(int j = 0; j < 3; ++j) a[j] = A[j];
            }
            else {
              Eigen::Matrix3d B = complement_transf(Psi);
              Eigen::Vector3d Psidot;
              Psidot << V[0], V[1], V[2];
              Eigen::Matrix3d Bdot = complement_transf_dot(Psi, Psidot);
              Eigen::Map<Eigen::Vector3d> PsiddotC(a), Psiddot(A);
              PsiddotC = B*Psiddot + Bdot*Psidot;
            }
#else
            std::cerr << " *** Warning: requested angular acceleration transformation from "
                      << " angularintype = " << angularintype << " to "
                      << " angularouttype = " << angularouttype << " is not supported\n";
#endif
          }
          break;
        default :
          std::cerr << " *** Warning: requested angular acceleration transformation from "
                    << " angularintype = " << angularintype << " to "
                    << " angularouttype = " << angularouttype << " is not supported\n";
          for(int j = 0; j < 3; ++j) a[j] = 0;
      }
    } break;

    default :
      std::cerr << " *** Warning: requested angular velocity transformation from "
                << " angularintype = " << angularintype << " to "
                << " angularouttype = " << angularouttype << " is not supported\n";
      for(int j = 0; j < 3; ++j) a[j] = 0;
  }
}
