#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/utilities.hpp>
#include <iostream>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

void tran_veloc(double rten[3][3], double rvec[3], double V[3], int angularintype,
                int angularouttype, bool rescalein, bool rescaleout, double v[3])
/*****************************************************************
 *
 *  Purpose:
 *     Transform the angular velocity 
 *
 *  Input:
 *     rten           : rotation tensor
 *     rvec           : Euler rotation vector (not necessarily rescaled)
 *     V              : input angular velocity
 *     angularintype  : type of the input angular velocity
 *                      0 -> spatial, 1 -> convected, 2 -> total
 *     angularouttype : type of the transformed angular velocity 
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
 *     v : transformed angular velocity
 *
 *******************************************************************/
{
  switch (angularintype) {

    case 0 : {
      switch(angularouttype) {
        case 0 : // spatial to spatial: v = V
          for(int j = 0; j < 3; ++j) v[j] = V[j];
          break;
        case 1 : // spatial to convected: v = R^T*V
          mat_mult_vec(rten, V, v, 1);
          break;
#ifdef USE_EIGEN3
        case 2 : { // spatial to total: v = T^{-T}*V
          Eigen::Vector3d Psi;
          if(rescaleout) {
            mat_to_vec(rten, Psi.data());
          }
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          }
          Eigen::Matrix3d T;
          tangential_transf(Psi, T);
          Eigen::Map<Eigen::Vector3d> omega(V), Psidot(v);
          Psidot = T.transpose().inverse()*omega;
        } break;
#endif
        default :
          std::cerr << " *** Warning: requested angular velocity transformation from " 
                    << " angularintype = " << angularintype << " to "
                    << " angularouttype = " << angularouttype << " is not supported\n";
          for(int j = 0; j < 3; ++j) v[j] = 0;
      }
    } break;

    case 1 : {
      switch(angularouttype) {
        case 0 : // convected to spatial: v = R*V
          mat_mult_vec(rten, V, v, 0);
          break;
        case 1 : // convected to convected: v = V
          for(int j = 0; j < 3; ++j) v[j] = V[j];
          break;
#ifdef USE_EIGEN3
        case 2 : { // convected to total: v = T^{-1}*V
          Eigen::Vector3d Psi;
          if(rescaleout) {
            mat_to_vec(rten, Psi.data());
          } 
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          } 
          Eigen::Matrix3d T;
          tangential_transf(Psi, T);
          Eigen::Map<Eigen::Vector3d> Omega(V), Psidot(v);
          Psidot = T.inverse()*Omega;
        } break;
#endif
        default :
          std::cerr << " *** Warning: requested angular velocity transformation from " 
                    << " angularintype = " << angularintype << " to "
                    << " angularouttype = " << angularouttype << " is not supported\n";
          for(int j = 0; j < 3; ++j) v[j] = 0;
      }
    } break;

    case 2 : {
      switch(angularouttype) {
#ifdef USE_EIGEN3
        case 0 : { // total to spatial: V = T^T*v
          Eigen::Vector3d Psi;
          if(rescalein) {
            mat_to_vec(rten, Psi.data());
          }
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          }
          Eigen::Map<Eigen::Vector3d> omega(v), Psidot(V);
          Eigen::Matrix3d T;
          tangential_transf(Psi, T);
          omega = T.transpose()*Psidot;
        } break;
        case 1 : { // total to convected: V = T*v
          Eigen::Vector3d Psi;
          if(rescalein) {
            mat_to_vec(rten, Psi.data());
          }
          else {
            Psi << rvec[0], rvec[1], rvec[2];
          }
          Eigen::Map<Eigen::Vector3d> Omega(v), Psidot(V);
          Eigen::Matrix3d T;
          tangential_transf(Psi, T);
          Omega = T*Psidot;
        } break;
#endif
        case 2 : { // total to total:
          if(rescalein == rescaleout) {
            for(int j = 0; j < 3; ++j) v[j] = V[j];
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
              for(int j = 0; j < 3; ++j) v[j] = V[j];
            }
            else {
              Eigen::Matrix3d B = complement_transf(Psi);
              Eigen::Map<Eigen::Vector3d> PsidotC(v), Psidot(V);
              PsidotC = B*Psidot;
            }
#else
            std::cerr << " *** Warning: requested angular velocity transformation from "
                      << " angularintype = " << angularintype << " to "
                      << " angularouttype = " << angularouttype << " is not supported\n";
#endif
          }
        }  break;
        default :
          std::cerr << " *** Warning: requested angular velocity transformation from "
                    << " angularintype = " << angularintype << " to "
                    << " angularouttype = " << angularouttype << " is not supported\n";
          for(int j = 0; j < 3; ++j) v[j] = 0;
      }
    } break;

    default :
      std::cerr << " *** Warning: requested angular velocity transformation from "
                << " angularintype = " << angularintype << " to "
                << " angularouttype = " << angularouttype << " is not supported\n";
      for(int j = 0; j < 3; ++j) v[j] = 0;
  }
}
