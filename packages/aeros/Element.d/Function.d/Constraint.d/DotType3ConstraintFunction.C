#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/DotType3ConstraintFunction.h>
#include <Element.d/Function.d/Constraint.d/exp-map.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for dot
// type 2 (variant 2) constraint function with double precision scalar
template<>
Eigen::Matrix<double,1,12>
Jacobian<double,DotType3ConstraintFunction>
::operator() (const Eigen::Matrix<double,12,1>& q, double)
{
  Eigen::Matrix<double,12,1> J;

  Eigen::Vector3d a0,b0,b0hat,c0,c0hat,a,bhat,chat,d1,d2;
  a0 << sconst(0), sconst(1), sconst(2);
  b0 << sconst(3), sconst(4), sconst(5);
  c0 << sconst(6), sconst(7), sconst(8);
  b0hat = b0.normalized();
  c0hat = c0.normalized();

  a = a0 + q.segment<3>(6) - q.segment<3>(0);

  if( (q.segment<3>(3).array() == 0).all() && (q.segment<3>(9).array() == 0).all() ) {

    J.segment<3>(0) = -b0hat - c0hat;
    J.segment<3>(3) = b0hat.cross(a);
    J.segment<3>(6) = b0hat + c0hat;
    J.segment<3>(9) = c0hat.cross(a);
    return J.transpose();
  }

  // rotation parameters
  Eigen::Vector3d v1 = q.segment<3>(3);
  Eigen::Vector3d v2 = q.segment<3>(9);

  // rotated axis
  Eigen::Quaternion<double> z1;
  z1.setFromOneVector(v1);
  bhat = z1.toRotationMatrix()*b0hat;

  Eigen::Quaternion<double> z2;
  z2.setFromOneVector(v2);
  chat = z2.toRotationMatrix()*c0hat;

  // partial derivatives of rotation matrices wrt rotation parameters
  double dRdvi_data[3][3];
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > dRdvi(&dRdvi_data[0][0]);

  // partial derivatives of constraint functions wrt x, y, z 
  for(int i = 0; i < 3; ++i) {
    J[0+i] = -bhat[i] - chat[i];
    J[6+i] =  bhat[i] + chat[i];
  }

  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(v1.data(), i, dRdvi_data);
    d1 = dRdvi*b0hat;
    Partial_R_Partial_EM3(v2.data(), i, dRdvi_data);
    d2 = dRdvi*c0hat;
    J[3+i] = d1.dot(a);
    J[9+i] = d2.dot(a);
  }

  return J.transpose();
}

// specializing the member function template of constraint hessian operator for dot
// type 2 (variant 2) constraint function with double precision scalar
template<> 
Eigen::Matrix<double,12,12>
Hessian<double,DotType3ConstraintFunction>
::operator() (const Eigen::Matrix<double,12,1>& q, double)
{
  Eigen::Matrix<double,12,12> H;
  H.setZero();

  Eigen::Vector3d a0,b0,c0,b0hat,c0hat,a,d1,d2;
  a0 << sconst(0), sconst(1), sconst(2);
  b0 << sconst(3), sconst(4), sconst(5);
  c0 << sconst(6), sconst(7), sconst(8);
  b0hat = b0.normalized();
  c0hat = c0.normalized();

  a = a0 + q.segment<3>(6) - q.segment<3>(0);

  if( (q.segment<3>(3).array() == 0).all() && (q.segment<3>(9).array() == 0).all() ) {

    Eigen::Matrix3d ab = a*b0hat.transpose();

    H(3,3) = -ab(1,1) - ab(2,2);
    H(4,4) = -ab(0,0) - ab(2,2);
    H(5,5) = -ab(0,0) - ab(1,1);
    H(3,4) = H(4,3) = 0.5*(ab(0,1) + ab(1,0));
    H(3,5) = H(5,3) = 0.5*(ab(0,2) + ab(2,0));
    H(4,5) = H(5,4) = 0.5*(ab(1,2) + ab(2,1));

    H(4,2) = H(2,4) = H(5,7) = H(7,5) = b0hat[0];
    H(4,8) = H(8,4) = H(5,1) = H(1,5) = -b0hat[0];

    H(3,8) = H(8,3) = H(5,0) = H(0,5) = b0hat[1];
    H(3,2) = H(2,3) = H(5,6) = H(6,5) = -b0hat[1];

    H(3,1) = H(1,3) = H(4,6) = H(6,4) = b0hat[2];
    H(3,7) = H(7,3) = H(4,0) = H(0,4) = -b0hat[2];

    Eigen::Matrix3d ac = a*c0hat.transpose();

    H(9,9)   = -ac(1,1) - ac(2,2);
    H(10,10) = -ac(0,0) - ac(2,2);
    H(11,11) = -ac(0,0) - ac(1,1);
    H(9,10)  = H(10,9)  = 0.5*(ac(0,1) + ac(1,0));
    H(9,11)  = H(11,9)  = 0.5*(ac(0,2) + ac(2,0));
    H(10,11) = H(11,10) = 0.5*(ac(1,2) + ac(2,1));

    H(10,2) = H(2,10) = H(11,7) = H(7,11) = c0hat[0];
    H(10,8) = H(8,10) = H(11,1) = H(1,11) = -c0hat[0];

    H(9,8) = H(8,9) = H(11,0) = H(0,11) = c0hat[1];
    H(9,2) = H(2,9) = H(11,6) = H(6,11) = -c0hat[1];

    H(9,1) = H(1,9) = H(10,6) = H(6,10) = c0hat[2];
    H(9,7) = H(7,9) = H(10,0) = H(0,10) = -c0hat[2];

    return H;
  }


  // rotation parameters
  Eigen::Vector3d v1 = q.segment<3>(3);
  Eigen::Vector3d v2 = q.segment<3>(9);

  // first and second partial derivatives of rotation matrices wrt rotation parameters
  double d2Rdvidvj_data[3][3], dRdvi_data[3][3];
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > d2Rdvidvj(&d2Rdvidvj_data[0][0]), dRdvi(&dRdvi_data[0][0]);

  for(int i=0; i<3; ++i) {
    for(int j=i; j<3; ++j) {
      Second_Partial_R_Partial_EM3(v1.data(), i, j, d2Rdvidvj_data);
      d1 = d2Rdvidvj*b0hat;
      H(3+i,3+j) = H(3+j,3+i) = d1.dot(a);

      Second_Partial_R_Partial_EM3(v2.data(), i, j, d2Rdvidvj_data);
      d2 = d2Rdvidvj*c0hat;
      H(9+i,9+j) = H(9+j,9+i) = d2.dot(a);
    }
  }
  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(v1.data(), i, dRdvi_data);
    d1 = dRdvi*b0hat;
    for(int j=0; j<3; ++j) {
      H(3+i,j)   = H(j,  3+i) = -d1[j];
      H(3+i,6+j) = H(6+j,3+i) =  d1[j];
    }

    Partial_R_Partial_EM3(v2.data(), i, dRdvi_data);
    d2 = dRdvi*c0hat;
    for(int j=0; j<3; ++j) {
      H(9+i,j)   = H(j,  9+i) = -d2[j];
      H(9+i,6+j) = H(6+j,9+i) =  d2[j];
    }
  }

  return H;
}

} // namespace Simo

#endif
