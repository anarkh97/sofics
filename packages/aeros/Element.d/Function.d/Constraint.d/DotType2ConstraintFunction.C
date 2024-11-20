#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/DotType2ConstraintFunction.h>
#include <Element.d/Function.d/Constraint.d/exp-map.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for dot
// type 2 constraint function with double precision scalar
template<>
Eigen::Matrix<double,1,9>
Jacobian<double,DotType2ConstraintFunction>
::operator() (const Eigen::Matrix<double,9,1>& q, double)
{
  Eigen::Matrix<double,9,1> J;

  Eigen::Vector3d a0,b0,b0hat,a,bhat,d1;
  a0 << sconst(0), sconst(1), sconst(2);
  b0 << sconst(3), sconst(4), sconst(5);
  b0hat = b0.normalized();
  int type = iconst(0);

  a = a0 + q.segment<3>(6) - q.segment<3>(0);

  if( (q.segment<3>(3).array() == 0).all() ) {

    J.segment<3>(0) = -b0hat;
    J.segment<3>(3) = b0hat.cross(a);
    J.segment<3>(6) = b0hat;
    if(type == 2) return -J.transpose(); else return J.transpose();
  }

  // rotation parameters
  Eigen::Vector3d v1 = q.segment<3>(3);

  // rotated axis
  Eigen::Quaternion<double> z1;
  z1.setFromOneVector(v1);
  bhat = z1.toRotationMatrix()*b0hat;

  // partial derivatives of rotation matrices wrt rotation parameters
  double dRdvi_data[3][3];
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > dRdvi(&dRdvi_data[0][0]);

  // partial derivatives of constraint functions wrt x, y, z 
  for(int i = 0; i < 3; ++i) {
    J[0+i] = -bhat[i];
    J[6+i] =  bhat[i];
  }

  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(v1.data(), i, dRdvi_data);
    d1 = dRdvi*b0hat;
    J[3+i] = d1.dot(a);
  }

  if(type == 2) return -J.transpose(); else return J.transpose();
}

// specializing the member function template of constraint hessian operator for dot
// type 2 constraint function with double precision scalar
template<>
Eigen::Matrix<double,9,9>
Hessian<double,DotType2ConstraintFunction>
::operator() (const Eigen::Matrix<double,9,1>& q, double)
{
  Eigen::Matrix<double,9,9> H;
  H.setZero();

  Eigen::Vector3d a0,b0,b0hat,a,d1;
  a0 << sconst(0), sconst(1), sconst(2);
  b0 << sconst(3), sconst(4), sconst(5);
  b0hat = b0.normalized();
  int type = iconst(0);

  a = a0 + q.segment<3>(6) - q.segment<3>(0);

  if( (q.segment<3>(3).array() == 0).all() ) {

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

    if(type == 2) return -H; return H;
  }

  // rotation parameters
  Eigen::Vector3d v1 = q.segment<3>(3);

  // first and second partial derivatives of rotation matrices wrt rotation parameters
  double d2Rdvidvj_data[3][3], dRdvi_data[3][3];
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > d2Rdvidvj(&d2Rdvidvj_data[0][0]), dRdvi(&dRdvi_data[0][0]);

  for(int i=0; i<3; ++i) {
    for(int j=i; j<3; ++j) {
      Second_Partial_R_Partial_EM3(v1.data(), i, j, d2Rdvidvj_data);
      d1 = d2Rdvidvj*b0hat;
      H(3+i,3+j) = H(3+j,3+i) = d1.dot(a);
    }
  }
  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(v1.data(), i, dRdvi_data);
    d1 = dRdvi*b0hat;
    for(int j=0; j<3; ++j) {
      H(3+i,j)   = H(j,  3+i) = -d1[j];
      H(3+i,6+j) = H(6+j,3+i) =  d1[j];
    }
  }

  if(type == 2) return -H; return H;
}

} // namespace Simo

#endif
