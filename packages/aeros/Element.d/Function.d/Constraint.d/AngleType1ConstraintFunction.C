#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/AngleType1ConstraintFunction.h>
#include <Element.d/Function.d/Constraint.d/exp-map.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for angle
// type 1 constraint function with double precision scalar
template<> 
Eigen::Matrix<double,1,6>
Jacobian<double,AngleType1ConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double)
{
  // return value
  Eigen::Matrix<double,6,1> J;

  Eigen::Vector3d a0,b0,a0hat,b0hat,ahat,bhat,d1,d2;
  a0 << sconst(0), sconst(1), sconst(2);
  b0 << sconst(3), sconst(4), sconst(5);
  a0hat = a0.normalized();
  b0hat = b0.normalized();
  int type = iconst(0);

  if( (q.array() == 0).all()) {

    double x = a0hat.dot(b0hat);     // cos(theta)
    double dacosdx = -1/sqrt(1-x*x); // derivative of arccos(x) w.r.t. x
    J.head<3>() = a0hat.cross(b0hat);
    J.tail<3>() = -J.head<3>();
    J *= dacosdx;
    if(type == 2) return -J.transpose(); else return J.transpose();
  }

  // rotation parameters
  Eigen::Vector3d v1 = q.segment<3>(0);
  Eigen::Vector3d v2 = q.segment<3>(3);

  // rotated axes
  Eigen::Quaternion<double> z1, z2;
  z1.setFromOneVector(v1); z2.setFromOneVector(v2);
  ahat = z1.toRotationMatrix()*a0hat;
  bhat = z2.toRotationMatrix()*b0hat;

  double x = ahat.dot(bhat);       // cos(theta)
  double dacosdx = -1/sqrt(1-x*x); // derivative of arccos(x) w.r.t. x

  // partial derivatives of rotation matrices wrt rotation parameters
  double dRdvi1_data[3][3], dRdvi2_data[3][3];
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > dRdvi1(&dRdvi1_data[0][0]), dRdvi2(&dRdvi2_data[0][0]);

  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(v1.data(), i, dRdvi1_data);
    Partial_R_Partial_EM3(v2.data(), i, dRdvi2_data);

    d1 = dRdvi1*a0hat;
    d2 = dRdvi2*b0hat;

    J[0+i] = dacosdx*(d1.dot(bhat));
    J[3+i] = dacosdx*(ahat.dot(d2));
  }

  if(type == 2) return -J.transpose(); else return J.transpose();
}

// specializing the member function template of constraint hessian operator for angle
// type 1 constraint function with double precision scalar
template<>
Eigen::Matrix<double,6,6>
Hessian<double,AngleType1ConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t)
{
  Eigen::Matrix<double,6,6> H;

  using std::sqrt;
  H.setZero();

  Eigen::Vector3d a0, b0, a0hat, b0hat, ahat, bhat, d1, d2;
  a0 << sconst(0), sconst(1), sconst(2);
  b0 << sconst(3), sconst(4), sconst(5);
  a0hat = a0.normalized();
  b0hat = b0.normalized();
  int type = iconst(0);

  if( (q.array() == 0).all()) {

    double x = a0hat.dot(b0hat);     // cos(theta)
    double dacosdx = -1/sqrt(1-x*x); // derivative of arccos(x) w.r.t. x
    double xdacosdx = x*dacosdx;     // d2acosdx2/(dacosdx^2), where d2acosdx2 is the second derivative of arccos(x) w.r.t. x

    Eigen::Matrix3d ab = dacosdx*a0hat*b0hat.transpose();

    H(0,0) = H(3,3) = -ab(1,1) - ab(2,2);
    H(1,1) = H(4,4) = -ab(0,0) - ab(2,2);
    H(2,2) = H(5,5) = -ab(0,0) - ab(1,1);

    H(0,1) = H(1,0) = H(3,4) = H(4,3) = 0.5*(ab(0,1) + ab(1,0));
    H(0,2) = H(2,0) = H(3,5) = H(5,3) = 0.5*(ab(0,2) + ab(2,0));
    H(1,2) = H(2,1) = H(4,5) = H(5,4) = 0.5*(ab(1,2) + ab(2,1));

    H(0,3) = H(3,0) =  ab(2,2) + ab(1,1);
    H(0,4) = H(4,0) = -ab(1,0);
    H(0,5) = H(5,0) = -ab(2,0);

    H(1,3) = H(3,1) = -ab(0,1);
    H(1,4) = H(4,1) =  ab(2,2) + ab(0,0);
    H(1,5) = H(5,1) = -ab(2,1);

    H(2,3) = H(3,2) = -ab(0,2);
    H(2,4) = H(4,2) = -ab(1,2);
    H(2,5) = H(5,2) =  ab(1,1) + ab(0,0);

    Eigen::Matrix<double,3,1> tmp1 = dacosdx*a0hat.cross(b0hat);     // tmp1 = J.head<3>() = -J.tail<3>()
    Eigen::Matrix<double,3,3> tmp2 = xdacosdx*tmp1*tmp1.transpose(); 
    Eigen::Matrix<double,6,6> tmp3; 
    tmp3 << tmp2, -tmp2,
           -tmp2,  tmp2; // now tmp3 = xdacosdx*J*J.transpose()

    H += tmp3;

    if(type == 2) return -H; else return H;
  }


  // rotation parameters
  Eigen::Vector3d v1 = q.segment<3>(0);
  Eigen::Vector3d v2 = q.segment<3>(3);

  // rotated axes
  Eigen::Quaternion<double> z1, z2;
  z1.setFromOneVector(v1); z2.setFromOneVector(v2);
  ahat = z1.toRotationMatrix()*a0hat;
  bhat = z2.toRotationMatrix()*b0hat;

  // evaluate the jacobian
  Eigen::Matrix<double,1,6> J;
  Jacobian<double,AngleType1ConstraintFunction> dfdq(sconst, iconst);
  J = dfdq(q, t);

  double x = ahat.dot(bhat);       // cos(theta)
  double dacosdx = -1/sqrt(1-x*x); // derivative of arccos(x) w.r.t. x
  double xdacosdx = x*dacosdx;     // d2acosdx2/(dacosdx^2), where d2acosdx2 is the second derivative of arccos(x) w.r.t. x

  // first and second partial derivatives of rotation matrices wrt rotation parameters
  double d2Rdvidvj_data[3][3], dRdvi_data[3][3], dRdvj_data[3][3];
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > d2Rdvidvj(&d2Rdvidvj_data[0][0]), dRdvi(&dRdvi_data[0][0]), dRdvj(&dRdvj_data[0][0]);

  for(int i=0; i<3; ++i) {
    for(int j=i; j<3; ++j) {
      Second_Partial_R_Partial_EM3(v1.data(), i, j, d2Rdvidvj_data);
      d1 = d2Rdvidvj*a0hat;

      Second_Partial_R_Partial_EM3(v2.data(), i, j, d2Rdvidvj_data);
      d2 = d2Rdvidvj*b0hat;

      H(i,j) = H(j,i) = dacosdx*d1.dot(bhat) + xdacosdx*(J[i]*J[j]);
      H(3+i,3+j) = H(3+j,3+i) = dacosdx*ahat.dot(d2) +  xdacosdx*(J[3+i]*J[3+j]);
    }
  }
  for(int i=0; i<3; ++i) {
    Partial_R_Partial_EM3(v1.data(), i, dRdvi_data);
    d1 = dRdvi*a0hat;

    for(int j=0; j<3; ++j) {
      Partial_R_Partial_EM3(v2.data(), j, dRdvj_data);
      d2 = dRdvj*b0hat;
      H(i,3+j) = H(3+j,i) = dacosdx*d1.dot(d2) + xdacosdx*(J[i]*J[3+j]);
    }
  }

  if(type == 2) return -H; else return H;
}

} // namespace Simo

#endif
