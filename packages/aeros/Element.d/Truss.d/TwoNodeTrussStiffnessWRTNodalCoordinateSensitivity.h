#ifdef USE_EIGEN3
#ifndef _TWONODETRUSSSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_
#define _TWONODETRUSSSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Truss.d/TrussElementTemplate.cpp>

template<typename Scalar>
class TwoNodeTrussStiffnessWRTNodalCoordinateSensitivity : public MatrixValuedFunction<6,6,6,Scalar,3,0,double>
{
  public:
    TrussElementTemplate<Scalar> ele;
    Eigen::Array<Scalar,2,1> globalx, globaly, globalz; // nodal coordinates
    Scalar E, A, preload; 

  public:
    TwoNodeTrussStiffnessWRTNodalCoordinateSensitivity(const Eigen::Array<double,3,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      E = sconst[0];
      A = sconst[1];
      preload = sconst[2];
    }

    Eigen::Matrix<Scalar,6,6> operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar)
    {
      Scalar x[2], y[2], z[2];
      x[0] = q[0];   x[1] = q[3];
      y[0] = q[1];   y[1] = q[4];
      z[0] = q[2];   z[1] = q[5];
      Eigen::Array<Scalar,6,6> estiff;
      ele.stiffnessMatrix(estiff.data(), (Scalar *) x, (Scalar *) y, (Scalar *) z, E, A, preload); 

      return estiff; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif


