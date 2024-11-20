#ifdef USE_EIGEN3
#ifndef _TWONODETRUSSSTRESSWRTNODALCOORDINATESENSITIVITY_H_
#define _TWONODETRUSSSTRESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Truss.d/TrussElementTemplate.cpp>

template<typename Scalar>
class TwoNodeTrussStressWRTNodalCoordinateSensitivity : public VectorValuedFunction<6,2,Scalar,11,0,double>
{
  public:
    TrussElementTemplate<Scalar> ele;
    Scalar E, A, W, Ta, preload; 
    Scalar elDisp[6];

  public:
    TwoNodeTrussStressWRTNodalCoordinateSensitivity(const Eigen::Array<double,11,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      E = sconst[0];
      A = sconst[1];
      W = sconst[2];
      Ta = sconst[3];
      preload = sconst[4];
      elDisp[0] = sconst[5];
      elDisp[1] = sconst[6];
      elDisp[2] = sconst[7];
      elDisp[3] = sconst[8];
      elDisp[4] = sconst[9];
      elDisp[5] = sconst[10];
    }

    Eigen::Matrix<Scalar,2,1> operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar)
    {
      Scalar x[2], y[2], z[2];
      x[0] = q[0];   x[1] = q[3];
      y[0] = q[1];   y[1] = q[4];
      z[0] = q[2];   z[1] = q[5];
      Eigen::Array<Scalar,2,1> stress;
      int strInd = 6;
      ele.vmst(stress.data(), (Scalar *) x, (Scalar *) y, (Scalar *) z, elDisp,
               E, A, W, Ta, preload, strInd); 

      return stress; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
