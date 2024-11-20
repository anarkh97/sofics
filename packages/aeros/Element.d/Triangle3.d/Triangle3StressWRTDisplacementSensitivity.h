#ifdef USE_EIGEN3
#ifndef _TRIANGLE3STRESSWRTDISPLACEMENTSENSITIVITY_H_
#define _TRIANGLE3STRESSWRTDISPLACEMENTSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Triangle3.d/Triangle3ElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class Triangle3StressWRTDisplacementSensitivity : public VectorValuedFunction<6,3,Scalar,8,1,double>
{
  public:
    Triangle3ElementTemplate<Scalar> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly; // nodal coordinates
    Scalar E, nu; // material properties

  public:
    Triangle3StressWRTDisplacementSensitivity(const Eigen::Array<double,8,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      E = sconst[6];
      nu = sconst[7]; 
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      Eigen::Matrix<Scalar,6,1> globalu = q;

      Eigen::Array<Scalar,3,7> elStress;
      Eigen::Array<Scalar,3,7> elStrain;

      ele.sands4(globalx.data(), globaly.data(), 
                 globalu.data(),elStress.data(),
                 elStrain.data(), E, nu);

      // return value:

      Eigen::Matrix<Scalar,3,1> v;
      v[0] = elStress(0,6); 
      v[1] = elStress(1,6);
      v[2] = elStress(2,6);

      return v; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif

