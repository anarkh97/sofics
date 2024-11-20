#ifdef USE_EIGEN3
#ifndef _SHEARPANELSTRESSWRTDISPLACEMENTSENSITIVITY_H_
#define _SHEARPANELSTRESSWRTDISPLACEMENTSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Shear.d/ShearPanelTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class ShearPanelStressWRTDisplacementSensitivity : public VectorValuedFunction<12,4,Scalar,14,0,double>
{
  public:
    ShearPanelTemplate<Scalar> ele;
    Eigen::Array<Scalar,4,1> globalx, globaly, globalz; // nodal coordinates
    Scalar E, G; // material properties

  public:
    ShearPanelStressWRTDisplacementSensitivity(const Eigen::Array<double,14,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      globalx = sconst.segment<4>(0).cast<Scalar>();
      globaly = sconst.segment<4>(4).cast<Scalar>();
      globalz = sconst.segment<4>(8).cast<Scalar>();
      E = sconst[12];
      G = sconst[13];
    }

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,12,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      Eigen::Matrix<Scalar,12,1> globalu = q;

      Eigen::Array<Scalar,6,4> elStress;
      Eigen::Array<Scalar,6,4> elStrain;
      Scalar vmssig, vmseps;

      ele.spstress(globalx.data(), globaly.data(), globalz.data(), globalu.data(), 
                   G, E, 0, 0, 
                   elStress.data(), elStrain.data(), vmssig, vmseps); 

      // return value:
      Eigen::Matrix<Scalar,4,1> v;
      v[3] = v[2] = v[1] = v[0] = vmssig; 

      return v; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
