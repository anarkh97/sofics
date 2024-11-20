#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTIFFNESSWRTTHICKNESSSENSITIVITY_H_
#define _SHELLELEMENTSTIFFNESSWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterial.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType1.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the stiffness matrix w.r.t the shell thickness

template<typename Scalar>
class ShellElementStiffnessWRTThicknessSensitivity : public MatrixValuedFunction<1,18,18,Scalar,68,3,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz; // nodal coordinates
    Eigen::Array<Scalar,18,1> globalu;                  // nodal displacements
    Scalar E, nu, rho, Ta, W; // material properties
    Eigen::Array<Scalar,9,1> cframe; // composite frame
    Eigen::Array<Scalar,42,1> coefs;
    Eigen::Array<Scalar,3,1> ndtemps;
    int type;
    int flag;
    int tflg;

  public:
    ShellElementStiffnessWRTThicknessSensitivity(const Eigen::Array<double,68,1>& sconst, const Eigen::Array<int,3,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      E = sconst[9];
      nu = sconst[10];
      rho = sconst[11];
      globalu.setZero();
      type = iconst[0];
      flag = iconst[1];
      tflg = iconst[2];
      if(type == 1) {
        cframe = sconst.segment<9>(12).cast<Scalar>();
        coefs = sconst.segment<42>(21).cast<Scalar>();
      }
      Ta = sconst[63]; 
      W = sconst[64];
      ndtemps = sconst.segment<3>(65).cast<Scalar>();
    }

    Eigen::Matrix<Scalar,18,18> operator() (const Eigen::Matrix<Scalar,1,1>& h, Scalar)
    {
      // inputs:
      // h = thickness
      ShellMaterial<Scalar> *gpmat;
      switch(type) {
        case 0 :
          gpmat = new ShellMaterialType0<Scalar>(E, h[0], nu, rho, Ta, W);
          break;
        case 1 : 
          gpmat = new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, h[0], Ta);
          break;
        default :
          std::cerr << " *** ERROR: ShellElementStiffnessWRTThicknessSensitivity is not defined for this case.\n";
          exit(-1);
      }

      // elm      <input>   Finite Element Number                           not actually used
      // nu       <input>   Poisson's Ratio (for an Isotropic Element)
      // globalX  <input>   X- Nodal Coordinates
      // globalY  <input>   Y- Nodal Coordinates
      // globalZ  <input>   Z- Nodal Coordinates
      // estiff   <output>  Element Stiffness 
      // type     <input>   Type of Constitutive Law (0, 1, 2, 3, or 4)     0 for linear elastic
      // flag     <input>   0: return local, 1: return global

      Eigen::Matrix<Scalar,18,18> estiff;
      ele.andesstf(0, estiff.data(), (Scalar*)NULL, nu, globalx.data(), globaly.data(), globalz.data(), globalu.data(),
                   type, gpmat, flag, tflg, ndtemps.data());
      delete gpmat;

      // return value:
      // element stiffness matrix
      return estiff; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif

