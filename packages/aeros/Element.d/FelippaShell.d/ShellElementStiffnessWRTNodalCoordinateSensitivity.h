#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>

// class template to facilitate computation of the sensitivities of the stiffness matrix w.r.t the nodal coordinates

template<typename Scalar>
class ShellElementStiffnessWRTNodalCoordinateSensitivity : public MatrixValuedFunction<9,18,18,Scalar,60,2,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,18,1> globalu; // nodal displacements
    Scalar E, nu, rho, h, Ta, W;       // material properties
    Eigen::Array<Scalar,9,1> cframe;   // composite frame
    Eigen::Array<Scalar,42,1> coefs;
    Eigen::Array<Scalar,3,1> ndtemps;
    int type;
    int tflg;

  public:
    ShellElementStiffnessWRTNodalCoordinateSensitivity(const Eigen::Array<double,60,1>& sconst, const Eigen::Array<int,2,1>& iconst)
    {
      globalu.setZero();
      E = sconst[0];
      nu = sconst[1];
      rho = sconst[2];
      h = sconst[3];
      type = iconst[0];
      tflg = iconst[1];
      if(type == 1 || type == 5) {
        cframe = sconst.segment<9>(4).cast<Scalar>();
        coefs = sconst.segment<42>(13).cast<Scalar>(); 
      }
      Ta = sconst[55];
      W = sconst[56];
      ndtemps = sconst.segment<3>(57).cast<Scalar>();
    }

    Eigen::Matrix<Scalar,18,18> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = Nodal coordinates [x0, y0, z0, x1, y1, z1, x2, y2, z2]
      Eigen::Matrix<Scalar,3,1> globalx, globaly, globalz; 
      globalx << q[0], q[3], q[6];
      globaly << q[1], q[4], q[7];
      globalz << q[2], q[5], q[8];

      ShellMaterial<Scalar> *gpmat;
      switch(type) {
        case 0 :
          gpmat = new ShellMaterialType0<Scalar>(E, h, nu, rho, Ta, W); 
          break;
        case 1 :
          gpmat = new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, h, Ta);
          break;
        case 5 :
          gpmat = new ShellMaterialType5<Scalar>(coefs.data(), cframe.data(), rho, h, Ta);
          break;
        default :
          std::cerr << " *** ERROR: ShellElementStiffnessWRTNodalCoordinateSensitivity is not defined for this case.\n";
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
      int flag = 1;
      Eigen::Matrix<Scalar,18,18> estiff;
      ele.andesstf(0, estiff.data(), (Scalar*)NULL, nu, globalx.data(), globaly.data(), globalz.data(),
                   globalu.data(), type, gpmat, flag, tflg, ndtemps.data());

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
