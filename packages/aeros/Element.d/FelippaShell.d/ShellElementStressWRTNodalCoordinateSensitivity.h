#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTRESSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTSTRESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal coordinates

template<typename Scalar>
class ShellElementStressWRTNodalCoordinateSensitivity : public VectorValuedFunction<9,3,Scalar,78,3,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,18,1> globalu; // element displacements
    Scalar E, nu, rho, h, Ta, W; // material properties
    Eigen::Array<Scalar,9,1> cframe; // composite frame
    Eigen::Array<Scalar,42,1> coefs;
    Eigen::Array<Scalar,3,1> ndtemps;
    int surface; // thru-thickness location at which stresses are to be evaluated
    int type;
    int sflg;

  public:
    ShellElementStressWRTNodalCoordinateSensitivity(const Eigen::Array<double,78,1>& sconst, const Eigen::Array<int,3,1>& iconst)
    {
      globalu = sconst.segment<18>(0).cast<Scalar>();
      E = sconst[18];
      nu = sconst[19];
      rho = sconst[20];
      h = sconst[21];
      surface = iconst[0];
      type = iconst[1];
      sflg = iconst[2];
      if(type == 1 || type == 5) {
        cframe = sconst.segment<9>(22).cast<Scalar>();
        coefs = sconst.segment<42>(31).cast<Scalar>(); 
      }
      Ta = sconst[73];
      W = sconst[74];
      ndtemps = sconst.segment<3>(75).cast<Scalar>();
    }
    
    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = nodal coordinates

      ShellMaterial<Scalar> *nmat;
      switch(type) {
        case 0 :
          nmat = new ShellMaterialType0<Scalar>(E, h, nu, rho, Ta, W);
          break;
        case 1 : 
          nmat = new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, h, Ta);
          break;
        case 5 :
          nmat = new ShellMaterialType5<Scalar>(coefs.data(), cframe.data(), rho, h, Ta);
          break;
        default :
          std::cerr << " *** ERROR: ShellElementStressWRTNodalCoordinateSensitivity is not defined for this case.\n";
          exit(-1);
      }
      Eigen::Array<Scalar,3,1> globalx; 
      Eigen::Array<Scalar,3,1> globaly; 
      Eigen::Array<Scalar,3,1> globalz; 
      globalx << q[0], q[3], q[6];
      globaly << q[1], q[4], q[7];
      globalz << q[2], q[5], q[8];

      // elm      <input>   Finite Element Number                           not actually used
      // maxstr   <input>   Maximum Number of Stresses                      7 (6 components of symmetric stress tensor + 1 von mises)
      // nu       <input>   Poisson's Ratio (for an Isotropic Element)
      // globalX  <input>   X- Nodal Coordinates
      // globalY  <input>   Y- Nodal Coordinates
      // globalZ  <input>   Z- Nodal Coordinates
      // stress   <output>  Stresses (Von Mises Stress) of the Element
      // ctyp     <input>   Type of Constitutive Law (0, 1, 2, 3, or 4)     0 for linear elastic
      // flag     <input>   0: stress, 1: strain, >= 2: internal variables
      // surface  <input>   1: upper, 2: median, 3: lower
      Eigen::Array<Scalar,7,3> stress;
      ele.andesvms(0, 7, nu, globalx.data(), globaly.data(), globalz.data(), globalu.data(),
                   stress.data(), type, nmat, 0, surface, sflg, ndtemps.data());
      delete nmat;

      // return value:
      // von mises stresses at nodes
      Eigen::Matrix<Scalar,3,1> v;
      v << stress(6,0), stress(6,1), stress(6,2);

      return v; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
