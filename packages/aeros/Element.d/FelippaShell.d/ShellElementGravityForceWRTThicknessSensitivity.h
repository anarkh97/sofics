#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTGRAVITYFORCEWRTTHICKNESSSENSITIVITY_H_
#define _SHELLELEMENTGRAVITYFORCEWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>

// class template to facilitate computation of the sensitivities of the gravity force w.r.t the shell thickness

template<typename Scalar>
class ShellElementGravityForceWRTThicknessSensitivity : public VectorValuedFunction<1,18,Scalar,14,1,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Matrix<Scalar,3,1> gravityAcceleration, globalx, globaly, globalz;
    Scalar rho; // mean density (mass per unit volume). For a composite this is the area density (mass per unit area) divided by the total thickness
    Scalar nsm; // non-structural mass (mass per unit area)
    int gravflg;

  public:
    ShellElementGravityForceWRTThicknessSensitivity(const Eigen::Array<double,14,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      gravityAcceleration = sconst.segment<3>(0).cast<Scalar>();
      globalx = sconst.segment<3>(3).cast<Scalar>();
      globaly = sconst.segment<3>(6).cast<Scalar>();
      globalz = sconst.segment<3>(9).cast<Scalar>();
      rho = sconst[12];
      nsm = sconst[13];
      gravflg = iconst[0];
    }

    Eigen::Matrix<Scalar,18,1> operator() (const Eigen::Matrix<Scalar,1,1>& q, Scalar)
    {
      // inputs:
      // q[0] = shell thickness

      Eigen::Matrix<Scalar,18,1> gravityForce;
      Scalar rhoh = nsm + rho*q[0];
      ele.andesgf(1, globalx.data(), globaly.data(), globalz.data(), gravityForce.data(),
                  gravityAcceleration.data(), gravflg, rhoh);

      // return value:
      return gravityForce; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
