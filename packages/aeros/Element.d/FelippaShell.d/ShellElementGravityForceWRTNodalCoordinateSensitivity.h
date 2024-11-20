#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>

// class template to facilitate computation of the sensitivities of the gravity force w.r.t the nodal coordinates

template<typename Scalar>
class ShellElementGravityForceWRTNodalCoordinateSensitivity : public VectorValuedFunction<9,18,Scalar,4,1,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Matrix<Scalar,3,1> gravityAcceleration;
    Scalar rhoh; // area density
    int gravflg;

  public:
    ShellElementGravityForceWRTNodalCoordinateSensitivity(const Eigen::Array<double,4,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      gravityAcceleration = sconst.segment<3>(0).cast<Scalar>();
      rhoh = sconst[3];
      gravflg = iconst[0];
    }

    Eigen::Matrix<Scalar,18,1> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = nodal coordinates
      
      Scalar x[3] = { q[0], q[3], q[6] };
      Scalar y[3] = { q[1], q[4], q[7] };
      Scalar z[3] = { q[2], q[5], q[8] };

      Eigen::Matrix<Scalar,18,1> gravityForce;

      ele.andesgf(1, x, y, z, gravityForce.data(),
                  gravityAcceleration.data(), gravflg, rhoh);

      // return value:
      return gravityForce; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
