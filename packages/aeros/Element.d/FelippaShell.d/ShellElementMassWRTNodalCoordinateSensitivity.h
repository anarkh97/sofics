#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTMASSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTMASSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>

// class template to facilitate computation of the sensitivities of the mass w.r.t the nodal coordinates

template<typename Scalar>
class ShellElementMassWRTNodalCoordinateSensitivity : public ScalarValuedFunction<9,Scalar,1,0,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Scalar rhoh; // area density

  public:
    ShellElementMassWRTNodalCoordinateSensitivity(const Eigen::Array<double,1,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      rhoh = sconst[0];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = nodal coordinates

      Scalar globalx[3] = { q[0], q[3], q[6] };
      Scalar globaly[3] = { q[1], q[4], q[7] };
      Scalar globalz[3] = { q[2], q[5], q[8] };

      Scalar mass;

      ele.andesms(1, globalx, globaly, globalz, (Scalar*) NULL, (Scalar*) NULL, mass, rhoh);
 
      return mass;
    }
};

#endif
#endif
