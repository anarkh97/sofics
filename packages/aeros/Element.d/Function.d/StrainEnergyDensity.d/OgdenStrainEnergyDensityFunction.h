#ifndef _OGDENSTRAINENERGYDENSITYFUNCTION_H_
#define _OGDENSTRAINENERGYDENSITYFUNCTION_H_

#include <limits>
#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

namespace Simo {

template<typename Scalar>
class OgdenStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    Eigen::Matrix<double,3,1> mu, alpha;

  public:
    OgdenStrainEnergyDensityFunction(double mu1, double mu2, double mu3, double alpha1, double alpha2, double alpha3)
    {
      mu << mu1, mu2, mu3;
      alpha << alpha1, alpha2, alpha3;
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const 
    {
      // inputs: deformation gradient
      // output: strain energy density
#ifdef HACK_PJSA_EIGEN3 
      using std::pow;

      Eigen::Matrix<Scalar,3,3> C = F.transpose()*F; // right Cauchy-Green strain
      Eigen::Matrix<Scalar,3,1> lambda; // principal stretches

      if(F.isIdentity(100*std::numeric_limits<double>::epsilon())) {
        // this is a workaround because eigenvalues are not differentiable when F is the identity.
        lambda = C.diagonal().cwiseSqrt();
      }
      else {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar,3,3> > es(C);
        lambda = es.eigenvalues().cwiseSqrt();
      }

      Scalar W = 0;
      for(int i = 0; i < 3; ++i) {
        if(mu[i] != 0) W += mu[i]/alpha[i]*(pow(lambda[0],alpha[i])+pow(lambda[1],alpha[i])+pow(lambda[2],alpha[i])-3);
      }
      return W;
#else
      return Scalar(0);
#endif
    }

};

} // namespace Simo

#endif
