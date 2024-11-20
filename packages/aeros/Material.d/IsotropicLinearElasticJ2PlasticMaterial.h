// Sriramajayam

/*
 * IsotropicLinearElasticJ2PlasticMaterial.h
 * DG++
 *
 * Created by Ramsharan Rangarajan on 09/20/2010.
 *
 * Copyright (c) 2006 Adrian Lew
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#ifndef ISOTROPICLINEARELASTICJ2PLASTICMATERIAL
#define ISOTROPICLINEARELASTICJ2PLASTICMATERIAL

//! \brief Class for material with isotropic linear elastic and J2
//! plastic material.


/** This is the class for a rate independent elastic-plastic material
    with isotropic linear elastic response and J2 plastic response
    with possible isotropic and/or kinematic hardening.

    The continuum model for the constitutive response of the material
    is described below.

    1. The infinitesimal strain \f$\varepsilon\f$ is computed as
    \f[\varepsilon = \frac{\nabla u + {\nabla u}^T}{2},\f]
    where \f$u\f$ is the displacement field.

    2. An additive decomposition of the strain into an elastic and
    plastic part is assumed. \f[ \varepsilon = \varepsilon^e +
    \varepsilon^p.\f]

    3. The Cauchy stress is computed using the elastic response of the
    material as
    \f[\sigma = {\bf C}:\varepsilon^e, \f]
    where \f${\bf C}\f$ is the tensor of elastic modulii.

    4. The yield function is given by
    \f[ f(\sigma, {\bf \sigma}^b,e_p) = \|{\bf s}-{\bf \sigma}^b\| - \sqrt{\frac{2}{3}}(\sigma_Y + Ke_p),\f]
    \f[\mbox{where} \quad {\bf s} = \mbox{dev}(\sigma ) = \sigma - \frac{\mbox{tr}(\sigma )}{3}{\bf 1}\f]
    is the deviatoric part of the Cauchy stress.
    \f$\sigma^b\f$ and \f$e^p\f$ are plastic internal variables.
    Read explanations for these and remaining terms below.

    a. Isotropic hardening: The yield stress of the material in the
    above flow rule is \f$(\sigma_Y + Ke^p)\f$, where \f$\sigma_Y\f$
    is the flow stress of the material in pure tension and \f$K\f$ is
    the isotropic hardening modulus. This strain hardening mechanism
    is understood to control the radius of the yield surface as a
    function of the effective plastic strain \f$e^p\f$, where the
    effective plastic strain rate is defined as \f$\dot{e}_p =
    \sqrt{\frac{2}{3}}\|\dot{\varepsilon}^p\|\f$. Therefore, the
    effective plastic strain is computed as
    \f[ e_p = e_p^0 + \int_0^t\sqrt{\frac{2}{3}}\|\dot{\varepsilon}^p\|, \f]
    where \f$e_p^0\f$ is the initial effective plastic strain.

    b. Kinematic hardening: The center of the yield surface is given
    by the back stress \f$\sigma^b\f$, which is a deviatoric
    tensor. The evolution of the back stress is governed by a flow
    rule which will be described next.

    5. An associative flow rule is assumed. Define \f$\xi := {\bf
    s}-\sigma^b.\f$ and \f${\bf n}:=\frac{\xi}{\|\xi\|}\f$.

    The evolution equation for the plastic strain is given by
    \f[ \dot{\varepsilon}^p = \lambda \frac{\partial f}{\partial \sigma} = \lambda \frac{\xi}{\|\xi\|} = \lambda {\bf n},\f]
    where \f$\lambda \f$ is the consistency parameter (to ensure plastic flow on the yield surface).

    The evolution of the effective plastic strain follows as
    \f[\dot{e}^p = \sqrt{\frac{2}{3}}\|\dot{\varepsilon}^p\| = \sqrt{\frac{2}{3}}\lambda.\f]

    The evolution of the back stress, also associative, is given by
    \f[\dot{\sigma}^b = -{\frac{2}{3}}H\frac{\partial f}{\partial {\bf \sigma}_b}
    = {\frac{2}{3}}H\frac{\partial f}{\partial \xi} = {\frac{2}{3}}H\lambda {\bf n},\f]
    where \f$H\f$ is a constant kinematic hardening modulus. The factor of \f$2/3\f$
    is quite arbitrary though standard in the literature.

    6. The Kuhn-Tucker loading/unloading conditions are given by
    \f$\lambda \geq 0\f$, \f$f \leq 0\f$ and \f$\lambda f=0\f$.

    7. The consitency condition to determine the parameter \f$\lambda
    \f$ is given by \f$\lambda \dot{f}=0\f$.

    Next, the numerical integration algorithm is described.

    Let the state of the system be known at time step \f$t_n\f$ be
    known as \f$\{e_p^n, \varepsilon^p_n, {\bf \sigma}^b_n\}\f$.

    Given the strain at the next time step \f$t_{n+1}\f$ as
    \f$\varepsilon_{n+1}\f$, the objective is to determine the new
    state of the system \f$\{e_p^{n+1},\varepsilon^p_{n+1}, {\bf
    \sigma}^b_{n+1}\} \f$ as well as the Cauchy stress at this state.

    With \f$\Delta \lambda = \Delta t \lambda_{n+1}\f$,
    implicit Euler integration of the flow rule gives
    \f[\varepsilon^p_{n+1} = \varepsilon^p_n + \Delta \lambda {\bf n}_{n+1},\f]
    \f[e^p_{n+1} = e^p_n + \sqrt{\frac{2}{3}}\Delta \lambda, \f]
    \f[ {\bf \sigma}^b_{n+1} = {\bf \sigma}^b_n + {\frac{2}{3}}\Delta \lambda H {\bf n}_{n+1}. \f]

    1. Evaluate a trial state by freezing the plastic flow.
    \f[\varepsilon^p_{trial} = \varepsilon^p_n, \f]
    \f[e^p_{trial} = e^p_n,\f]
    \f[\varepsilon^e_{trial} = \varepsilon_{n+1}-\varepsilon^p_{trial}, \f]
    \f[\sigma^b_{trial} = \sigma^b_{n}, \f]
    \f[\sigma_{trial} = {\bf C}:\varepsilon^e_{trial}, \f]
    \f[{\bf s}_{trial} = \mbox{dev}(\sigma_{trial}),\f]
    \f[\xi_{trial} = {\bf s}_{trial}-{\bf \sigma}^b_{trial}, \f]
    \f[{\bf n}_{trial} = \frac{\xi_{trial}}{\|\xi_{trial}\|}.\f]

    2. Evaluate the yield function at the trial state as
    \f[f_{trial} = f(\sigma_{trial}, \sigma^b_{trial}, e^p_{trial}). \f]

    3. If \f$f_{trial}\leq 0\f$, then this strain increment step is
    purely elastic. Set
    \f[\varepsilon^p_{n+1}=\varepsilon^p_{trial},e^p_{n+1}=e^p_n,~ {\bf \sigma}^b_{n+1}=\sigma^b_{trial},
    ~\sigma = \sigma_{trial},~{\bf C}^{ep} = {\bf C}.\f]

    4. If not, this step is elastoplastic and the final state of the
    material is different from the trial state.

    To determine the final state, we need to compute the consistency
    parameter. The following calculations demonstrate how.

    \f[{\bf s}_{n+1} = \mbox{dev}(\sigma_{n+1}) = \mbox{dev}(\sigma_{trial}-{\bf C}:(\varepsilon^p_{n+1}-\varepsilon^p_n)) = {\bf s}_{trial}-2\mu
    \Delta\lambda {\bf n}_{n+1}.\f]

    \f[{\bf \sigma}^b_{n+1} = {\bf \sigma}^b_n + {\frac{2}{3}}\Delta \lambda H {\bf n}_{n+1}
    = {\bf \sigma}^b_{trial} + {\frac{2}{3}}\Delta \lambda H {\bf n}_{n+1}.\f]

    \f[\mbox{Therefore,} \quad \xi_{n+1} = {\bf s}_{n+1}-{\bf \sigma}^b_{n+1} = \xi_{trial}-\Delta \lambda(2\mu+{\frac{2}{3}}H){\bf n}_{n+1}.\f]

    This implies that \f$\xi_{n+1}\f$ and \f$\xi_{trial}\f$ are
    collinear, i.e., \f${\bf n}_{n+1} = {\bf n}_{trial}\f$.  Hence
    \f${\bf n}_{n+1}\f$ is known.

    From the above equation for \f$\xi_{n+1}\f$, we get
    \f[\|\xi_{n+1}\| = \xi_{n+1}\cdot {\bf n}_{n+1} = \|\xi_{trial}\|-\Delta \lambda(2\mu+{\frac{2}{3}}H).\f]

    At the new state of the material, we should have
    \f[f_{n+1} = f(\sigma_{n+1}, {\bf \sigma}^b_{n+1}, e^p_{n+1})=0.\f]

    \f[\mbox{i.e.,}~f_{n+1} = \|\xi_{n+1}\| - \sqrt{\frac{2}{3}}\left(\sigma_Y+K(e^p_n+\sqrt{\frac{2}{3}}\Delta \lambda) \right) = 0. \f]

    \f[\mbox{i.e.,}~\|\xi_{trial}\|-\Delta \lambda(2\mu+{\frac{2}{3}}H)- \sqrt{\frac{2}{3}}\left(\sigma_Y+K\left(e^p_n+\sqrt{\frac{2}{3}}\Delta \lambda\right)\right)=0.\f]

    \f[\mbox{i.e.,}~\|\xi_{trial}\|-\sqrt{\frac{2}{3}}\left(\sigma_Y+Ke^p_n\right) -\Delta \lambda \left(2\mu + {\frac{2}{3}}\left(H+K\right)\right) = 0.\f]

    \f[\mbox{i.e.,}~f_{trial}- \Delta \lambda \left(2\mu + \frac{2}{3}\left(H+K\right)\right) = 0.\f]

    \f[\Rightarrow ~ \Delta \lambda = \frac{f_{trial}}{2\mu +\frac{2}{3}\left(H +K\right)}. \f]

    Determine the new state of the material as
    \f[\varepsilon^p_{n+1} = \varepsilon^p_n + \Delta \lambda {\bf n}_{n+1}, \f]
    \f[e^p_{n+1} = e^p_n + \sqrt{\frac{2}{3}}\Delta \lambda,\f]
    \f[{\bf \sigma}^b_{n+1} = {\bf \sigma}^b_n + \frac{2}{3}\Delta \lambda H {\bf n}_{n+1}, \f]
    \f[ \sigma = {\bf C}:(\varepsilon_{n+1}-\varepsilon^p_{n+1}).\f]

    Finally, the consistent elasto-plastic modulii can be computed as

    \f[{\bf C}^{ep} =  \kappa {\bf 1}\otimes {\bf 1} + 2\mu \theta_{n+1} \left({\bf I}-\frac{1}{3}{\bf 1}\otimes {\bf 1}\right)
    -2\mu \overline{\theta}_{n+1} {\bf n}_{n+1}\otimes {\bf n}_{n+1},\f]

    \f[\mbox{where}~ \theta_{n+1} = 1- \frac{2\mu\Delta \lambda}{\|\xi_{trial}\|}, \f]
    \f[\mbox{and}~ \overline{\theta}_{n+1} = \frac{1}{1+\frac{K+H}{3\mu}}-1+\theta_{n+1}.\f]


    THE MAIN ASSUMPTIONS MADE IN THE CLASS ARE:

    a. The problem is strain driven. That is, the state of the
    material is updated based on given strain increments.

    b. The material has no initial plastic strain.

    c. The material has no initial back stress.

    d. The state of the material IS UPDATED during constitutive
    response calculations in explicitly specified otherwise.

    e. The tangents returned are the consistent elasto-plastic modulii.

    f. The material is rate independent. Therefore, there is no
    explicit need to specify the time step anywhere.
 */


#include "Material.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <limits>

class IsotropicLinearElasticJ2PlasticMaterial
{
 public:

  //! Constructor
  //! \param iLambda Lamé constant for linear elastic material
  //! \param iMu Lamé constant for linear elastic material
  //! \param iSigmaY Flow stress from uniaxial tension test
  //! \param iK Isotropic hardening modulus, zero by default
  //! \param iH Kinematic hardening modulus, zero by default
  //! \param iTol Tolerance for convergence of nonlinear solve, 1.0e-6 by default
  //! \param ik Strength coefficient for generalized power-law hardening, zero by default
  //! \param in Hardening exponent for generalized power-law hardening, zero by default
  IsotropicLinearElasticJ2PlasticMaterial(double iLambda, double iMu,
                                          double iSigmaY, double iK=0., double iH=0.,
                                          double iTol = 1.0e-6,
                                          double iequivEPSplasticF = std::numeric_limits<double>::infinity(),
                                          double ik=0., double in=0.);

  //! Destructor
  ~IsotropicLinearElasticJ2PlasticMaterial();

  //! Copy constructor
  //! \param Mat Input. Object to be copied.
  IsotropicLinearElasticJ2PlasticMaterial(const IsotropicLinearElasticJ2PlasticMaterial &Mat);

  //! Cloning mechanism
  IsotropicLinearElasticJ2PlasticMaterial * Clone() const;

  //! Compute the elasto-plastic response of material
  //! Returns true if calculations went well.
  //! \param Fnp1 Input. Deformation gradient at new time step.
  //! \param CauchyStress Output. Computed Cauchy stress.
  //! \param Cep Output. Elasto-plastic tangent, computed if requested.
  //! \param UpdateFlag. Input. State of material is updated if set to true. Defaulted to true.
  //! \param dt Input. Time increment.
  bool ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
                                                std::vector<double> *CauchyStress,
                                                std::vector<double> *Cep = 0,
                                                const bool UpdateFlag = true,
                                                double dt = 0.);

  //! Return the plastic strain in the material
  std::vector<double> GetMaterialPlasticStrain() const;

  //! Return the equivalent plastic strain in the material
  double GetMaterialEquivalentPlasticStrain() const;

  //! Return the back stress in the material
  std::vector<double> GetMaterialBackStress() const;

  //! Return the isotropic hardening modulus
  double GetIsotropicHardeningModulus() const;

  //! Return the kinematic hardening modulus
  double GetKinematicHardeningModulus() const;

  //! Return the flow stress from a uniaxial tension test
  double GetYieldStressFromTensionTest() const;

  //! Return bulk modulus of material
  double GetBulkModulus() const;

  //! Return shear modulus of material
  double GetShearModulus() const;

  //! Returns dissipated energy in material
  double GetDissipatedEnergy() const;

  //! Returns the equivalent plastic strain at failure
  double GetEquivalentPlasticStrainAtFailure() const;

  //! Returns the tolerance for convergence of nonlinear solve
  double GetTolerance() const;

  //! Set the plastic strain in the material
  //! \param EPSplastic Input. Plastic strain.
  void SetMaterialPlasticStrain(const std::vector<double> &EPSplastic);

  //! Set the equivalent plastic strain in the material
  //! \param equivEPSplastic Input. Equivalent plastic strain.
  void SetMaterialEquivalentPlasticStrain(double equivEPSplastic);

  //! Set the back stress in the material
  //! \param BackStress Input. Back stress.
  void SetMaterialBackStress(const std::vector<double> &BackStress);

  //! Set the stress and strain at a point in the stress-strain curve
  //! \param iEqPlasticStrain Input. Equivalent plastic straint at a point.
  //! \param iYieldStress Input. Yield stress at a point.
  void SetExperimentalCurveData(double iEqPlasticStrain, double iYieldStress);

  //! Set the data at a point in the yield stress scale factor vs. effective plastic strain rate curve
  //! \param iEqPlasticStrainRate Input. Equivalent plastic strain rate at a point.
  //! \param iYieldStressScaleFactor Input. Yield stress at a point.
  void SetExperimentalCurveData2(double iEqPlasticStrainRate, double iYieldStressScaleFactor);

  //! Checks if state of the material lies on or within yield surface.
  //! \param CS Input. Cauchy stress.
  //! \param TOL Input. Tolerance for comparing yield function with zero. Defaulted to 1e-6.
  //! \param dt Input. Time increment. Defaulted to 0.
  //! Returns true if f<=TOL and false otherwise.
  bool CheckMaterialState(const std::vector<double> &CS, const double TOL=1.e-6,
                          const double dt = 0.) const;

 private:

  //! Compute the linear elastic response
  //! Returns true if calculations went well.
  //! \param EPS Input. Symmetric infinitesimal strain tensor, NOT DEFORMATION GRADIENT.
  //! \param CS Output. Cauchy stress.
  //! \param C Output. Elastic modulii.
  bool ComputeElasticConstitutiveResponse(const std::vector<double> &EPS,
                                          std::vector<double> *CS,
                                          std::vector<double> *C=0) const;

  //! Returns the yield stress using the generalized power law
  //! \param eqP Input. Equivalent plastic strain.
  //! \param K Output. Isotropic hardening modulus.
  double GetYieldStressUsingGeneralizedPowerLaw(const double eqP, double &K) const;

  //! Returns the yield stress by interpolating the experimental stress-strain curve
  //! \param eqP Input. Equivalent plastic strain.
  //! \param K Output. Isotropic hardening modulus.
  double GetYieldStressUsingExperimentalCurve(const double eqP, double &K) const;

  //! Returns the yield stress scale factor by interpolating the experimental curve
  //! \param eqPdot Input. Equivalent plastic strain time derivative.
  //! \param R Output. Derivative of the scale factor w.r.t. eqPdot.
  double GetScaleFactorUsingExperimentalCurve(const double eqPdot, double &R) const;

  //! Evaluate the yield function
  //! \param CauchyStress Input. Cauchy stress.
  //! \param SigmaB Input. Back stress.
  //! \param eqP Input. Equivalent plastic strain increment.
  //! \param K Output. Isotropic hardening modulus.
  //! \param dt Input. Time increment.
  double EvaluateYieldFunction(const double * CauchyStress,
                               const double * SigmaB,
                               const double DeltaEqP,
                               double &K,
                               const double dt) const;

  //! Compute the deviatoric part of a 3x3 tensor
  //! \param T Input. Should have size 9.
  std::vector<double> Deviatoric(const double * T) const;

  //! Compute the norm of a 3x3 tensor
  //! \param T Input. Should have size 9.
  double Norm(const double * T) const;

  //! Bulk modulus
  double Kappa;

  //! Shear modulus
  double Mu;

  //! Flow stress in uniaxial tension for linear hardening / constant term for generalized power-law hardening
  double SigmaY;

  //! Isotropic hardening modulus for linear hardening / linear term for generalized power-law hardening
  double K;

  //! Kinematic hardening modulus
  double H;

  //! Strength coefficient for generalized power-law hardening
  double k;

  //! Hardening exponent for generalized power-law hardening
  double n;

  //! Tolerance for convergence of nonlinear solve
  double Tol;

  //! Equivalent plastic strain at failure
  double equivEPSplasticF;

  //! Linear elastic material
  IsotropicLinearElastic * ILE;

  //! Plastic strain in material
  std::vector<double> EPSplastic;

  //! Equivalent plastic strain in material
  double equivEPSplastic;

  //! Back stress of material
  std::vector<double> BackStress;

  //! Plastic strain values in experimental stress-strain curve
  std::vector<double> ExpEqPlasticStrain;

  //! Yield stress values in experimental stress-strain curve
  std::vector<double> ExpYieldStress;

  //! Plastic strain rate values in experimental stress-strain scaling curve
  std::vector<double> ExpEqPlasticStrainRate;

  //! Yield stress scaling values in experimental stress-strain scaling curve
  std::vector<double> ExpYieldStressScale;
};


#endif
