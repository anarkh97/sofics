// Sriramajayam

/*
 * IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h
 * DG++
 *
 * Created by Ramsharan Rangarajan on 11/07/2010.
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


#ifndef ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESSMATERIAL
#define ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESSMATERIAL

#include <Material.d/ElastoPlasticPlaneStressMaterial.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

//! \brief Class for computing the constitutive response of a
//! isotropic linear elastic and J2 plastic material under a plane stress.

/**
   Since the material is assumed to be in a state of plane stress, it
   is convenient to adopt vector notations for stresses and strains
   instead of tensor notations. Namely,

   \f[\text{Any stress}~\sigma =
   \begin{bmatrix}
   \sigma_{11} \\ \sigma_{22} \\ \sigma_{12}
   \end{bmatrix},\f]

   \f[\text{Any strain}~\varepsilon =
   \begin{bmatrix}
   \epsilon_{11}\\ \epsilon_{22} \\ 2\epsilon_{12}
   \end{bmatrix}.\f]

   The continuum model is summarized below.

   1. Symmetric strain computed from displacements

   \f[ \varepsilon =
   \begin{bmatrix}u_{1,1}\\ u_{2,2} \\ u_{1,2}+u_{2,1}\end{bmatrix},
   \quad {\bf u} = \text{displacement field}.\f]

   2. Addititive decomposition of strain into elastic and plastic parts

   \f[ \varepsilon = \varepsilon^e + \varepsilon^p, \quad
   \quad \varepsilon^e = \text{elastic strain}, ~
   \varepsilon^p = \text{plastic strain}.\f]

   3. Calculation of Cauchy stress

   \f[ \sigma = {\bf C}\varepsilon^e, \f]
   \f[\text{where}~{\bf C} = \frac{E}{1-\nu^2}
   \begin{bmatrix}
   1 && \nu && 0 \\
   \nu && 1 && 0 \\
   0 && 0 && \frac{1-\nu}{2}
   \end{bmatrix}. \f]

   4. Kinematic hardening mechanism:

   \f[ \xi = \sigma-\sigma^b, \quad \sigma^b = \text{back stress}.  \f]

   5. Associative flow rule for plastic strain

   \f[\dot{\varepsilon}^p = \lambda {\bf P} \xi, \f]
   \f[ \text{where} ~ {\bf P} = \frac{1}{3} \begin{bmatrix}
   2 && -1 && 0 \\
   -1 && 2 && 0 \\
   0 && 0 && 6
   \end{bmatrix}\f]

   The matrix \f${\bf P}\f$ helps relate the plastic strain to the
   deviatoric stress.

   6. Associative flow rule for back stress

   \f[\dot{\sigma}^b = \frac{2}{3}\lambda H \xi, \f]

   7. Associative flow rule for equivalent plastic strain

   \f[ \dot{e}^p = \sqrt{\frac{2}{3}}\lambda \sqrt{\xi^T{\bf P} \xi}.\f]

   8. J2 Plastic flow with isotropic hardening

   \f[ f:= \sqrt{\xi^T{\bf P}\xi} - \sqrt{\frac{2}{3}}(\sigma_Y + Ke^p).\f]

   We have used
   \f[ E = \text{Young's modulus}, \f]
   \f[ \nu = \text{Poisson's ratio}, \f]
   \f[ H = \text{Kinematic hardening modulus}, \f]
   \f[ K = \text{Isotropic hardening modulus}, \f]
   \f[ \sigma_Y = \text{Yield stress in 1D}.\f]


   Next, the numerical algorithm is described.

   It is assumed that the state of the matrial is known at time
   \f$t_n\f$ as \f$\{\varepsilon^p_n, \sigma^b_n, e^p_n\}\f$.

   Given a new total strain \f$\varepsilon_{n+1}\f$ at time
   \f$t_{n+1}\f$, the objective is to find the new state of the
   system.

   Let \f$\Delta\lambda = (t_{n+1}-t_n)\lambda_{n+1}\f$. Implicit
   integration of the flow rules give

   \f[ \varepsilon^p_{n+1} = \varepsilon^p_n + \Delta\lambda {\bf
   P}\xi_{n+1},\f]

   \f[ \sigma^b_{n+1} = \sigma^b_n + \frac{2}{3}\Delta\lambda H \xi_{n+1}, \f]

   \f[ e^p_{n+1} = e^p_n + \sqrt{\frac{2}{3}}\Delta\lambda \sqrt{\xi_{n+1}^T{\bf P}\xi_{n+1}}. \f]

   1. Evaluate trial elastic state by freezing plastic flow.

   \f[ \varepsilon^e_{trial} = \varepsilon_{n+1}-\varepsilon^p_{n}, \f]

   \f[ \sigma_{trial} = {\bf C}\varepsilon^e_{trial}, \f]

   \f[ \sigma^b_{trial} = \sigma^b_n, \f]

   \f[ \xi_{trial} = \sigma_{trial}-\sigma^b_{trial}, \f]

   \f[ e^p_{trial} = e^p_n, \f]

   \f[ f_{trial} = \sqrt{\xi_{trial}^T{\bf P}\xi_{trial}} - \sqrt{\frac{2}{3}}(\sigma_Y+Ke^p_{trial}).\f]

   2. If \f$f_{trial}\leq 0\f$, then this step is purely
   elastic. There is no need to update plastic variables of the
   material. Return
   \f[ \sigma_{n+1} = \sigma_{trial}.\f]

   3. Otherwise, this step is elastoplastic. The consistency parameter
   \f$\Delta \lambda\f$ is non-zero and needs to be determined. The
   following calculations demonstrate how.

   \f[\sigma_{n+1}
   = {\bf C}(\varepsilon_{n+1}-\varepsilon^p_{n+1})
   = \sigma_{trial} - {\bf C}(\varepsilon^p_{n+1}-\varepsilon^p_{n})
   = \sigma_{trial} - \Delta\lambda{\bf C}{\bf P}\xi_{n+1}.
   \f]

   \f[\sigma^b_{n+1} = \sigma^b_{trial} + \frac{2}{3}\Delta\lambda H \xi_{n+1}.\f]

   \f[\text{So,} ~ \xi_{n+1}
   = \sigma_{n+1}-\sigma^b_{n+1}
   = \xi_{trial} - \left(\frac{2}{3}H\Delta\lambda {\bf I}+\Delta\lambda{\bf C}{\bf P}\right)\xi_{n+1}.\f]

   \f[\Rightarrow~
   \left(\left(1+\frac{2}{3}H\Delta\lambda\right){\bf I} + \Delta\lambda{\bf C}{\bf P}\right)\xi_{n+1} = \xi_{trial}.\f]

   \f[\text{i.e.,}~\xi_{n+1} := {\bf M}(\Delta\lambda)\xi_{trial}.\f]

   The yield condition at time \f$t_{n+1}\f$ has to be satisfied. So,

   \f[ \sqrt{\xi_{n+1}^T{\bf P}\xi_{n+1}} -
   \sqrt{\frac{2}{3}}\left(\sigma_Y+K\left(e^p_n+\sqrt{\frac{2}{3}}\Delta\lambda\bar{f}_{n+1}\right)\right) = 0.\f]

   The consistency parameter is computed by numerically solving the above
   equation after substituting the expression for \f$\xi_{n+1}\f$ in
   terms of \f$\Delta\lambda\f$.

   The bisection method is used here. This is mainly for
   simplicity. Note that it is always possible to bracket the yeild
   function. This is becuase \f$f>0\f$ at the the trial state
   corresponding to \f$\Delta\lambda=0\f$.

   On the otherhand, if the entire increment in strain is plastic,
   then \f$f\f$ has to be negative. Otherwise, there is no hope for a
   solution while theoretically, it is known that a unique solution
   exists. To compute this extreme, we can set

   \f[ e^p_{n+1} = e^p_n + \sqrt{\frac{2}{3}}\|\varepsilon_{n+1}-\varepsilon^p_n \|. \f]

   \f[ \text{Then,}~ \Delta\lambda_{max} = \frac{3}{2}\frac{e^p_{n+1}-e^p_n}{\sigma_Y+Ke^p_{n+1}}.\f]


   The main ASSUMPTIONS are:

   1. The problem is strain driven, i.e., the state of the system is
   updated given a new strain.

   2. There is no initial plastic strain.

   3. There is no initial backstress.

   4. The state of the material IS UPDATED unless explicitly specified
   otherwise.

   5. Since there is no rate dependence, there is no need to specify
   the time step anywhere.

*/

class IsotropicLinearElasticJ2PlasticPlaneStressMaterial : public ElastoPlasticPlaneStressMaterial
{
 public:

  //! Constructor
  //! \param iLambda Lamé constant for elastic response
  //! \param iMu Lamé constant for elastic response
  //! \param iSigmaY Yield stress in 1D
  //! \param iK Isotropic hardening modulus
  //! \param iH Kinematic hardening modulus
  //! \param iTol Tolerance for convergence of nonlinear solve
  //! \param iequivEPSplasticF Equivalent plastic strain at failure
  IsotropicLinearElasticJ2PlasticPlaneStressMaterial(double iLambda, double iMu,
                                                     double iSigmaY, double iK = 0.,
                                                     double iH = 0., double iTol = 1.0e-6,
                                                     double iequivEPSplasticF = std::numeric_limits<double>::infinity());

  //! Destructor
  ~IsotropicLinearElasticJ2PlasticPlaneStressMaterial();

  //! Copy constructor
  //! \param MatObj Object to be copied
  IsotropicLinearElasticJ2PlasticPlaneStressMaterial(const IsotropicLinearElasticJ2PlasticPlaneStressMaterial &MatObj);

  //! Cloning
  IsotropicLinearElasticJ2PlasticPlaneStressMaterial * Clone() const;

  //! Compute the elastoplastic constitutive response.
  //! Returns true if calculations went well and false otherwise.
  //! \param Fnp1 Input. Deformation gradient at new state of material. Size 9x1.
  //! \param CauchyStress Output. Has size 9x1.
  //! \param Cep Output. Algorithmic elastoplastic tangent. If requested, has size 81x1.
  //! \param UpdateFlag Input. Material state updated if true. Note that by default, material state is updated.
  //! \param dt Input. Time increment.
  bool ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
                                                std::vector<double> * CauchyStress,
                                                std::vector<double> * Cep = 0,
                                                const bool UpdateFlag = true,
                                                const double dt = 0.);

  //! Returns the plastic strain in material (3x1 vector)
  const std::vector<double> & GetMaterialPlasticStrain() const;

  //! Returns equivalent plastic strain in material
  double GetMaterialEquivalentPlasticStrain() const;

  //! Returns back stress in material (3x1 vector)
  const std::vector<double> & GetMaterialBackStress() const;

  //! Returns the Isotropic hardening modulus
  double GetIsotropicHardeningModulus() const;

  //! Returns the kinematic hardening modulus
  double GetKinematicHardeningModulus() const;

  //! Returns the yield stress in 1D tension test
  double GetYieldStressFromTensionTest() const;

  //! Returns the bulk modulus of material
  double GetBulkModulus() const;

  //! Returns shear modulus of material
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

  //! Set the (x,y) values at a point in the yield stress vs. effective plastic strain curve
  //! \param iEqPlasticStrain Input. Equivalent plastic strain at a point.
  //! \param iYieldStress Input. Yield stress at a point.
  void SetExperimentalCurveData(double iEqPlasticStrain, double iYieldStress);

  //! Set the (x,y) values at a point in the yield stress scale factor vs. effective plastic strain rate curve
  //! \param iEqPlasticStrainRate Input. Equivalent plastic strain rate at a point.
  //! \param iYieldStressScaleFactor Input. Yield stress at a point.
  void SetExperimentalCurveData2(double iEqPlasticStrainRate, double iYieldStressScaleFactor);

  //! Print the internal variables
  void Print();

  //! Checks if the state of the material lies within the yield surface.
  //! \param CS Input. Cauchy stress 9x1 vector.
  //! \param TOL Input. Tolerance to use for check.
  //! \param dt Input. Time increment.
  //! The tolerance is non-dimensional. The check performed is
  //! \f[\frac{f}{\sigma_Y}<TOL~\Rightarrow~\text{material state OK}. \f]
  bool CheckMaterialState(const std::vector<double> &CS, const double TOL = 1.e-6,
                          const double dt = 0.) const;

 protected:

  //! Compute the elastic constitutive response of material
  //! Returns true if calculations went well and false otherwise
  //! \param EPS Input. Elastic strain 3x1.
  //! \param CS Output. Computed Cauchy stress 3x1.
  //! \param C Output. Elastic modulii 3x3.
  bool ComputeElasticConstitutiveResponse(const std::vector<double> &EPS,
                                          std::vector<double> *CS,
                                          std::vector<double> *C=0) const;

  //! Returns the yield stress by interpolating the experimental stress-strain curve
  //! \param eqP Input. Equivalent plastic strain.
  //! \param K Output. Isotropic hardening modulus.
  double GetYieldStressUsingExperimentalCurve(const double eqP, double &K) const;

  //! Returns the yield stress scale factor by interpolating the experimental curve
  //! \param eqPdot Input. Equivalent plastic strain time derivative.
  //! \param R Output. Derivative of the scale factor w.r.t. eqPdot.
  double GetScaleFactorUsingExperimentalCurve(const double eqPdot, double &R) const;

  //! Evaluates the yield function
  //! \param Xi Input. \f$\xi = \sigma-\sigma^b\f$. Size 3x1.
  //! \param DeltaEqP Input. Equivalent plastic strain increment.
  //! \param K Output. Isotropic hardening modulus.
  //! \param dt Input. Time increment.
  double EvaluateYieldFunction(const double * Xi, const double DeltaEqP, double &K, const double dt) const;

  //! Evaluates the norm of the deviatoric part of \f$\xi = \sigma-\sigma^b\f$.
  //! \param Xi Input. \f$\xi = \sigma-\sigma^b\f$, size 3x1.
  double ComputeJ2(const double *Xi) const;

  //! Compute \f$\xi_{n+1}\f$ given \f$\Delta \lambda\f$ and \f$\xi_{trial}\f$.
  //! Returns true if calculation was completed.
  //! \param Xitrial Input. \f$\xi_{trial}\f$, size 3x1
  //! \param lambda Input. \f$\Delta\lambda\f$.
  //! \param Xi Output. \f$\xi_{n+1}\f$ such that
  //! \f[\left(\left(1+\frac{2}{3}H\Delta\lambda\right){\bf I} + \Delta\lambda{\bf C}{\bf P}\right)\xi_{n+1} = \xi_{trial}.\f]
  bool ComputeXiGivenConsistencyParameter(const double * Xitrial, const double lambda, double * Xi) const;

 private:
  //! Young's modulus
  double E;

  //! Poisson ratio
  double nu;

  //! Yield stress in 1D
  double SigmaY;

  //! Isotropic hardening modulus
  double K;

  //! Kinematic hardening modulus
  double H;

  //! Tolerance for convergence of nonlinear solve
  double Tol;

  //! Equivalent plastic strain at failure
  double equivEPSplasticF;

  //! Plastic strain
  std::vector<double> EPSplastic;

  //! Back stress
  std::vector<double> BackStress;

  //! Equivalent plastic strain
  double equivEPSplastic;

  //! Pre-computed tensor coefficients
  double t00, t01, t22;

  //! Plastic strain values in experimental stress-strain curve
  std::vector<double> ExpEqPlasticStrain;

  //! Yield stress values in experimental stress-strain curve
  std::vector<double> ExpYieldStress;

  //! Plastic strain rate values in experimental stress-strain scaling curve
  std::vector<double> ExpEqPlasticStrainRate;

  //! Yield stress scaling values in experimental stress-strain scaling curve
  std::vector<double> ExpYieldStressScale;
};

// Return plastic strain
inline const std::vector<double>& IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetMaterialPlasticStrain() const
{
  return EPSplastic;
}

// Return equivalent plastic strain
inline double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetMaterialEquivalentPlasticStrain() const
{ return equivEPSplastic; }

// Return back stress
inline const std::vector<double>& IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetMaterialBackStress() const
{
  return BackStress;
}


#endif
