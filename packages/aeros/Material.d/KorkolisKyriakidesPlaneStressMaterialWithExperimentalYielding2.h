// Sriramajayam

/*
 * KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2.h
 * DG++
 *
 * Created by Ramsharan Rangarajan on 12/03/2010.
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


#ifndef KORKOLISKYRIAKIDESPLANESTRESSMATERIALWITHEXPERIMENTALYIELDING2
#define KORKOLISKYRIAKIDESPLANESTRESSMATERIALWITHEXPERIMENTALYIELDING2

#include <Material.d/ElastoPlasticPlaneStressMaterial.h>
#include <iostream>
#include <vector>
#include <cmath>


class KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2 : public ElastoPlasticPlaneStressMaterial
{
 public:

  //! Constructor
  //! \param iLambda Lame constant for elastic response (SI units)
  //! \param iMu Lame constant for elastic response (SI units)
  KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2(double iLambda, double iMu, double iTol1 = 1.0e-6, double iTol2 = 1.0e-6);

  //! Destructor
  virtual ~KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2();

  //! Copy constructor
  //! \param MatObj Object to be copied
  KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2
    (const KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2 &MatObj);

  //! Cloning
  virtual KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2 * Clone() const;

  //! Compute the elastoplastic constitutive response.
  //! Returns true if calculations went well and false otherwise.
  //! \param Fnp1 Input. Deformation gradient at new state of material. Size 9x1.
  //! \param CauchyStress Output. Has size 9x1.
  //! \param Cep Output. Algorithmic elastoplastic tangent. If requested, has size 81x1
  //! \param UpdateFlag Input. Material state updated if true. Set to true by default.
  //! \param dt Input. Time increment.
  virtual bool ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
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

  //! Returns the bulk modulus of material
  double GetBulkModulus() const;

  //! Returns shear modulus of material
  double GetShearModulus() const;

  //! Returns dissipated energy in material
  double GetDissipatedEnergy() const;

  //! Set the plastic strain in the material
  void SetMaterialPlasticStrain(const std::vector<double> &EPSplastic);

  //! Set the equivalent plastic strain in the material
  void SetMaterialEquivalentPlasticStrain(double equivEPSplastic);

  //! Set the back stress in the material
  void SetMaterialBackStress(const std::vector<double> &BackStress);

  //! Checks if the state of the material lies within the yield surface.
  //! \param CS Input. Cauchy stress 9x1 vector
  //! \param TOL Input. Tolerance to use for check
  //! The tolerance is non-dimensional. The check performed is
  //! \f[\frac{f}{\sigma_Y}<TOL~\Rightarrow~\text{material state OK}. \f]
  virtual bool CheckMaterialState(const std::vector<double> &CS, const double TOL = 1.e-6) const;

 protected:

  // Compute the elastic constitutive response of material
  //! Returns true if calculations went well and false otherwise
  //! \param EPS Input. Elastic strain, 3x1
  //! \param CS Output. Computed cauchy stress, 3x1
  //! \param C Output. Elastic modulii 3x3
  virtual bool ComputeElasticConstitutiveResponse(const std::vector<double> &EPS,
                                                  std::vector<double> *CS,
                                                  std::vector<double> *C=0) const;

  //! Evaluates the yield function
  //! \param CS Input.  Cauchy stress. Size 3x1.
  //! \param eqP input. Equivalent plastic strain.
  virtual double EvaluateYieldFunction(const double * CS, const double eqP) const;

 private:
  //! Evaluate derivative of yield function with respect to
  //! Cauchy stress given \f$\sigma, e^p\f$.
  //! The derivative is computed numerically with a central difference formula.
  //! The perturbation equals \f$1.e-6 \times \sigma_Y\f$.
  //! \param CS Input. Cauchy stress, size 3x1.
  //! \param eqP Input. \f$e^p\f$.
  std::vector<double> EvaluateDerivativeOfYieldFunction(const double * CS, const double eqP) const;

  //! Computes the two linear transformations of stress involved in
  //! this model.
  //! \param CS Input. Cauchy stress, size 3x1.
  //! \param L1CS Output. First transformation of \f$\sigma\f$, size 3x1.
  //! \param L2CS Output. Second transformation of \f$\sigma\f$, size 3x1.
  void ComputeLinearTransformationsOfStress(const double * CS,
                                            double * L1CS, double * L2CS) const;


  //! Returns the yield stress by interpolating the experimental stress-strain curve
  //! \param eqP Equivalent plastic strain
  double GetYieldStressUsingExperimentalCurve(const double eqP) const;

  //! Associative flow rule for \f$\sigma\f$, given
  //! \f$\Delta \lambda\f$ and \f${\bf N}\f$.
  //! \param CStrial Input. \f$\sigma_{trial}\f$, size 3x1
  //! \param lambda Input. Value of \f$\Delta \lambda\f$.
  //! \param N Input. Value of \f${\bf N}\f$, size 3x1
  //! \param CS output. Computed value of \f$\xi\f$.
  void ComputeCauchyStressGivenConsistencyParameterAndDirection(const double * CStrial, const double lambda,
                                                                const double * N, double * CS) const;

  //! Associative flow rule to compute equivalent plastic strain
  //! \f$e^p\f$ given \f$\delta \lambda\f$ and \f${\bf N}\f$.
  //! \param lambda Input. Value of \f$\Delta \lambda\f$.
  //! \param N Input. Value of \f${\bf N}\f$, size 3x1
  double ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(const double lambda,
                                                                             const double * N) const;

  //! Finds \f$\Delta\lambda\f$ such that \f$f=0\f$, given a direction \f${\bf N}\f$
  //! for the evolution of strain.
  //! Returns true if solution was found and false otherwise.
  //! \param CStrial Input. \f$\sigma_{trial}\f$, size 3x1
  //! \param N Input. Value of \f${\bf N}\f$, size 3x1
  //! \param lambda_L Input. Guess for \f$\Delta\lambda\f$ for which \f$f\leq 0\f$.
  //! \param lambda_R Input. Guess for \f$\Delta\lambda\f$ for which \f$f\geq 0\f$.
  //! \param lambda Output. Solved value of \f$\Delta\lambda\f$.
  bool ComputeConsistencyParameterGivenDirection(const double * CStrial, const double * N,
                                                 double lambda_L, double lambda_R, double &lambda) const;

  //! Computes the two eigenvalues of a symmetric 2x2 tensr given in vector form.
  //! \param X Input. Vector form of symmetric tensor, size 3x1
  //! \param A Output. First eigenvalue
  //! \param B Output. Second eigenvalue
  void ComputeEigenvalues(const double * X, double &A, double &B) const;

  //! Computes the tensor norm of a deviatoric stress tensor
  //! given in vector form
  //! \param S Input. Symmetric deviator stress tensor, size 3x1.
  double DeviatoricStressNorm(const double * S) const;

  //! Computes the tensor norm of a deviatoric strain tensor
  //! given in vector form
  //! \param S Input. Symmetric deviator strain tensor, size 3x1.
  double DeviatoricStrainNorm(const double * S) const;

  //! Young's modulus
  double E;

  //! Poisson ratio
  double nu;

  //! Tolerances for convergence of nonlinear solve
  double Tol1, Tol2;

  //! Plastic strain
  std::vector<double> EPSplastic;

  //! Back stress
  std::vector<double> BackStress;

  //! Equivalent plastic strain
  double equivEPSplastic;

  //! Plastic strain values in experimental stress-strain curve
  static const double ExpEqPlasticStrain[];

  //! Yield stress values in experimental stress-strain curve
  static const double ExpYieldStress[];

};

#endif
