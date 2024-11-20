// Sriramajayam

/*
 * KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding.cpp
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


#include "KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding.h"

// Constructor
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding(double iLambda, double iMu, double iTol1, double iTol2)
{
  // Youngs modulus
  E = iMu*(3.*iLambda + 2.*iMu)/(iLambda + iMu);

  // Poisson ratio
  nu = 0.5*iLambda/(iLambda+iMu);

  // Tolerances for convergence of nonlinear solve
  Tol1 = (iTol1 > 0) ? iTol1 : 1.0e-6;
  Tol2 = (iTol2 > 0) ? iTol2 : 1.0e-6;

  // Zero initial plastic strain and equivalent plastic strain
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<3; i++) {
    EPSplastic.push_back( 0. );
    BackStress.push_back( 0. );
  }
  equivEPSplastic = 0.;
}

// Destructor
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
~KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding()
{}

// Copy constructor
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding
(const KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding & Mat)
  : E(Mat.E), nu(Mat.nu), Tol1(Mat.Tol1), Tol2(Mat.Tol2)
{
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<3; i++) {
    EPSplastic.push_back( Mat.EPSplastic[i] );
    BackStress.push_back( Mat.BackStress[i] );
  }
  equivEPSplastic = Mat.equivEPSplastic;
}

// Cloning
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding *
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::Clone() const
{ return new KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding(*this); }


// Return plastic strain
const std::vector<double> & KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
GetMaterialPlasticStrain() const
{
/* PJSA
  std::vector<double> EP(9,0.);
  EP[0] = EPSplastic[0];
  EP[4] = EPSplastic[1];
  EP[8] = -(EP[0]+EP[4]);
  EP[1] = EP[3] = 0.5*EPSplastic[2];
  EP[2] = EP[5] = EP[6] = EP[7] = 0.;
  return EP;
*/
  return EPSplastic;
}

// Return back stress
const std::vector<double> & KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
GetMaterialBackStress() const
{
  return BackStress;
}

// Return equivalent plastic strain
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
GetMaterialEquivalentPlasticStrain() const
{ return equivEPSplastic; }

// Return Bulk modulus
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
GetBulkModulus() const
{ return E/(3.*(1.-2.*nu)); }

// Return shear modulus
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
GetShearModulus() const
{ return 0.5*E/(1.+nu); }

// Return dissipated energy
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
GetDissipatedEnergy() const
{
  std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::GetDissipatedEnergy()- "
           <<"not implemented.\n";
  return 0;
}

// Set the plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
SetMaterialPlasticStrain(const std::vector<double> &iEPSplastic)
{ for(int i = 0; i < 3; ++i) EPSplastic[i] = iEPSplastic[i]; }

// Set the plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
SetMaterialBackStress(const std::vector<double> &)
{}

// Set the equivalent plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
SetMaterialEquivalentPlasticStrain(double iEquivEPSplastic)
{ equivEPSplastic = iEquivEPSplastic; }

// Check if state of material lies within yield surface
bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
CheckMaterialState(const std::vector<double> &CS, const double TOL) const
{
  // Cauchy stress in vector form
  double CSvec[3] = { CS[0], CS[4], CS[1] };

  // Evaluate yield function
  double Fval = EvaluateYieldFunction(CSvec, equivEPSplastic);

  // Check for tolerance after non-dimensionalizing Fval
  double SigmaYScale = ExpYieldStress[0];
  if( Fval/SigmaYScale < TOL )
    return true;
  else
    return false;
}

// Compute the elastic constitutive response
bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ComputeElasticConstitutiveResponse(const std::vector<double> &EPS,
                                   std::vector<double> *CS,
                                   std::vector<double> *C) const
{
  if( int(CS->size())<3 )
    CS->resize(3);

  double v = E/(1.-nu*nu);
  double Cmat[3][3] = {{v*1., v*nu, 0.},
                       {v*nu, v, 0.},
                       {0., 0., v*(1.-nu)/2.}};

  for(int i=0; i<3; i++)
    {
      (*CS)[i] = 0.;
      for(int j=0; j<3; j++)
        (*CS)[i] += Cmat[i][j]*EPS[j];
    }

  if( C!=0 )
    {
      if( int(C->size())<9 )
        C->resize( 9 );
      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
          (*C)[3*i+j] = Cmat[i][j];
    }

  return true;
}

// Compute elastoplastic response of material
bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
                                         std::vector<double> *CauchyStress,
                                         std::vector<double> *Cep,
                                         const bool UpdateFlag,
                                         const double dt)
{
  // Notify that elastoplastic tangents are not computed.
  if( Cep!=0 )
    {
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::ComputeElastoPlasticConstitutiveResponse()- "
               <<"Elastoplastic tangents not implemented.\n";
      return false;
    }

  // Resize outputs if required
  if( int(CauchyStress->size())<9 )
    CauchyStress->resize(9);

  // Vector form of symmetric strain
  double EPS[3] = { Fnp1[0]-1., Fnp1[4]-1., Fnp1[1]+Fnp1[3] };

  // Compute trial elastic state by freezing plastic flow
  std::vector<double> EPSelas(3,0.);
  for(int i=0; i<3; i++)
    EPSelas[i] = EPS[i]-EPSplastic[i];

  // Compute trial Cauchy stress
  std::vector<double> CStrial(3,0.);
  if( !ComputeElasticConstitutiveResponse(EPSelas, &CStrial) )
    {
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::ComputeElastoPlasticConstitutiveResponse()- "
               <<"Could not compute elastic constitutive response at trial state.\n";
      return false;
    }

  // Evaluate yield function at trial state
  double Ftrial = EvaluateYieldFunction(&CStrial[0], equivEPSplastic);

  // Use some tolerance for checking yield function value
  double TOL = ExpYieldStress[0]*Tol2; // ORIG: ExpYieldStress[0]*1.e-6

  if( Ftrial<0 /*ORIG: TOL*/ )
    {
      // This step is purely elastic
      // No need to update plastic variables
      // Convert Cauchy stress to tensor form
      (*CauchyStress)[0] = CStrial[0];
      (*CauchyStress)[4] = CStrial[1];
      (*CauchyStress)[1] = (*CauchyStress)[3] = CStrial[2];
      (*CauchyStress)[2] = (*CauchyStress)[5] = (*CauchyStress)[6] = (*CauchyStress)[7] = (*CauchyStress)[8] = 0.;

      return true;
    }

  else
    {
      // This step is elasto-plastic.
      // Need to compute consistency parameter, called lambda here.
      // Need to compute direction of increment, called N here.

      // Initial guess for N = df/d(sigma) at (CStrial, eqPtrial)
      std::vector<double> Ntrial = EvaluateDerivativeOfYieldFunction(&CStrial[0], equivEPSplastic);
      double normNtrial = DeviatoricStressNorm(&Ntrial[0]);
      if( normNtrial < 1.e-6 )
        {
          std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::ComputeElastoPlasticConstitutiveResponse()- "
                   <<"Normal to yield surface has close to zero magnitude.\n";
          return false;
        }

      // Maximum permitted value for lambda- useful in solution
      // Assume that the entire increment in strain is plastic.
      double lambdaMax = DeviatoricStrainNorm(&EPSelas[0])/normNtrial;

      // Iterate to find lambda and N such that
      // f(lambda, N) = 0 and df/d(sigma) at (lambda, N) = N.

      // Solution (lambda, N)
      double lambda = 0.;
      std::vector<double> N(3,0.);

      // Initialize N to Ntrial
      N = Ntrial;

      // Max. number of iterations
      int MaxIt = 1000;

      // Flag to check for convergence
      bool CONVERGED = false;

      for(int it=0; it<MaxIt && !CONVERGED; it++)
        {
          // Find lambda such that f = 0.
          if( !ComputeConsistencyParameterGivenDirection(&CStrial[0], &N[0], lambdaMax, 0., lambda) )
            {
              std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::ComputeElastoPlasticConstitutiveResponse()- "
                       <<"Could not compute consistency parameter for given direction.\n";
              return false;
            }

          // Update guess for Cauchy stress
          double CSUpdate[3] = {0.,0.,0.};
          ComputeCauchyStressGivenConsistencyParameterAndDirection(&CStrial[0], lambda, &N[0], CSUpdate);

          // Update guess for equivalent plastic strain
          double eqPupdate = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda, &N[0]);

          // Update guess for normal direction
          N = EvaluateDerivativeOfYieldFunction(CSUpdate, eqPupdate);

          // Update Cauchy stress and eqP with this normal direction
          ComputeCauchyStressGivenConsistencyParameterAndDirection(&CStrial[0], lambda, &N[0], CSUpdate);
          eqPupdate = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda, &N[0]);

          // Check if this updated state is a solution
          double Fval = EvaluateYieldFunction(CSUpdate, eqPupdate);
          if(  fabs(Fval) < TOL )
            // This (lambda, N) is a solution
            CONVERGED = true;
        }

      // Check if solution converged
      if( !CONVERGED )
        {
          std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::ComputeElastoPlasticConstitutiveResponse()- "
                   <<"Could not solve for new state of system.\n";
          return false;
        }

      // Compute actual elastic strain = total strain - plastic strain
      for(int i=0; i<3; i++)
        EPSelas[i] = EPS[i]-EPSplastic[i] - lambda*N[i];

      // Compute Cauchy stress
      std::vector<double> CS(3,0.);
      ComputeElasticConstitutiveResponse(EPSelas, &CS);

      // Convert to tensor form
      (*CauchyStress)[0] = CS[0];
      (*CauchyStress)[4] = CS[1];
      (*CauchyStress)[1] = (*CauchyStress)[3] = CS[2];
      (*CauchyStress)[2] = (*CauchyStress)[5] = (*CauchyStress)[6] = (*CauchyStress)[7] = (*CauchyStress)[8] = 0.;

      // Update state of system unles specified otherwise
      if( UpdateFlag )
        {
          // Update plastic strain
          for(int i=0; i<3; i++)
            EPSplastic[i] += lambda*N[i];

          // Update equivalent plastic strain
          equivEPSplastic += sqrt(2./3.)*lambda*DeviatoricStressNorm(&N[0]);
        }
    }
  return true;
}

// Evaluate yield function
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
EvaluateYieldFunction(const double * CS, const double eqP) const
{
  // Compute two linear transformations of Cauchy stress
  double L1CS[]  = {0.,0.,0.};
  double L2CS[] = {0.,0.,0.};
  ComputeLinearTransformationsOfStress(CS, L1CS, L2CS);

  // A scale for yield stress
  double SigmaY = ExpYieldStress[0];

  // Compute L1CS/SigmaY and L2CS/SigmaY (avoid large numbers)
  double L1hat[3], L2hat[3];
  for(int i=0; i<3; i++)
    {
      L1hat[i] = L1CS[i]/SigmaY;
      L2hat[i] = L2CS[i]/SigmaY;
    }

  // Compute eigenvalues of L1hat
  double a=0., b=0.;
  ComputeEigenvalues(L1hat, a, b);

  // Compute eigenvalues of L2hat
  double A=0., B=0.;
  ComputeEigenvalues(L2hat, A, B);

  // Compute effective stress
  double SigmaHat8 = (pow(a-b,8.)+pow(2.*A+B,8.)+pow(A+2.*B,8.))/2.;
  double EffStress = SigmaY*pow(SigmaHat8,0.125);

  // Compute effective yield stress
  double EffYieldStress = GetYieldStressUsingExperimentalCurve(eqP);

  return EffStress - EffYieldStress;
}

// Evaluate derivative of yield function with respect to Cauchy stress
std::vector<double> KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
EvaluateDerivativeOfYieldFunction(const double * CS, const double eqP) const
{
  // Derivative wrt Cauchy stress
  std::vector<double> N(3,0.);

  // Scale for yield stress
  double SigmaY = ExpYieldStress[0];

  // How much to perturb
  double deltaSigma = (1.e-6)*SigmaY;

  double MyCS[] = {CS[0], CS[1], CS[2]};

  for(int i=0; i<3; i++)
    {
      // Plus direction
      MyCS[i] = CS[i]+deltaSigma;
      double Fplus = EvaluateYieldFunction(MyCS, eqP);

      // Minus direction
      MyCS[i] = CS[i]-deltaSigma; // ORIG: CS[i]-2.*deltaSigma;
      double Fminus = EvaluateYieldFunction(MyCS, eqP);

      // Compute derivative
      N[i] = (Fplus-Fminus)/(2.*deltaSigma);

      // Restore component of CS.
      MyCS[i] = CS[i];
    }
  return N;
}

// Linear transformations of stress involved in this model
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ComputeLinearTransformationsOfStress(const double * CS, double * L1CS, double * L2CS) const
{
  // Model parameters
  double a1 = 0.78, a2 = 1.15, a3 = 0.85, a4 = 0.89, a5 = 1.06, a6 = 1.03, a7 = 1.0, a8 = 1.0;

  // Compute L1CS
  L1CS[0] = (a1/3.)*(2.*CS[0]-CS[1]);
  L1CS[1] = (a2/3.)*(2.*CS[1]-CS[0]);
  L1CS[2] = a7*CS[2];

  // Entries for transformation L2
  double Lvec[] = { (8.*a5 - 2.*a3 - 2.*a6 + 2.*a4)/9.,
                    (4.*a6 - 4.*a4 - 4.*a5 + 1.*a3)/9.,
                    (4.*a3 - 4.*a5 - 4.*a4 + 1.*a6)/9.,
                    (8.*a4 - 2.*a6 - 2.*a3 + 2.*a5)/9.,
                    a8 };

  // Compute L2CS
  L2CS[0] = Lvec[0]*CS[0]+Lvec[1]*CS[1];
  L2CS[1] = Lvec[2]*CS[0]+Lvec[3]*CS[1];
  L2CS[2] = Lvec[4]*CS[2];
}

// Compute yield stress by interpolating experimental stress-strain curve
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
GetYieldStressUsingExperimentalCurve(const double eqP) const
{
  double SigmaY = 0.;
  if( eqP<=ExpEqPlasticStrain[0] )
    SigmaY = ExpYieldStress[0];
  else if ( eqP>=ExpEqPlasticStrain[30] )
    SigmaY = ExpYieldStress[30];
  else
    {
      for(int i=0; i<30; i++)
        if( eqP>=ExpEqPlasticStrain[i] && eqP<=ExpEqPlasticStrain[i+1] )
          {
            double lambda = (ExpEqPlasticStrain[i+1]-eqP)/(ExpEqPlasticStrain[i+1]-ExpEqPlasticStrain[i]);
            SigmaY = lambda*ExpYieldStress[i] + (1.-lambda)*ExpYieldStress[i+1];
          }
    }
  return SigmaY;
}

// Compute \f$\sigma\f$ given \f$\sigma_{trial}, \lambda, {\bf N} \f$.
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ComputeCauchyStressGivenConsistencyParameterAndDirection(const double * CStrial, const double lambda,
                                                         const double * N, double * CS) const
{
  // Cauchy stress = CStrial - lambda*C*N.

  // Elastic modulii.
  double v = E/(1.-nu*nu);
  double Cmat[3][3] = {{v*1., v*nu, 0.},
                       {v*nu, v, 0.},
                       {0., 0., v*(1.-nu)/2.}};

  for(int i=0; i<3; i++)
    {
      CS[i] = CStrial[i];
      for(int j=0; j<3; j++)
        CS[i] -= lambda*Cmat[i][j]*N[j];
    }
}

// Compute equivalent plastic strain given lambda and N.
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(const double lambda, const double * N) const
{
  return equivEPSplastic + sqrt(2./3.)*lambda*DeviatoricStressNorm(N);
}

bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ComputeConsistencyParameterGivenDirection(const double * CStrial, const double * N,
                                          double lambda_L, double lambda_R,
                                          double & lambda) const
{
  // Assume that f>=0 for lambda_L and f<=0 at lambda_R.

  double CS[3] = {0.,0.,0.};
  double eqP = 0.;

  // Compute yield function at lambda_R
  ComputeCauchyStressGivenConsistencyParameterAndDirection(CStrial, lambda_R, N, CS);
  eqP = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda_R, N);
  double F_R = EvaluateYieldFunction(CS, eqP);

  // Compute yield function at lambda_L
  ComputeCauchyStressGivenConsistencyParameterAndDirection(CStrial, lambda_L, N, CS);
  eqP = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda_L, N);
  double F_L = EvaluateYieldFunction(CS, eqP);

  // Check that the bisection method can be applied
  if( F_L>0. || F_R<0. )
    {
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::ComputeConsistencyParameterGivenDirection()- "
               <<"Could not bracket function for bisection method.\n";
      return false;
    }

  // Maximum number of iterations
  int MaxIt = 1000;

  // Flag for convergence
  bool CONVERGED = false;

 // Use some tolerance for checking convergence
  double TOL = Tol1*ExpYieldStress[0]; // ORIG: 1e-6*ExpYieldStress[0]

  // Solution = lambda

  // Iterate
  for(int i=0; i<MaxIt && !CONVERGED; i++)
    {
      // New guess for solution
      lambda = 0.5*(lambda_L + lambda_R);
      ComputeCauchyStressGivenConsistencyParameterAndDirection(CStrial, lambda, N, CS);
      eqP = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda, N);

      // Evaluate yield function at this point
      double F = EvaluateYieldFunction(CS, eqP);

      // Check if root has been found
      if(std::abs(F) < TOL || (lambda_L-lambda_R)/2 < 0)  //ORIG: if( fabs(F)<TOL )
        CONVERGED = true;

      // Otherwise update bracket
      else
        {
          if(!((F < 0 && F_R < 0) || (F > 0 && F_R > 0))) //ORIG: if( F<0. )
            {
              lambda_L = lambda;
              F_L = F;
            }
          else
            {
              lambda_R = lambda;
              F_R = F;
            }
        }
    }

  // Check if solution converged
  if( !CONVERGED )
    {
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::ComputeConsistencyParameterGivenDirection()- "
               <<"Could not compute consistency parameter for given direction.\n";
      return false;
    }
  else
    return true;
}

// Compute eigenvalues
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ComputeEigenvalues(const double * X, double &A, double &B) const
{
  // Sum of roots
  double sor = X[0]+X[1];

  // Square root of discriminant
  double disc = sqrt( pow(X[0]-X[1],2.)+4.*pow(X[2],2.) );

  A = 0.5*(sor + disc);
  B = 0.5*(sor - disc);
}

// Compute the norm of a deviatoric stress tensor
// given in vector form [Sxx, Syy, Sxy]
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
DeviatoricStressNorm(const double * S) const
{
  return sqrt( 2.*(pow(S[0],2.)+pow(S[1],2.)+pow(S[2],2.)+S[0]*S[1]) );
}

// Compute the norm of a deviatoric strain tensor
// given in vector form [Sxx, Syy, 2Sxy]
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
DeviatoricStrainNorm(const double * S) const
{
  return sqrt( 2.*(pow(S[0],2.)+pow(S[1],2.)+S[0]*S[1]) + 0.5*pow(S[2],2.) );
}

// The 6.89475908677537e6 factor converts ksi to Pa
//#define USE_SI_UNITS
#ifdef USE_SI_UNITS
const double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ExpYieldStress[31] = { 39.126*6.89475908677537e6,
                       40.855*6.89475908677537e6,
                       42.015*6.89475908677537e6,
                       42.671*6.89475908677537e6,
                       43.292*6.89475908677537e6,
                       44.198*6.89475908677537e6,
                       44.806*6.89475908677537e6,
                       45.048*6.89475908677537e6,
                       45.631*6.89475908677537e6,
                       46.427*6.89475908677537e6,
                       47.087*6.89475908677537e6,
                       47.635*6.89475908677537e6,
                       48.458*6.89475908677537e6,
                       49.095*6.89475908677537e6,
                       49.848*6.89475908677537e6,
                       50.185*6.89475908677537e6,
                       50.659*6.89475908677537e6,
                       50.987*6.89475908677537e6,
                       51.291*6.89475908677537e6,
                       51.534*6.89475908677537e6,
                       51.850*6.89475908677537e6,
                       52.200*6.89475908677537e6,
                       52.500*6.89475908677537e6,
                       52.720*6.89475908677537e6,
                       52.900*6.89475908677537e6,
                       53.050*6.89475908677537e6,
                       53.170*6.89475908677537e6,
                       53.270*6.89475908677537e6,
                       53.370*6.89475908677537e6,
                       53.450*6.89475908677537e6,
                       53.580*6.89475908677537e6 };
#else
// these units are psi
const double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ExpYieldStress[31] = { 39.126e3,
                       40.855e3,
                       42.015e3,
                       42.671e3,
                       43.292e3,
                       44.198e3,
                       44.806e3,
                       45.048e3,
                       45.631e3,
                       46.427e3,
                       47.087e3,
                       47.635e3,
                       48.458e3,
                       49.095e3,
                       49.848e3,
                       50.185e3,
                       50.659e3,
                       50.987e3,
                       51.291e3,
                       51.534e3,
                       51.850e3,
                       52.200e3,
                       52.500e3,
                       52.720e3,
                       52.900e3,
                       53.050e3,
                       53.170e3,
                       53.270e3,
                       53.370e3,
                       53.450e3,
                       53.580e3 };
#endif

const double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding::
ExpEqPlasticStrain[31] = { 0.0000e-2,
                           0.0187e-2,
                           0.0513e-2,
                           0.0885e-2,
                           0.1409e-2,
                           0.2665e-2,
                           0.5801e-2,
                           0.7106e-2,
                           1.1467e-2,
                           1.8328e-2,
                           2.4497e-2,
                           3.0011e-2,
                           3.9337e-2,
                           4.8031e-2,
                           5.9797e-2,
                           6.6361e-2,
                           7.6791e-2,
                           8.6381e-2,
                           9.7461e-2,
                           10.9280e-2,
                           13.0000e-2,
                           16.0000e-2,
                           20.0000e-2,
                           24.0000e-2,
                           28.0000e-2,
                           32.0000e-2,
                           36.0000e-2,
                           40.0000e-2,
                           45.0000e-2,
                           50.0000e-2,
                           60.0000e-2 };
