// Sriramajayam

/*
 * KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2.cpp
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


#include "KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2.h"

// Constructor
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2(double iLambda, double iMu, double iTol1, double iTol2)
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
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
~KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2()
{}

// Copy constructor
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2
(const KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2 & Mat)
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
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2 *
KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::Clone() const
{ return new KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2(*this); }

// Return plastic strain
const std::vector<double> & KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
const std::vector<double> & KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
GetMaterialBackStress() const
{
  return BackStress;
}

// Return equivalent plastic strain
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
GetMaterialEquivalentPlasticStrain() const
{ return equivEPSplastic; }

// Return Bulk modulus
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
GetBulkModulus() const
{ return E/(3.*(1.-2.*nu)); }

// Return shear modulus
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
GetShearModulus() const
{ return 0.5*E/(1.+nu); }

// Return dissipated energy
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
GetDissipatedEnergy() const
{
  std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::GetDissipatedEnergy()- "
           <<"not implemented.\n";
  return 0;
}

// Set the plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
SetMaterialPlasticStrain(const std::vector<double> &iEPSplastic)
{ for(int i = 0; i < 3; ++i) EPSplastic[i] = iEPSplastic[i]; }

// Set the plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
SetMaterialBackStress(const std::vector<double> &)
{}

// Set the equivalent plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
SetMaterialEquivalentPlasticStrain(double iEquivEPSplastic)
{ equivEPSplastic = iEquivEPSplastic; }

// Check if state of material lies within yield surface
bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
                                         std::vector<double> *CauchyStress,
                                         std::vector<double> *Cep,
                                         const bool UpdateFlag,
                                         const double dt)
{
  // Notify that elastoplastic tangents are not computed.
  if( Cep!=0 )
    {
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::ComputeElastoPlasticConstitutiveResponse()- "
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
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::ComputeElastoPlasticConstitutiveResponse()- "
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
          std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::ComputeElastoPlasticConstitutiveResponse()- "
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
              std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::ComputeElastoPlasticConstitutiveResponse()- "
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
          std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::ComputeElastoPlasticConstitutiveResponse()- "
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
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
std::vector<double> KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
GetYieldStressUsingExperimentalCurve(const double eqP) const
{
  double SigmaY = 0.;
  if( eqP<=ExpEqPlasticStrain[0] )
    SigmaY = ExpYieldStress[0];
  else if ( eqP>=ExpEqPlasticStrain[53] )
    SigmaY = ExpYieldStress[53];
  else
    {
      for(int i=0; i<53; i++)
        if( eqP>=ExpEqPlasticStrain[i] && eqP<=ExpEqPlasticStrain[i+1] )
          {
            double lambda = (ExpEqPlasticStrain[i+1]-eqP)/(ExpEqPlasticStrain[i+1]-ExpEqPlasticStrain[i]);
            SigmaY = lambda*ExpYieldStress[i] + (1.-lambda)*ExpYieldStress[i+1];
          }
    }
  return SigmaY;
}

// Compute \f$\sigma\f$ given \f$\sigma_{trial}, \lambda, {\bf N} \f$.
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(const double lambda, const double * N) const
{
  return equivEPSplastic + sqrt(2./3.)*lambda*DeviatoricStressNorm(N);
}

bool KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::ComputeConsistencyParameterGivenDirection()- "
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
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::ComputeConsistencyParameterGivenDirection()- "
               <<"Could not compute consistency parameter for given direction.\n";
      return false;
    }
  else
    return true;
}

// Compute eigenvalues
void KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
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
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
DeviatoricStressNorm(const double * S) const
{
  return sqrt( 2.*(pow(S[0],2.)+pow(S[1],2.)+pow(S[2],2.)+S[0]*S[1]) );
}

// Compute the norm of a deviatoric strain tensor
// given in vector form [Sxx, Syy, 2Sxy]
double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
DeviatoricStrainNorm(const double * S) const
{
  return sqrt( 2.*(pow(S[0],2.)+pow(S[1],2.)+S[0]*S[1]) + 0.5*pow(S[2],2.) );
}

// The 6.89475908677537e3 factor converts psi to Pa
//#define USE_SI_UNITS
#ifdef USE_SI_UNITS
const double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
ExpYieldStress[54] = { 19540.0*6.89475908677537e3,
                       25490.0*6.89475908677537e3,
                       28480.0*6.89475908677537e3,
                       31310.0*6.89475908677537e3,
                       34180.0*6.89475908677537e3,
                       37040.0*6.89475908677537e3,
                       39720.0*6.89475908677537e3,
                       41340.0*6.89475908677537e3,
                       42360.0*6.89475908677537e3,
                       43040.0*6.89475908677537e3,
                       43490.0*6.89475908677537e3,
                       43870.0*6.89475908677537e3,
                       44140.0*6.89475908677537e3,
                       44310.0*6.89475908677537e3,
                       44620.0*6.89475908677537e3,
                       44880.0*6.89475908677537e3,
                       45030.0*6.89475908677537e3,
                       45370.0*6.89475908677537e3,
                       45520.0*6.89475908677537e3,
                       45580.0*6.89475908677537e3,
                       45730.0*6.89475908677537e3,
                       45880.0*6.89475908677537e3,
                       46120.0*6.89475908677537e3,
                       46280.0*6.89475908677537e3,
                       46380.0*6.89475908677537e3,
                       46480.0*6.89475908677537e3,
                       46620.0*6.89475908677537e3,
                       46700.0*6.89475908677537e3,
                       46770.0*6.89475908677537e3,
                       46850.0*6.89475908677537e3,
                       46910.0*6.89475908677537e3,
                       47000.0*6.89475908677537e3,
                       47080.0*6.89475908677537e3,
                       47150.0*6.89475908677537e3,
                       47230.0*6.89475908677537e3,
                       47300.0*6.89475908677537e3,
                       47370.0*6.89475908677537e3,
                       47430.0*6.89475908677537e3,
                       47480.0*6.89475908677537e3,
                       47520.0*6.89475908677537e3,
                       47570.0*6.89475908677537e3,
                       47680.0*6.89475908677537e3,
                       47760.0*6.89475908677537e3,
                       47840.0*6.89475908677537e3,
                       47870.0*6.89475908677537e3,
                       47910.0*6.89475908677537e3,
                       47970.0*6.89475908677537e3,
                       47990.0*6.89475908677537e3,
                       48010.0*6.89475908677537e3,
                       48010.0*6.89475908677537e3,
                       48000.0*6.89475908677537e3,
                       48020.0*6.89475908677537e3,
                       48060.0*6.89475908677537e3,
                       48061.0*6.89475908677537e3 };
#else
// these units are psi
const double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
ExpYieldStress[54] = { 19540.0,
                       25490.0,
                       28480.0,
                       31310.0,
                       34180.0,
                       37040.0,
                       39720.0,
                       41340.0,
                       42360.0,
                       43040.0,
                       43490.0,
                       43870.0,
                       44140.0,
                       44310.0,
                       44620.0,
                       44880.0,
                       45030.0,
                       45370.0,
                       45520.0,
                       45580.0,
                       45730.0,
                       45880.0,
                       46120.0,
                       46280.0,
                       46380.0,
                       46480.0,
                       46620.0,
                       46700.0,
                       46770.0,
                       46850.0,
                       46910.0,
                       47000.0,
                       47080.0,
                       47150.0,
                       47230.0,
                       47300.0,
                       47370.0,
                       47430.0,
                       47480.0,
                       47520.0,
                       47570.0,
                       47680.0,
                       47760.0,
                       47840.0,
                       47870.0,
                       47910.0,
                       47970.0,
                       47990.0,
                       48010.0,
                       48010.0,
                       48000.0,
                       48020.0,
                       48060.0,
                       48061.0 };
#endif

const double KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2::
ExpEqPlasticStrain[54] = { 0.0,
                           0.000584,
                           0.000881,
                           0.001162,
                           0.001453,
                           0.001771,
                           0.00216,
                           0.002594,
                           0.003028,
                           0.003426,
                           0.003776,
                           0.004092,
                           0.004378,
                           0.004625,
                           0.005038,
                           0.005383,
                           0.005688,
                           0.006384,
                           0.006821,
                           0.00708,
                           0.007673,
                           0.00852,
                           0.01021,
                           0.01189,
                           0.01353,
                           0.0152,
                           0.0175,
                           0.01893,
                           0.02063,
                           0.02234,
                           0.02406,
                           0.02575,
                           0.02743,
                           0.02913,
                           0.0308,
                           0.03251,
                           0.03424,
                           0.03597,
                           0.03768,
                           0.03927,
                           0.0411,
                           0.04541,
                           0.04898,
                           0.05233,
                           0.05408,
                           0.0567,
                           0.05839,
                           0.05924,
                           0.06095,
                           0.06181,
                           0.06268,
                           0.06355,
                           0.06442,
                           1.0 };
