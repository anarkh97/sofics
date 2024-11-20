// Sriramajayam

/*
 * KorkolisKyriakidesPlaneStressMaterial.cpp
 * DG++
 *
 * Created by Ramsharan Rangarajan on 11/29/2010.
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


#include "KorkolisKyriakidesPlaneStressMaterial.h"

// Constructor
KorkolisKyriakidesPlaneStressMaterial::
KorkolisKyriakidesPlaneStressMaterial(double iLambda, double iMu,
                                      double iSigmaY, double iK,
                                      double iH, double iTol1, double iTol2)
{
  // Youngs modulus
  E = iMu*(3.*iLambda + 2.*iMu)/(iLambda + iMu);

  // Poisson ratio
  nu = 0.5*iLambda/(iLambda+iMu);

  // Constants related to plasticity
  SigmaY = iSigmaY;
  K = iK;
  H = iH;

  // Tolerances for convergence of nonlinear solve
  Tol1 = (iTol1 > 0) ? iTol1 : 1.0e-6;
  Tol2 = (iTol2 > 0) ? iTol2 : 1.0e-6;

  // Zero initial plastic strain and backstress
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<3; i++)
    {
      EPSplastic.push_back( 0. );
      BackStress.push_back( 0. );
    }
  equivEPSplastic = 0.;
}

// Destructor
KorkolisKyriakidesPlaneStressMaterial::~KorkolisKyriakidesPlaneStressMaterial()
{}

// Copy constructor
KorkolisKyriakidesPlaneStressMaterial::
KorkolisKyriakidesPlaneStressMaterial(const KorkolisKyriakidesPlaneStressMaterial &Mat)
  : E(Mat.E), nu(Mat.nu), SigmaY(Mat.SigmaY), K(Mat.K), H(Mat.H), Tol1(Mat.Tol1), Tol2(Mat.Tol2)
{
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<3; i++)
    {
      EPSplastic.push_back( Mat.EPSplastic[i] );
      BackStress.push_back( Mat.BackStress[i] );
    }
  equivEPSplastic = Mat.equivEPSplastic;
}

// Cloning
KorkolisKyriakidesPlaneStressMaterial *
KorkolisKyriakidesPlaneStressMaterial::Clone() const
{ return new KorkolisKyriakidesPlaneStressMaterial(*this); }

// Return plastic strain
const std::vector<double> & KorkolisKyriakidesPlaneStressMaterial::
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

// Return equivalent plastic strain
double KorkolisKyriakidesPlaneStressMaterial::
GetMaterialEquivalentPlasticStrain() const
{ return equivEPSplastic; }


// Return back stress
const std::vector<double> & KorkolisKyriakidesPlaneStressMaterial::
GetMaterialBackStress() const
{
/* PJSA
  std::vector<double> BS(9,0.);
  BS[0] = BackStress[0];
  BS[4] = BackStress[1];
  BS[8] = -(BS[0]+BS[4]);
  BS[1] = BS[3] = BackStress[2];
  BS[2] = BS[5] = BS[6] = BS[7] = 0.;
  return BS;
*/
  return BackStress;
}

// Return isotropic hardening modulus
double KorkolisKyriakidesPlaneStressMaterial::
GetIsotropicHardeningModulus() const
{ return K; }

// Return kinematic hardening modulus
double KorkolisKyriakidesPlaneStressMaterial::
GetKinematicHardeningModulus() const
{ return H; }

// Return yield stress for 1D tension test
double KorkolisKyriakidesPlaneStressMaterial::
GetYieldStressFromTensionTest() const
{ return SigmaY; }

// Return Bulk modulus
double KorkolisKyriakidesPlaneStressMaterial::
GetBulkModulus() const
{ return E/(3.*(1.-2.*nu)); }

// Return shear modulus
double KorkolisKyriakidesPlaneStressMaterial::
GetShearModulus() const
{ return 0.5*E/(1.+nu); }

// Return dissipated energy
double KorkolisKyriakidesPlaneStressMaterial::
GetDissipatedEnergy() const
{ return SigmaY*equivEPSplastic; }

// Set the plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterial::
SetMaterialPlasticStrain(const std::vector<double> &iEPSplastic)
{ for(int i = 0; i < 3; ++i) EPSplastic[i] = iEPSplastic[i]; }

// Set the equivalent plastic strain in the material
void KorkolisKyriakidesPlaneStressMaterial::
SetMaterialEquivalentPlasticStrain(double iEquivEPSplastic)
{ equivEPSplastic = iEquivEPSplastic; }

// Set the back stress in the material
void KorkolisKyriakidesPlaneStressMaterial::
SetMaterialBackStress(const std::vector<double> &iBackStress)
{ for(int i = 0; i < 3; ++i) BackStress[i] = iBackStress[i]; }

// Print all of the internal variables
void KorkolisKyriakidesPlaneStressMaterial::
Print()
{
  std::cerr << "Plastic Strain = " << EPSplastic[0] << " " << EPSplastic[1] << " " << EPSplastic[2] << std::endl;
  std::cerr << "Back Stress = " << BackStress[0] << " " << BackStress[1] << " " << BackStress[2] << std::endl;
  std::cerr << "Equivalent Plastic Strain = " << equivEPSplastic << std::endl;
}

// Compute the elastic constitutive response
bool KorkolisKyriakidesPlaneStressMaterial::
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

// Linear transformations of stress involved in this model
void KorkolisKyriakidesPlaneStressMaterial::
ComputeLinearTransformationsOfStress(const double * Xi, double * L1Xi, double * L2Xi) const
{
  // Model parameters
  double a1 = 0.78, a2 = 1.15, a3 = 0.85, a4 = 0.89, a5 = 1.06, a6 = 1.03, a7 = 1.0, a8 = 1.0;

  // Compute L1Xi
  L1Xi[0] = (a1/3.)*(2.*Xi[0]-Xi[1]);
  L1Xi[1] = (a2/3.)*(2.*Xi[1]-Xi[0]);
  L1Xi[2] = a7*Xi[2];

  // Entries for transformation L2
  double Lvec[] = { (8.*a5 - 2.*a3 - 2.*a6 + 2.*a4)/9.,
                    (4.*a6 - 4.*a4 - 4.*a5 + 1.*a3)/9.,
                    (4.*a3 - 4.*a5 - 4.*a4 + 1.*a6)/9.,
                    (8.*a4 - 2.*a6 - 2.*a3 + 2.*a5)/9.,
                    a8 };

  // Compute L2Xi
  L2Xi[0] = Lvec[0]*Xi[0]+Lvec[1]*Xi[1];
  L2Xi[1] = Lvec[2]*Xi[0]+Lvec[3]*Xi[1];
  L2Xi[2] = Lvec[4]*Xi[2];
}

// Compute eigenvalues
void KorkolisKyriakidesPlaneStressMaterial::
ComputeEigenvalues(const double * X, double &A, double &B) const
{
  // Sum of roots
  double sor = X[0]+X[1];

  // Square root of discriminant
  double disc = sqrt( pow(X[0]-X[1],2.)+4.*pow(X[2],2.) );

  A = 0.5*(sor + disc);
  B = 0.5*(sor - disc);
}

// Evaluate yield function
double KorkolisKyriakidesPlaneStressMaterial::
EvaluateYieldFunction(const double * Xi, const double eqP) const
{
  // Compute two linear transformations of Xi
  double L1Xi[]  = {0.,0.,0.};
  double L2Xi[] = {0.,0.,0.};
  ComputeLinearTransformationsOfStress(Xi, L1Xi, L2Xi);

  // Compute L1Xi/SigmaY and L2Xi/SigmaY (avoid large numbers)
  double L1hat[3], L2hat[3];
  for(int i=0; i<3; i++)
    {
      L1hat[i] = L1Xi[i]/SigmaY;
      L2hat[i] = L2Xi[i]/SigmaY;
    }

  // Compute eigenvalues of L1hat
  double a=0., b=0.;
  ComputeEigenvalues(L1hat, a, b);

  // Compute eigenvalues of L2hat
  double A=0., B=0.;
  ComputeEigenvalues(L2hat, A, B);

  double SigmaHat8 = (pow(a-b,8.)+pow(2.*A+B,8.)+pow(A+2.*B,8.))/2.;
  double EffStress = SigmaY*pow(SigmaHat8,0.125);

  double Rad = SigmaY+K*eqP;

  return EffStress - Rad;
}

// Check if state of material lies within yield surface
bool KorkolisKyriakidesPlaneStressMaterial::
CheckMaterialState(const std::vector<double> &CS, const double TOL) const
{
  // Compute Xi = Cauchy stress - Back stress (in vector form)
  double Xi[3] = { CS[0]-BackStress[0],
                   CS[4]-BackStress[1],
                   CS[1]-BackStress[2] };

  // Evaluate yield function
  double Fval = EvaluateYieldFunction(Xi, equivEPSplastic);

  // Check for tolerance after non-dimensionalizing Fval
  if( Fval/SigmaY < TOL )
    return true;
  else
    return false;
}

// Evaluate derivative of yield function with respect to Cauchy stress
std::vector<double> KorkolisKyriakidesPlaneStressMaterial::
EvaluateDerivativeOfYieldFunction(const double * Xi, const double eqP) const
{
  // Derivative wrt Cauchy stress = Derivative wrt Xi.
  std::vector<double> N(3,0.);

  // How much to perturb
  double deltaSigma = (1.e-6)*SigmaY;

  double MyXi[] = {Xi[0], Xi[1], Xi[2]};

  for(int i=0; i<3; i++)
    {
      // Plus direction
      MyXi[i] = Xi[i]+deltaSigma;
      double Fplus = EvaluateYieldFunction(MyXi, eqP);

      // Minus direction
      MyXi[i] = Xi[i]-deltaSigma; // ORIG: Xi[i]-2.*deltaSigma
      double Fminus = EvaluateYieldFunction(MyXi, eqP);

      // Compute derivative
      N[i] = (Fplus-Fminus)/(2.*deltaSigma);

      // Restore component of Xi.
      MyXi[i] = Xi[i];
    }
  return N;
}

// Compute the norm of a deviatoric stress tensor
// given in vector form [Sxx, Syy, Sxy]
double KorkolisKyriakidesPlaneStressMaterial::
DeviatoricStressNorm(const double * S) const
{
  return sqrt( 2.*(pow(S[0],2.)+pow(S[1],2.)+pow(S[2],2.)+S[0]*S[1]) );
}

// Compute the norm of a deviatoric strain tensor
// given in vector form [Sxx, Syy, 2Sxy]
double KorkolisKyriakidesPlaneStressMaterial::
DeviatoricStrainNorm(const double * S) const
{
  return sqrt( 2.*(pow(S[0],2.)+pow(S[1],2.)+S[0]*S[1]) + 0.5*pow(S[2],2.) );
}

// Compute Xi given Xitrial, lambda and N.
void KorkolisKyriakidesPlaneStressMaterial::
ComputeXiGivenConsistencyParameterAndDirection(const double * Xitrial, const double lambda,
                                               const double * N, double * Xi) const
{
  // Cauchy stress = CStrial - lambda*C*N.
  // Back stress   = BStrial + (2/3)*lambda*H*N
  // Xi = Cauchy Stress - Back Stress

  // Elastic modulii.
  double v = E/(1.-nu*nu);
  double Cmat[3][3] = {{v*1., v*nu, 0.},
                       {v*nu, v, 0.},
                       {0., 0., v*(1.-nu)/2.}};



  for(int i=0; i<3; i++)
    {
      Xi[i] = Xitrial[i] + (2./3.)*lambda*H*N[i];
      for(int j=0; j<3; j++)
        Xi[i] -= lambda*Cmat[i][j]*N[j];
    }
}

// Compute equivalent plastic strain given lambda and N.
double KorkolisKyriakidesPlaneStressMaterial::
ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(const double lambda, const double * N) const
{
  return equivEPSplastic + sqrt(2./3.)*lambda*DeviatoricStressNorm(N);
}

// Compute elastoplastic response of material
bool KorkolisKyriakidesPlaneStressMaterial::
ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
                                         std::vector<double> *CauchyStress,
                                         std::vector<double> *Cep,
                                         const bool UpdateFlag,
                                         const double dt)
{
  // Notify that elastoplastic tangents are not computed.
  if( Cep!=0 )
    {
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
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
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
               <<"Could not compute elastic constitutive response at trial state.\n";
      return false;
    }

  // Trial Xi
  double Xitrial[] = { CStrial[0]-BackStress[0],
                       CStrial[1]-BackStress[1],
                       CStrial[2]-BackStress[2] };

  // Evaluate yield function at trial state
  double Ftrial = EvaluateYieldFunction(Xitrial, equivEPSplastic);

  // Use some tolerance for checking yield function value
  double TOL = SigmaY*Tol2; // ORIG: SigmaY*1.e-6

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

      // Initial guess for N = df/d(sigma) at (Xitrial, eqPtrial)
      std::vector<double> Ntrial = EvaluateDerivativeOfYieldFunction(Xitrial, equivEPSplastic);
      double normNtrial = DeviatoricStressNorm(&Ntrial[0]);
      if( normNtrial < 1.e-6 )
        {
          std::cerr<<"\nKorkolisKyriakidesPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
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
          if( !ComputeConsistencyParameterGivenDirection(Xitrial, &N[0], lambdaMax, 0., lambda) )
            {
              std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                       <<"Could not compute consistency parameter for given direction.\n";
              return false;
            }

          // Update guess for Xi
          double XiUpdate[3] = {0.,0.,0.};
          ComputeXiGivenConsistencyParameterAndDirection(Xitrial, lambda, &N[0], XiUpdate);

          // Update guess for equivalent plastic strain
          double eqPupdate = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda, &N[0]);

          // Update guess for normal direction
          N = EvaluateDerivativeOfYieldFunction(XiUpdate, eqPupdate);

          // Update Xi and eqP with this normal direction
          ComputeXiGivenConsistencyParameterAndDirection(Xitrial, lambda, &N[0], XiUpdate);
          eqPupdate = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda, &N[0]);

          // Check if this updated state is a solution
          double Fval = EvaluateYieldFunction(XiUpdate, eqPupdate);
          if(  fabs(Fval) < TOL )
            // This (lambda, N) is a solution
            CONVERGED = true;
        }

      // Check if solution converged
      if( !CONVERGED )
        {
          std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
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
          for(int i=0; i<3; i++)
            {
              // Update plastic strain
              EPSplastic[i] += lambda*N[i];

              // Update backstress
              BackStress[i] += (2./3.)*lambda*H*N[i];
            }

          // Update equivalent plastic strain
          equivEPSplastic += sqrt(2./3.)*lambda*DeviatoricStressNorm(&N[0]);
        }
    }
  return true;
}

bool KorkolisKyriakidesPlaneStressMaterial::
ComputeConsistencyParameterGivenDirection(const double * Xitrial, const double * N,
                                          double lambda_L, double lambda_R,
                                          double & lambda) const
{
  // Assume that f>=0 for lambda_L and f<=0 at lambda_R.

  double Xi[3] = {0.,0.,0.};
  double eqP = 0.;

  // Compute yield function at lambda_R
  ComputeXiGivenConsistencyParameterAndDirection(Xitrial, lambda_R, N, Xi);
  eqP = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda_R, N);
  double F_R = EvaluateYieldFunction(Xi, eqP);

  // Compute yield function at lambda_L
  ComputeXiGivenConsistencyParameterAndDirection(Xitrial, lambda_L, N, Xi);
  eqP = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda_L, N);
  double F_L = EvaluateYieldFunction(Xi, eqP);

  // Check that the bisection method can be applied
  if( F_L>0. || F_R<0. )
    {
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterial::ComputeConsistencyParameterGivenDirection()- "
               <<"Could not bracket function for bisection method.\n";
      return false;
    }

  // Maximum number of iterations
  int MaxIt = 1000;

  // Flag for convergence
  bool CONVERGED = false;

 // Use some tolerance for checking convergence
  double TOL = Tol1*SigmaY; // ORIG: 1.e-6*SigmaY

  // Solution = lambda

  // Iterate
  for(int i=0; i<MaxIt && !CONVERGED; i++)
    {
      // New guess for solution
      lambda = 0.5*(lambda_L + lambda_R);
      ComputeXiGivenConsistencyParameterAndDirection(Xitrial, lambda, N, Xi);
      eqP = ComputeEquivalentPlasticStrainGivenConsistencyParameterAndDirection(lambda, N);

      // Evaluate yield function at this point
      double F = EvaluateYieldFunction(Xi, eqP);

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
      std::cerr<<"\n KorkolisKyriakidesPlaneStressMaterial::ComputeConsistencyParameterGivenDirection()- "
               <<"Could not compute consistency parameter for given direction.\n";
      return false;
    }
  else
    return true;
}

