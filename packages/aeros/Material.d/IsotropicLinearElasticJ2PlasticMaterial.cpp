// Sriramajayam

/*
 * IsotropicLinearElasticJ2PlasticMaterial.cpp
 * DG++
 *
 * Created by Ramsharan Rangarajan on 09/21/2010.
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


#include "IsotropicLinearElasticJ2PlasticMaterial.h"


// Constructor
IsotropicLinearElasticJ2PlasticMaterial::
IsotropicLinearElasticJ2PlasticMaterial(double iLambda, double iMu,
                                        double iSigmaY, double iK, double iH,
                                        double iTol, double iequivEPSplasticF,
                                        double ik, double in)
  :Mu(iMu), SigmaY(iSigmaY), K(iK), H(iH), Tol(iTol), equivEPSplasticF(iequivEPSplasticF), k(ik), n(in)
{
  // Compute bulk modulus
  Kappa = iLambda + (2./3.)*iMu;

  // Create material for elastic response
  ILE = new IsotropicLinearElastic(iLambda, iMu);

  // Initialize plastic strain to zero.
  EPSplastic.clear();
  for(int i=0; i<9; i++)
    EPSplastic.push_back( 0. );

  // Initialize equivalent plastic strain to zero
  equivEPSplastic = 0.;

  // Initialize back stress to zero
  BackStress.clear();
  for(int i=0; i<9; i++)
    BackStress.push_back( 0. );
}

// Destructor
IsotropicLinearElasticJ2PlasticMaterial::
~IsotropicLinearElasticJ2PlasticMaterial()
{ delete ILE; }

// Copy constructor
IsotropicLinearElasticJ2PlasticMaterial::
IsotropicLinearElasticJ2PlasticMaterial(const IsotropicLinearElasticJ2PlasticMaterial &Mat)
  : Kappa(Mat.Kappa), Mu(Mat.Mu), SigmaY(Mat.SigmaY), K(Mat.K), H(Mat.H), Tol(Mat.Tol), equivEPSplasticF(Mat.equivEPSplasticF),
    k(Mat.k), n(Mat.n), ExpEqPlasticStrain(Mat.ExpEqPlasticStrain), ExpYieldStress(Mat.ExpYieldStress),
    ExpEqPlasticStrainRate(Mat.ExpEqPlasticStrainRate), ExpYieldStressScale(Mat.ExpYieldStressScale)
{
  ILE = new IsotropicLinearElastic(*Mat.ILE);
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<9; i++)
    {
      EPSplastic.push_back( Mat.EPSplastic[i] );
      BackStress.push_back( Mat.BackStress[i] );
    }
  equivEPSplastic = Mat.equivEPSplastic;
}

// Cloning
IsotropicLinearElasticJ2PlasticMaterial *
IsotropicLinearElasticJ2PlasticMaterial::Clone() const
{ return new IsotropicLinearElasticJ2PlasticMaterial(*this); }

// Return plastic strain in material
std::vector<double> IsotropicLinearElasticJ2PlasticMaterial::GetMaterialPlasticStrain() const
{ return EPSplastic; }

// Resturn the equivalent plastic strain in material
double IsotropicLinearElasticJ2PlasticMaterial::GetMaterialEquivalentPlasticStrain() const
{ return equivEPSplastic; }

// Return back stress in material
std::vector<double> IsotropicLinearElasticJ2PlasticMaterial::GetMaterialBackStress() const
{ return BackStress; }

// Return isotropic hardening modulus
double IsotropicLinearElasticJ2PlasticMaterial::GetIsotropicHardeningModulus() const
{ return K; }

// Return kinematic hardening modulus
double IsotropicLinearElasticJ2PlasticMaterial::GetKinematicHardeningModulus() const
{ return H; }

// Return flow stress from uniaxial tension test
double IsotropicLinearElasticJ2PlasticMaterial::GetYieldStressFromTensionTest() const
{ return SigmaY; }

// Return bulk modulus
double IsotropicLinearElasticJ2PlasticMaterial::GetBulkModulus() const
{ return Kappa; }

// Return shear modulus of material
double IsotropicLinearElasticJ2PlasticMaterial::GetShearModulus() const
{ return Mu; }

// Return dissipated energy
double IsotropicLinearElasticJ2PlasticMaterial::GetDissipatedEnergy() const
{
  if(SigmaY < 0 || ExpEqPlasticStrainRate.size() > 0) {
    std::cerr << " *** WARNING: IsotropicLinearElasticJ2PlasticMaterial::GetDissipatedEnergy is not implemented for nonlinear isotropic hardening.\n";
    return 0;
  }
  return (SigmaY + 0.5*K*equivEPSplastic)*equivEPSplastic;
}

// Return equivalent plastic strain at failure
double IsotropicLinearElasticJ2PlasticMaterial::GetEquivalentPlasticStrainAtFailure() const
{ return equivEPSplasticF; }

// Return tolerance for convergence of nonlinear solve
double IsotropicLinearElasticJ2PlasticMaterial::GetTolerance() const
{ return Tol; }

// Set the plastic strain in the material
void IsotropicLinearElasticJ2PlasticMaterial::SetMaterialPlasticStrain(const std::vector<double> &iEPSplastic)
{ for(int i = 0; i < 9; ++i) EPSplastic[i] = iEPSplastic[i]; }

// Set the equivalent plastic strain in the material
void IsotropicLinearElasticJ2PlasticMaterial::SetMaterialEquivalentPlasticStrain(double iEquivEPSplastic)
{ equivEPSplastic = iEquivEPSplastic; }

// Set the back stress in the material
void IsotropicLinearElasticJ2PlasticMaterial::SetMaterialBackStress(const std::vector<double> &iBackStress)
{ for(int i = 0; i < 9; ++i) BackStress[i] = iBackStress[i]; }

// Compute linear elastic response
bool IsotropicLinearElasticJ2PlasticMaterial::
ComputeElasticConstitutiveResponse(const std::vector<double> &EPS,
                                   std::vector<double> *CS,
                                   std::vector<double> *C) const
{
  // Resize of necessary
  if( int(CS->size())<9 )
    CS->resize( 9 );

  if( C )
    if( int(C->size())<81 )
      C->resize( 81 );

  // Convert strain to deformation gradient
  std::vector<double> F = EPS;
  F[0] += 1.;
  F[4] += 1.;
  F[8] += 1.;

  return ILE->GetConstitutiveResponse(&F, CS, C);
}

// Compute deviatoric part of tensor
std::vector<double> IsotropicLinearElasticJ2PlasticMaterial::Deviatoric(const double * T) const
{
  // 1/3 of trace of T
  double TRby3 = (T[0] + T[4] + T[8])/3.;

  // Identity
  double I[] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};

  // Compute deviatoric part of T as T-(1/3)trace(T)*I.
  std::vector<double> devT(9);
  for(int i=0; i<9; i++)
    devT[i] = T[i]-TRby3*I[i];

  return devT;
}

// Compute norm of tensor
double IsotropicLinearElasticJ2PlasticMaterial::Norm(const double *T) const
{
  double norm2 = 0.;
  for(int i=0; i<9; i++)
    norm2 += pow(T[i],2.);
  return sqrt(norm2);
}

// Compute yield stress and isotropic hardening modulus by interpolating experimental stress-strain curve
double IsotropicLinearElasticJ2PlasticMaterial::
GetYieldStressUsingExperimentalCurve(const double eqP, double &K) const
{
  double SigmaY = 0.;
  const int N = std::min(ExpEqPlasticStrain.size(), ExpYieldStress.size());
  if( eqP<=ExpEqPlasticStrain[0] ) {
    SigmaY = ExpYieldStress[0];
    K = (ExpYieldStress[1]-ExpYieldStress[0])/(ExpEqPlasticStrain[1]-ExpEqPlasticStrain[0]);
  }
  else if ( eqP>=ExpEqPlasticStrain[N-1] ) {
    K = (ExpYieldStress[N-1]-ExpYieldStress[N-2])/(ExpEqPlasticStrain[N-1]-ExpEqPlasticStrain[N-2]);
    SigmaY = ExpYieldStress[N-1] + K*(eqP-ExpEqPlasticStrain[N-1]);
  }
  else
    {
      for(int i=0; i<N-1; i++)
        if( eqP>=ExpEqPlasticStrain[i] && eqP<ExpEqPlasticStrain[i+1] )
          {
            double lambda = (ExpEqPlasticStrain[i+1]-eqP)/(ExpEqPlasticStrain[i+1]-ExpEqPlasticStrain[i]);
            SigmaY = lambda*ExpYieldStress[i] + (1.-lambda)*ExpYieldStress[i+1];
            K = (ExpYieldStress[i+1]-ExpYieldStress[i])/(ExpEqPlasticStrain[i+1]-ExpEqPlasticStrain[i]);
            break;
          }
    }
  return SigmaY;
}

// Compute yield stress scaling factor and its derivate w.r.t eqPdot by interpolating experimental stress-strain curve
double IsotropicLinearElasticJ2PlasticMaterial::
GetScaleFactorUsingExperimentalCurve(const double eqPdot, double &R) const
{
  double ScaleF = 0.;
  const int N = std::min(ExpEqPlasticStrainRate.size(), ExpYieldStressScale.size());
  if( eqPdot<=ExpEqPlasticStrainRate[0] ) {
    ScaleF = ExpYieldStressScale[0];
    R = (ExpYieldStressScale[1]-ExpYieldStressScale[0])/(ExpEqPlasticStrainRate[1]-ExpEqPlasticStrainRate[0]);
  }
  else if ( eqPdot>=ExpEqPlasticStrainRate[N-1] ) {
    R = (ExpYieldStressScale[N-1]-ExpYieldStressScale[N-2])/(ExpEqPlasticStrainRate[N-1]-ExpEqPlasticStrainRate[N-2]);
    ScaleF = ExpYieldStressScale[N-1] + R*(eqPdot-ExpEqPlasticStrainRate[N-1]);
  }
  else
    {
      for(int i=0; i<N-1; i++)
        if( eqPdot>=ExpEqPlasticStrainRate[i] && eqPdot<ExpEqPlasticStrainRate[i+1] )
          {
            double lambda = (ExpEqPlasticStrainRate[i+1]-eqPdot)/(ExpEqPlasticStrainRate[i+1]-ExpEqPlasticStrainRate[i]);
            ScaleF = lambda*ExpYieldStressScale[i] + (1.-lambda)*ExpYieldStressScale[i+1];
            R = (ExpYieldStressScale[i+1]-ExpYieldStressScale[i])/(ExpEqPlasticStrainRate[i+1]-ExpEqPlasticStrainRate[i]);
            break;
          }
    }
  return ScaleF;
}

// Compute yield stress and isotropic hardening modulus by generalized power law
double IsotropicLinearElasticJ2PlasticMaterial::
GetYieldStressUsingGeneralizedPowerLaw(const double eqP, double &K) const
{
  using std::pow;
  const double &a = IsotropicLinearElasticJ2PlasticMaterial::SigmaY;
  const double &c = IsotropicLinearElasticJ2PlasticMaterial::K;
  double E = 9*Kappa*Mu/(3*Kappa+Mu);
  double eps0 = pow((E-c)/k,1./(n-1));
  double SigmaY = a + k*pow(eqP+eps0,n) + c*eqP;
  K = c + n*k*pow(eqP+eps0,n-1);

  return SigmaY;
}

// Evaluate yield function
double IsotropicLinearElasticJ2PlasticMaterial::
EvaluateYieldFunction(const double * CauchyStress,
                      const double * SigmaB,
                      const double DeltaEqP,
                      double &K,
                      const double dt) const
{
  // linear hardening: f = norm(xi) - sqrt(2/3) * (SigmaY + K*eqP),
  // or generalized power law hardening: f = norm(xi) - sqrt(2/3) * (a + k*(eqP+eps0)^n + c*eqP),
  // or piecewise linear hardening: f = norm(xi) - sqrt(2/3) * GetYieldStressUsingExperimentalCurve(eqP).
  // where
  // xi = dev(Cauchy stress)-SigmaB
  // DeltaEqP = equivalent plastic strain increment

  // Equivalent plastic strain
  double eqP = equivEPSplastic + DeltaEqP;

  // Deviatoric part of CauchyStress
  std::vector<double> S = Deviatoric(&CauchyStress[0]);

  // xi = S-SigmaB
  double xi[9];
  for(int i=0; i<9; i++)
    xi[i] = S[i]-SigmaB[i];

  // Norm of xi
  double normXI = Norm(xi);

  // Evaluate radius of yield surface
  double YSrad;
  if(SigmaY >= 0) {
    if(k == 0) { // linear isotropic hardening
      K = GetIsotropicHardeningModulus();
      YSrad = sqrt(2./3.)*(SigmaY+K*eqP);
    }
    else { // generalized power law isotropic hardening
      YSrad = sqrt(2./3.)*GetYieldStressUsingGeneralizedPowerLaw(eqP, K);
    }
  }
  else { // piecewise linear isotropic hardening
    YSrad = sqrt(2./3.)*GetYieldStressUsingExperimentalCurve(eqP, K);
  }

  // Scale yield stress due to strain-rate dependency
  if(ExpEqPlasticStrainRate.size() > 0 && dt > 0) {
    double R;
    double ScaleF = GetScaleFactorUsingExperimentalCurve(DeltaEqP/dt, R);
    K = sqrt(3./2.)*YSrad*R/dt + ScaleF*K;
    YSrad *= ScaleF;
  }

  return normXI - YSrad;
}

// Compute elasto-plastic constitutive response
bool IsotropicLinearElasticJ2PlasticMaterial::
ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
                                         std::vector<double> *CauchyStress,
                                         std::vector<double> *Cep,
                                         const bool UpdateFlag,
                                         const double dt)
{
  // Resize output Cauchy stress
  if( int(CauchyStress->size())<9 )
    CauchyStress->resize(9);

  // Resize output tangents if requested
  if( Cep )
    if( int(Cep->size())<81 )
      Cep->resize( 81 );

  // Check for failure
  if(equivEPSplastic >= equivEPSplasticF) {
    if( Cep ) for(int i=0; i<81; i++) (*Cep)[i] = 0;
    for(int i=0; i<9; i++) (*CauchyStress)[i] = 0;
    return true;
  }

  // Elastic modulii
  std::vector<double> * Ce = 0;
  if( Cep )
    Ce = new std::vector<double>(81);

  // Identity tensor
  double I[9] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};

  // Compute symmetric infinitesimal strain tensor
  std::vector<double> EPS(9);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      EPS[3*i+j] = 0.5*(Fnp1[3*i+j]+Fnp1[3*j+i])-I[3*i+j];

  // EVALUATE TRIAL STATE BY FREEZING PLASTIC FLOW

  // Trial elastic strain
  std::vector<double> epsEtrial(9);
  for(int i=0; i<9; i++)
    epsEtrial[i] = EPS[i]-EPSplastic[i];

  // Trial Cauchy stress
  std::vector<double> CStrial(9);
  if( !ComputeElasticConstitutiveResponse(epsEtrial, &CStrial, Ce) )
    {
      std::cerr<<"\n IsotropicLinearElasticJ2PlasticMaterial::ComputeElastoPlasticConstitutiveResponse()- "
               <<"Could not compute elastic response of material. \n";
      return false;
    }

  // Evaluate yield function and isotropic hardening modulus at trial state
  double K;
  double Ftrial = EvaluateYieldFunction(&CStrial[0], &BackStress[0], 0., K, dt);


  // COMPUTE NEW STATE OF MATERIAL

  if( Ftrial<= 0.)
    {
      // This step is purely elastic.
      // Final state of the material is the trial state.

      for(int i=0; i<9; i++)
        (*CauchyStress)[i] = CStrial[i];

      if( Cep )
        for(int i=0; i<81; i++)
          (*Cep)[i] = (*Ce)[i];

      // No need to update plastic internal variables
    }

  else
    {
      // This step is elasto-plastic.
      // The trial state is not the actual state of the material

      // Compute consistency parameter
      double dLambda = Ftrial/(2.*Mu + (2./3.)*(H+K));

      // Strial
      std::vector<double> Strial = Deviatoric(&CStrial[0]);

      // XItrial = Strial-BackStress
      std::vector<double> XItrial(9);
      for(int i=0; i<9; i++)
        XItrial[i] = Strial[i]-BackStress[i];

      double normXItrial = Norm(&XItrial[0]);
      if( normXItrial<1.e-8 )
        {
          std::cerr<<"\n IsotropicLinearElasticJ2PlasticMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                   <<"Norm of trial stress close to zero. \n";
          return false;
        }

      // Unit tensor XItrial/norm(XItrial) = XI_(n+1)/norm(XI_(n+1)).
      double N[9];
      for(int i=0; i<9; i++)
        N[i] = XItrial[i]/normXItrial;

      // New elastic strain and back stress
      std::vector<double> NewEPSelastic(9), NewBackStress(9);

      int nItMax = 20;
      for(int j=0; j<nItMax; ++j) { // Newton iterations to solve for dLambda, if required

        for(int i=0; i<9; i++)
            NewEPSelastic[i] = EPS[i] - EPSplastic[i] - dLambda*N[i];

        // New Cauchy Stress
        if( !ComputeElasticConstitutiveResponse(NewEPSelastic, CauchyStress) )
           {
            std::cerr<<"\n IsotropicLinearElasticJ2PlasticMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                     <<"Could not compute elastic response.\n";
            return false;
          }
 
        if(SigmaY >= 0 && k == 0) break; // for linear isotropic hardening dLambda has a closed-form solution
        else {
          for(int i=0; i<9; i++)
            NewBackStress[i] = BackStress[i] + (2./3.)*H*dLambda*N[i];
          double F = EvaluateYieldFunction(&(*CauchyStress)[0], &NewBackStress[0], sqrt(2./3.)*dLambda, K, dt);
          if(std::abs(F) <= Tol*(SigmaY >= 0 ? SigmaY : ExpYieldStress[0]))
            break;
          else
            dLambda += F/(2.*Mu + (2./3.)*(H+K));
        }
        if( j+1 == nItMax )
        {
          std::cerr<<"\n IsotropicLinearElasticJ2PlasticMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                   <<"Iterations to compute consistency parameter did not converge.\n";
          return false;
        }
      }

      // Check for failure
      if(equivEPSplasticF < std::numeric_limits<double>::infinity() &&
         equivEPSplastic+sqrt(2./3.)*sqrt(2./3.)*dLambda >= equivEPSplasticF)
        {
          if( Cep ) for(int i=0; i<81; i++) (*Cep)[i] = 0;
          for(int i=0; i<9; i++) (*CauchyStress)[i] = 0;
        }

      else if( Cep )
        {
          // Evaluate the consistent elasto-plastic modulii
          double theta    = 1. -2.*Mu*dLambda/normXItrial;
          double thetabar = 1./(1.+(K+H)/(3.*Mu)) + theta -1.;
        
          for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
              for(int k=0; k<3; k++)
                for(int L=0; L<3; L++)
                  (*Cep)[27*i+9*j+3*k+L] =
                    Kappa*I[3*i+j]*I[3*k+L] +
                    2.*Mu*theta*( 0.5*(I[3*i+k]*I[3*j+L] + I[3*i+L]*I[3*j+k]) - (1./3.)*I[3*i+j]*I[3*k+L] ) -
                    2.*Mu*thetabar*N[3*i+j]*N[3*k+L];
        }

      // Update plastic strain and backstress if requested
      if( UpdateFlag==true )
        {
          equivEPSplastic += sqrt(2./3.)*dLambda;
          for(int i=0; i<9; i++)
            {
              EPSplastic[i] += dLambda*N[i];
              BackStress[i] += (2./3.)*H*dLambda*N[i];
            }
        }
    }

  if( Ce )
    delete Ce;

  return true;
}

// Checks if state of the material lies on or within yield surface
bool IsotropicLinearElasticJ2PlasticMaterial::
CheckMaterialState(const std::vector<double> &CS, const double TOL, const double dt) const
{
  double K;
  double f = EvaluateYieldFunction(&CS[0], &BackStress[0], 0., K, 0.);
  if( (SigmaY >= 0 && f < TOL*SigmaY) || (SigmaY < 0 && f < TOL*ExpYieldStress[0]) )
    return true;
  else
    return false;
}

// Set the (x,y) values at a point in the yield stress vs. effective plastic strain curve
void IsotropicLinearElasticJ2PlasticMaterial::
SetExperimentalCurveData(double iEqPlasticStrain, double iYieldStress)
{
  ExpEqPlasticStrain.push_back(iEqPlasticStrain);
  ExpYieldStress.push_back(iYieldStress);
}

// Set the (x,y) values at a point in the yield stress scale factor vs. effective plastic strain rate curve
void IsotropicLinearElasticJ2PlasticMaterial::
SetExperimentalCurveData2(double iEqPlasticStrainRate, double iYieldStressScale)
{
  ExpEqPlasticStrainRate.push_back(iEqPlasticStrainRate);
  ExpYieldStressScale.push_back(iYieldStressScale);
}

