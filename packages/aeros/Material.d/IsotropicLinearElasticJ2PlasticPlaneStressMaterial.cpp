// Sriramajayam

/*
 * IsotropicLinearElasticJ2PlasticPlaneStressMaterial.cpp
 * DG++
 *
 * Created by Ramsharan Rangarajan on 11/08/2010.
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

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
#include "IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h"


// Constructor
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
IsotropicLinearElasticJ2PlasticPlaneStressMaterial(double iLambda, double iMu,
                                                   double iSigmaY, double iK,
                                                   double iH, double iTol,
                                                   double iequivEPSplasticF)
{
  // Youngs modulus
  E = iMu*(3.*iLambda + 2.*iMu)/(iLambda + iMu);

  // Poisson ratio
  nu = 0.5*iLambda/(iLambda+iMu);

  // Constants related to plasticity
  SigmaY = iSigmaY;
  K = iK;
  H = iH;
  Tol = (iTol > 0) ? iTol : 1.0e-6;
  equivEPSplasticF = iequivEPSplasticF;

  // Zero initial plastic strain and backstress
  EPSplastic.clear();
  BackStress.clear();
  for(int i=0; i<3; i++)
    {
      EPSplastic.push_back( 0. );
      BackStress.push_back( 0. );
    }
  equivEPSplastic = 0.;

  // Pre-compute tensor coefficients
  t00 = (2./3.)*H + E/(3.*(1.-nu*nu))*(2.-nu);
  t01 = E/(3.*(1.-nu*nu))*(2.*nu-1.);
  t22 = (2./3.)*H + E/(1.-nu*nu)*(1.-nu);
}

// Destructor
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::~IsotropicLinearElasticJ2PlasticPlaneStressMaterial()
{}

// Copy constructor
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
IsotropicLinearElasticJ2PlasticPlaneStressMaterial(const IsotropicLinearElasticJ2PlasticPlaneStressMaterial &Mat)
  :E(Mat.E), nu(Mat.nu), SigmaY(Mat.SigmaY), K(Mat.K), H(Mat.H), Tol(Mat.Tol), t00(Mat.t00), t01(Mat.t01), t22(Mat.t22),
   equivEPSplasticF(Mat.equivEPSplasticF), ExpEqPlasticStrain(Mat.ExpEqPlasticStrain), ExpYieldStress(Mat.ExpYieldStress),
   ExpEqPlasticStrainRate(Mat.ExpEqPlasticStrainRate), ExpYieldStressScale(Mat.ExpYieldStressScale)
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
IsotropicLinearElasticJ2PlasticPlaneStressMaterial *
IsotropicLinearElasticJ2PlasticPlaneStressMaterial::Clone() const
{ return new IsotropicLinearElasticJ2PlasticPlaneStressMaterial(*this); }

// Return isotropic hardening modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetIsotropicHardeningModulus() const
{ return K; }

// Return kinematic hardening modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetKinematicHardeningModulus() const
{ return H; }

// Return yield stress for 1D tension test
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetYieldStressFromTensionTest() const
{ return SigmaY; }

// Return Bulk modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetBulkModulus() const
{ return E/(3.*(1.-2.*nu)); }

// Return shear modulus
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetShearModulus() const
{ return 0.5*E/(1.+nu); }

// Return dissipated energy
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetDissipatedEnergy() const
{ 
  if(SigmaY < 0 || ExpEqPlasticStrainRate.size() > 0) {
    std::cerr << " *** WARNING: IsotropicLinearElasticJ2PlasticPlaneStressMaterial::GetDissipatedEnergy is not implemented for nonlinear isotropic hardening.\n";
    return 0;
  }
  return (SigmaY + 0.5*K*equivEPSplastic)*equivEPSplastic;
}

// Return equivalent plastic strain at failure
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetEquivalentPlasticStrainAtFailure() const
{ return equivEPSplasticF; }

// Return tolerance for convergence of nonlinear solve
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
GetTolerance() const
{ return Tol; }

// Set the plastic strain in the material
void IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
SetMaterialPlasticStrain(const std::vector<double> &iEPSplastic)
{ for(int i = 0; i < 3; ++i) EPSplastic[i] = iEPSplastic[i]; }

// Set the equivalent plastic strain in the material
void IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
SetMaterialEquivalentPlasticStrain(double iEquivEPSplastic)
{ equivEPSplastic = iEquivEPSplastic; }

// Set the back stress in the material
void IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
SetMaterialBackStress(const std::vector<double> &iBackStress)
{ for(int i = 0; i < 3; ++i) BackStress[i] = iBackStress[i]; }

// Print all of the internal variables
void IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
Print()
{
  std::cerr << "Plastic Strain = " << EPSplastic[0] << " " << EPSplastic[1] << " " << EPSplastic[2] << std::endl;
  std::cerr << "Back Stress = " << BackStress[0] << " " << BackStress[1] << " " << BackStress[2] << std::endl;
  std::cerr << "Equivalent Plastic Strain = " << equivEPSplastic << std::endl;
}

// Compute the elastic constitutive response
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
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

// Evaluate norm of deviatoric part of \f$\sigma-\sigma^b\f$
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeJ2(const double *Xi) const
{
/* ORIG:
  double P[3][3] = {{2./3.,  -1./3., 0.},
                    {-1./3., 2./3.,  0.},
                    {0., 0., 2.}};

  // Norm of deviatoric part of Xi = Transpose(Xi) * P * Xi.
  double normXi2 = 0.;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      normXi2 += Xi[i]*P[i][j]*Xi[j];

  return sqrt(normXi2);
*/
  return sqrt(2/3.*(Xi[0]*Xi[0] - Xi[0]*Xi[1] + Xi[1]*Xi[1]) + 2*Xi[2]*Xi[2]);
}

// Compute yield stress and isotropic hardening modulus by interpolating experimental stress-strain curve
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
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
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
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

// Evaluate yield function
double IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
EvaluateYieldFunction(const double * Xi, const double DeltaEqP, double &K, const double dt) const
{
  double eqP = equivEPSplastic + DeltaEqP;

  // Norm of Xi
  double J2 = ComputeJ2(Xi);

  // Evaluate radius of yield surface
  double YSrad;
  if(SigmaY >= 0) {
    K = GetIsotropicHardeningModulus();
    YSrad = sqrt(2./3.) * (SigmaY + K*eqP);
  }
  else {
    YSrad = sqrt(2./3.) * GetYieldStressUsingExperimentalCurve(eqP, K);
  }

  // Scale yield stress due to strain-rate dependency
  if(ExpEqPlasticStrainRate.size() > 0 && dt > 0) {
    double R;
    double ScaleF = GetScaleFactorUsingExperimentalCurve(DeltaEqP/dt, R);
    K = sqrt(3./2.)*YSrad*R/dt + ScaleF*K;
    YSrad *= ScaleF;
  }

  return J2 - YSrad;
}

// Check if state of material lies within yield surface
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
CheckMaterialState(const std::vector<double> &CS,
                   const double TOL,
                   const double dt) const
{
  // Compute Xi = Cauchy stress - Back stress (in vector form)
  double Xi[3] = { CS[0]-BackStress[0],
                   CS[4]-BackStress[1],
                   CS[1]-BackStress[2] };

  // Evaluate yield function
  double K;
  double Fval = EvaluateYieldFunction(Xi, 0, K, dt);

  // Check for tolerance after non-dimensionalizing Fval.
  if( (SigmaY >= 0 && Fval < TOL*SigmaY) || (SigmaY < 0 && Fval < TOL*ExpYieldStress[0]) )
    return true;
  else
    return false;
}

// Compute elastoplastic response of material
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1,
                                         std::vector<double> *CauchyStress,
                                         std::vector<double> *Cep,
                                         const bool UpdateFlag,
                                         const double dt)
{
#ifndef USE_EIGEN3
  // Notify that elastoplastic tangents are not computed.
  if( Cep!=0 )
    {
      std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
               <<"Elastoplastic tangents implementation requires Eigen3 library.\n";
      return false;
    }
#endif
  // Resize output tangents if requested
  if( Cep )
    if( int(Cep->size())<9 )
      Cep->resize( 9 );

  // Resize outputs if required
  if( int(CauchyStress->size())<9 )
    CauchyStress->resize(9);

  // Check for failure
  if(equivEPSplastic >= equivEPSplasticF) {
    if( Cep ) for(int i=0; i<9; i++) (*Cep)[i] = 0;
    for(int i=0; i<9; i++) (*CauchyStress)[i] = 0;
    return true;
  }

  // Elastic modulii
  std::vector<double> * Ce = 0;
  if( Cep )
    Ce = new std::vector<double>(9);

  // Vector form of symmetric strain
  double EPS[3] = { Fnp1[0]-1., Fnp1[4]-1., Fnp1[1]+Fnp1[3] };

  // Compute trial elastic state by freezing plastic flow
  std::vector<double> EPSelas(3,0.);
  for(int i=0; i<3; i++)
    EPSelas[i] = EPS[i]-EPSplastic[i];

  // Compute trial Cauchy stress
  std::vector<double> CStrial(3,0.);

  if( !ComputeElasticConstitutiveResponse(EPSelas, &CStrial, Ce) )
    {
      std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
               <<"Could not compute elastic response at trial state.\n";
      return false;
    }

  // Trial Xi
  double Xitrial[] = { CStrial[0]-BackStress[0],
                       CStrial[1]-BackStress[1],
                       CStrial[2]-BackStress[2] };

  // Evaluate yield function at trial state
  double K;
  double Ftrial = EvaluateYieldFunction(Xitrial, 0, K, dt);

  // Use some tolerance for checking yield function value
  // Note: I observe a relationship between TOL and the nltol under NONLINEAR
  // looks like TOL should be at least an order of magnitude smaller than nltol
  double TOL = (SigmaY >= 0) ? SigmaY*Tol : ExpYieldStress[0]*Tol;

  if( Ftrial<0 /*ORIG: TOL*/ )
    {
      // This step is purely elastic.
      // No need to update plastic variables.

      // Convert Cauchy stress to tensor form
      (*CauchyStress)[0] = CStrial[0];
      (*CauchyStress)[4] = CStrial[1];
      (*CauchyStress)[1] = (*CauchyStress)[3] = CStrial[2];
      (*CauchyStress)[2] = (*CauchyStress)[5] = (*CauchyStress)[6] = (*CauchyStress)[7] = (*CauchyStress)[8] = 0.;

      if( Cep )
        for(int i=0; i<9; i++)
          (*Cep)[i] = (*Ce)[i];

    }

  else
    {
      // This step is elasto-plastic.
      // Need to compute consistency parameter, called lambda here.
      // Use bisection method to compute lambda.

      // Max. number of iteration = 1000.
      int nItMax = 1000;
      // Tolerance for checking convergence if f<TOL
      bool CONVERGED = false;

      // lambda = 0 is the trial state and Ftrial is surely > 0.
      double lambda_R = 0.;
      double F_R = Ftrial;

      // The other extreme is when there is plenty of plastic strain.
      double epMAX = equivEPSplastic + sqrt(2./3.)*sqrt( pow(EPS[0]-EPSplastic[0],2.) +
                                                         pow(EPS[1]-EPSplastic[1],2.) +
                                                         0.5*pow(EPS[2]-EPSplastic[2],2.) );  // due the factor of 2 in the strains.

      // Step size for lambda to probe where F < 0.
      double KMAX;
      double lambdaStep = 1.5*(epMAX-equivEPSplastic)/((SigmaY >= 0) ? SigmaY+K*epMAX : GetYieldStressUsingExperimentalCurve(epMAX,KMAX));

      // For computing Xi as a function of lambda.
      double Xi[] = {0.,0.,0.};

      // Proble for lambda_L that corresponds to F < 0.
      double lambda_L = 0.;
      double F_L = 0.;
      bool NegFfound = false;
      int counter = 0;
      while( NegFfound==false )
        {
          // Increment counter
          counter++;
          // Increment lamba by one step
          lambda_L += lambdaStep;
        
          // Evaluate corresponding value of Xi
          if( !ComputeXiGivenConsistencyParameter(Xitrial, lambda_L, Xi) )
            {
              std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                       <<"Could not evaluate stress state for some guess of consistency parameter.\n";
              return false;
            }
        
          // Evaluate yield function
          F_L = EvaluateYieldFunction(Xi, sqrt(2./3.)*lambda_L*ComputeJ2(Xi), K, dt);
        
          if( F_L<=0. )
            // lambda corresponding to negative value of F has been found.
            NegFfound = true;
          else
            {
              // F is still positive.
              // This means F = 0 lies to the left of lambda_L.
              // So shift the right side of the bracket to this point.
              lambda_R = lambda_L;
              F_R      = F_L;
            }
        
          // Sound a warning if lambda is orders of magnitude away from lambdaStep.
          if( counter%100 == 0) {
            std::cerr<<"\n WARNING: "
                     <<"IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                     <<"Consistency parameter exceeds "<<counter<<" times initial estimate.\n";
          }
          // Upper limit on how far to probe- avoid an infinite loop.
          if( counter > 1000 )
            return false;
        }

      // Check that F has been bracketed.
      if( F_L>0. || F_R<0. )
        {
          std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                   <<"Could not bracket yield function to solve for consistency parameter.\n";
          return false;
        }

      // Required solution lies between lambda_L and lambda_R
      double lambda = 0.;

      // Iterate
      for(int it=0; it<nItMax && !CONVERGED; it++)
        {
          lambda = (lambda_L + lambda_R)/2.;
          if( !ComputeXiGivenConsistencyParameter(Xitrial, lambda, Xi) )
            {
              std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                       <<"Could not evaluate stress state for some guess of consistency parameter.\n";
              return false;
            }
          double F = EvaluateYieldFunction(Xi, sqrt(2./3.)*lambda*ComputeJ2(Xi), K, dt);
        
          if(std::abs(F) < TOL || (lambda_L-lambda_R)/2 < 0)  //ORIG: if( std::abs(F) < TOL )
            CONVERGED = true;
          else
            {
              if(!((F < 0 && F_R < 0) || (F > 0 && F_R > 0))) //ORIG: if( F<0. )
                {
                  lambda_L = lambda;
                  F_L      = F;
                }
              else
                {
                  lambda_R = lambda;
                  F_R      = F;
                }
            }
        }

      if( !CONVERGED )
        {
          std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                   <<"Iterations to compute consistency parameter did not converge.\n";
          return false;
        }

      // Compute P*Xi
      double PXi[] = { (2./3.)*Xi[0] - (1./3.)*Xi[1],
                       (2./3.)*Xi[1] - (1./3.)*Xi[0],
                       2.*Xi[2] };

      // Compute new elastic strain
      std::vector<double> EPSelastic(3);
      for(int i=0; i<3; i++)
        EPSelastic[i] = EPS[i] - EPSplastic[i] - lambda*PXi[i];

      // Compute new Cauchy stress
      std::vector<double> CS(3);
      if( !ComputeElasticConstitutiveResponse(EPSelastic, &CS, 0) )
        {
          std::cerr<<"\n IsotropicLinearElasticJ2PlasticPlaneStressMaterial::ComputeElastoPlasticConstitutiveResponse()- "
                   <<"Could not compute elastic response.\n";
          return false;
        }

      // Convert vector to tensor form
      (*CauchyStress)[0] = CS[0];
      (*CauchyStress)[4] = CS[1];
      (*CauchyStress)[1] = (*CauchyStress)[3] = CS[2];
      (*CauchyStress)[2] = (*CauchyStress)[5] = (*CauchyStress)[6] = (*CauchyStress)[7] = (*CauchyStress)[8] = 0.;

      // Check for failure
      if(equivEPSplasticF < std::numeric_limits<double>::infinity() &&
         equivEPSplastic+sqrt(2./3.)*lambda*ComputeJ2(Xi) >= equivEPSplasticF)
        {
          if( Cep ) for(int i=0; i<9; i++) (*Cep)[i] = 0;
          for(int i=0; i<9; i++) (*CauchyStress)[i] = 0;
        }

      else if( Cep )
        {
#ifdef USE_EIGEN3
          // Evaluate the consistent elasto-plastic modulii
          using Eigen::Matrix3d;
          using Eigen::Vector3d;

          // Compute A*P*xi where A is the modified (algorithmic tangent modulus)
          Matrix3d C; C << (*Ce)[0], (*Ce)[1], (*Ce)[2],
                           (*Ce)[3], (*Ce)[4], (*Ce)[5],
                           (*Ce)[6], (*Ce)[7], (*Ce)[8];
          Matrix3d P; P << 2/3.,  -1./3, 0,
                           -1/3., 2/3.,  0,
                           0,     0,     2;
          Matrix3d A = (C.inverse() + lambda/(1 + 2/3.*lambda*H)*P).inverse();
          Vector3d xi; xi << Xi[0], Xi[1], Xi[2];
          Vector3d Pxi = P*xi;
          Vector3d APxi = A*Pxi;

          // Compute N = A*P*xi / sqrt(xi^T*P*A*P*Xi)
          double s2 = Pxi.dot(APxi);
          Vector3d N = APxi/sqrt(s2);

          // Compute beta
          double theta1 = 1 + 2/3.*H*lambda;
          double theta2 = 1 - 2/3.*K*lambda;
          double fbar2 = xi.dot(Pxi);
          double beta = 2/3. * theta1/theta2 * fbar2 * (K*theta1 + H*theta2) / s2;

          for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
              (*Cep)[3*i+j] = A(i,j) - N[i]*N[j]/(1+beta);
#endif
        }

      // If requested, update state of material
      if( UpdateFlag==true )
        {
          for(int i=0; i<3; i++)
            {
              // Update plastic strain
              EPSplastic[i] += lambda*PXi[i];        
              // Update back stress
              BackStress[i] += (2./3.)*H*lambda*Xi[i];
            }
          // Update equivalent plastic strain
          equivEPSplastic += sqrt(2./3.)*lambda*ComputeJ2(Xi);
        }

    }

  if( Ce )
    delete Ce;

  return true;
}

// Given a value of consistency parameter and trial Xi, compute new Xi
bool IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
ComputeXiGivenConsistencyParameter(const double * Xitrial, const double lambda, double * Xi) const
{
/* ORIG:
  double v = lambda*E/(3.*(1.-nu*nu));

  double a = 1.+(2./3.)*lambda*H + v*(2.-nu);
  double b = v*(2.*nu-1.);
  double c = 1.+(2./3.)*lambda*H + v*3.*(1.-nu);

  // M*Xi = Xitrial,
  // where M = (1+(2/3)*lambda*H) + lambda*C*P
  // M = [a b 0]
  //     [b a 0]
  //     [0 0 c]

  double a2mb2 = pow(a,2.)-pow(b,2.);
*/
  double a = 1.+lambda*t00;
  double b =    lambda*t01;
  double c = 1.+lambda*t22;

  double a2mb2 = a*a - b*b;

  // Check that M is invertible.
  if( fabs(c*a2mb2) < 1.e-6 )
    return false;
  else
    {
      Xi[0] = (a*Xitrial[0]-b*Xitrial[1])/a2mb2;
      Xi[1] = (a*Xitrial[1]-b*Xitrial[0])/a2mb2;
      Xi[2] = Xitrial[2]/c;
      return true;
    }
}

// Set the (x,y) values at a point in the yield stress vs. effective plastic strain curve
void IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
SetExperimentalCurveData(double iEqPlasticStrain, double iYieldStress)
{
  ExpEqPlasticStrain.push_back(iEqPlasticStrain);
  ExpYieldStress.push_back(iYieldStress);
}

// Set the (x,y) values at a point in the yield stress scale factor vs. effective plastic strain rate curve
void IsotropicLinearElasticJ2PlasticPlaneStressMaterial::
SetExperimentalCurveData2(double iEqPlasticStrainRate, double iYieldStressScale)
{
  ExpEqPlasticStrainRate.push_back(iEqPlasticStrainRate);
  ExpYieldStressScale.push_back(iYieldStressScale);
}

