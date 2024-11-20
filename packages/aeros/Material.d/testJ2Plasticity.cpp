// Sriramajayam

/*
 * testJ2Plasticity.cpp
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
#include <fstream>

// Print material constants and internal variables of material.
// Also provide what these values should be.
void PrintMaterialDetails(const IsotropicLinearElasticJ2PlasticMaterial &Mat, 
			  const double Kappa, const double Mu, 
			  const double SigmaY, const double K, const double H, 
			  const double equivEPSplastic,
			  const std::vector<double> EPSplastic, 
			  const std::vector<double> BackStress);


// Load and unload material.
void PerformLoadUnloadTest(IsotropicLinearElasticJ2PlasticMaterial &Mat, char * filename);

int main()
{
  // Material constants (all in MPa)
  double E = 150.e3;
  double Mu = 60.e3;
  double SigmaY = 300.;
  double K = 0.; 
  double H = 0.; 
  double Lambda = Mu*(E-2.*Mu)/(3.*Mu-E);
  double Kappa = Lambda + (2./3.)*Mu;
  
  std::vector<double> ZeroVec(9,0.);
  for(int i=0; i<9; i++)
    ZeroVec[i] = 0.;
  
  // Create elasto-plastic material
  IsotropicLinearElasticJ2PlasticMaterial Mat(Lambda, Mu, SigmaY, K, H);
  
  // Print details of material
  std::cout<<"\n\nDetails of material:";
  PrintMaterialDetails(Mat, Kappa, Mu, SigmaY, K, H, 0., ZeroVec, ZeroVec);

  // Test copy constructor
  IsotropicLinearElasticJ2PlasticMaterial Copy(Mat);
  std::cout<<"\n\nDetails of copy of material (testing copy constructor):";
  PrintMaterialDetails(Copy, Kappa, Mu, SigmaY, K, H, 0., ZeroVec, ZeroVec);
  
  // Test cloning 
  IsotropicLinearElasticJ2PlasticMaterial * CL = Mat.Clone();
  std::cout<<"\n\nDetails of clone of material (testing cloning):";
  PrintMaterialDetails(*CL, Kappa, Mu, SigmaY, K, H, 0., ZeroVec, ZeroVec);
  
  // Perform a loading-unloading test
  std::cout<<"\n\nSimulating loading-unloading tests.";
  
  // Without any hardening
  IsotropicLinearElasticJ2PlasticMaterial Mat0(Lambda, Mu, SigmaY,0.,0.); 
  PerformLoadUnloadTest(Mat0, (char*)"NoHardening.dat");
  
  // With isotropic hardening
  IsotropicLinearElasticJ2PlasticMaterial MatI(Lambda, Mu, SigmaY, 100.e3, 0.);
  PerformLoadUnloadTest(MatI, (char*)"IsoHardening.dat");
  
  // With kinematic hardening
  IsotropicLinearElasticJ2PlasticMaterial MatK(Lambda, Mu, SigmaY, 0., 100.e3);
  PerformLoadUnloadTest(MatK, (char*)"KinHardening.dat");  
  
  // With isotropic and kinematic hardening
  IsotropicLinearElasticJ2PlasticMaterial MatIK(Lambda, Mu, SigmaY, 100.e3, 100.e3); 
  PerformLoadUnloadTest(MatIK, (char*)"IsoKinHardening.dat");
    
  std::cout<<"\n\nChecks: ";
  std::cout<<"\nInitial yield strain should equal "<<SigmaY/(2.*Mu);
  std::cout<<"\nInitial yield stress should equal "<<SigmaY;
  std::cout<<"\nSlope of initial loading curve should equal "<<2.*Mu;
  
  std::cout<<"..done.\n\n";
 }




void PrintMaterialDetails(const IsotropicLinearElasticJ2PlasticMaterial &Mat,
			  const double Kappa, const double Mu, 
			  const double SigmaY, const double K, 
			  const double H, const double equivEPSplastic,
			  const std::vector<double> EPSplastic, 
			  const std::vector<double> BackStress)
{
  std::cout<<"\n\nBulk modulus = "<<Mat.GetBulkModulus()
	   <<"\nShould read as "<<Kappa;
  
  std::cout<<"\n\nShear modulus = "<<Mat.GetShearModulus()
	   <<"\nShould read as "<<Mu;
  
  std::cout<<"\n\nIsotropic hardening modulus = "<<Mat.GetIsotropicHardeningModulus()
	   <<"\nShould read as "<<K;
  
  std::cout<<"\n\nKinematic hardening modulus = "<<Mat.GetKinematicHardeningModulus()
	   <<"\nShould read as "<<H;
  
  std::cout<<"\n\nYield stress from uniaxial tension test = "<<Mat.GetYieldStressFromTensionTest()
	   <<"\nShould read as "<<SigmaY;
  
  std::cout<<"\n\nPlastic strain in material = (";
  for(int i=0; i<9; i++)
    std::cout<<" "<<Mat.GetMaterialPlasticStrain()[i];
  std::cout<<").\nShould read as (";
  for(int i=0; i<9; i++)
    std::cout<<" "<<EPSplastic[i];
  std::cout<<")";
  
  std::cout<<"\n\nEquivalent plastic strain in material = "<<Mat.GetMaterialEquivalentPlasticStrain()
	   <<"\nShould read as "<<equivEPSplastic;
  
  std::cout<<"\n\nBack stress in material = (";
  for(int i=0; i<9; i++)
    std::cout<<" "<<Mat.GetMaterialBackStress()[i];
  std::cout<<").\nShould read as (";
  for(int i=0; i<9; i++)
    std::cout<<" "<<BackStress[i];
  std::cout<<")";
  
  std::cout<<"\n";
}


void PerformLoadUnloadTest(IsotropicLinearElasticJ2PlasticMaterial &Mat, char * filename)
{
  // Total strain at the end of test
  double TotStrain = 5.*1.e-3;
  
  // Number of loading steps
  int nSteps = 40;
  
  // Strain increment history
  std::vector<double> FxxHistory(4*nSteps+1);
 
  // Initially no strain
  FxxHistory[0] = 1.;

  // Load upto +Total strain
  for(int n=1; n<nSteps; n++)
    FxxHistory[n] =  FxxHistory[n-1]+TotStrain/double(nSteps);
  
  // Unload to -Total strain
  for(int n=nSteps; n<3*nSteps; n++)
    FxxHistory[n] = FxxHistory[n-1]-TotStrain/double(nSteps);
  
  // Load to zero strain
  for(int n=3*nSteps; n<=4*nSteps; n++)
    FxxHistory[n] = FxxHistory[n-1]+TotStrain/double(nSteps);

  // Write data to file
  std::fstream plot;
  plot.open(filename, std::ios::out);
  plot<<"# %Strain vs Sigma_xx-Sigma_yy.";
  
  // Load-unload material incrementally
  std::vector<double> F(9);
  for(int i=0; i<9; i++)
    F[i] = 0.;
  F[0] = F[4] = F[8] = 1.;
  for(unsigned int step=0; step<FxxHistory.size(); step++)
    {
      // Deformation gradient at this loading step
      F[0] = FxxHistory[step];
      
      // Compute new Cauchy stress
      std::vector<double> CauchyStress(9);
      Mat.ComputeElastoPlasticConstitutiveResponse(F, &CauchyStress);
      plot <<"\n"<<(FxxHistory[step]-1.)*100.
	   <<"\t"<<CauchyStress[0]-CauchyStress[4]
	   <<"\t"<<Mat.GetMaterialEquivalentPlasticStrain();
      
      // Check if material state is acceptable
      if( !Mat.CheckMaterialState(CauchyStress) )
	std::cerr<<"\nMaterial state UNACCEPTABLE :( \n";
    }
  plot.close();
}
  


