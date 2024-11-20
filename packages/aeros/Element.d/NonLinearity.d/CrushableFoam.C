#include <Element.d/NonLinearity.d/CrushableFoam.h>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <cstddef>
#include <cfloat>

#ifdef USE_EIGEN3
void
transfromToSecondPiolaStess(const Eigen::Ref<const Eigen::Matrix3d> F, Tensor_d0s2_Ss12 &stress)
{
  // Converts cauchy stress to second piola-kirchoff.
  // Conversion required as GaussIntgElement uses second pk stress.
  Eigen::Matrix3d sigma;
  Eigen::Matrix3d S;
  stress.assignTo(sigma);

  double J = F.determinant();
  if(J==0) {
    fprintf(stderr, "***ERROR: Deformation gradient has zero determinant.\n");
    exit(-1);
  }

  Eigen::Matrix3d Finv = F.inverse();
  S = J*Finv*sigma*(Finv.transpose());
  stress = S;
}
#endif

void
CrushableFoam::getStress(Tensor *stress, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: CrushableFoam::getStress is not implemented\n";
}

void 
CrushableFoam::getElasticity(Tensor *_tm) const
{
  Tensor_d0s4_Ss12s34 &tm = static_cast<Tensor_d0s4_Ss12s34 &>(*_tm);
  tm.setZero();

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = (nu != 0) ? lambda/nu : E;

  (tm)[0][0] = lambdadivnu*(1-nu);
  (tm)[1][1] = lambdadivnu*(1-2*nu)/2;
  (tm)[2][2] = lambdadivnu*(1-2*nu)/2;
  (tm)[3][3] = lambdadivnu*(1-nu);
  (tm)[4][4] = lambdadivnu*(1-2*nu)/2;
  (tm)[5][5] = lambdadivnu*(1-nu);
  (tm)[0][3] = lambdadivnu*nu;
  (tm)[3][0] = lambdadivnu*nu;
  (tm)[0][5] = lambdadivnu*nu;
  (tm)[5][0] = lambdadivnu*nu;
  (tm)[3][5] = lambdadivnu*nu;
  (tm)[5][3] = lambdadivnu*nu;
}

void 
CrushableFoam::getTangentMaterial(Tensor *tm, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: CrushableFoam::getTangentMaterial is not implemented\n";
}

void 
CrushableFoam::getStressAndTangentMaterial(Tensor *stess, Tensor *tm, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: CrushableFoam::getStressAndTangentMaterial is not implemented\n";
}

void 
CrushableFoam::updateStates(Tensor &en, Tensor &enp, double *state, double temp)
{
  std::cerr << "WARNING: CrushableFoam::updateStates is not implemented\n";
}

void
CrushableFoam::integrate(Tensor *_stress, Tensor *_tm, Tensor &_en, Tensor  &_enp,
                                 double *staten, double *statenp, double temp, Tensor *cache, double dt) const
{
  //////////////////////////////////////////////////////////////////////////////
  /// A.Reyes et. al. Implementation of Constitutive Model for Aluminum Foam ///
  //////////////////////////////////////////////////////////////////////////////
  using std::sqrt;
  using std::pow;

#ifdef USE_EIGEN3
  Tensor_d0s4_Ss12s34 &tm = static_cast<Tensor_d0s4_Ss12s34 &>(*_tm);
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress);
  
  Eigen::Matrix<double,3,3,Eigen::RowMajor> Fnp;
  Eigen::Matrix<double,3,3,Eigen::RowMajor> L;
  try
  {
    Fnp = dynamic_cast<Tensor_d1s2_full &>(*cache)[0].matrix();
  }
  catch(const std::exception& e)
  {
    std::cerr << "***ERROR: Deformation tensor is required in stress evaluation routine for CrushableFoam materials.\n";
    exit(-1);
  }
  try {
    L = dynamic_cast<Tensor_d1s2_full &>(*cache)[1].matrix();
  }
  catch(const std::exception& e) {
    std::cerr << "***ERROR: Velocity gradient tensor is required in stress evaluation routine for CrushableFoam materials.\n";
    exit(-1);
  }
  
  // check for failure
  if(statenp != 0 && statenp[12] >= epsF) {
    stress.setZero();
    return;
  }

  // subtract thermal strain
  double e0 = (temp-Tref)*alpha;
  enp[0] -= e0;
  enp[3] -= e0;
  enp[5] -= e0;

  if(statenp == 0) {
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

    tm[0][0] = lambdadivnu*(1-nu);
    tm[1][1] = lambdadivnu*(1-2*nu)/2;
    tm[2][2] = lambdadivnu*(1-2*nu)/2;
    tm[3][3] = lambdadivnu*(1-nu);
    tm[4][4] = lambdadivnu*(1-2*nu)/2;
    tm[5][5] = lambdadivnu*(1-nu);
    tm[0][3] = lambdadivnu*nu;
    tm[3][0] = lambdadivnu*nu;
    tm[0][5] = lambdadivnu*nu;
    tm[5][0] = lambdadivnu*nu;
    tm[3][5] = lambdadivnu*nu;
    tm[5][3] = lambdadivnu*nu;

    stress = tm||enp;

    // convert 2nd Piola Kirchoff stress to Cauchy stress
    double J = Fnp.determinant();
    if(J==0) {
      fprintf(stderr, "***ERROR: Deformation gradient has zero determinant.\n");
      exit(-1);
    }

    Eigen::Matrix3d S;
    stress.assignTo(S);

    Eigen::Matrix3d cauchy = (1./J)*Fnp*S*Fnp.transpose();
    statenp[6]  = cauchy(0,0); statenp[7]  = cauchy(0,1);  statenp[8]  = cauchy(0,2);
                               statenp[9]  = cauchy(1,1); statenp[10]  = cauchy(1,2);
    	                                                  statenp[11]  = cauchy(2,2);
    
  }
  else {
    //state: from 0 to 5, plastic rate of deformation; from 6 to 11, cauchy stress at t_n; 12 equivalent plastic strain

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double mu = E/(2*(1+nu));
    double bulk = lambda + (2./3)*mu;
    double v, s; // value and slope of yield stress vs. effective plastic strain curve

    Eigen::Matrix3d stressn;
    Eigen::Matrix3d strialnp;
    Eigen::Matrix3d strialdevnp;
    Eigen::Matrix3d D;
    Eigen::Matrix3d W;

    stressn << staten[6],  staten[7],  staten[8],
               staten[7],  staten[9], staten[10],
               staten[8], staten[10], staten[11];
 
    D = 0.5*(L + L.transpose()); // rate of deformation tensor
    W = 0.5*(L - L.transpose()); // spin tensor

    // assuming D is completely elastic
    strialnp = stressn + dt*(W*stressn + stressn*(W.transpose())) 
             + dt*(lambda*D.trace()*Eigen::Matrix3d::Identity() + 2*mu*D); 

    strialdevnp = strialnp - (strialnp.trace()/3)*Eigen::Matrix3d::Identity();

    double strialmnp = (1./3)*strialnp.trace();
    double strialVMnp = sqrt((3./2)*(strialdevnp*strialdevnp).trace());
    double stressFactor = 1 + pow(alphaDF/3, 2);

    double sigmatrialnp = sqrt(1./stressFactor)*sqrt(pow(strialVMnp, 2) + pow(alphaDF*strialmnp, 2));

    double trialyieldnp;
    if(!ysst) {
      // Yield Criterion from Deshpande and Fleck -- behavior beyond densification not defined
      bool hasReachedDensification = std::abs(staten[12]-epsD) < 1e-6;
      trialyieldnp = (hasReachedDensification) ? -DBL_MAX : 
                     sigmatrialnp - (sigP+gamma*(staten[12]/epsD)-alpha2*log(1-pow(staten[12]/epsD, beta)));
    }
    else {
      ysst->getValAndSlopeAlt2(staten[12], &v, &s);
      trialyieldnp = sigmatrialnp - v;
    }

    if (trialyieldnp <= 0) {
 
      //fprintf(stderr, "je suis dans la yield surface\n");
      getElasticity(_tm);

      stress = strialnp;
      for (int i=0; i<6; ++i) {
        statenp[i] = staten[i];
        statenp[i+6] = stress[i];
      }

    }
    else {

      //fprintf(stderr, "je suis en dehors de la yield surface\n");

      int i,j,k,l;
      double A1 = (3*mu)/stressFactor;
      double A2 = (bulk*pow(alphaDF, 2))/stressFactor;
      double dsigmahat; // derivative of hardeining fuction and sigmaHat wrt deltahateps
      double sigmnp, sigVMnp, compf;

      double deltahateps = 0; 
      if(ysst || dt > 0) {
        // Newton-Raphson algorithm for solution of the return mapping equation
        for (i=0; i<20; ++i) { 
          if(ysst) ysst->getValAndSlopeAlt2(staten[12]+deltahateps, &v, &s);
          else {
            compf = (staten[12] + deltahateps)/epsD;
            v = sigP + gamma*compf - alpha2*log(1 - pow(compf, beta));
            s = gamma/epsD + (alpha2*beta/epsD)/(pow(compf, 1- beta) - compf);
          }

	  sigmnp = strialmnp*v/(v+A2*deltahateps); // mean stress at n+1
          sigVMnp = strialVMnp*v/(v+A1*deltahateps); // von mises stress at n+1
          sigmatrialnp = sqrt((1./stressFactor)*(pow(sigVMnp,2)+pow(alphaDF*sigmnp,2)));
          trialyieldnp = sigmatrialnp - v;

          if(isnan(trialyieldnp)) { // nan means foam has densified; behavior beyond this point not defined
	    // AN: might want to implement volumetric failure instead
	    deltahateps =  epsD - staten[12];
	    break;
          }
          if(std::abs(trialyieldnp) <= tol*(ysst ? ysst->getVal(0) : sigP)) break;

          dsigmahat = s*sigmatrialnp/v - (1/(stressFactor*sigmatrialnp))*(pow(sigVMnp, 2)*(s + A1)/(v + A1*deltahateps) + 
                      pow(alphaDF*sigmnp, 2)*(s + A2)/(v + A2*deltahateps));

          deltahateps += trialyieldnp/(s - dsigmahat);
        }
      }
   
      if(i==20)
        fprintf(stderr, " ***WARNING: Newton iteration for plastic multiplier did not converge.\n");

      Tensor_d0s2_Ss12 normalnp;
      normalnp = (3./2)*(1./strialVMnp)*strialdevnp;
      double deltaepsv = (pow(alphaDF, 2)*sigmnp*deltahateps)/(stressFactor*v);
      double deltaepsVM = (sigVMnp*deltahateps)/(stressFactor*v);

      statenp[12] = staten[12] + deltahateps;

      stress = strialnp;
      for (i=0; i<6; ++i) {
        statenp[i] = staten[i] + (deltaepsVM/dt)*normalnp[i];
        stress[i] -= 2*mu*deltaepsVM*normalnp[i]; 
 
 	if(i==0 || i==3 || i==5) {
          statenp[i] += (1./3)*(deltaepsv/dt);
          stress[i] -= bulk*deltaepsv;
        }

        statenp[6+i] = stress[i];
      }

    }

    // TODO: Tangent modulus (tm) not updated ... matrix not needed for explicit dynamic analysis.

    transfromToSecondPiolaStess(Fnp, stress);

  }

#else
  std::cerr << "***ERROR: CrushableFoam::Integrate requires Eigen.\n";
  exit(-1);
#endif

}

void
CrushableFoam::integrate(Tensor *_stress, Tensor &_en, Tensor  &_enp,
                                 double *staten, double *statenp, double temp, Tensor *cache, double dt) const
{
  //////////////////////////////////////////////////////////////////////////////
  /// A.Reyes et. al. Implementation of Constitutive Model for Aluminum Foam ///
  //////////////////////////////////////////////////////////////////////////////

  using std::sqrt;
  using std::pow;

#ifdef USE_EIGEN3
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > Fnp = static_cast<Tensor_d1s2_full &>(*cache)[0].matrix();
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > L = static_cast<Tensor_d1s2_full &>(*cache)[1].matrix();

  // check for failure
  if(statenp != 0 && statenp[12] >= epsF) {
    stress.setZero();
    return;
  }

  // subtract thermal strain
  double e0 = (temp-Tref)*alpha;
  enp[0] -= e0;
  enp[3] -= e0;
  enp[5] -= e0;

  if(statenp == 0) {

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

    Tensor_d0s4_Ss12s34 tm;

    tm[0][0] = lambdadivnu*(1-nu);
    tm[1][1] = lambdadivnu*(1-2*nu)/2;
    tm[2][2] = lambdadivnu*(1-2*nu)/2;
    tm[3][3] = lambdadivnu*(1-nu);
    tm[4][4] = lambdadivnu*(1-2*nu)/2;
    tm[5][5] = lambdadivnu*(1-nu);
    tm[0][3] = lambdadivnu*nu;
    tm[3][0] = lambdadivnu*nu;
    tm[0][5] = lambdadivnu*nu;
    tm[5][0] = lambdadivnu*nu;
    tm[3][5] = lambdadivnu*nu;
    tm[5][3] = lambdadivnu*nu;

    stress = tm||enp;

    // convert 2nd Piola Kirchoff stress to Cauchy stress
    double J = Fnp.determinant();
    if(J==0) {
      fprintf(stderr, "***ERROR: Deformation gradient has zero determinant.\n");
      exit(-1);
    }

    Eigen::Matrix3d S;
    stress.assignTo(S);

    Eigen::Matrix3d cauchy = (1./J)*Fnp*S*Fnp.transpose();
    statenp[6]  = cauchy(0,0); statenp[7]  = cauchy(0,1);  statenp[8]  = cauchy(0,2);
                               statenp[9]  = cauchy(1,1); statenp[10]  = cauchy(1,2);
    							                                        statenp[11]  = cauchy(2,2);

  }
  else {
    //state: from 0 to 5, plastic rate of deformation; from 6 to 11, cauchy stress at t_n; 12 equivalent plastic strain

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double mu = E/(2*(1+nu));
    double bulk = lambda + (2./3)*mu;
    double v, s; // value and slope of yield stress vs. effective plastic strain curve

    // NOTE: might need to use E_foam and nu_foam instead of base aluminum properties for bulk and shear modulus.
    Eigen::Matrix3d stressn;
    Eigen::Matrix3d strialnp;
    Eigen::Matrix3d strialdevnp;
    Eigen::Matrix3d D;
    Eigen::Matrix3d W;

    stressn << staten[6],  staten[7],  staten[8],
               staten[7],  staten[9], staten[10],
               staten[8], staten[10], staten[11];
 
    D = 0.5*(L + L.transpose()); // rate of deformation tensor
    W = 0.5*(L - L.transpose()); // spin tensor

    // assuming D is completely elastic
    strialnp = stressn + dt*(W*stressn + stressn*(W.transpose())) 
             + dt*(lambda*D.trace()*Eigen::Matrix3d::Identity() + 2*mu*D); 

    strialdevnp = strialnp - (strialnp.trace()/3)*Eigen::Matrix3d::Identity();

    double strialmnp = (1./3)*strialnp.trace();
    double strialVMnp = sqrt((3./2)*(strialdevnp*strialdevnp).trace());
    double stressFactor = 1 + pow(alphaDF/3, 2);

    double sigmatrialnp = sqrt(1./stressFactor)*sqrt(pow(strialVMnp, 2) + pow(alphaDF*strialmnp, 2));

    double trialyieldnp;
    if(!ysst) {
      // Yield Criterion from Deshpande and Fleck -- behavior beyond densification not defined
      bool hasReachedDensification = std::abs(staten[12]-epsD) < 1e-6;
      trialyieldnp = (hasReachedDensification) ? -DBL_MAX : 
                     sigmatrialnp - (sigP+gamma*(staten[12]/epsD)-alpha2*log(1-pow(staten[12]/epsD, beta)));
    }
    else {
      ysst->getValAndSlopeAlt2(staten[12], &v, &s);
      trialyieldnp = sigmatrialnp - v;
    }

    if (trialyieldnp <= 0) {
 
      //fprintf(stderr, "je suis dans la yield surface\n");

      stress = strialnp;
      for (int i=0; i<6; ++i) {
        statenp[i] = staten[i];
        statenp[i+6] = stress[i];
      }

    }
    else {

      //fprintf(stderr, "je suis en dehors de la yield surface\n");

      int i,j,k,l;
      double A1 = (3*mu)/stressFactor;
      double A2 = (bulk*pow(alphaDF, 2))/stressFactor;
      double dsigmahat; // derivative of hardeining fuction and sigmaHat wrt deltahateps
      double sigmnp, sigVMnp, compf;

      double deltahateps = 0;
      if(ysst || dt > 0) {
        // Newton-Raphson algorithm for solution of the return mapping equation
        for (i=0; i<20; ++i) { 
          if(ysst) ysst->getValAndSlopeAlt2(staten[12]+deltahateps, &v, &s);
          else {
            compf = (staten[12] + deltahateps)/epsD;
            v = sigP + gamma*compf - alpha2*log(1 - pow(compf, beta));
            s = gamma/epsD + (alpha2*beta/epsD)/(pow(compf, 1- beta) - compf);
          }

	        sigmnp = strialmnp*v/(v+A2*deltahateps); // mean stress at n+1
          sigVMnp = strialVMnp*v/(v+A1*deltahateps); // von mises stress at n+1
          sigmatrialnp = sqrt((1./stressFactor)*(pow(sigVMnp,2)+pow(alphaDF*sigmnp,2)));
          trialyieldnp = sigmatrialnp - v;

          if(isnan(trialyieldnp)) { // nan means foam has densified; behavior beyond this point not defined
	    // AN: might want to implement volumetric failure instead
	    deltahateps =  epsD - staten[12];
	    break;
	  }
          if(std::abs(trialyieldnp) <= tol*(ysst ? ysst->getVal(0) : sigP)) break;

          dsigmahat = s*sigmatrialnp/v - (1/(stressFactor*sigmatrialnp))*(pow(sigVMnp, 2)*(s + A1)/(v + A1*deltahateps) + 
                      pow(alphaDF*sigmnp, 2)*(s + A2)/(v + A2*deltahateps));

          deltahateps += trialyieldnp/(s - dsigmahat);
        }
   
        if(i==20)
          fprintf(stderr, " ***WARNING: Newton iteration for plastic multiplier did not converge.\n");

        Tensor_d0s2_Ss12 normalnp;
        normalnp = (3./2)*(1./strialVMnp)*strialdevnp; 
        double deltaepsv = (pow(alphaDF, 2)*sigmnp*deltahateps)/(stressFactor*v);
        double deltaepsVM = (sigVMnp*deltahateps)/(stressFactor*v);

        statenp[12] = staten[12] + deltahateps;

        stress = strialnp;
        for (i=0; i<6; ++i) {
          statenp[i] = staten[i] + (deltaepsVM/dt)*normalnp[i];
          stress[i] -= 2*mu*deltaepsVM*normalnp[i]; 
  
          if(i==0 || i==3 || i==5) {
            statenp[i] += (1./3)*(deltaepsv/dt);
            stress[i] -= bulk*deltaepsv;
          }

          statenp[6+i] = stress[i];
        }

      }
      else if (dt==0) { // used for post processing --> dt=0;
        for (int i=0; i<6; ++i)
          stress[i] = statenp[6+i]; // should contain updated cauchy stress (t_n+1)
      }

    }

    transfromToSecondPiolaStess(Fnp, stress);

  }

#else
  std::cerr << "***ERROR: CrushableFoam::Integrate requires Eigen.\n";
  exit(-1);
#endif

}

void 
CrushableFoam::initStates(double *st)
{
  for (int i=0; i<13; ++i)
    st[i] = 0;
}

extern GreenLagrangeStrain greenLagrangeStrain;

StrainEvaluator *
CrushableFoam::getStrainEvaluator() const
{
  return &greenLagrangeStrain;
} 

bool
CrushableFoam::getBackStress(double *statenp, Tensor *_backstress)
{
  std::cerr << "WARNING: CrushableFoam::getBackStress is not implemented\n";
  return true;
}

bool
CrushableFoam::getPlasticStrain(double *statenp, Tensor *_plasticstrain)
{
  // Returns plastic part of rate of deformation tensor (D).

  Tensor_d0s2_Ss12 *plasticstrain = static_cast<Tensor_d0s2_Ss12 *>(_plasticstrain);
  for (int i=0; i<6; ++i) {
    (*plasticstrain)[i] = statenp[i];
  }

  return true;
}

double
CrushableFoam::getStrainEnergyDensity(Tensor &_enp, double *statenp, double temp)
{
/*
  std::cerr << " *** WARNING: CrushableFoam::getStrainEnergyDensity is not implemented.\n";
  return 0;
*/

  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 eplastnp;
  Tensor_d0s2_Ss12 eelastnp;

  eplastnp.buildTensorOf(statenp);
  eelastnp = enp-eplastnp;

  double lambda = E*nu/((1+nu)*(1-2*nu));
  double mu = E/(2*(1+nu));

  double I1 = eelastnp.getTrace();
  return lambda/2*I1*I1 + mu*eelastnp.innerProduct();
}

double
CrushableFoam::getDissipatedEnergy(double *statenp)
{
  // fprintf(stdout, "In CrushableFoam::getDissipatedEnergy()\n");
  if(ysst) {
    std::cerr << " *** WARNING: CrushableFoam::getDissipatedEnergy is not implemented for the YSST card.\n";
    return 0;
  }
  else {
    double compactionFactor = statenp[12]/epsD;
    double v = sigP + gamma*compactionFactor - alpha2*log(1 - pow(compactionFactor, beta));
    return statenp[12]*v;
  }
}

inline void
CrushableFoam::print(std::ostream &out) const
{
  out << "CrushableFoam " << rhoF << " " << E << " " << nu << " " << sigP << " " << alpha2 << " " << gamma << " ";
  out << beta << " " << epsD;
}

void
CrushableFoam::print2(std::ostream &out) const
{
  if(epsF != std::numeric_limits<double>::infinity() || tol != 1e-6) {
    out << " " << epsF << " " << tol;
  }
}
