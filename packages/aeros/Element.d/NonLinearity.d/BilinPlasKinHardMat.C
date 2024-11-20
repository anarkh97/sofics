#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <cstddef>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

template <class T> inline double delta(T a, T b) { return (a==b) ? 1.0 : 0.0; }

template<int e>
ElasPlasKinHardMat<e>::ElasPlasKinHardMat(StructProp *p)
{
  E = p->E;
  nu = p->nu;
  Ep = p->Ep;         // Tangent modulus from the strain-stress curve
  sigE = p->sigE;     // Yield equivalent stress
  theta = 0;
  Tref = p->Ta;
  alpha = p->W;
  epsF = std::numeric_limits<double>::infinity();
  ysst = NULL;
  yssrt = NULL;
}

template<int e>
void
ElasPlasKinHardMat<e>::getStress(Tensor *stress, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: ElasPlasKinHardMat<e>::getStress is not implemented\n";
}

template<int e>
void 
ElasPlasKinHardMat<e>::getElasticity(Tensor *_tm) const
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

template<int e>
void 
ElasPlasKinHardMat<e>::getTangentMaterial(Tensor *tm, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: ElasPlasKinHardMat<e>::getTangentMaterial is not implemented\n";
}

template<int e>
void 
ElasPlasKinHardMat<e>::getStressAndTangentMaterial(Tensor *stess, Tensor *tm, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: ElasPlasKinHardMat<e>::getStressAndTangentMaterial is not implemented\n";
}

template<int e>
void 
ElasPlasKinHardMat<e>::updateStates(Tensor &en, Tensor &enp, double *state, double temp)
{
  std::cerr << "WARNING: ElasPlasKinHardMat<e>::updateStates is not implemented\n";
}

template<int e>
void
ElasPlasKinHardMat<e>::integrate(Tensor *_stress, Tensor *_tm, Tensor &_en, Tensor  &_enp,
                                 double *staten, double *statenp, double temp, Tensor *, double dt) const
{
  //////////////////////////////////////////////////////////////////////////////
  /// Simo and Hughes - Computational Inelasticity - Springer -1998- (p:124) ///
  //////////////////////////////////////////////////////////////////////////////

  // theta == 0 corresponds to Kinematic hardening and theta == 1 to isotropic 
  // hardening. Note: this is now a member variable with default value 0

  Tensor_d0s4_Ss12s34 &tm = static_cast<Tensor_d0s4_Ss12s34 &>(*_tm);
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress);

  // check for failure
  if(statenp != 0 && statenp[12] >= epsF) {
    stress.setZero();
    tm.setZero();
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
  }
  else {
    //state: from 0 to 5, plastic strain; from 6 to 11, center of the yield surface; 12 equivalent plastic strain

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double mu = E/(2*(1+nu));
    double bulk = lambda + (2./3)*mu;
    double Hprime = (E != Ep) ? (E*Ep)/(E-Ep) : std::numeric_limits<double>::max();
    double Hk = (1-theta)*Hprime; // kinematic hardening modulus
    double v, s; // value and slope of yield stress vs. effective plastic strain curve
    double v2, s2; // value and slope of yield stress scaling factor vs. effective plastic strain rate curve

    Tensor_d0s2_Ss12 betan;
    Tensor_d0s2_Ss12 edevnp;
    Tensor_d0s2_Ss12 eplastn;
    Tensor_d0s2_Ss12 temp;
    Tensor_d0s2_Ss12 eplastdevn;
    Tensor_d0s2_Ss12 strialnp;
    Tensor_d0s2_Ss12 xitrialnp;

    eplastn.buildTensorOf(staten);
    betan.buildTensorOf(staten+6);

    enp.getDeviation(edevnp);
    eplastn.getDeviation(eplastdevn);

    strialnp = 2*mu*(edevnp - eplastdevn);
    xitrialnp = strialnp - betan;

    double xitrialnpnorm = sqrt(xitrialnp.innerProduct());

    double trialyieldnp;
    if(!ysst) trialyieldnp = xitrialnpnorm - sqrt(2./3)*(sigE+theta*Hprime*staten[12]); // Yield Criterion (Simo & Hughes eq. 3.3.6)
    else {
      ysst->getValAndSlopeAlt2(staten[12], &v, &s);
      trialyieldnp = xitrialnpnorm - sqrt(2./3)*v;
      Hprime = Hk + s;
    }

    if (trialyieldnp <= 0) {
 
      //fprintf(stderr, "je suis dans la yield surface\n");

      for (int i=0; i<13; ++i) {
        statenp[i] = staten[i] ;
      }

      getElasticity(_tm);

      temp = (enp - eplastn);
      stress = tm || temp;
    }
    else {

      //fprintf(stderr, "je suis en dehors de la yield surface\n");

      int i,j,k,l;

      double plastmult = trialyieldnp/(2*mu+(2./3)*Hprime); 
      if(ysst || (yssrt && dt > 0)) {
        // Newton-Raphson algorithm for solution of the return mapping equation
        for (i=0; i<20; ++i) { 
          if(ysst) ysst->getValAndSlopeAlt2(staten[12]+sqrt(2./3)*plastmult, &v, &s);
          else {
            v = sigE+theta*Hprime*(staten[12]+sqrt(2./3)*plastmult);
            s = theta*Hprime;
          }
          if(yssrt && dt > 0) {
            yssrt->getValAndSlopeAlt2(sqrt(2./3)*plastmult/dt, &v2, &s2);
            s = v*s2/dt + v2*s;
            v *= v2;
          }
          trialyieldnp = xitrialnpnorm - (2*mu + 2./3*Hk)*plastmult - sqrt(2./3)*v;
          Hprime = Hk + s;
          if(std::abs(trialyieldnp) <= tol*(ysst ? ysst->getVal(0) : sigE)) break;
          plastmult += trialyieldnp/(2*mu + (2./3)*Hprime);
        }
      }
    
      Tensor_d0s2_Ss12 normalnp = (1./xitrialnpnorm)*xitrialnp;
      double thetanp = 1 - 2*mu*plastmult/xitrialnpnorm;
      double thetaprimenp = 1/(1+Hprime/(3*mu)) - (1 - thetanp);

      statenp[12] = staten[12] + sqrt(2./3)*plastmult;

      for (i=0; i<6; ++i) {
        statenp[i] = staten[i] + plastmult*normalnp[i];
        statenp[6+i] = staten[6+i] + sqrt(2./3)*Hk*(statenp[12]-staten[12])*normalnp[i];   
      }

      stress = ((bulk*3*(enp - edevnp) + strialnp) - (2*mu*plastmult*normalnp));

      for (i=0; i<3; ++i)
        for (j=i; j<3; ++j)
          for (k=0; k<3; ++k)
            for (l=k; l<3; ++l)
              tm[i*(5-i)/2+j][k*(5-k)/2+l] = (bulk*delta(i,j)*delta(k,l))
               +(2*mu*thetanp*((delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))/2.-delta(i,j)*delta(k,l)/3.))
               - 2*mu*thetaprimenp*normalnp[i*(5-i)/2+j]*normalnp[k*(5-k)/2+l];
    }
  }
}

template<int e>
void
ElasPlasKinHardMat<e>::integrate(Tensor *_stress, Tensor &_en, Tensor  &_enp,
                                 double *staten, double *statenp, double temp, Tensor *, double dt) const
{
  //////////////////////////////////////////////////////////////////////////////
  /// Simo and Hughes - Computational Inelasticity - Springer -1998- (p:124) ///
  //////////////////////////////////////////////////////////////////////////////

  // theta == 0 corresponds to Kinematic hardening and theta == 1 to isotropic 
  // hardening. Note: this is now a member variable with default value 0

  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress);

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
  }
  else {
    //state: from 0 to 5, plastic strain; from 6 to 11, center of the yield surface; 12 equivalent plastic strain

    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double mu = E/(2*(1+nu));
    double bulk = lambda + (2./3)*mu;
    double Hprime = (E != Ep) ? (E*Ep)/(E-Ep) : std::numeric_limits<double>::max();
    double Hk = (1-theta)*Hprime; // kinematic hardening modulus
    double v, s; // value and slope of yield stress vs. effective plastic strain curve
    double v2, s2; // value and slope of yield stress scaling factor vs. effective plastic strain rate curve

    Tensor_d0s2_Ss12 betan;
    Tensor_d0s2_Ss12 edevnp;
    Tensor_d0s2_Ss12 eplastn;
    Tensor_d0s2_Ss12 temp;
    Tensor_d0s2_Ss12 eplastdevn;
    Tensor_d0s2_Ss12 strialnp;
    Tensor_d0s2_Ss12 xitrialnp;

    eplastn.buildTensorOf(staten);
    betan.buildTensorOf(staten+6);

    enp.getDeviation(edevnp);
    eplastn.getDeviation(eplastdevn);

    strialnp = 2*mu*(edevnp - eplastdevn);
    xitrialnp = strialnp - betan;

    double xitrialnpnorm = sqrt(xitrialnp.innerProduct());

    double trialyieldnp;
    if(!ysst) trialyieldnp = xitrialnpnorm - sqrt(2./3)*(sigE+theta*Hprime*staten[12]); // Yield Criterion (Simo & Hughes eq. 3.3.6)
    else {
      ysst->getValAndSlopeAlt2(staten[12], &v, &s);
      trialyieldnp = xitrialnpnorm - sqrt(2./3)*v;
      Hprime = Hk + s;
    }

    if (trialyieldnp <= 0) {
 
      //fprintf(stderr, "je suis dans la yield surface\n");

      for (int i=0; i<13; ++i) {
        statenp[i] = staten[i] ;
      }

      Tensor_d0s4_Ss12s34 tm;
      getElasticity(&tm);

      temp = (enp - eplastn);
      stress = tm || temp;
    }
    else {

      //fprintf(stderr, "je suis en dehors de la yield surface\n");

      int i,j,k,l;

      double plastmult = trialyieldnp/(2*mu+(2./3)*Hprime); 
      if(ysst || (yssrt && dt > 0)) {
        // Newton-Raphson algorithm for solution of the return mapping equation
        for (i=0; i<20; ++i) {
          if(ysst) ysst->getValAndSlopeAlt2(staten[12]+sqrt(2./3)*plastmult, &v, &s);
          else {
            v = sigE+theta*Hprime*(staten[12]+sqrt(2./3)*plastmult);
            s = theta*Hprime;
          }
          if(yssrt && dt > 0) {
            yssrt->getValAndSlopeAlt2(sqrt(2./3)*plastmult/dt, &v2, &s2);
            s = v*s2/dt + v2*s;
            v *= v2;
          }
          trialyieldnp = xitrialnpnorm - (2*mu + 2./3*Hk)*plastmult - sqrt(2./3)*v;
          Hprime = Hk + s;
          if(std::abs(trialyieldnp) <= tol*(ysst ? ysst->getVal(0) : sigE)) break;
          plastmult += trialyieldnp/(2*mu + (2./3)*Hprime);
        }
      }
    
      Tensor_d0s2_Ss12 normalnp = (1./xitrialnpnorm)*xitrialnp;

      statenp[12] = staten[12] + sqrt(2./3)*plastmult;

      for (i=0; i<6; ++i) {
        statenp[i] = staten[i] + plastmult*normalnp[i];
        statenp[6+i] = staten[6+i] + sqrt(2./3)*Hk*(statenp[12]-staten[12])*normalnp[i];   
      }

      stress = ((bulk*3*(enp - edevnp) + strialnp) - (2*mu*plastmult*normalnp));
    }
  }
}

template<int e>
void 
ElasPlasKinHardMat<e>::initStates(double *st)
{
  for (int i=0; i<13; ++i)
    st[i] = 0;
}

extern LinearStrain linearStrain;
extern GreenLagrangeStrain greenLagrangeStrain;
extern LogarithmicStrain logarithmicStrain;

template<int e>
StrainEvaluator *
ElasPlasKinHardMat<e>::getStrainEvaluator() const
{
  switch(e) {
    case 0: return &linearStrain; break;
    case 1: return &greenLagrangeStrain; break;
    case 2: return &logarithmicStrain; break;
  }
  return NULL;
} 

template<int e>
bool
ElasPlasKinHardMat<e>::getBackStress(double *statenp, Tensor *_backstress)
{
  Tensor_d0s2_Ss12 * backstress = static_cast<Tensor_d0s2_Ss12 *>(_backstress);
  for (int i=0; i<6; ++i) {
    (*backstress)[i] = statenp[6+i];
  }

  return true;
}

template<int e>
bool
ElasPlasKinHardMat<e>::getPlasticStrain(double *statenp, Tensor *_plasticstrain)
{
  Tensor_d0s2_Ss12 * plasticstrain = static_cast<Tensor_d0s2_Ss12 *>(_plasticstrain);
  for (int i=0; i<6; ++i) {
    (*plasticstrain)[i] = statenp[i];
  }

  return true;
}

template<int e>
double
ElasPlasKinHardMat<e>::getStrainEnergyDensity(Tensor &_enp, double *statenp, double temp)
{
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

template<int e>
double
ElasPlasKinHardMat<e>::getDissipatedEnergy(double *statenp)
{
  if(ysst || yssrt) {
    std::cerr << " *** WARNING: ElasPlasKinHardMat::getDissipatedEnergy is not implemented for nonlinear isotropic hardening.\n";
    return 0;
  }
  else if(theta == 0) { // kinematic hardening
    return statenp[12]*sigE;
  }
  else { // linear isotropic hardening
    double Hprime = (E != Ep) ? (E*Ep)/(E-Ep) : std::numeric_limits<double>::max(); 
    return statenp[12]*(sigE+0.5*theta*Hprime*statenp[12]);
  }
}

template<>
inline void
ElasPlasKinHardMat<0>::print(std::ostream &out) const
{
  out << "BilinearPlastic " << rho << " " << E << " " << nu << " " << Ep << " " << sigE << " " << theta << " " << Tref << " " << alpha;
}

template<>
inline void 
ElasPlasKinHardMat<1>::print(std::ostream &out) const 
{
  out << "FiniteStrainPlastic " << rho << " " << E << " " << nu << " " << Ep << " " << sigE << " " << theta << " " << Tref << " " << alpha;
}

template<>
inline void 
ElasPlasKinHardMat<2>::print(std::ostream &out) const 
{
  out << "LogStrainPlastic " << rho << " " << E << " " << nu << " " << Ep << " " << sigE << " " << theta << " " << Tref << " " << alpha;
}

template<int e>
void
ElasPlasKinHardMat<e>::print2(std::ostream &out) const
{
  if(epsF != std::numeric_limits<double>::infinity() || tol != 1e-6 || yssrtid != 0) {
    out << " " << epsF << " " << tol << " " << yssrtid;
  }
}
