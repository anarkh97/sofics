#include <Element.d/NonLinearity.d/SimoPlasticMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <cstddef>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#endif

SimoPlasticMat::SimoPlasticMat(double _rho, double _E, double _nu,
                               double _Ep, double _sigE, double _theta,
                               double _tol, int _yssrtid)
{
  rho = _rho;
  E = _E;
  nu = _nu;
  Ep = _Ep;
  sigE = _sigE;
  theta = _theta;
  tol = _tol;
  yssrtid = _yssrtid;
  ysst = NULL;
  yssrt = NULL;
}

int
SimoPlasticMat::getNumStates() const
{
  // the internal variables are : the plastic Green-Lagrange strain tensor Cp-I (6 doubles),
  // the center of the yield surface in principal deviatoric stress space (3 doubles),
  // the equivalent plastic strain (1 double),
  return 10;
}

void
SimoPlasticMat::getStress(Tensor *stress, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: SimoPlasticMat::getStress is not implemented\n";
}

void 
SimoPlasticMat::getElasticity(Tensor *_tm) const
{
  std::cerr << "WARNING: SimoPlasticMat::getElasticity is not implemented\n";
}

void 
SimoPlasticMat::getTangentMaterial(Tensor *tm, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: SimoPlasticMat::getTangentMaterial is not implemented\n";
}

void 
SimoPlasticMat::getStressAndTangentMaterial(Tensor *stess, Tensor *tm, Tensor &strain, double *state, double temp)
{
  std::cerr << "WARNING: SimoPlasticMat::getStressAndTangentMaterial is not implemented\n";
}

void 
SimoPlasticMat::updateStates(Tensor &en, Tensor &enp, double *state, double temp)
{
  std::cerr << "WARNING: SimoPlasticMat::updateStates is not implemented\n";
}

void
SimoPlasticMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &_en, Tensor  &_enp,
                          double *staten, double *statenp, double temp,
                          Tensor *_cache, double dt) const
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Array3d> epsetrialnp(&static_cast<Tensor_d0s2_Ss12_diag &>(_enp)[0]);
  Eigen::Map<Eigen::Array3d> betanp(&static_cast<Tensor_d0s2_Ss12_diag &>(*_stress)[0]);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > a(&static_cast<Tensor_d0s4_Ss12s34_diag &>(*_tm)[0][0]);
  Tensor_d1s2_full & cache = static_cast<Tensor_d1s2_full &>(*_cache);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > Fnp = cache[0].matrix(); // deformation gradient
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > n   = cache[1].matrix(); // eigenvectors of benptrial

  // back stress
  Eigen::Map<Eigen::Array3d> xibarn(&staten[6]), xibarnp(&statenp[6]);
  // consistency parameter
  double &gamman = staten[9], &gammanp = statenp[9];
  // plastic right Cauchy-Green strain tensor at t_n
  Eigen::Matrix3d Cpn;
  Cpn << staten[0]+1, staten[1],   staten[2],
         staten[1],   staten[3]+1, staten[4],
         staten[2],   staten[4],   staten[5]+1;

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double mu = E/(2*(1+nu)); // shear modulus
  double bulk = lambda + (2./3)*mu;
  double Hprime = (E != Ep) ? (E*Ep)/(E-Ep) : std::numeric_limits<double>::max();
  double Hk = (1-theta)*Hprime; // kinematic hardening modulus
  double v, s; // value and slope of yield stress vs. effective plastic strain curve
  double v2, s2; // value and slope of yield stress scaling factor vs. effective plastic strain rate curve

  // compute trial deviatoric elastic logarithmic principal stretches and deviatoric principal Kirchhoff stresses
  Eigen::Array3d epsebartrialnp = epsetrialnp - 1/3.*epsetrialnp.sum();
  Eigen::Array3d betabartrialnp = 2*mu*epsebartrialnp;

  // evaluate the yield function
  Eigen::Array3d zetabartrialnp = betabartrialnp - xibarn;
  double zetabartrialnpnorm = sqrt((zetabartrialnp*zetabartrialnp).sum());

  double trialyieldnp;
  if(!ysst) trialyieldnp = zetabartrialnpnorm - sqrt(2./3)*(sigE+theta*Hprime*staten[9]); // Yield Criterion (Simo & Hughes eq. 3.3.6)
  else {
    ysst->getValAndSlopeAlt2(staten[9], &v, &s);
    trialyieldnp = zetabartrialnpnorm - sqrt(2./3)*v;
    Hprime = Hk + s;
  }

  if (trialyieldnp <= 0) {
 
    //fprintf(stderr, "je suis dans la yield surface\n");

    for (int i=0; i<10; ++i) statenp[i] = staten[i];

    // compute principal Kirchhoff stress
    betanp = bulk*3*(epsetrialnp - epsebartrialnp) + betabartrialnp;

    // compute tangent modulus a = ∂β/∂ε
    a = bulk*Eigen::Matrix3d::Ones() + 2*mu*(Eigen::Matrix3d::Identity()-1/3.*Eigen::Matrix3d::Ones());
  }
  else {

    //fprintf(stderr, "je suis en dehors de la yield surface\n");

    double plastmult = trialyieldnp/(2*mu+(2./3)*Hprime); 
    if(ysst || (yssrt && dt > 0)) {
      // Newton-Raphson algorithm for solution of the return mapping equation
      for (int i=0; i<20; ++i) {
        if(ysst) ysst->getValAndSlopeAlt2(staten[9]+sqrt(2./3)*plastmult, &v, &s);
        else {
          v = sigE+theta*Hprime*(staten[9]+sqrt(2./3)*plastmult);
          s = theta*Hprime;
        }
        if(yssrt && dt > 0) {
          yssrt->getValAndSlopeAlt2(sqrt(2./3)*plastmult/dt, &v2, &s2);
          s = v*s2/dt + v2*s;
          v *= v2;
        }
        trialyieldnp = zetabartrialnpnorm - (2*mu + 2./3*Hk)*plastmult - sqrt(2./3)*v;
        Hprime = Hk + s;
        if(std::abs(trialyieldnp) <= tol*(ysst ? ysst->getVal(0) : sigE)) break;
        plastmult += trialyieldnp/(2*mu + (2./3)*Hprime);
      }
    }
   
    Eigen::Array3d normalnp = (1./zetabartrialnpnorm)*zetabartrialnp;

    // update equivalent plastic strain
    gammanp = gamman + sqrt(2./3)*plastmult;

    // update elastic logarithmic principal stretches
    Eigen::Array3d epsenp = epsetrialnp - plastmult*normalnp;

    // update backstress
    xibarnp = xibarn + sqrt(2./3)*Hk*(gammanp-gamman)*normalnp;

    // update plastic right Cauchy-Green strain tensor Cpnp = Fnp^T*benp^{-1}*Fnp where benp = n*exp(2*epsenp)*n^T
    Eigen::Matrix3d Cpnp = Fnp.transpose()*n*((1/(2*epsenp).exp()).matrix().asDiagonal())*n.transpose()*Fnp;
    statenp[0] = Cpnp(0,0)-1;
    statenp[1] = Cpnp(0,1);
    statenp[2] = Cpnp(0,2);
    statenp[3] = Cpnp(1,1)-1;
    statenp[4] = Cpnp(1,2);
    statenp[5] = Cpnp(2,2)-1;

    // compute principal Kirchhoff stress
    betanp = bulk*3*(epsetrialnp - epsebartrialnp) + betabartrialnp - (2*mu*plastmult*normalnp);

    // compute the tangent modulus a = ∂β/∂ε
    double thetanp = 1 - 2*mu*plastmult/zetabartrialnpnorm;
    double thetaprimenp = 1/(1+Hprime/(3*mu)) - (1 - thetanp);
    a = bulk*Eigen::Matrix3d::Ones()
      + 2*mu*thetanp*(Eigen::Matrix3d::Identity()-1/3.*Eigen::Matrix3d::Ones())
      - 2*mu*thetaprimenp*normalnp.matrix()*normalnp.matrix().transpose();
  }
#endif
}

void
SimoPlasticMat::integrate(Tensor *_stress, Tensor &, Tensor  &_enp,
                          double *staten, double *statenp, double temp,
                          Tensor *_cache, double dt) const
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Array3d> epsetrialnp(&static_cast<Tensor_d0s2_Ss12_diag &>(_enp)[0]);
  Eigen::Map<Eigen::Array3d> betanp(&static_cast<Tensor_d0s2_Ss12_diag &>(*_stress)[0]);
  Tensor_d1s2_full & cache = static_cast<Tensor_d1s2_full &>(*_cache);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > Fnp = cache[0].matrix(); // deformation gradient
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > n   = cache[1].matrix(); // eigenvectors of benptrial

  // back stress
  Eigen::Map<Eigen::Array3d> xibarn(&staten[6]), xibarnp(&statenp[6]);
  // consistency parameter
  double &gamman = staten[9], &gammanp = statenp[9];
  // plastic right Cauchy-Green strain tensor at t_n
  Eigen::Matrix3d Cpn;
  Cpn << staten[0]+1, staten[1],   staten[2],
         staten[1],   staten[3]+1, staten[4],
         staten[2],   staten[4],   staten[5]+1;

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double mu = E/(2*(1+nu)); // shear modulus
  double bulk = lambda + (2./3)*mu;
  double Hprime = (E != Ep) ? (E*Ep)/(E-Ep) : std::numeric_limits<double>::max();
  double Hk = (1-theta)*Hprime; // kinematic hardening modulus
  double v, s; // value and slope of yield stress vs. effective plastic strain curve
  double v2, s2; // value and slope of yield stress scaling factor vs. effective plastic strain rate curve

  // compute trial deviatoric elastic logarithmic principal stretches and deviatoric principal Kirchhoff stresses
  Eigen::Array3d epsebartrialnp = epsetrialnp - 1/3.*epsetrialnp.sum();
  Eigen::Array3d betabartrialnp = 2*mu*epsebartrialnp;

  // evaluate the yield function
  Eigen::Array3d zetabartrialnp = betabartrialnp - xibarn;
  double zetabartrialnpnorm = sqrt((zetabartrialnp*zetabartrialnp).sum());

  double trialyieldnp;
  if(!ysst) trialyieldnp = zetabartrialnpnorm - sqrt(2./3)*(sigE+theta*Hprime*staten[9]); // Yield Criterion (Simo & Hughes eq. 3.3.6)
  else {
    ysst->getValAndSlopeAlt2(staten[9], &v, &s);
    trialyieldnp = zetabartrialnpnorm - sqrt(2./3)*v;
    Hprime = Hk + s;
  }

  if (trialyieldnp <= 0) {
 
    //fprintf(stderr, "je suis dans la yield surface\n");

    for (int i=0; i<10; ++i) statenp[i] = staten[i];

    // compute principal Kirchhoff stress
    betanp = bulk*3*(epsetrialnp - epsebartrialnp) + betabartrialnp;
  }
  else {

    //fprintf(stderr, "je suis en dehors de la yield surface\n");

    double plastmult = trialyieldnp/(2*mu+(2./3)*Hprime); 
    if(ysst || (yssrt && dt > 0)) {
      // Newton-Raphson algorithm for solution of the return mapping equation
      for (int i=0; i<20; ++i) {
        if(ysst) ysst->getValAndSlopeAlt2(staten[9]+sqrt(2./3)*plastmult, &v, &s);
        else {
          v = sigE+theta*Hprime*(staten[9]+sqrt(2./3)*plastmult);
          s = theta*Hprime;
        }
        if(yssrt && dt > 0) {
          yssrt->getValAndSlopeAlt2(sqrt(2./3)*plastmult/dt, &v2, &s2);
          s = v*s2/dt + v2*s;
          v *= v2;
        }
        trialyieldnp = zetabartrialnpnorm - (2*mu + 2./3*Hk)*plastmult - sqrt(2./3)*v;
        Hprime = Hk + s;
        if(std::abs(trialyieldnp) <= tol*(ysst ? ysst->getVal(0) : sigE)) break;
        plastmult += trialyieldnp/(2*mu + (2./3)*Hprime);
      }
    }
   
    Eigen::Array3d normalnp = (1./zetabartrialnpnorm)*zetabartrialnp;

    // update equivalent plastic strain
    gammanp = gamman + sqrt(2./3)*plastmult;

    // update elastic logarithmic principal stretches
    Eigen::Array3d epsenp = epsetrialnp - plastmult*normalnp;

    // update backstress
    xibarnp = xibarn + sqrt(2./3)*Hk*(gammanp-gamman)*normalnp;

    // update plastic right Cauchy-Green strain tensor Cpnp = Fnp^T*benp^{-1}*Fnp where benp = n*exp(2*epsenp)*n^T
    Eigen::Matrix3d Cpnp = Fnp.transpose()*n*((1/(2*epsenp).exp()).matrix().asDiagonal())*n.transpose()*Fnp;
    statenp[0] = Cpnp(0,0)-1;
    statenp[1] = Cpnp(0,1);
    statenp[2] = Cpnp(0,2);
    statenp[3] = Cpnp(1,1)-1;
    statenp[4] = Cpnp(1,2);
    statenp[5] = Cpnp(2,2)-1;

    // compute principal Kirchhoff stress
    betanp = bulk*3*(epsetrialnp - epsebartrialnp) + betabartrialnp - (2*mu*plastmult*normalnp);
  }
#endif
}

void 
SimoPlasticMat::initStates(double *st)
{
  for(int i=0; i<getNumStates(); ++i) st[i] = 0;
}

extern ElasticLogarithmicPrincipalStretches elasticLogarithmicPrincipalStretches;

StrainEvaluator *
SimoPlasticMat::getStrainEvaluator() const
{
  return &elasticLogarithmicPrincipalStretches;
} 

bool
SimoPlasticMat::getBackStress(double *statenp, Tensor *_backstress)
{
  Tensor_d0s2_Ss12 * backstress = static_cast<Tensor_d0s2_Ss12 *>(_backstress);
  (*backstress)[0] = statenp[6];
  (*backstress)[1] = 0;
  (*backstress)[2] = 0;
  (*backstress)[3] = statenp[7];
  (*backstress)[4] = 0;
  (*backstress)[3] = statenp[8];

  return true;
}

bool
SimoPlasticMat::getPlasticStrain(double *statenp, Tensor *_plasticstrain)
{
  Tensor_d0s2_Ss12 * plasticstrain = static_cast<Tensor_d0s2_Ss12 *>(_plasticstrain);
  for (int i=0; i<6; ++i) {
    (*plasticstrain)[i] = statenp[i];
  }

  return true;
}

double
SimoPlasticMat::getEquivPlasticStrain(double *statenp)
{
  return statenp[9];
}

double
SimoPlasticMat::getStrainEnergyDensity(Tensor &_enp, double *statenp, double temp)
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_enp);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double mu = E/(2*(1+nu));                  // shear modulus
  double trace = eps[0]+eps[1]+eps[2];

  return 0.5*lambda*trace*trace + mu*(eps[0]*eps[0]+eps[1]*eps[1]+eps[2]*eps[2]);
}

double
SimoPlasticMat::getDissipatedEnergy(double *statenp)
{
  if(ysst || yssrt) {
    std::cerr << " *** WARNING: ElasPlasKinHardMat::getDissipatedEnergy is not implemented for nonlinear isotropic hardening.\n";
    return 0;
  }
  else if(theta == 0) { // kinematic hardening
    return statenp[9]*sigE;
  }
  else { // linear isotropic hardening
    double Hprime = (E != Ep) ? (E*Ep)/(E-Ep) : std::numeric_limits<double>::max();
    return statenp[9]*(sigE+0.5*theta*Hprime*statenp[9]);
  }
}

void
SimoPlasticMat::print(std::ostream &out) const
{
  out << "SimoPlastic " << rho << " " << E << " " << nu << " " << Ep << " " << sigE << " " << theta;
}

