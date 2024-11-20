#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/SimoElasticMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

SimoElasticMat::SimoElasticMat(double _rho, double _E, double _nu)
{
  rho = _rho;
  E = _E;
  nu = _nu;
}

NLMaterial *
SimoElasticMat::clone() const
{
  return new SimoElasticMat(*this);
}

void
SimoElasticMat::getStress(Tensor *_stress, Tensor &_strain, double*, double)
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double M = (nu==0) ? E : lambda*(1-nu)/nu; // P-wave modulus
  
  beta[0] = M*eps[0] + lambda*(eps[1]+eps[2]);
  beta[1] = M*eps[1] + lambda*(eps[0]+eps[2]);
  beta[2] = M*eps[2] + lambda*(eps[0]+eps[1]);
}

void 
SimoElasticMat::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "SimoElasticMat::getTangentMaterial is not implemented\n"; exit(-1);
}

void 
SimoElasticMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "SimoElasticMat::getStressAndTangentMaterial is not implemented\n"; exit(-1);
}

void 
SimoElasticMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_strain,
                          double *, double *, double, Tensor *, double) const
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);
  Tensor_d0s4_Ss12s34_diag *tm = static_cast<Tensor_d0s4_Ss12s34_diag *>(_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double M = (nu==0) ? E : lambda*(1-nu)/nu; // P-wave modulus

  (*tm)[0][0] = (*tm)[1][1] = (*tm)[2][2] = M;
  (*tm)[0][1] = (*tm)[1][0] = (*tm)[0][2] = (*tm)[2][0] = (*tm)[1][2] = (*tm)[2][1] = lambda;

  beta[0] = M*eps[0] + lambda*(eps[1]+eps[2]);
  beta[1] = M*eps[1] + lambda*(eps[0]+eps[2]);
  beta[2] = M*eps[2] + lambda*(eps[0]+eps[1]);
}

void
SimoElasticMat::integrate(Tensor *_stress, Tensor &, Tensor &_strain,
                          double *, double *, double, Tensor *, double) const
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double M = (nu==0) ? E : lambda*(1-nu)/nu; // P-wave modulus

  beta[0] = M*eps[0] + lambda*(eps[1]+eps[2]);
  beta[1] = M*eps[1] + lambda*(eps[0]+eps[2]);
  beta[2] = M*eps[2] + lambda*(eps[0]+eps[1]);
}

double
SimoElasticMat::getStrainEnergyDensity(Tensor &_strain, double *, double)
{
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu)); // Lame's 1st parameter
  double mu = E/(2*(1+nu));                  // shear modulus
  double trace = eps[0]+eps[1]+eps[2];

  return 0.5*lambda*trace*trace + mu*(eps[0]*eps[0]+eps[1]*eps[1]+eps[2]*eps[2]);
}

void
SimoElasticMat::print(std::ostream &out) const
{
  out << "SimoElastic " << rho << " " << E << " " << nu;
}

void
SimoElasticMat::getMaterialConstants(std::vector<double> &c)
{
  c.resize(2);
  c[0] = E;
  c[1] = nu;
}

extern LogarithmicPrincipalStretches logarithmicPrincipalStretches;

StrainEvaluator *
SimoElasticMat::getStrainEvaluator() const
{
  return &logarithmicPrincipalStretches;
}

