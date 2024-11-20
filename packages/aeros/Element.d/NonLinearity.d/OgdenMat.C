#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

template<int _m, int _n>
OgdenMat::OgdenMat(double _rho, double (&_mu)[_m], double (&_alpha)[_m], double (&_K)[_n])
{
  rho = _rho;
  m = std::min(_m,9);
  for(int i=0; i<m; ++i) { mu[i] = _mu[i]; alpha[i] = _alpha[i]; }
  n = std::min(_n,9);
  for(int i=0; i<n; ++i) K[i] = _K[i];
}

NLMaterial *
OgdenMat::clone() const
{
  return new OgdenMat(*this);
}

void
OgdenMat::getStress(Tensor *_stress, Tensor &_strain, double*, double)
{
  using std::pow;
  using std::exp;
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);

  double trace = eps[0]+eps[1]+eps[2];
  double J = exp(trace);
  double dUdJ = 0;
  for(int i=0; i<n; ++i) dUdJ += K[i]*(i+1)*pow(J-1,2*i+1);
  beta[0] = beta[1] = beta[2] = J*dUdJ;

  double lambdabar[3] = { exp(eps[0]-trace/3), exp(eps[1]-trace/3), exp(eps[2]-trace/3) };
  for(int i=0; i<m; ++i) {
    double w[3] = { mu[i]*pow(lambdabar[0],alpha[i]), mu[i]*pow(lambdabar[1],alpha[i]), mu[i]*pow(lambdabar[2],alpha[i]) };
    double toto = (w[0]+w[1]+w[2])/3;
    beta[0] += w[0]-toto;
    beta[1] += w[1]-toto;
    beta[2] += w[2]-toto;
  }
}

void 
OgdenMat::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "OgdenMat::getTangentMaterial is not implemented\n"; exit(-1);
}

void 
OgdenMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "OgdenMat::getStressAndTangentMaterial is not implemented\n"; exit(-1);
}

void 
OgdenMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_strain,
                    double *, double *, double, Tensor *, double) const
{
  using std::pow;
  using std::exp;
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);
  Tensor_d0s4_Ss12s34_diag *tm = static_cast<Tensor_d0s4_Ss12s34_diag *>(_tm);

  double trace = eps[0]+eps[1]+eps[2];
  double J = exp(trace);
  double dUdJ = 0, d2UdJ2 = 0; 
  for(int i=0; i<n; ++i) {
    dUdJ += K[i]*(i+1)*pow(J-1,2*i+1);
    d2UdJ2 += K[i]*(i+1)*(2*i+1)*pow(J-1,2*i);
  }
  beta[0] = beta[1] = beta[2] = J*dUdJ;
  (*tm)[0][0] = (*tm)[1][1] = (*tm)[2][2] = (*tm)[0][1] = (*tm)[0][2] = (*tm)[1][2] = d2UdJ2*J*J + J*dUdJ;

  double lambdabar[3] = { exp(eps[0]-trace/3), exp(eps[1]-trace/3), exp(eps[2]-trace/3) };
  for(int i=0; i<m; ++i) {
    double w[3] = { mu[i]*pow(lambdabar[0],alpha[i]), mu[i]*pow(lambdabar[1],alpha[i]), mu[i]*pow(lambdabar[2],alpha[i]) };
    double toto = (w[0]+w[1]+w[2])/3;
    beta[0] += w[0]-toto;
    beta[1] += w[1]-toto;
    beta[2] += w[2]-toto;
    (*tm)[0][0] += alpha[i]*(w[0]+toto)/3;
    (*tm)[1][1] += alpha[i]*(w[1]+toto)/3;
    (*tm)[2][2] += alpha[i]*(w[2]+toto)/3;
    (*tm)[0][1] += alpha[i]*(-w[0]-w[1]+toto)/3;
    (*tm)[0][2] += alpha[i]*(-w[0]-w[2]+toto)/3;
    (*tm)[1][2] += alpha[i]*(-w[1]-w[2]+toto)/3;
  }

  (*tm)[1][0] = (*tm)[0][1];
  (*tm)[2][0] = (*tm)[0][2];
  (*tm)[2][1] = (*tm)[1][2];
}

void
OgdenMat::integrate(Tensor *_stress, Tensor &, Tensor &_strain,
                    double *, double *, double, Tensor *, double) const
{
  using std::pow;
  using std::exp;
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);
  Tensor_d0s2_Ss12_diag &beta = static_cast<Tensor_d0s2_Ss12_diag &>(*_stress);

  double trace = eps[0]+eps[1]+eps[2];
  double J = exp(trace);
  double dUdJ = 0;
  for(int i=0; i<n; ++i) dUdJ += K[i]*(i+1)*pow(J-1,2*i+1);
  beta[0] = beta[1] = beta[2] = J*dUdJ;

  double lambdabar[3] = { exp(eps[0]-trace/3), exp(eps[1]-trace/3), exp(eps[2]-trace/3) };
  for(int i=0; i<m; ++i) {
    double w[3] = { mu[i]*pow(lambdabar[0],alpha[i]), mu[i]*pow(lambdabar[1],alpha[i]), mu[i]*pow(lambdabar[2],alpha[i]) };
    double toto = (w[0]+w[1]+w[2])/3;
    beta[0] += w[0]-toto;
    beta[1] += w[1]-toto;
    beta[2] += w[2]-toto;
  }
}

double
OgdenMat::getStrainEnergyDensity(Tensor &_strain, double *, double)
{
  using std::exp;
  using std::pow;
  Tensor_d0s2_Ss12_diag &eps = static_cast<Tensor_d0s2_Ss12_diag &>(_strain);

  double trace = eps[0]+eps[1]+eps[2];
  double J = exp(trace);
  double U = 0;
  for(int i=0; i<n; ++i) U += K[i]/2*pow(J-1,2*i+2);

  double W = 0;
  double lambdabar[3] = { exp(eps[0]-trace/3), exp(eps[1]-trace/3), exp(eps[2]-trace/3) };
  for(int i=0; i<m; ++i) W += mu[i]/alpha[i]*(pow(lambdabar[0],alpha[i])+pow(lambdabar[1],alpha[i])+pow(lambdabar[2],alpha[i])-3);

  return W+U;
}

void
OgdenMat::print(std::ostream &out) const
{
  out << "Ogden " << rho;
  for(int i=0; i<m; ++i) out << " " << mu[i];
  for(int i=0; i<m; ++i) out << " " << alpha[i];
  for(int i=0; i<n; ++i) out << " " << K[i];
}

void
OgdenMat::getMaterialConstants(std::vector<double> &c)
{
  c.resize(2*m);
  for(int i=0; i<m; ++i) { c[i] = mu[i]; c[m+i] = alpha[i]; }
}

extern LogarithmicPrincipalStretches logarithmicPrincipalStretches;

StrainEvaluator *
OgdenMat::getStrainEvaluator() const
{
  return &logarithmicPrincipalStretches;
}

template OgdenMat::OgdenMat(double, double (&)[1], double (&)[1], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[2], double (&)[2], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[3], double (&)[3], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[4], double (&)[4], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[5], double (&)[5], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[6], double (&)[6], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[7], double (&)[7], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[8], double (&)[8], double (&)[1]);
template OgdenMat::OgdenMat(double, double (&)[9], double (&)[9], double (&)[1]);

template OgdenMat::OgdenMat(double, double (&)[2], double (&)[2], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[3], double (&)[3], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[4], double (&)[4], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[5], double (&)[5], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[6], double (&)[6], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[7], double (&)[7], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[8], double (&)[8], double (&)[2]);
template OgdenMat::OgdenMat(double, double (&)[9], double (&)[9], double (&)[2]);
