#include <cmath>
#include <limits>
#include <Element.d/NonLinearity.d/NeoHookeanMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

static double sym_matlib_inverse(double *A, double *Ainv)
{
  double det = A[0]*(A[3]*A[5]-A[4]*A[4])
              -A[1]*(A[1]*A[5]-A[4]*A[2])
              +A[2]*(A[1]*A[4]-A[3]*A[2]);

  if (fabs(det) < std::numeric_limits<double>::epsilon()) return 0.e0;

  double detinv = 1./det;
  Ainv[0] = detinv*( A[3]*A[5]-A[4]*A[4]);
  Ainv[1] = detinv*(-A[1]*A[5]+A[2]*A[4]);
  Ainv[2] = detinv*( A[1]*A[4]-A[2]*A[3]);
  Ainv[3] = detinv*( A[0]*A[5]-A[2]*A[2]);
  Ainv[4] = detinv*(-A[0]*A[4]+A[2]*A[1]);
  Ainv[5] = detinv*( A[0]*A[3]-A[1]*A[1]);
  
  return det;
}


NeoHookeanMat::NeoHookeanMat(double _rho, double E, double nu)
{
  rho = _rho;
  lambda = E*nu/((1+nu)*(1-2*nu));
  mu = E/(2*(1+nu));
}

NLMaterial *
NeoHookeanMat::clone() const
{
  return new NeoHookeanMat(*this);
}

void
NeoHookeanMat::getStress(Tensor *_stress, Tensor &_strain, double*, double)
{
  Tensor_d0s2_Ss12 &strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);

  using std::log;
  // compute right Cauchy-Green tensor C = F^t*F = 2*E+I
  double C[6] = { 1+2*strain[0],  2*strain[1],   2*strain[2],
                                1+2*strain[3],   2*strain[4],
                                               1+2*strain[5] };

  // compute PK2 stresses
  double Cinv[6];
  double detC = sym_matlib_inverse(C,Cinv);

  if(detC < std::numeric_limits<double>::epsilon()) {
    std::cerr << "NeoHookeanMat::getStress: close to negative jacobian\n";
    return;
  }

  double defVol = 0.5*log(detC);
  double p = lambda*defVol;
  double coef = p-mu;

  if(_stress) {
    Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress);

    stress[0] = coef*Cinv[0] + mu;
    stress[1] = coef*Cinv[1];
    stress[2] = coef*Cinv[2];
    stress[3] = coef*Cinv[3] + mu;
    stress[4] = coef*Cinv[4];
    stress[5] = coef*Cinv[5] + mu;
  }
}

void 
NeoHookeanMat::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "NeoHookeanMat::getTangentMaterial is not implemented\n"; exit(-1);
}

void 
NeoHookeanMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double)
{
  std::cerr << "NeoHookeanMat::getStressAndTangentMaterial is not implemented\n"; exit(-1);
}

void 
NeoHookeanMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_strain,
                         double *, double *, double, Tensor *, double) const
{
  Tensor_d0s2_Ss12 &strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);

  using std::log;
  // compute right Cauchy-Green tensor C = F^t*F = 2*E+I
  double C[6] = { 1+2*strain[0],  2*strain[1],   2*strain[2],
                                1+2*strain[3],   2*strain[4],
                                               1+2*strain[5] };

  // compute PK2 stresses and derivatives wrt E
  double Cinv[6];
  double detC = sym_matlib_inverse(C,Cinv);

  if(detC < std::numeric_limits<double>::epsilon()) {
    std::cerr << "NeoHookeanMat::integrate: close to negative jacobian\n";
    return;
  }

  double defVol = 0.5*log(detC);
  double p = lambda*defVol;
  double coef = p-mu;

  if(_stress) {
    Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress);

    stress[0] = coef*Cinv[0] + mu;
    stress[1] = coef*Cinv[1];
    stress[2] = coef*Cinv[2];
    stress[3] = coef*Cinv[3] + mu;
    stress[4] = coef*Cinv[4];
    stress[5] = coef*Cinv[5] + mu;
  }

  if(_tm) {
    Tensor_d0s4_Ss12s34 &tm = static_cast<Tensor_d0s4_Ss12s34 &>(*_tm);

    double Cinv2[21] = { Cinv[0]*Cinv[0], Cinv[0]*Cinv[1], Cinv[0]*Cinv[2], Cinv[0]*Cinv[3], Cinv[0]*Cinv[4], Cinv[0]*Cinv[5],
                                          Cinv[1]*Cinv[1], Cinv[1]*Cinv[2], Cinv[1]*Cinv[3], Cinv[1]*Cinv[4], Cinv[1]*Cinv[5],
                                                           Cinv[2]*Cinv[2], Cinv[2]*Cinv[3], Cinv[2]*Cinv[4], Cinv[2]*Cinv[5],
                                                                            Cinv[3]*Cinv[3], Cinv[3]*Cinv[4], Cinv[3]*Cinv[5],
                                                                                             Cinv[4]*Cinv[4], Cinv[4]*Cinv[5],
                                                                                                              Cinv[5]*Cinv[5]};

    double lamin1 = lambda-coef;
    double lamin2 = lambda-2*coef;
    double twoc = 2*coef;
    tm[0][0] = lamin2*Cinv2[0];
    tm[0][1] = tm[1][0] = lamin2*Cinv2[1];
    tm[0][2] = tm[2][0] = lamin2*Cinv2[2];
    tm[0][3] = tm[3][0] = lambda*Cinv2[3] - twoc*Cinv2[6];
    tm[0][4] = tm[4][0] = lambda*Cinv2[4] - twoc*Cinv2[7];
    tm[0][5] = tm[5][0] = lambda*Cinv2[5] - twoc*Cinv2[11];

    tm[1][1] = lamin1*Cinv2[6] - coef*Cinv2[3];
    tm[1][2] = tm[2][1] = lamin1*Cinv2[7] - coef*Cinv2[4];
    tm[1][3] = tm[3][1] = lamin2*Cinv2[8];
    tm[1][4] = tm[4][1] = lamin1*Cinv2[9] - coef*Cinv2[12];
    tm[1][5] = tm[5][1] = lambda*Cinv2[10] - twoc*Cinv2[13];

    tm[2][2] = lamin1*Cinv2[11] - coef*Cinv2[5 ];
    tm[2][3] = tm[3][2] = lambda*Cinv2[12] - twoc*Cinv2[9];
    tm[2][4] = tm[4][2] = lamin1*Cinv2[13] - coef*Cinv2[10];
    tm[2][5] = tm[5][2] = lamin2*Cinv2[14];

    tm[3][3] = lamin2*Cinv2[15];
    tm[3][4] = tm[4][3] = lamin2*Cinv2[16];
    tm[3][5] = tm[5][3] = lambda*Cinv2[17] - twoc*Cinv2[18];

    tm[4][4] = lamin1*Cinv2[18] - coef*Cinv2[17];
    tm[4][5] = tm[5][4] = lamin2*Cinv2[19];

    tm[5][5] = lamin2*Cinv2[20];
  }
}

void
NeoHookeanMat::integrate(Tensor *_stress, Tensor &, Tensor &_strain,
                         double *, double *, double, Tensor *, double) const
{
  Tensor_d0s2_Ss12 &strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);

  using std::log;
  // compute right Cauchy-Green tensor C = F^t*F = 2*E+I
  double C[6] = { 1+2*strain[0],  2*strain[1],   2*strain[2],
                                1+2*strain[3],   2*strain[4],
                                               1+2*strain[5] };

  // compute PK2 stresses
  double Cinv[6];
  double detC = sym_matlib_inverse(C,&Cinv[0]);

  if(detC < std::numeric_limits<double>::epsilon()) {
    std::cerr << "NeoHookeanMat::integrate: close to negative jacobian\n";
    return;
  }

  double defVol = 0.5*log(detC);
  double p = lambda*defVol;
  double coef = p-mu;

  if(_stress) {
    Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(*_stress);

    stress[0] = coef*Cinv[0] + mu;
    stress[1] = coef*Cinv[1];
    stress[2] = coef*Cinv[2];
    stress[3] = coef*Cinv[3] + mu;
    stress[4] = coef*Cinv[4];
    stress[5] = coef*Cinv[5] + mu;
  }
}

double
NeoHookeanMat::getStrainEnergyDensity(Tensor &_strain, double *, double)
{
  Tensor_d0s2_Ss12 &strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);

  using std::log;
  // compute right Cauchy-Green tensor C = F^t*F = 2*E+I
  double C[6] = { 1+2*strain[0],  2*strain[1],   2*strain[2],
                                1+2*strain[3],   2*strain[4],
                                               1+2*strain[5] };

  double detC = C[0]*(C[3]*C[5]-C[4]*C[4])
               -C[1]*(C[1]*C[5]-C[4]*C[2])
               +C[2]*(C[1]*C[4]-C[3]*C[2]);

  if(detC < std::numeric_limits<double>::epsilon()) {
    std::cerr << "NeoHookeanMat::getStrainEnergyDensity: close to negative jacobian\n";
    return 0;
  }

  double defVol = 0.5*log(detC);
  return mu*(strain[0]+strain[3]+strain[5]-defVol) + lambda/2*defVol*defVol;
}

void
NeoHookeanMat::getMaterialConstants(std::vector<double> &c)
{
  c.resize(2);
  c[0] = lambda;
  c[1] = mu;
}

extern GreenLagrangeStrain greenLagrangeStrain;

StrainEvaluator *
NeoHookeanMat::getStrainEvaluator() const
{
  return &greenLagrangeStrain;
}

void
NeoHookeanMat::print(std::ostream &out) const
{
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  out << "NeoHookean " << rho << " " << E << " " << nu;
}
