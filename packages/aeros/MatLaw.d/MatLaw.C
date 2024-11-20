#include "MatLaw.h"

MatLaw::MatLaw(double _rho, double _E, double _nu, double _Tref, double _alpha)
{
  rho = _rho;
  nu = _nu;
  E = _E;
  Tref = _Tref;
  alpha = _alpha;
}

void
MatLaw::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  Tensor_d0s4_Ss12s34 tm;
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  getTangentMaterial(&tm, _strain, state, temp);

  double e0 = (temp-Tref)*alpha;
  strain[0] -= e0;
  strain[3] -= e0;
  strain[5] -= e0;

  (*stress) = tm||strain;
}

void 
MatLaw::getTangentMaterial(Tensor *_tm, Tensor &, double*, double temp)
{
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = (nu != 0) ? lambda/nu : E;

  (*tm)[0][0] = lambdadivnu*(1-nu);
  (*tm)[1][1] = lambdadivnu*(1-2*nu)/2;
  (*tm)[2][2] = lambdadivnu*(1-2*nu)/2;
  (*tm)[3][3] = lambdadivnu*(1-nu);
  (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
  (*tm)[5][5] = lambdadivnu*(1-nu);
  (*tm)[0][3] = lambdadivnu*nu;
  (*tm)[3][0] = lambdadivnu*nu;
  (*tm)[0][5] = lambdadivnu*nu;
  (*tm)[5][0] = lambdadivnu*nu;
  (*tm)[3][5] = lambdadivnu*nu;
  (*tm)[5][3] = lambdadivnu*nu;
}

void 
MatLaw::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
{
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = (nu != 0) ? lambda/nu : E;

  (*tm)[0][0] = lambdadivnu*(1-nu);
  (*tm)[1][1] = lambdadivnu*(1-2*nu)/2;
  (*tm)[2][2] = lambdadivnu*(1-2*nu)/2;
  (*tm)[3][3] = lambdadivnu*(1-nu);
  (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
  (*tm)[5][5] = lambdadivnu*(1-nu);
  (*tm)[0][3] = lambdadivnu*nu;
  (*tm)[3][0] = lambdadivnu*nu;
  (*tm)[0][5] = lambdadivnu*nu;
  (*tm)[5][0] = lambdadivnu*nu;
  (*tm)[3][5] = lambdadivnu*nu;
  (*tm)[5][3] = lambdadivnu*nu;

  double e0 = (temp-Tref)*alpha;
  strain[0] -= e0;
  strain[3] -= e0;
  strain[5] -= e0;

  (*stress) =(*tm)||strain;
}

void 
MatLaw::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                  double *, double *, double temp)
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = (nu != 0) ? lambda/nu : E;

  (*tm)[0][0] = lambdadivnu*(1-nu);
  (*tm)[1][1] = lambdadivnu*(1-2*nu)/2;
  (*tm)[2][2] = lambdadivnu*(1-2*nu)/2;
  (*tm)[3][3] = lambdadivnu*(1-nu);
  (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
  (*tm)[5][5] = lambdadivnu*(1-nu);
  (*tm)[0][3] = lambdadivnu*nu;
  (*tm)[3][0] = lambdadivnu*nu;
  (*tm)[0][5] = lambdadivnu*nu;
  (*tm)[5][0] = lambdadivnu*nu;
  (*tm)[3][5] = lambdadivnu*nu;
  (*tm)[5][3] = lambdadivnu*nu;

  double e0 = (temp-Tref)*alpha;
  enp[0] -= e0;
  enp[3] -= e0;
  enp[5] -= e0;

  (*stress) = (*tm)||enp;
}

void
MatLaw::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                  double *, double *, double temp)
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);

  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double lambdadivnu = (nu != 0) ? lambda/nu : E;

  Tensor_d0s4_Ss12s34 *tm = new Tensor_d0s4_Ss12s34();

  (*tm)[0][0] = lambdadivnu*(1-nu);
  (*tm)[1][1] = lambdadivnu*(1-2*nu)/2;
  (*tm)[2][2] = lambdadivnu*(1-2*nu)/2;
  (*tm)[3][3] = lambdadivnu*(1-nu);
  (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
  (*tm)[5][5] = lambdadivnu*(1-nu);
  (*tm)[0][3] = lambdadivnu*nu;
  (*tm)[3][0] = lambdadivnu*nu;
  (*tm)[0][5] = lambdadivnu*nu;
  (*tm)[5][0] = lambdadivnu*nu;
  (*tm)[3][5] = lambdadivnu*nu;
  (*tm)[5][3] = lambdadivnu*nu;

  double e0 = (temp-Tref)*alpha;
  enp[0] -= e0;
  enp[3] -= e0;
  enp[5] -= e0;

  (*stress) = (*tm)||enp;

  delete tm;
}

extern LinearStrain linearStrain;

StrainEvaluator *
MatLaw::getStrainEvaluator() const
{
  return &linearStrain;
}

