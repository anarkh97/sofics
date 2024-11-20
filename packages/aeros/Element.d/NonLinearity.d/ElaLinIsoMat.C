#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Utils.d/NodeSpaceArray.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

ElaLinIsoMat::ElaLinIsoMat(StructProp *p)
{
  E = p->E;
  nu = p->nu;
  rho = p->rho;
  Tref = p->Ta;
  alphas[0] = alphas[3] = alphas[5] = alpha = p->W;
  alphas[1] = alphas[2] = alphas[4] = 0;
  m_tm = 0;
  ymtt = p->ymtt;
  ctett = p->ctett;
}

ElaLinIsoMat::ElaLinIsoMat(double _rho)
{
  rho = _rho;
  nu = 0;
  E = 0;
  Tref = 0;
  alphas[0] = alphas[3] = alphas[5] = alpha = 0;
  alphas[1] = alphas[2] = alphas[4] = 0;
  m_tm = 0;
  ymtt = NULL;
  ctett = NULL;
}

ElaLinIsoMat::ElaLinIsoMat(double _rho, double _E, double _nu)
{
  rho = _rho;
  nu = _nu;
  E = _E;
  Tref = 0;
  alphas[0] = alphas[3] = alphas[5] = alpha = 0;
  alphas[1] = alphas[2] = alphas[4] = 0;
  m_tm = 0;
  ymtt = NULL;
  ctett = NULL;
}

ElaLinIsoMat::ElaLinIsoMat(double _rho, double _E, double _nu, double _Tref, double _alpha)
{
  rho = _rho;
  nu = _nu;
  E = _E;
  Tref = _Tref;
  alphas[0] = alphas[3] = alphas[5] = alpha = _alpha;
  alphas[1] = alphas[2] = alphas[4] = 0;
  m_tm = 0;
  ymtt = NULL;
  ctett = NULL;
}

ElaLinIsoMat::ElaLinIsoMat(double _rho, double C[6][6], double _Tref, double _alphas[6])
{
  rho = _rho;
  Tref = _Tref;
  setTangentMaterial(C);
  setThermalExpansionCoef(_alphas);
  ymtt = NULL;
  ctett = NULL;
}

ElaLinIsoMat::~ElaLinIsoMat()
{
  if(m_tm) delete m_tm;
}

NLMaterial *
ElaLinIsoMat::clone() const
{
  return new ElaLinIsoMat(*this);
}

void
ElaLinIsoMat::setTangentMaterial(double C[6][6])
{
  m_tm = new Tensor_d0s4_Ss12s34();
  int index_map[6] = { 0,3,5,1,4,2 };
  for(int i=0; i<6; ++i)
    for(int j=0; j<6; ++j)
      (*m_tm)[i][j] = C[index_map[i]][index_map[j]];
}

void
ElaLinIsoMat::setThermalExpansionCoef(double _alphas[6])
{
  int index_map[6] = { 0,3,5,1,4,2 };
  for(int i=0; i<6; ++i) alphas[i] = _alphas[index_map[i]];
}

void
ElaLinIsoMat::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  Tensor_d0s4_Ss12s34 tm;
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  getTangentMaterial(&tm, _strain, state, temp);

  // subtract thermal strain
  if(m_tm) {
    for(int i=0; i<6; ++i) strain[i] -= (temp-Tref)*alphas[i];
  }
  else {
    double alpha = (ctett) ? ctett->getValAlt(temp) : ElaLinIsoMat::alpha;
    double e0 = (temp-Tref)*alpha;
    strain[0] -= e0;
    strain[3] -= e0;
    strain[5] -= e0;
  }

  (*stress) = tm||strain;
}

void 
ElaLinIsoMat::getTangentMaterial(Tensor *_tm, Tensor &, double*, double temp)
{
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  if(m_tm) *tm = *m_tm;
  else {
    double E = (ymtt) ? ymtt->getValAlt(temp) : ElaLinIsoMat::E;
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

    (*tm)[0][0] = (*tm)[3][3] = (*tm)[5][5] = lambdadivnu*(1-nu);
    (*tm)[1][1] = (*tm)[2][2] = (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
    (*tm)[0][3] = (*tm)[3][0] = (*tm)[0][5] = (*tm)[5][0] = (*tm)[3][5] = (*tm)[5][3] = lambda;
  }
}

void 
ElaLinIsoMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
{
  Tensor_d0s2_Ss12 & strain = static_cast<Tensor_d0s2_Ss12 &>(_strain);
  Tensor_d0s2_Ss12 * stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 * tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  if(m_tm) {
    *tm = *m_tm;

    // subtract thermal strain
    for(int i=0; i<6; ++i) strain[i] -= (temp-Tref)*alphas[i];
  }
  else {
    double E = (ymtt) ? ymtt->getValAlt(temp) : ElaLinIsoMat::E;
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

    (*tm)[0][0] = (*tm)[3][3] = (*tm)[5][5] = lambdadivnu*(1-nu);
    (*tm)[1][1] = (*tm)[2][2] = (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
    (*tm)[0][3] = (*tm)[3][0] = (*tm)[0][5] = (*tm)[5][0] = (*tm)[3][5] = (*tm)[5][3] = lambda;

    // subtract thermal strain
    double alpha = (ctett) ? ctett->getValAlt(temp) : ElaLinIsoMat::alpha;
    double e0 = (temp-Tref)*alpha;
    strain[0] -= e0;
    strain[3] -= e0;
    strain[5] -= e0;
  }

  (*stress) =(*tm)||strain;
}

void 
ElaLinIsoMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                        double *, double *, double temp, Tensor *, double) const
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);

  if(m_tm) {
    *tm = *m_tm;

    // subtract thermal strain
    for(int i=0; i<6; ++i) enp[i] -= (temp-Tref)*alphas[i];
  }
  else {
    double E = (ymtt) ? ymtt->getValAlt(temp) : ElaLinIsoMat::E;
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

    (*tm)[0][0] = (*tm)[3][3] = (*tm)[5][5] = lambdadivnu*(1-nu);
    (*tm)[1][1] = (*tm)[2][2] = (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
    (*tm)[0][3] = (*tm)[3][0] = (*tm)[0][5] = (*tm)[5][0] = (*tm)[3][5] = (*tm)[5][3] = lambda;

    // subtract thermal strain
    double alpha = (ctett) ? ctett->getValAlt(temp) : ElaLinIsoMat::alpha;
    double e0 = (temp-Tref)*alpha;
    enp[0] -= e0;
    enp[3] -= e0;
    enp[5] -= e0;
  }

  (*stress) = (*tm)||enp;
}

void
ElaLinIsoMat::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                        double *, double *, double temp, Tensor *, double) const
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);

  if(m_tm) {
    // subtract thermal strain
    for(int i=0; i<6; ++i) enp[i] -= (temp-Tref)*alphas[i];

    (*stress) = (*m_tm)||enp;
  }
  else {
    double E = (ymtt) ? ymtt->getValAlt(temp) : ElaLinIsoMat::E;
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    double lambdadivnu = (nu != 0) ? lambda/nu : E;

    Tensor_d0s4_Ss12s34 *tm = new Tensor_d0s4_Ss12s34();

    (*tm)[0][0] = (*tm)[3][3] = (*tm)[5][5] = lambdadivnu*(1-nu);
    (*tm)[1][1] = (*tm)[2][2] = (*tm)[4][4] = lambdadivnu*(1-2*nu)/2;
    (*tm)[0][3] = (*tm)[3][0] = (*tm)[0][5] = (*tm)[5][0] = (*tm)[3][5] = (*tm)[5][3] = lambda;

    // subtract thermal strain
    double alpha = (ctett) ? ctett->getValAlt(temp) : ElaLinIsoMat::alpha;
    double e0 = (temp-Tref)*alpha;
    enp[0] -= e0;
    enp[3] -= e0;
    enp[5] -= e0;

    (*stress) = (*tm)||enp;

    delete tm;
  }
}

double
ElaLinIsoMat::getStrainEnergyDensity(Tensor &_enp, double *, double temp)
{
  Tensor_d0s2_Ss12 &enp = static_cast<Tensor_d0s2_Ss12 &>(_enp);

  if(m_tm) {
    // subtract thermal strain
    for(int i=0; i<6; ++i) enp[i] -= (temp-Tref)*alphas[i];

    Tensor_d0s2_Ss12 stress;
    stress = (*m_tm)||enp;
    return 0.5*(enp[0]*stress[0] + enp[3]*stress[3] + enp[5]*stress[5]
           + 2.0*(enp[1]*stress[1] + enp[2]*stress[2] + enp[4]*stress[4]));
  }
  else {
    double E = (ymtt) ? ymtt->getValAlt(temp) : ElaLinIsoMat::E;
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2*(1+nu));

    // subtract thermal strain
    double alpha = (ctett) ? ctett->getValAlt(temp) : ElaLinIsoMat::alpha;
    double e0 = (temp-Tref)*alpha;
    enp[0] -= e0;
    enp[3] -= e0;
    enp[5] -= e0;

    double I1 = enp.getTrace();
    return lambda/2*I1*I1 + mu*enp.innerProduct();
  }
}

void
ElaLinIsoMat::getMaterialConstants(std::vector<double> &c)
{
  c.resize(2);
  c[0] = E*nu/((1+nu)*(1-2*nu)); // Lamé's first parameter
  c[1] = E/(2*(1+nu));           // Lamé's second parameter (shear modulus)
}

extern LinearStrain linearStrain;

StrainEvaluator *
ElaLinIsoMat::getStrainEvaluator() const
{
  return &linearStrain;
}

extern GreenLagrangeStrain greenLagrangeStrain;

StrainEvaluator *
StVenantKirchhoffMat::getStrainEvaluator() const
{
  return &greenLagrangeStrain;
}

NLMaterial * 
StVenantKirchhoffMat::clone() const
{
  return new StVenantKirchhoffMat(*this);
}

extern LogarithmicStrain logarithmicStrain;

StrainEvaluator *
HenckyMat::getStrainEvaluator() const
{
  return &logarithmicStrain;
}

NLMaterial *
HenckyMat::clone() const
{
  return new HenckyMat(*this);
}

