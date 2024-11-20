#include <Element.d/NonLinearity.d/2DMat.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Math.d/TTensor.h>

ElaLinIsoMat2D::ElaLinIsoMat2D(StructProp *p)
{
  E = p->E;
  nu = p->nu;
  t = p->eh;
  rho = p->rho;
  Tref = p->Ta;
  alpha = p->W;
}

ElaLinIsoMat2D::ElaLinIsoMat2D(double _rho, double _E, double _nu, double _t,
                               double _Tref, double _alpha)
{
  rho = _rho; nu = _nu; E = _E; t = _t; Tref = _Tref; alpha = _alpha;
}

void
ElaLinIsoMat2D::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  SymTensor<SymTensor<double,2>,2> tm;
  SymTensor<double,2> & strain = static_cast<SymTensor<double,2> &>(_strain);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);

  double e0 = alpha*(temp-Tref);
  strain[0] -= e0;
  strain[1] -= e0;

  getTangentMaterial(&tm, _strain, state, temp);
  (*stress) = tm||strain;
}

void 
ElaLinIsoMat2D::getTangentMaterial(Tensor *_tm, Tensor &, double*, double temp)
{
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2>  *>(_tm);

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;
}

void 
ElaLinIsoMat2D::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
{
  SymTensor<double,2> & strain = static_cast<SymTensor<double,2> &>(_strain);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2> *>(_tm);

  double e0 = alpha*(temp-Tref);
  strain[0] -= e0;
  strain[1] -= e0;

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;

  (*stress) =(*tm)||strain;
}

void 
ElaLinIsoMat2D::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                          double *staten, double *statenp, double temp,
                          Tensor *, double) const
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2> *>(_tm);

  double e0 = alpha*(temp-Tref);
  enp[0] -= e0;
  enp[1] -= e0;

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;

  (*stress) = (*tm)||enp;
}

void
ElaLinIsoMat2D::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                          double *staten, double *statenp, double temp,
                          Tensor *, double) const
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = new SymTensor<SymTensor<double,2>,2>();

  double e0 = alpha*(temp-Tref);
  enp[0] -= e0;
  enp[1] -= e0;

  double E11 = t*E/(1-nu*nu);
  double E12 = nu*E11;
  double E33 = t*E/(1+nu);
  (*tm)[0][0] = E11;
  (*tm)[1][1] = E11;
  (*tm)[2][2] = E33/2;
  (*tm)[0][1] = (*tm)[1][0] = E12;
  (*tm)[0][2] = (*tm)[2][0] = 0;
  (*tm)[1][2] = (*tm)[2][1] = 0;

  (*stress) = (*tm)||enp;

  delete tm;
}

extern LinearStrain2D<9> linStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
ElaLinIsoMat2D::getGenStrainEvaluator()
{
  return &linStrain2D;
}

extern GLStrain2D<9> glStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
StVenantKirchhoffMat2D::getGenStrainEvaluator()
{
  return &glStrain2D;
}

