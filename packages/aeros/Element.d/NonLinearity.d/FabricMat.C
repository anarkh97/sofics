#include <Element.d/NonLinearity.d/FabricMat.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Math.d/TTensor.h>

FabricMat::FabricMat(StructProp *p)
{
  t = p->eh;
  rho = p->rho;
  Tref = p->Ta;
  alpha = p->W;
  strain_measure = GREEN_LAGRANGE;
}

FabricMat::FabricMat(double _rho, int _x_ymst_id, int _y_ymst_id, double _Gxy, double _nuxy, double _nuyx,
                     double _t, double _Tref, double _alpha, StrainMeasure _strain_measure)
{
  rho = _rho; x_ymst_id = _x_ymst_id; y_ymst_id = _y_ymst_id; Gxy = _Gxy; nuxy = _nuxy; nuyx = _nuyx;
  t = _t; Tref = _Tref; alpha = _alpha; strain_measure = _strain_measure; x_ymst = y_ymst = 0;
}

void
FabricMat::setEDProps(MFTTData *_ymst)
{
  if(x_ymst_id != 0 && _ymst && _ymst->getID() == x_ymst_id) x_ymst = _ymst;
  if(y_ymst_id != 0 && _ymst && _ymst->getID() == y_ymst_id) y_ymst = _ymst;
}

void
FabricMat::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  SymTensor<SymTensor<double,2>,2> tm;
  SymTensor<double,2> & strain = static_cast<SymTensor<double,2> &>(_strain);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);

  double e0 = alpha*(temp-Tref);
  strain[0] -= e0;
  strain[1] -= e0;

  double Ex = (x_ymst) ? x_ymst->getValAlt(strain[0]) : 0;
  double Ey = (y_ymst) ? y_ymst->getValAlt(strain[1]) : 0;

  (*stress)[0] = t/(1-nuxy*nuyx)*(Ex*strain[0] + nuyx*Ex*strain[1]);
  (*stress)[1] = t/(1-nuxy*nuyx)*(nuxy*Ey*strain[0] + Ey*strain[1]);
  (*stress)[2] = t*2*Gxy*strain[2];
}

void 
FabricMat::getTangentMaterial(Tensor *, Tensor &, double*, double)
{
  std::cerr << "WARNING: FabricMat::getTangentMaterial is not implemented\n";
}

void 
FabricMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
{
  std::cerr << "WARNING: FabricMat::getStressAndTangentMaterial is not implemented\n";
}

void 
FabricMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                     double *staten, double *statenp, double temp,
                     Tensor *, double) const
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2> *>(_tm);

  double e0 = alpha*(temp-Tref);
  enp[0] -= e0;
  enp[1] -= e0;

  double Ex = (x_ymst) ? x_ymst->getValAlt(enp[0]) : 0;
  double Ey = (y_ymst) ? y_ymst->getValAlt(enp[1]) : 0;

  (*stress)[0] = t/(1-nuxy*nuyx)*(Ex*enp[0] + nuyx*Ex*enp[1]);
  (*stress)[1] = t/(1-nuxy*nuyx)*(nuxy*Ey*enp[0] + Ey*enp[1]);
  (*stress)[2] = t*2*Gxy*enp[2];

  (*tm)[0][0] = t/(1-nuxy*nuyx)*Ex;
  (*tm)[0][1] = t/(1-nuxy*nuyx)*nuyx*Ex;
  (*tm)[1][0] = t/(1-nuxy*nuyx)*nuxy*Ey;
  (*tm)[1][1] = t/(1-nuxy*nuyx)*Ey;
  (*tm)[2][2] = t*Gxy;
  (*tm)[0][2] = (*tm)[2][0] = (*tm)[1][2] = (*tm)[2][1] = 0;
}

void
FabricMat::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                     double *staten, double *statenp, double temp,
                     Tensor *, double) const
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);

  double e0 = alpha*(temp-Tref);
  enp[0] -= e0;
  enp[1] -= e0;

  double Ex = (x_ymst) ? x_ymst->getValAlt(enp[0]) : 0;
  double Ey = (y_ymst) ? y_ymst->getValAlt(enp[1]) : 0;

  (*stress)[0] = t/(1-nuxy*nuyx)*(Ex*enp[0] + nuyx*Ex*enp[1]);
  (*stress)[1] = t/(1-nuxy*nuyx)*(nuxy*Ey*enp[0] + Ey*enp[1]);
  (*stress)[2] = t*2*Gxy*enp[2];
}

extern LinearStrain2D<9> linStrain2D;
extern GLStrain2D<9> glStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
FabricMat::getGenStrainEvaluator()
{
  switch(strain_measure) {
    case INFINTESIMAL: return &linStrain2D;
    default:
    case GREEN_LAGRANGE: return &glStrain2D;
  }
}

