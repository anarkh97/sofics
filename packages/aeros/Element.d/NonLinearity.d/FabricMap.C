#include <Element.d/NonLinearity.d/FabricMap.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Math.d/TTensor.h>

//reference:
// Thomas Borrvall, Curtis Ehle, and Troy Stratton, "A Fabric Material Model with Stress Map Functionality in LS-DYNA"
// 10th European LS-DYNA Conference 2015, Würzburg, Germany

FabricMap::FabricMap(StructProp *p)
{
  t = p->eh;
  rho = p->rho;
  Tref = p->Ta;
  alpha = p->W;
  strain_measure = GREEN_LAGRANGE;
}

FabricMap::FabricMap(double _rho, int _pxx_map_id, int _pyy_map_id, int _sxy_map_id,
                     double _t, double _Tref, double _alpha, StrainMeasure _strain_measure)
{
  rho = _rho; pxx_map_id = _pxx_map_id; pyy_map_id = _pyy_map_id; sxy_map_id = _sxy_map_id;
  t = _t; Tref = _Tref; alpha = _alpha; strain_measure = _strain_measure; pxx_map = pyy_map = 0; sxy_map = 0;
}

void
FabricMap::setS1DProps(MFTTData *_ss1dt)
{
  if(sxy_map_id != 0 && _ss1dt && _ss1dt->getID() == sxy_map_id) sxy_map = _ss1dt;
}

void
FabricMap::setS2DProps(SS2DTData *_ss2dt)
{
  if(pxx_map_id != 0 && _ss2dt && _ss2dt->getID() == pxx_map_id) pxx_map = _ss2dt;
  if(pyy_map_id != 0 && _ss2dt && _ss2dt->getID() == pyy_map_id) pyy_map = _ss2dt;
}  

void
FabricMap::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  SymTensor<SymTensor<double,2>,2> tm;
  SymTensor<double,2> & strain = static_cast<SymTensor<double,2> &>(_strain);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);

  double e0 = alpha*(temp-Tref);
  strain[0] -= e0;
  strain[1] -= e0;

  double ex = sqrt(1+2*strain[0]) - 1, 
         ey = sqrt(1+2*strain[1]) - 1; // engineering strains

  if(pxx_map) {
    if(!pxx_map->engineeringFlag)
      (*stress)[0] = t*pxx_map->getValAlt(strain[0], strain[1]);
    else
      (*stress)[0] = t*pxx_map->getValAlt(ex, ey)/(1+ex);
  }
  else {
    (*stress)[0] = 0;
  }

  if(pyy_map) {
    if(!pyy_map->engineeringFlag)
      (*stress)[1] = t*pyy_map->getValAlt(strain[1], strain[0]);
    else
      (*stress)[1] = t*pyy_map->getValAlt(ey, ex)/(1+ey);
  }
  else {
    (*stress)[1] = 0;
  }

  (*stress)[2] = (sxy_map) ? t*sxy_map->getValAlt(strain[2]) : 0;
}

void 
FabricMap::getTangentMaterial(Tensor *, Tensor &, double*, double)
{
  std::cerr << "WARNING: FabricMap::getTangentMaterial is not implemented\n";
}

void 
FabricMap::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
{
  std::cerr << "WARNING: FabricMap::getStressAndTangentMaterial is not implemented\n";
}

void 
FabricMap::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                     double *staten, double *statenp, double temp,
                     Tensor *, double) const
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);
  SymTensor<SymTensor<double,2>,2> * tm = static_cast<SymTensor<SymTensor<double,2>,2> *>(_tm);

  double e0 = alpha*(temp-Tref);
  enp[0] -= e0;
  enp[1] -= e0;

  double ex = sqrt(1+2*enp[0]) - 1,
         ey = sqrt(1+2*enp[1]) - 1; // engineering strains

  if(pxx_map) {
    if(!pxx_map->engineeringFlag) {
      double S, dSdE[2];
      pxx_map->getValAndSlopeAlt(enp[0], enp[1], &S, &(dSdE[0]), &(dSdE[1]));
      (*stress)[0] = t*S;
      (*tm)[0][0] = t*dSdE[0];
      (*tm)[0][1] = t*dSdE[1];
    }
    else {
      double P, dPde[2];
      pxx_map->getValAndSlopeAlt(ex, ey, &P, &(dPde[0]), &(dPde[1]));
      (*stress)[0] = t*P/(1+ex);
      (*tm)[0][0] = t*1/sqrt(1+2*enp[0])*(dPde[0] - (*stress)[0])/(1+ex);
      (*tm)[0][1] = t*1/sqrt(1+2*enp[1])*dPde[1]/(1+ex);
    }
  }
  else {
    (*stress)[0] = 0;
    (*tm)[0][0] = 0;
    (*tm)[0][1] = 0;
  }

  if(pyy_map) {  
    if(!pyy_map->engineeringFlag) {
      double S, dSdE[2];
      pyy_map->getValAndSlopeAlt(enp[1], enp[0], &S, &(dSdE[1]), &(dSdE[0]));
      (*stress)[1] = t*S;
      (*tm)[1][0] = t*dSdE[0];
      (*tm)[1][1] = t*dSdE[1];
    }
    else {
      double P, dPde[2];
      pyy_map->getValAndSlopeAlt(ey, ex, &P, &(dPde[1]), &(dPde[0]));
      (*stress)[1] = t*P/(1+ey);
      (*tm)[1][0] = t*1/sqrt(1+2*enp[0])*dPde[0]/(1+ey);
      (*tm)[1][1] = t*1/sqrt(1+2*enp[1])*(dPde[1]-(*stress)[1])/(1+ey);
    }
  }
  else {
    (*stress)[1] = 0;
    (*tm)[1][0] = 0;
    (*tm)[1][1] = 0; 
  }

  if(sxy_map) { 
    double S, dSdE;
    sxy_map->getValAndSlopeAlt(enp[2], &S, &dSdE);
    (*stress)[2] = t*S;
    (*tm)[2][2] = t*dSdE/2;
  } 
  else { 
    (*stress)[2] = 0;
    (*tm)[2][2] = 0;
  }
  (*tm)[0][2] = (*tm)[2][0] = (*tm)[1][2] = (*tm)[2][1] = 0;
}

void
FabricMap::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                     double *staten, double *statenp, double temp,
                     Tensor *, double) const
{
  SymTensor<double,2> & enp = static_cast<SymTensor<double,2> &>(_enp);
  SymTensor<double,2> * stress = static_cast<SymTensor<double,2> *>(_stress);

  double e0 = alpha*(temp-Tref);
  enp[0] -= e0;
  enp[1] -= e0;

  double ex = sqrt(1+2*enp[0]) - 1,
         ey = sqrt(1+2*enp[1]) - 1; // engineering strains

  if(pxx_map) {
    if(!pxx_map->engineeringFlag)
      (*stress)[0] = t*pxx_map->getValAlt(enp[0], enp[1]);
    else
      (*stress)[0] = t*pxx_map->getValAlt(ex, ey)/(1+ex);
  }
  else {
    (*stress)[0] = 0;
  }

  if(pyy_map) {
    if(!pyy_map->engineeringFlag)
      (*stress)[1] = t*pyy_map->getValAlt(enp[1], enp[0]);
    else
      (*stress)[1] = t*pyy_map->getValAlt(ey, ex)/(1+ey);
  }
  else {
    (*stress)[1] = 0;
  }

  (*stress)[2] = (sxy_map) ? t*sxy_map->getValAlt(enp[2]) : 0;
}

extern LinearStrain2D<9> linStrain2D;
extern GLStrain2D<9> glStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
FabricMap::getGenStrainEvaluator()
{
  switch(strain_measure) {
    case INFINTESIMAL: return &linStrain2D;
    default:
    case GREEN_LAGRANGE: return &glStrain2D;
  }
}

