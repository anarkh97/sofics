#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

template<typename Material>
MaterialWrapper<Material>::MaterialWrapper(Material *_mat)
{
  mat = _mat->Clone();
}

template<typename Material>
int
MaterialWrapper<Material>::getNumStates() const
{
  return 0;
}

template<>
inline int
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getNumStates() const
{
  return 19;
}

template<typename Material>
void
MaterialWrapper<Material>::getStress(Tensor *_stress, Tensor &_strain, double*, double temp)
{
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, NULL);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(state) {
    // set the internal variables    
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(state[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(state[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(state[18]);
  }
  
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);
    
  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];
  
  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, NULL);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  delete clone;
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStress is not implemented\n";
}

template<typename Material>
void 
MaterialWrapper<Material>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double temp)
{
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // first elasticity tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, &ltangents);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double* state, double temp)
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(state) {
    // set the internal variables    
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(state[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(state[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(state[18]);
  }
  
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);
    
  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor
  std::vector<double> ltangents; // consistent tangent modulus

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];
  
  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, &ltangents);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];

  delete clone;
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double* state, double temp)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getTangentMaterial is not implemented\n";
}

template<typename Material>
void 
MaterialWrapper<Material>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
{
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // first elasticity tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, &ltangents);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double* state, double temp)
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(state) {
    // set the internal variables
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(state[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(state[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(state[18]);
  }

  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor
  std::vector<double> ltangents; // consistent tangent modulus

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, &ltangents);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];

  delete clone;
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double* state, double temp)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStressAndTangentMaterial is not implemented\n";
}

template<typename Material>
void 
MaterialWrapper<Material>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                     double *, double *statenp, double, Tensor *, double) const
{
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // first elasticity tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, &ltangents);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];
}

template<typename Material>
void
MaterialWrapper<Material>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                                     double *, double *statenp, double, Tensor *, double) const
{
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, NULL);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                                                    double *staten, double *statenp, double,
                                                                    Tensor *, double dt) const
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(staten) {
    // set the internal variables
    // these should be the converged values at the beginning of the time step 
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(staten[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(staten[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(staten[18]);
  }

  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor
  std::vector<double> ltangents; // consistent tangent modulus

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, &ltangents, true, dt);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];

  if(statenp) {
    // get the updated internal variables
    PlasticStrain = clone->GetMaterialPlasticStrain();
    for (int i = 0; i < 9; ++i) statenp[i] = PlasticStrain[i];
    BackStress = clone->GetMaterialBackStress();
    for (int i = 0; i < 9; ++i) statenp[9+i] = BackStress[i];
    statenp[18] = clone->GetMaterialEquivalentPlasticStrain();
  }

  delete clone;
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                                                               double *staten, double *statenp, double,
                                                                               Tensor *, double) const
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate is not implemented\n";
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                                                                    double *staten, double *statenp, double,
                                                                    Tensor *, double dt) const
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(staten) {
    // set the internal variables
    // these should be the converged values at the beginning of the time step 
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(staten[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(staten[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(staten[18]);
  }

  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, NULL, true, dt);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  if(statenp) {
    // get the updated internal variables
    PlasticStrain = clone->GetMaterialPlasticStrain();
    for (int i = 0; i < 9; ++i) statenp[i] = PlasticStrain[i];
    BackStress = clone->GetMaterialBackStress();
    for (int i = 0; i < 9; ++i) statenp[9+i] = BackStress[i];
    statenp[18] = clone->GetMaterialEquivalentPlasticStrain();
  }

  delete clone;
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                                                                               double *staten, double *statenp, double,
                                                                               Tensor *, double) const
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate is not implemented\n";
}

template<typename Material>
void
MaterialWrapper<Material>::initStates(double *state)
{
  for(int i = 0; i < getNumStates(); ++i) state[i] = 0;
}

template<typename Material>
double 
MaterialWrapper<Material>::getDensity()
{ 
  return mat->GetDensityInReference();
} 

template<>
inline double
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getDensity()
{
  return 0.0;
}

template<>
inline double
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getDensity()
{
  return 0.0;
}

extern DeformationGradient deformationGradient;

template<typename Material>
StrainEvaluator *
MaterialWrapper<Material>::getStrainEvaluator() const
{
  return &deformationGradient;
}

template<typename Material>
double
MaterialWrapper<Material>::getEquivPlasticStrain(double *statenp)
{ 
  return 0;
}

template<>
inline double
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getEquivPlasticStrain(double *statenp)
{
  return statenp[18];
}

template<typename Material>
double
MaterialWrapper<Material>::getStrainEnergyDensity(Tensor &_enp, double *, double)
{
  std::cerr << "WARNING: MaterialWrapper<Material>::getStrainEnergyDensity is not implemented\n";
  return 0.0;
}

template<typename Material>
void
MaterialWrapper<Material>::setSDProps(MFTTData *ysst)
{
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::setSDProps(MFTTData *ysst)
{
  double SigmaY = mat->GetYieldStressFromTensionTest();
  if(SigmaY < 0 && ysst && ysst->getID() == -int(SigmaY)) {
    for(int i=0; i<ysst->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData(ysst->getT(i), ysst->getV(i));
    }
  }
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::setSDProps(MFTTData *ysst)
{
  double SigmaY = mat->GetYieldStressFromTensionTest();
  if(SigmaY < 0 && ysst && ysst->getID() == -int(SigmaY)) {
    for(int i=0; i<ysst->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData(ysst->getT(i), ysst->getV(i));
    }
  }
}

template<typename Material>
void
MaterialWrapper<Material>::setSRDProps(MFTTData *yssrt)
{
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::setSRDProps(MFTTData *yssrt)
{
  if(yssrtid > 0 && yssrt && yssrt->getID() == yssrtid) {
    for(int i=0; i<yssrt->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData2(yssrt->getT(i), yssrt->getV(i));
    }
  }
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::setSRDProps(MFTTData *yssrt)
{
  if(yssrtid > 0 && yssrt && yssrt->getID() == yssrtid) {
    for(int i=0; i<yssrt->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData2(yssrt->getT(i), yssrt->getV(i));
    }
  }
}

template<>
inline void
MaterialWrapper<IsotropicLinearElastic>::print(std::ostream &out) const
{
  double rho = mat->GetDensityInReference();
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  out << "IsotropicLinearElastic " << rho << " " << E << " " << nu;
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::print(std::ostream &out) const
{
  double rho = 0.0;
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  double sigmaY = mat->GetYieldStressFromTensionTest();
  double K = mat->GetIsotropicHardeningModulus();
  double H = mat->GetKinematicHardeningModulus();
  double Tol = mat->GetTolerance();
  double epsF = mat->GetEquivalentPlasticStrainAtFailure();
  out << "IsotropicLinearElasticJ2Plastic " << rho << " " << " " << E << " " << nu << " " << sigmaY << " " << K << " " << H
      << " " << Tol << " " << epsF << " " << yssrtid;
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::print(std::ostream &out) const
{
  double rho = 0.0;
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  double sigmaY = mat->GetYieldStressFromTensionTest();
  double K = mat->GetIsotropicHardeningModulus();
  double H = mat->GetKinematicHardeningModulus();
  double Tol = mat->GetTolerance();
  double epsF = mat->GetEquivalentPlasticStrainAtFailure();
  out << "IsotropicLinearElasticJ2PlasticPlaneStress " << rho << " " << E << " " << nu << " " << sigmaY << " " << K << " " << H
      << " " << Tol << " " << epsF << " " << yssrtid;
}

template<typename Material>
void
MaterialWrapper<Material>::getMaterialConstants(std::vector<double> &c)
{
  std::cerr << "material law does not implement getMaterialConstants function\n";
}

