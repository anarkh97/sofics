#include <Element.d/NonLinearity.d/PlaneStressMat.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/BrittleFractureTB.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Element.d/NonLinearity.d/PronyViscoElastic.h>
#include <Math.d/TTensor.h>
#include <Utils.d/NodeSpaceArray.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double _t)
: BaseMaterial(p1), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double _t)
: BaseMaterial(p1, p2), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double _t)
: BaseMaterial(p1, p2, p3), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double _t)
: BaseMaterial(p1, p2, p3, p4), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double _t)
: BaseMaterial(p1, p2, p3, p4, p5), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                                             double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                                             double p9, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                                             double p9, double p10, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, 
                                             double p9, double p10, double p11, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, 
                                             double p9, double p10, double p11, double p12, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                                             double p9, double p10, double p11, double p12, double p13, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                                             double p9, double p10, double p11, double p12, double p13, double p14, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                                             double p9, double p10, double p11, double p12, double p13, double p14, double p15, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15), t(_t)
{
}

template<typename BaseMaterial>
PlaneStressMat<BaseMaterial>::PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                                             double p9, double p10, double p11, double p12, double p13, double p14, double p15,
                                             double p16, double _t)
: BaseMaterial(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16), t(_t)
{
}

template<typename BaseMaterial>
NLMaterial *
PlaneStressMat<BaseMaterial>::clone() const
{
  return new PlaneStressMat<BaseMaterial>(*this);
}

template<typename BaseMaterial>
int
PlaneStressMat<BaseMaterial>::getNumStates() const
{
  return BaseMaterial::getNumStates()+3;
}

template<typename BaseMaterial>
void
PlaneStressMat<BaseMaterial>::initStates(double *state)
{
  BaseMaterial::initStates(state);
  const int i = BaseMaterial::getNumStates();
  state[i] = state[i+1] = state[i+2] = 0;
}

template<typename BaseMaterial>
void
PlaneStressMat<BaseMaterial>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  std::cerr << "PlaneStressMat<BaseMaterial>::getStress is not implemented\n";
}

template<typename BaseMaterial>
void 
PlaneStressMat<BaseMaterial>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                        double *staten, double *statenp, double temp,
                                        Tensor *cache, double dt) const
{
#ifdef USE_EIGEN3
  // Reference:
  // Klinkel, S., & Govindjee, S. (2013). Using finite strain 3D‐material models in beam and shell elements. Engineering Computations.

  SymTensor<double,2> &enp = static_cast<SymTensor<double,2> &>(_enp);       // [exx eyy exy]
  SymTensor<double,2> *stress = static_cast<SymTensor<double,2> *>(_stress); // [sxx syy sxy]
  SymTensor<SymTensor<double,2>,2> *tm = static_cast<SymTensor<SymTensor<double,2>,2> *>(_tm);

  Tensor_d0s2_Ss12 stress3d;    // [sxx, sxy, sxz, syy, syz, szz]
  Tensor_d0s2_Ss12 en3d, enp3d; // [exx, exy, exz, eyy, eyz, ezz]
  Tensor_d0s4_Ss12s34 tm3d;
  Tensor *cache3d = BaseMaterial::getStrainEvaluator()->getCacheInstance();

  Eigen::Map<Eigen::Matrix<double,6,1> > E(&enp3d[0]), S(&stress3d[0]);
  Eigen::Map<Eigen::Matrix<double,6,6,Eigen::RowMajor> > C(&tm3d[0][0]);
  Eigen::Map<Eigen::Vector3d> Em(&enp[0]), Sm(&(*stress)[0]);
  Eigen::Vector3d Ez, Sz; // [ezz exz eyz], [szz sxz syz]
  Eigen::Matrix3d Czz;
  double SzNorm0;

  if(staten) Ez = Eigen::Map<Eigen::Vector3d>(staten+BaseMaterial::getNumStates());
  else Ez.setZero();

  const int maxit = 5;
  const double tol = 1e-5, tolabs = 1e-10;

  for(int i=0; i<maxit; ++i) {
    E << Em[0], Em[2], 0.5*Ez[1], Em[1], 0.5*Ez[2], Ez[0];
    BaseMaterial::integrate(&stress3d, &tm3d, en3d, enp3d, staten, statenp, temp, cache3d, dt);
    Sz << S[5], S[2], S[4];
    if(i==0) SzNorm0 = Sz.norm();
    else if(Sz.norm() < std::max(tolabs, tol*SzNorm0)) break;
    Czz << tm3d[5][5], tm3d[5][2], tm3d[5][4],
           tm3d[2][5], tm3d[2][2], tm3d[2][4],
           tm3d[4][5], tm3d[4][2], tm3d[4][4];
    Ez -= Czz.fullPivLu().solve(Sz);
  }

  Sm << t*S[0], t*S[3], t*S[1];

  Eigen::Matrix3d Cmm, Cmz, Czm;
  Cmm << tm3d[0][0], tm3d[0][3], tm3d[0][1],
         tm3d[3][0], tm3d[3][3], tm3d[3][1],
         tm3d[1][0], tm3d[1][3], tm3d[1][1];
  Cmz << tm3d[0][5], tm3d[0][2], tm3d[0][4],
         tm3d[3][5], tm3d[3][2], tm3d[3][4],
         tm3d[1][5], tm3d[1][2], tm3d[1][4];
  Czm << tm3d[5][0], tm3d[5][3], tm3d[5][1],
         tm3d[2][0], tm3d[2][3], tm3d[2][1],
         tm3d[4][0], tm3d[4][3], tm3d[4][1];
  Eigen::Matrix3d Cpsc = Cmm - Cmz*Czz.fullPivLu().solve(Czm);
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) (*tm)[i][j] = t*Cpsc(i,j);

  if(statenp) Eigen::Map<Eigen::Vector3d>(statenp+BaseMaterial::getNumStates()) = Ez;

  if(cache) delete cache;
#else
  std::cerr << "PlaneStressMat<BaseMaterial>::integrate is not implemented\n";
#endif
}

template<typename BaseMaterial>
void
PlaneStressMat<BaseMaterial>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                                        double *staten, double *statenp, double temp,
                                        Tensor *cache, double dt) const
{
#ifdef USE_EIGEN3
  // Reference:
  // Klinkel, S., & Govindjee, S. (2013). Using finite strain 3D‐material models in beam and shell elements. Engineering Computations.

  SymTensor<double,2> &enp = static_cast<SymTensor<double,2> &>(_enp);       // [exx eyy exy]
  SymTensor<double,2> *stress = static_cast<SymTensor<double,2> *>(_stress); // [sxx syy sxy]

  Tensor_d0s2_Ss12 stress3d;    // [sxx, sxy, sxz, syy, syz, szz]
  Tensor_d0s2_Ss12 en3d, enp3d; // [exx, exy, exz, eyy, eyz, ezz]
  Tensor_d0s4_Ss12s34 tm3d;
  Tensor *cache3d = BaseMaterial::getStrainEvaluator()->getCacheInstance();

  Eigen::Map<Eigen::Matrix<double,6,1> > E(&enp3d[0]), S(&stress3d[0]);
  Eigen::Map<Eigen::Matrix<double,6,6,Eigen::RowMajor> > C(&tm3d[0][0]);
  Eigen::Map<Eigen::Vector3d> Em(&enp[0]), Sm(&(*stress)[0]);
  Eigen::Vector3d Ez, Sz; // [ezz exz eyz], [szz sxz syz]
  Eigen::Matrix3d Czz;
  double SzNorm0;

  if(staten) Ez = Eigen::Map<Eigen::Vector3d>(staten+BaseMaterial::getNumStates());
  else Ez.setZero();

  const int maxit = 5;
  const double tol = 1e-5, tolabs = 1e-10;

  for(int i=0; i<maxit; ++i) {
    E << Em[0], Em[2], 0.5*Ez[1], Em[1], 0.5*Ez[2], Ez[0];
    BaseMaterial::integrate(&stress3d, &tm3d, en3d, enp3d, staten, statenp, temp, cache3d, dt);
    Sz << S[5], S[2], S[4];
    if(i==0) SzNorm0 = Sz.norm();
    else if(Sz.norm() < std::max(tolabs, tol*SzNorm0)) break;
    Czz << tm3d[5][5], tm3d[5][2], tm3d[5][4],
           tm3d[2][5], tm3d[2][2], tm3d[2][4],
           tm3d[4][5], tm3d[4][2], tm3d[4][4];
    Ez -= Czz.fullPivLu().solve(Sz);
  }

  Sm << t*S[0], t*S[3], t*S[1];

  if(statenp) Eigen::Map<Eigen::Vector3d>(statenp+BaseMaterial::getNumStates()) = Ez;

  if(cache3d) delete cache3d;
#else
  std::cerr << "PlaneStressMat<BaseMaterial>::integrate is not implemented\n";
#endif
}

template<typename BaseMaterial>
void
PlaneStressMat<BaseMaterial>::print(std::ostream &out) const
{
  out << "PlaneStress";
  BaseMaterial::print(out);
  out << " " << t;
}

extern LinearStrain2D<9> linStrain2D;

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<ElaLinIsoMat>::getGenStrainEvaluator()
{
  return &linStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >::getGenStrainEvaluator()
{
  return &linStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<ElasPlasKinHardMat<0> >::getGenStrainEvaluator()
{
  return &linStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >::getGenStrainEvaluator()
{
  return &linStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >::getGenStrainEvaluator()
{ 
  return &linStrain2D;
}

extern GLStrain2D<9> glStrain2D;

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<StVenantKirchhoffMat>::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >::getGenStrainEvaluator()
{ 
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<ElasPlasKinHardMat<1> >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<NeoHookeanMat>::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<NeoHookeanMat> >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<PronyViscoElastic<NeoHookeanMat> >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> > >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<MooneyRivlinMat>::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<MooneyRivlinMat> >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<PronyViscoElastic<MooneyRivlinMat> >::getGenStrainEvaluator()
{
  return &glStrain2D;
}

template<>
inline
GenStrainEvaluator<TwoDTensorTypes<9> > *
PlaneStressMat<BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> > >::getGenStrainEvaluator()
{
  return &glStrain2D;
}
