#ifndef _PRONYVISCOELASTIC_H_
#define _PRONYVISCOELASTIC_H_

#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NeoHookeanMat.h>
#include <Element.d/NonLinearity.d/MooneyRivlinMat.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Element.d/NonLinearity.d/FabricMap.h>
#include <Element.d/NonLinearity.d/FabricMat.h>
#include <Math.d/TTensor.h>

class Tensor;


struct stress_policy_green_langrage
{
    // Stress tensor is symmetric and needs 6 elements to store
    typedef Tensor_d0s4_Ss12s34 d0s4_S;
    typedef Tensor_d0s2_Ss12 d0s2_S;
    const static std::size_t stride = 6;
};

struct stress_policy_stretches
{
    // Stress tensor is diagonal and needs 3 elements to store (For example as used in the OgdenMat)
    typedef Tensor_d0s4_Ss12s34_diag d0s4_S;
    typedef Tensor_d0s2_Ss12_diag    d0s2_S;
    const static std::size_t stride = 3;
};

struct stress_policy_2d
{
    typedef SymTensor<SymTensor<double,2>,2> d0s4_S;
    typedef SymTensor<double,2> d0s2_S;
    const static std::size_t stride = 3;
};

template<typename Material, class tensor_policy = stress_policy_green_langrage>
class PronyViscoElastic : public Material
{
  public:
    PronyViscoElastic(double p1,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6, double p7,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, int i1, int i2, int i3, double p2, double p3, double p4, FabricMap::StrainMeasure i4,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3); // for FabricMap
    PronyViscoElastic(double p1, int i1, int i2, double p2, double p3, double p4, double p5, double p6, double p7, FabricMat::StrainMeasure i3,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3); // for FabricMat

    int getNumStates() const override;

    void initStates(double *) override;

    void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

    void print(std::ostream &out) const override;

    NLMaterial * clone() const override;

  private:
    double ginf;
    double g1, tau1;
    double g2, tau2;
    double g3, tau3;
};


#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/PronyViscoElastic.C>
#endif

#endif
