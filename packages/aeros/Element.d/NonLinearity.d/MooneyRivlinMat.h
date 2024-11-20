#ifndef _MOONEYRIVLINMAT_H_
#define _MOONEYRIVLINMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class MooneyRivlinMat : public NLMaterial
{
  protected:
    double rho; // density
    double mu1, mu2, kappa, c00, c01, c10, M; // material properties

  public:
    MooneyRivlinMat(double _rho, double _mu1, double _mu2, double _kappa);

    int getNumStates() const override { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp) override;

    void getElasticity(Tensor *tm) const override {};

    void updateStates(Tensor &en, Tensor &enp, double *state, double temp) override {};

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp) override;
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp,
                   Tensor *cache, double dt=0) const override;

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp,
                   Tensor *cache, double dt=0) const override;

    void initStates(double *) override {};

    double getDensity() override { return rho; }

    StrainEvaluator * getStrainEvaluator() const override;

    double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp) override;

    void print(std::ostream &out) const override;

    NLMaterial * clone() const override;

    void getMaterialConstants(std::vector<double> &c) override;
};

#endif
