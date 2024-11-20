#ifndef _OGDENMAT_H_
#define _OGDENMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class OgdenMat : public NLMaterial
{
  protected:
    // isotropic material properties
    double rho; // density
    double mu[9], alpha[9]; // material properties characterizing distortional response
    int m, n; // number of terms in the Ogden series
    double K[9]; // material properties characterizing volumetric response

  public:
    template<int _m, int _n>
    OgdenMat(double _rho, double (&_mu)[_m], double (&_alpha)[_m], double (&_K)[_n]);

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
