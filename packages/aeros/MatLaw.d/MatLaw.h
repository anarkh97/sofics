#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <iostream>
#include <cstdlib>

// This is an example demonstrating how to define a material law to use with READ sub-command under MATLAW in Aero-S

class MatLaw : public NLMaterial
{
    double rho, E, nu, alpha, Tref;

  public:
    MatLaw(double _rho, double _E, double _nu, double _Tref, double _alpha);
    int getNumStates() const override { return 0; }
    void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp) override;
    void getElasticity(Tensor *tm) const override {}
    void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;
    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp) override;
    void updateStates(Tensor& en, Tensor& enp, double *state, double temp) {}
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp) const override;
    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp) const override;
    void initStates(double *) override {}
    StrainEvaluator * getStrainEvaluator() const override;
};

extern "C" NLMaterial * materialLoader(int nval, double *v)
{
  if(nval < 3) { std::cerr << " ERROR: materialLoader is passed too few arguments\n"; exit(-1); }
  double rho = v[0];
  double E = v[1];
  double nu = v[2];
  double Tref = (nval >= 4) ? v[3] : 0;
  double alpha = (nval >= 5) ? v[4] : 0;
  return new MatLaw(rho, E, nu, Tref, alpha);
}

