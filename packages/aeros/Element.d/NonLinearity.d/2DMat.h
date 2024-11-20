#ifndef _2DMAT_H_
#define _2DMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

class ElaLinIsoMat2D : public NLMaterial
{
   protected:
     double rho, E, nu, t, Tref, alpha;

   public:
     ElaLinIsoMat2D(StructProp *p);
     ElaLinIsoMat2D(double _rho, double _E, double _nu, double _t, double _Tref, double _alpha);

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

     GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator() override;

     double getDensity() override { return rho; }

     double getThickness() const override { return t; }

     double getReferenceTemperature() override { return Tref; }
};

// same equation as ElaLinIsoMat2D but with different Green-Lagrange strain evaluator
// (also known as St. Venant-Kirchhoff hyperelastic material
class StVenantKirchhoffMat2D : public ElaLinIsoMat2D
{
  public:
    StVenantKirchhoffMat2D(StructProp *p) : ElaLinIsoMat2D(p) {}
    StVenantKirchhoffMat2D(double rho, double E, double nu, double t, double Tref, double alpha) : ElaLinIsoMat2D(rho, E, nu, t, Tref, alpha) {}

    GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();
};


#endif
