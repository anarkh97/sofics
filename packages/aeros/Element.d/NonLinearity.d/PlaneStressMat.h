#ifndef _PLANESTRESSMAT_H_
#define _PLANESTRESSMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

template<class BaseMaterial>
class PlaneStressMat : public BaseMaterial
{
     double t;

   public:
     PlaneStressMat(double p1, double t);
     PlaneStressMat(double p1, double p2, double t);
     PlaneStressMat(double p1, double p2, double p3, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10,
                    double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10, 
                    double p11, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10, 
                    double p11, double p12, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10,
                    double p11, double p12, double p13, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10,
                    double p11, double p12, double p13, double p14, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10,
                    double p11, double p12, double p13, double p14, double p15, double t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9, double p10,
                    double p11, double p12, double p13, double p14, double p15, double p16, double t);

     int getNumStates() const override;

     void initStates(double *) override;

     void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;

     void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp,
                    Tensor *cache, double dt=0) const override;

     void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp,
                    Tensor *cache, double dt=0) const override;

     void print(std::ostream &out) const override;

     NLMaterial * clone() const override;

     GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator() override;

     double getThickness() const override { return t; }
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/PlaneStressMat.C>
#endif

#endif
