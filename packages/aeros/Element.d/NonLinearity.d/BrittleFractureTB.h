#ifndef _BRITTLEFRACTURETB_H_
#define _BRITTLEFRACTURETB_H_

class Tensor;

template<typename BaseMaterial>
class BrittleFractureTB : public BaseMaterial
{
    double maxprs;
    double exponent; 
    double Kf; //stress impulse

  public:
    BrittleFractureTB(StructProp *p) : BaseMaterial(p) {}
    BrittleFractureTB(double p1, 
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double p10, double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double p10, double p11, double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double p10, double p11, double p12, double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double p10, double p11, double p12, double p13, double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double p10, double p11, double p12, double p13, double p14, double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double p10, double p11, double p12, double p13, double p14, double p15, double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9,
                      double p10, double p11, double p12, double p13, double p14, double p15, double p16, double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}

    int getNumStates() const override;

    void initStates(double *) override;

    void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, Tensor *cache, double dt) const override;

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, Tensor *cache, double dt) const override;

    double getDamage(double *statenp) const override;

    void print2(std::ostream &out) const override;

    NLMaterial * clone() const override;
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/BrittleFractureTB.C>
#endif

#endif
