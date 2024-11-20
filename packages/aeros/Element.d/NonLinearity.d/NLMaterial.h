#ifndef _NLMATERIAL_H_
#define _NLMATERIAL_H_

#include <stdexcept>
#include <iostream>
#include <vector>

class StrainEvaluator;
class Tensor;
class Tensor_d0s2_Ss12;
template <typename Tensor> class GenStrainEvaluator;
template <int n> class TwoDTensorTypes;
class MFTTData;
class SS2DTData;

class NLMaterial
{
   public:
     NLMaterial() {}  

     virtual ~NLMaterial() {}

     virtual int getNumStates() const = 0;

     virtual void getTangentMaterial(Tensor *tm, Tensor &strain, double *state, double temp) = 0;

     virtual void getElasticity(Tensor *tm) const = 0;

     virtual void getStress(Tensor *stress, Tensor &strain, double *state, double temp) = 0; // returns conjugate stress

     virtual void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *state, double temp) = 0;

     virtual void updateStates(Tensor& en, Tensor& enp, double *state, double temp) = 0;

     virtual void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                            double *staten, double *statenp, double temp,
                            Tensor *cache, double dt=0) const = 0;

     virtual void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                            double *staten, double *statenp, double temp,
                            Tensor *cache, double dt=0) const = 0;

     virtual void initStates(double *) = 0;

     virtual double getDensity() { return 0; }

     virtual StrainEvaluator * getStrainEvaluator() const { return NULL; } // return default strain evaluator

     virtual GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator() { return NULL; }

     virtual double getEquivPlasticStrain(double *statenp) { return 0; }

     virtual bool getBackStress(double *statenp, Tensor *backstress) { return false; }

     virtual bool getPlasticStrain(double *statenp, Tensor *plasticstrain) { return false; }

     virtual double getDamage(double *statenp) const { return 0; }

     virtual double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp) {
       std::cerr << "material law does not implement getStrainEnergyDensity function\n";
       return 0;
     }

     virtual double getDissipatedEnergy(double *statenp) { return 0; }

     virtual double getThickness() const { return 0; }

     virtual double getReferenceTemperature() { return 0; }

     virtual double getPosdefifyTol() { return -1; }

     virtual void print(std::ostream &out) const {
       throw std::range_error("material law does not implement print function");
     }

     virtual void print2(std::ostream &out) const {}

     virtual NLMaterial * clone() const {
       std::cerr << "material law does not implement clone function\n";
       return 0;
     }

     virtual void setTangentMaterial(double C[6][6]) {
       std::cerr << "material law does not implement setTangentMaterial function\n";
     }

     virtual void setThermalExpansionCoef(double alpha[6]) {
       std::cerr << "material law does not implement setThermalExpansionCoef function\n";
     }

     virtual void setTDProps(MFTTData *ymtt, MFTTData *ctett) {}
     virtual void setSDProps(MFTTData *ysst) {}
     virtual void setSRDProps(MFTTData *yssrt) {}
     virtual void setS1DProps(MFTTData *ss1dt) {}
     virtual void setS2DProps(SS2DTData *ss2dt) {}
     virtual void setEDProps(MFTTData *ymst) {}

     virtual void getMaterialConstants(std::vector<double> &c) {
       std::cerr << "material law does not implement getMaterialConstants function\n";
     }
};

#endif
