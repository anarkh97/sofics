#ifndef _FABRICMAP_H_
#define _FABRICMAP_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

class FabricMap : public NLMaterial
{
   public:
     enum StrainMeasure { INFINTESIMAL = 0, GREEN_LAGRANGE = 1 };

protected:
	int pxx_map_id, pyy_map_id, sxy_map_id;
	double rho, t, Tref, alpha;
     enum StrainMeasure strain_measure;
	SS2DTData *pxx_map, *pyy_map;
	MFTTData *sxy_map;

public:
     FabricMap(StructProp *p);
     FabricMap(double _rho, int _pxx_map_id, int _pyy_map_id, int _sxy_map_id, double _t, double _Tref, double _alpha,
               StrainMeasure _strain_measure);

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

	void setS1DProps(MFTTData *ss1dt) override;
	void setS2DProps(SS2DTData *ss2dt) override;
};

#endif
