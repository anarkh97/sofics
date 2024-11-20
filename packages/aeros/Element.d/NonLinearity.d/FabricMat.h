#ifndef _FABRICMAT_H_
#define _FABRICMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

class FabricMat : public NLMaterial
{
public:
	enum StrainMeasure { INFINTESIMAL = 0, GREEN_LAGRANGE = 1 };

protected:
	int x_ymst_id, y_ymst_id;
	double rho, Gxy, nuxy, nuyx, t, Tref, alpha;
	enum StrainMeasure strain_measure;
	MFTTData *x_ymst, *y_ymst;

public:
	FabricMat(StructProp *p);
	FabricMat(double _rho, int _x_ymst_id, int _y_ymst_id, double _Gxy, double _nuxy, double _nyx, double _t,
	          double _Tref, double _alpha, StrainMeasure _strain_measure);

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

	void setEDProps(MFTTData *ymst) override;
};

#endif
