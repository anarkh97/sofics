#ifndef _NEOHOOKEANMAT_H_
#define _NEOHOOKEANMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class NeoHookeanMat : public NLMaterial
{
protected:
	double rho; // density
	double lambda, mu; // material properties

public:
	NeoHookeanMat(double _rho, double _E, double _nu);

	int getNumStates() const override { return 0; }

	void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;

	void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp) override;

	void getElasticity(Tensor *tm) const override {};

	void updateStates(Tensor &en, Tensor &enp, double *state, double temp) override {};

	void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp) override;

	void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
	               double *staten, double *statenp, double temp,
	               Tensor *cache, double dt) const override;

	void integrate(Tensor *stress, Tensor &en, Tensor &enp,
	               double *staten, double *statenp, double temp,
	               Tensor *cache, double dt) const override;

	void initStates(double *) override {};

	double getDensity() override { return rho; }

	StrainEvaluator * getStrainEvaluator() const override;

	double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp) override;

	void print(std::ostream &out) const override;

	NLMaterial * clone() const override;

	void getMaterialConstants(std::vector<double> &c) override;
};

#endif
