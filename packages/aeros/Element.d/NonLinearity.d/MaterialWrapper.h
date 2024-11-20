#ifndef _MATERIALWRAPPER_H_
#define _MATERIALWRAPPER_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Material.d/Material.h>
#include <Material.d/IsotropicLinearElasticJ2PlasticMaterial.h>
#include <Material.d/IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h>

class Tensor;

// This is a wrapper class to interface Adrian Lew's material library with FEM
template<typename Material>
class MaterialWrapper : public NLMaterial
{
protected:
	Material *mat;
	double lambda;
	double mu;
	int yssrtid;
	double posdefifyTol;

public:
	MaterialWrapper(Material*);
	MaterialWrapper(double*);
	~MaterialWrapper() { delete mat; }

	int getNumStates() const override;

	void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;

	void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp) override;

	void getElasticity(Tensor *tm) const override {}

	void updateStates(Tensor &en, Tensor &enp, double *state, double temp) override {}

	void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp) override;

	void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
	               double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

	void integrate(Tensor *stress, Tensor &en, Tensor &enp,
	               double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

	void initStates(double *) override;

	double getDensity() override;

	StrainEvaluator * getStrainEvaluator() const override;

	double getEquivPlasticStrain(double *statenp) override;

	double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp) override;

	double getPosdefifyTol() override { return posdefifyTol; }

	Material* getMaterial() { return mat; }

	void print(std::ostream &out) const override;

	void setSDProps(MFTTData *ysst) override;
	void setSRDProps(MFTTData *yssrt) override;

	void getMaterialConstants(std::vector<double> &c) override;
};

template<>
inline
MaterialWrapper<IsotropicLinearElastic>::MaterialWrapper(double *params)
{
	double rho    = params[0];
	double E      = params[1];
	double nu     = params[2];
	lambda = E*nu/((1.+nu)*(1.-2.*nu));
	mu     = E/(2.*(1.+nu));
	mat = new IsotropicLinearElastic(lambda,mu,rho);
	posdefifyTol = -1;
}

template<>
inline
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::MaterialWrapper(double *params)
{
	double rho    = params[0];
	double E      = params[1];
	double nu     = params[2];
	double sigmaY = params[3];
	double K      = params[4];
	double H      = params[5];
	double Tol    = params[6];
	double epsF   = (params[7] <= 0) ? std::numeric_limits<double>::infinity() : params[7];
	yssrtid   = int(params[8]);
	lambda = E*nu/((1.+nu)*(1.-2.*nu));
	mu     = E/(2.*(1.+nu));
	mat = new IsotropicLinearElasticJ2PlasticMaterial(lambda,mu,sigmaY,K,H,Tol,epsF);
	posdefifyTol = -1;
}

template<>
inline
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::MaterialWrapper(double *params)
{
	double rho    = params[0];
	double E      = params[1];
	double nu     = params[2];
	double sigmaY = params[3];
	double K      = params[4];
	double H      = params[5];
	double Tol    = params[6];
	double epsF   = (params[7] <= 0) ? std::numeric_limits<double>::infinity() : params[7];
	yssrtid   = int(params[8]);
	lambda = E*nu/((1.+nu)*(1.-2.*nu));
	mu     = E/(2.*(1.+nu));
	mat = new IsotropicLinearElasticJ2PlasticPlaneStressMaterial(lambda,mu,sigmaY,K,H,Tol,epsF);
	posdefifyTol = -1;
}

#ifdef _TEMPLATE_FIX_
#include <Element.d/NonLinearity.d/MaterialWrapper.C>
#endif

#endif
