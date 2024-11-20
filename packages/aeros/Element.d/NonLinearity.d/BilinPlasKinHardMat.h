#ifndef _BILINPLASKINHARDMAT_H_
#define _BILINPLASKINHARDMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Utils.d/MFTT.h>
#include <limits>

class StructProp;

template<int e>
class ElasPlasKinHardMat : public NLMaterial
{
protected:
	// Ep is the tangent modulus from the uniaxial strain-stress curve (Et for DYNA3D Material type 3)
	// theta is hardening parameter which specifies and arbitrary combination of isotropic and kinematic hardening (beta for DYNA3D Material type 3)
	// theta = 0 (default) is purely kinematic hardening, while theta = 1 is purely isotropic hardening
	double rho, E, nu, Ep, sigE, theta;
	// alpha is the thermal expansion coefficient and Tref is the reference temperature
	double alpha, Tref;
	// epsF is the equivalent plastic strain at failure
	double epsF;
	// strain dependent material properties
	MFTTData *ysst;
	// tolerance for convergence of nonlinear solve
	double tol;
	// rate dependent material properties
	int yssrtid;
	MFTTData *yssrt;

public:
	ElasPlasKinHardMat(StructProp *p);
	ElasPlasKinHardMat(double _rho, double _E, double _nu, double _Ep, double _sigE, double _theta = 0,
					   double _Tref = 0, double _alpha = 0, double _epsF = std::numeric_limits<double>::infinity(),
					   double _tol = 1e-6, int _yssrtid = 0)
	{ rho = _rho; E = _E; nu = _nu; Ep = _Ep; sigE = _sigE; theta = _theta; Tref = _Tref; alpha = _alpha; epsF = _epsF;
		tol = _tol; yssrtid = _yssrtid; ysst = NULL; yssrt = NULL; }

	void getStress(Tensor *stress, Tensor &strain, double *, double temp) override;

	void getTangentMaterial(Tensor *tm, Tensor &strain, double *, double temp) override;

	void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *, double temp) override;

	void getElasticity(Tensor *tm) const override;

	void updateStates(Tensor &en, Tensor &enp, double *state, double temp) override;

	int getNumStates() const override { return 13; } // the internal variables are : the plastic strain (6 doubles),
	// the center of the yield surface in sigma space (6 doubles),
	// and the equivalent plastic strain (1 double)

	void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
				   double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

	void integrate(Tensor *stress, Tensor &en, Tensor &enp,
				   double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

	void initStates(double *) override;

	double getDensity() override { return rho; }

	double getReferenceTemperature() override { return Tref; }

	StrainEvaluator * getStrainEvaluator() const override;

	double getEquivPlasticStrain(double *statenp) override { return statenp[12]; }

	bool getBackStress(double *statenp, Tensor *backstress) override;

	bool getPlasticStrain(double *statenp, Tensor *plasticstrain) override;

	double getDamage(double *statenp) const override { return (statenp[12] >= epsF) ? 1 : 0; }

	double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp = 0.0) override;

	double getDissipatedEnergy(double *statenp) override;

	void print(std::ostream &out) const override;

	void print2(std::ostream &out) const override;

	void setSDProps(MFTTData *_ysst) override { if(sigE < 0 && _ysst && _ysst->getID() == -int(sigE)) ysst = _ysst; }

	void setSRDProps(MFTTData *_yssrt) override { if(yssrtid > 0 && _yssrt && _yssrt->getID() == yssrtid) yssrt = _yssrt; }
};

typedef ElasPlasKinHardMat<0> BilinPlasKinHardMat;
typedef ElasPlasKinHardMat<1> FiniteStrainPlasKinHardMat;
typedef ElasPlasKinHardMat<2> LogStrainPlasKinHardMat;

#endif
