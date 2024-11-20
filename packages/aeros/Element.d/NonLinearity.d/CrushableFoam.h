#ifndef _CRUSHABLEFOAM_H_
#define _CRUSHABLEFOAM_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Utils.d/MFTT.h>
#include <limits>
#include <cmath>

class CrushableFoam : public NLMaterial
{
protected:
	// base material properties -- not used.
	//double rho, Ep, sigE, theta;
	double E, nu;
    	// foam properties
    	double rhoF, sigP, alpha2, beta, gamma, epsD, alphaDF;
	// alpha is the thermal expansion coefficient and Tref is the reference temperature
	double alpha = 0, Tref = 0;
	// epsF is the equivalent plastic strain at failure
	double epsF;
	// strain dependent material properties
	MFTTData *ysst;
	// tolerance for convergence of nonlinear solve
	double tol;

public:
	CrushableFoam(double _rhoF, double _EF, double _nuF, double _sigP, double _alpha2,
				  double _gamma, double _beta, double _epsD, double _alphaDF = sqrt(9./2), 
                  double _tol = 1e-6, double _epsF=std::numeric_limits<double>::infinity())
	{ 
		rhoF = _rhoF; E = _EF; nu = _nuF; sigP = _sigP; alpha2 = _alpha2; gamma = _gamma;
		beta = _beta; epsD = _epsD; alphaDF = _alphaDF; tol = _tol; epsF = _epsF; ysst = NULL;
	}

	void getStress(Tensor *stress, Tensor &strain, double *, double temp) override;

	void getTangentMaterial(Tensor *tm, Tensor &strain, double *, double temp) override;

	void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *, double temp) override;

	void getElasticity(Tensor *tm) const override;

	void updateStates(Tensor &en, Tensor &enp, double *state, double temp) override;

	// the internal variables are : the plastic rate of deformation (6 doubles),
	// the cauchy stress at time t_n (6 doubles),
	// and the equivalent plastic strain (1 double)
	int getNumStates() const override { return 13; }

	void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
		       double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

	void integrate(Tensor *stress, Tensor &en, Tensor &enp,
		       double *staten, double *statenp, double temp, Tensor *cache, double dt=0) const override;

	void initStates(double *) override;

	double getDensity() override { return rhoF; }

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

	void setSDProps(MFTTData *_ysst) override { if(sigP < 0 && _ysst && _ysst->getID() == -int(sigP)) ysst = _ysst; }

	void setSRDProps(MFTTData *_yssrt) override {  
		fprintf(stderr, "***Error: CrushableFoam::setSRDProps() not implemented.\n");
		exit(-1);
    }
};

#endif
