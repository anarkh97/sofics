#ifndef _EXPMAT_H_
#define _EXPMAT_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>

#include <iostream>
#include <iterator>

class MFTTData;

class ExpMat : public NLMaterial
{
public:
	int optctv; // constitutive law (1 for hypoelastic, 3 for elasto viscoplastic, 5 for j2 elasto plastic)
	double ematpro[20]; // material properties: Young's modulus, Poisson's ratio, mass density, etc.
	int optcor0; // warping correction (0: off, 1: on), default = 1
	// reference: Belytschko, Wong and Chiang, CMAME, 1992, vol. 96, pp. 93-107
	//            "Advances in one point quadrature shell elements" (see eq. 30)
	int optcor1; // transverse shear projection (0: off, 1: on), default = 0
	// reference: Belytschko, Wong and Chiang, CMAME, 1992, vol. 96, pp. 93-107
	//            "Advances in one point quadrature shell elements" (see eq. 35)
	int optprj;  // rotation projection (0: off, 1: drilling only, 2: coupled drilling and rigid body modes), default = 1
	// reference: Belytschko and Leviathan, CMAME, 1994, vol. 115, pp. 277-286
	//            "Projection schemes for one-point quadrature shell elements"
	int opthgc;  // hourglass control (perturbation type) (0: off, 1: on), default = 1
	// reference: Chiang, M.S. thesis, Northwestern Univ., 1992
	double prmhgc[3]; // hourglass control parameters, default = [ 2.5e-3, 2.5e-3, 2.5e-3 ]
	int ngqpt2;  // gq rule for through thickness (1-12,14,16,24), default = 3
	MFTTData *ysst; // yield stress vs. effective plastic strain
	int yssrtid;
	MFTTData *yssrt; // yield stress scaling factor vs. effective plastic strain rate

	ExpMat(int _optctv, double e1, double e2, double e3, double e4=0, double e5=0, double e6=0,
		   double e7=0, double e8=0, double e9=0, double e10=0, double e11=0, double e12=0, double e13=0,
		   double e14=0, double e15=0, double e16=0, double e17=0, double e18=0, double e19=0.833)
	{ optctv = _optctv;
	  ematpro[0] = e1; ematpro[1] = e2; ematpro[2] = e3; ematpro[3] = e4; ematpro[4] = e5;
	  ematpro[5] = e6; ematpro[6] = e7; ematpro[7] = e8; ematpro[8] = e9; ematpro[9] = e10;
	  ematpro[10] = e11; ematpro[11] = e12; ematpro[12] = e13; ematpro[13] = e14;
	  ematpro[14] = e15; ematpro[15] = e16; ematpro[16] = e17; ematpro[17] = e18;
	  ematpro[18] = e19; // shear correction factor
	  optcor0 = 1;
	  optcor1 = 0;
	  optprj = 1;
	  opthgc = 1;
	  prmhgc[0] = prmhgc[1] = prmhgc[2] = 2.5e-3;
	  ngqpt2 = 3;
	  ysst = NULL;
	  yssrtid = 0;
	  yssrt = NULL;
	}

	int getNumStates() const override { return 0; }

	void getStress(Tensor *stress, Tensor &strain, double*, double) override
	{ std::cerr << "ExpMat::getStress is not implemented\n"; }

	void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double) override
	{ std::cerr << "ExpMat::getTangentMaterial is not implemented\n"; }

	void getElasticity(Tensor *tm) const override
	{ std::cerr << "ExpMat::getElasticity is not implemented\n"; }

	void updateStates(Tensor &en, Tensor &enp, double *state, double) override {}

	void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double) override
	{ std::cerr << "ExpMat::getStressAndTangentMaterial is not implemented\n"; }

	void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
				   double *staten, double *statenp, double, Tensor *cache, double=0) const override
	{ std::cerr << "ExpMat::integrate is not implemented\n"; }

	void integrate(Tensor *stress, Tensor &en, Tensor &enp,
				   double *staten, double *statenp, double, Tensor *cache, double=0) const override
	{ std::cerr << "ExpMat::integrate is not implemented\n"; }

	void initStates(double *) override {}

	double getDensity() override { return ematpro[2]; }

	StrainEvaluator * getStrainEvaluator() const override
	{ std::cerr << "ExpMat::getStrainEvaluator is not implemented\n"; return NULL; }

	void print(std::ostream &out) const override {
	  std::string type;
	  switch (optctv) {
		case 1: type = "HypoElastic"; break;
		case 5: type = "J2Plasticity"; break;
		case 6: type = "KK1"; break;
		case 7: type = "KK2"; break;
		default: throw std::range_error("Unknown material law type");
	  }
	  out << type << " ";
	  std::copy(&ematpro[0], &ematpro[8], std::ostream_iterator<double>(out, " "));
	  out << yssrtid << " " << optcor0 << " " << optcor1 << " " << optprj << " " << opthgc << " "
		  << prmhgc[0] << " " << prmhgc[1] << " " << prmhgc[2] << " " << ngqpt2 << " "
		  << ematpro[18];
	}

	void setSDProps(MFTTData *_ysst) override { if(optctv == 5 && _ysst && _ysst->getID() == -int(ematpro[3])) ysst = _ysst; }
	void setSRDProps(MFTTData *_yssrt) override { if(optctv == 5 && yssrtid > 0 &&  _yssrt && _yssrt->getID() == yssrtid) yssrt = _yssrt; }
};

#endif
