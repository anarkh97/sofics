#ifndef _DEMLE2D_H_
#define _DEMLE2D_H_

#include <cmath>
#include <Element.d/DEM.d/DEMElement.h>

class DGMLE2d_LM: public DEMLM {
public:
	virtual void init() {};
	virtual void ldir(int,double[2],complex<double>*)=0;
};


class DGMLE2d_1_LM: public DGMLE2d_LM {
public:
	int nDofs() const override { return 2; }
	int type() const override { return 201; }
	void ldir(int,double[2],complex<double>*) override;
};

class DGMLE2d_4_LM: public DGMLE2d_LM {
public:
	int nDofs() const override { return 4; }
	int type() const override { return 202; }
	void ldir(int,double[2],complex<double>*) override;
};

class DGMLE2d: public DEMElement {
public:
	int ndir;
	virtual void dir(int,complex<double>*) {};
	int o;
	DGMLE2d(int _o, int* nodenums);
	int defaultLMType() const override { return 201; }
	bool dgmFlag() const override { return true; }
	bool condensedFlag() const override  {
		for(int i=0;i<nFaces();i++) if (bc[i]==3) return false;
		return true;
//   return false;
	}

	bool storeMatricesFlag() override { return true; }
	int nPolynomialDofs() const override { return 0; }
	int nEnrichmentDofs() const override { return ndir; }
	int nGeomNodes() const override { return o*o; }
	int nFaces() const override { return 4; }
	int nFaceCorners(int fi) const override { return 2; }
	int *faceCorners(int fi) override {
		int *fc = new int[2];
		if (fi==1) { fc[0] = nn[0]; fc[1] = nn[o-1]; }
		else if (fi==2) { fc[0] = nn[o-1]; fc[1] = nn[o*o-1]; }
		else if (fi==3) { fc[0] = nn[o*o-1]; fc[1] = nn[o*(o-1)]; }
		else if (fi==4) { fc[0] = nn[o*(o-1)]; fc[1] = nn[0]; }
		return fc;
	}
	int polyDofType() const override { return DofSet::Xdisp|DofSet::Ydisp; }
	int polyDofsPerNode() const override { return 2; }

	virtual void getRef(double *xyz,double *xy);
	void createM(complex<double>*)  override;
	void createRHS(complex<double>*) override;
	void createSol(double *xyz, complex<double>*,
	               complex<double>*) override;

	void enrichmentF(double *x, complex<double> *f) override;
	void polynomialF(double *x, double *f) override;
};


class DGMLE2d_4: public DGMLE2d {
public:
	DGMLE2d_4(int _o, int* nodenums);

	int getElementType() const override { return 1200; }
	void dir(int,complex<double>*) override;
	int defaultLMType() const override { return 201; }
};

class DGMLE2d_16: public DGMLE2d {
public:
	DGMLE2d_16(int _o, int* nodenums);

	int getElementType() const override { return 1201; }
	void dir(int,complex<double>*) override;
	int defaultLMType() const override { return 202; }
};


class DEMLE2d: public DGMLE2d {
public:
	DEMLE2d(int _o, int* nodenums);
	bool dgmFlag() const override { return false; }
	int nPolynomialDofs() const  override { return 2*o*o; }

	void createM(complex<double>*)  override;
	void createRHS(complex<double>*) override;
};

class DEMLE2d_4: public DEMLE2d {
public:
	DEMLE2d_4(int _o, int* nodenums);

	int getElementType() const override { return 1220; }
	void dir(int,complex<double>*) override;
	int defaultLMType() const override { return 201; }
};

#endif

