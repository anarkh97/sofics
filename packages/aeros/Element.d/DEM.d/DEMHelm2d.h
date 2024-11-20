#ifndef _DEMHELM2D_H_
#define _DEMHELM2D_H_

#include <cmath>
#include <Element.d/DEM.d/DEMElement.h>

class DGMHelm2d_LM: public DEMLM {
public:
	void init() override {};
	virtual complex<double> coef(int)=0;
};


class DGMHelm2d_1_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 1; }
	int type() const override { return 1; }
	complex<double> coef(int) override;
};

class DGMHelm2d_2_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 2; }
	int type() const override { return 2; }
	complex<double> coef(int) override;
};

class DGMHelm2d_4_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 4; }
	int type() const override { return 3; }
	complex<double> coef(int) override;
};

class DGMHelm2d_8_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 8; }
	int type() const override { return 4; }
	complex<double> coef(int) override;
};

class DGMHelm2d_Eva_LM: public DEMLM {
public:
	virtual complex<double> ldir(int,double[2])=0;
};

class DGMHelm2d_Eva1_2_LM: public DGMHelm2d_Eva_LM {
	int ndofs;
	double gamma;
	double normal[2];
public:
	DGMHelm2d_Eva1_2_LM() {
		// needs to be initialized
		ndofs = 0;
	}
	void init() override;
	int nDofs() const override { return ndofs; }
	int type() const override { return 5; }
	complex<double> ldir(int,double[2]) override;
};


class DGMHelm2d: public DEMElement {
public:
	int ndir;
	virtual complex<double> dir(int);
	int o;
	DGMHelm2d(int _o, int* nodenums);
	int defaultLMType() const override { return 1; }
	bool dgmFlag() const override { return true; }
	bool condensedFlag() const override  {
		for(int i=0;i<nFaces();i++) if (bc[i]==3) return false;
// return false;
		return true;
	}
// bool condensedFlag() const override  { return false; }
	bool storeMatricesFlag() override { return true; }
// bool storeMatricesFlag() override { return false; }
	int nPolynomialDofs() const override { return 0; }
	int nEnrichmentDofs() const override { return ndir; }
	int nGeomNodes() const override { return (o>0)?o*o:((-o)*(-o+1))/2; }
	int nFaces() const override { return (o>0)?4:3; }
	int nFaceCorners(int fi) const override { return 2; }
	int *faceCorners(int fi) override {
		int *fc = new int[2];
		if (o>0) {
			if (fi==1) { fc[0] = nn[0]; fc[1] = nn[o-1]; }
			else if (fi==2) { fc[0] = nn[o-1]; fc[1] = nn[o*o-1]; }
			else if (fi==3) { fc[0] = nn[o*o-1]; fc[1] = nn[o*(o-1)]; }
			else if (fi==4) { fc[0] = nn[o*(o-1)]; fc[1] = nn[0]; }
			return fc;
		} else {
			if (fi==3) { fc[0] = nn[0]; fc[1] = nn[-o-1]; }
			else if (fi==1) { fc[0] = nn[-o-1]; fc[1] = nn[((-o)*(-o+1))/2-1]; }
			else if (fi==2) { fc[0] = nn[((-o)*(-o+1))/2-1]; fc[1] = nn[0]; }
			return fc;
		}
	}
	int polyDofType() const override { return DofSet::Helm; }
	int polyDofsPerNode() const override { return 1; }

	virtual void getRef(double *xyz,double *xy);
	void createM(complex<double>*) override;
	void interfMatrix(int fi, DEMElement*,complex<double>*) override;
	void createRHS(complex<double>*) override;
	void createSol(double *xyz, complex<double>*,
				   complex<double>*) override;
};


class DGMHelm2d_4: public DGMHelm2d {
public:
	DGMHelm2d_4(int _o, int* nodenums);

	int getElementType() const override { return 1100; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 1; }
};

class DGMHelm2d_4t: public DGMHelm2d {
public:
	DGMHelm2d_4t(int _o, int* nodenums);

	int getElementType() const override { return 1110; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 1; }
};

class DGMHelm2d_8: public DGMHelm2d {
public:
	DGMHelm2d_8(int _o, int* nodenums);

	int getElementType() const override { return 1101; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 2; }
};

class DGMHelm2d_8t: public DGMHelm2d {
public:
	DGMHelm2d_8t(int _o, int* nodenums);

	int getElementType() const override { return 1111; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 2; }
};

class DGMHelm2d_16: public DGMHelm2d {
public:
	DGMHelm2d_16(int _o, int* nodenums);

	int getElementType() const override { return 1102; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 3; }
};

class DGMHelm2d_32: public DGMHelm2d {
public:
	DGMHelm2d_32(int _o, int* nodenums);

	int getElementType() const override { return 1103; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 4; }
};


class DGMHelm2d_Eva: public DGMHelm2d {
public:
	double normal[2];
	double gamma;
	DGMHelm2d_Eva(int _o, int* nodenums);
};

class DGMHelm2d_Eva2_8: public DGMHelm2d_Eva {
public:
	DGMHelm2d_Eva2_8(int _o, int* nodenums);

	int getElementType() const override { return 1104; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 5; }
};


class DEMHelm2d: public DGMHelm2d {
public:
	DEMHelm2d(int _o, int* nodenums);
	bool dgmFlag() const override { return false; }
	int nPolynomialDofs() const  override { return (o>0)?o*o:((-o)*(-o+1))/2; }

	void createM(complex<double>*) override;
	void interfMatrix(int fi, DEMElement*,complex<double>*) override;
	void createRHS(complex<double>*) override;
};


class DEMHelm2d_4: public DEMHelm2d {
public:
	DEMHelm2d_4(int _o, int* nodenums);

	int getElementType() const override { return 1120; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 1; }
};

class DEMHelm2d_4t: public DEMHelm2d {
public:
	DEMHelm2d_4t(int _o, int* nodenums);

	int getElementType() const override { return 1130; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 1; }
};

class DEMHelm2d_8: public DEMHelm2d {
public:
	DEMHelm2d_8(int _o, int* nodenums);

	int getElementType() const override { return 1121; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 2; }
};

class DEMHelm2d_8t: public DEMHelm2d {
public:
	DEMHelm2d_8t(int _o, int* nodenums);

	int getElementType() const override { return 1131; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 2; }
};

class DEMHelm2d_16: public DEMHelm2d {
public:
	DEMHelm2d_16(int _o, int* nodenums);

	int getElementType() const override { return 1122; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 3; }
};

class DEMHelm2d_32: public DEMHelm2d {
public:
	DEMHelm2d_32(int _o, int* nodenums);

	int getElementType() const override { return 1123; }
	complex<double> dir(int) override;
	int defaultLMType() const override { return 4; }
};

#endif

