#ifndef _DEMHELM3D_H_
#define _DEMHELM3D_H_

#include <cmath>
#ifndef SANDIA
#include <Element.d/DEM.d/DEMElement.h>
#include <Element.d/Helm.d/IntegFunction.h>
#else
#include "DEMElement.h"
#endif

class DGMHelm3d_LM: public DEMLM {
public:
    virtual void init() {};
    virtual void ldir(int,double*,double*,complex<double>*)=0;
};


class DGMHelm3d_1_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 1; }
    int type() const override { return 51; }
    void ldir(int,double*,double*,complex<double>*) override;
};

class DGMHelm3d_4_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 4; }
    int type() const override { return 52; }
    void ldir(int,double*,double*,complex<double>*) override;
};

class DGMHelm3d_8_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 8; }
    int type() const override { return 53; }
    void ldir(int,double*,double*,complex<double>*) override;
};

class DGMHelm3d_12_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 12; }
    int type() const override { return 54; }
    void ldir(int,double*,double*,complex<double>*) override;
};




class DGMHelm3d: public DEMElement {
protected:
public:
    void HelmDGMEMatrices3d(double *xyz,
                            int ndir, complex<double> *dirs,
                            int *nldirs, complex<double> *ldirs,
                            double kappa, int *sflags, double *xsc,
                            double *xc,
                            complex<double> *K, complex<double> *L);
    void HelmDGMEMatricesExactFace3d(double *xyz,
                                     int ndir, complex<double> *dirs,
                                     int nldir, complex<double> *ldirs,
                                     double kappa,  int sflag,
                                     double *xsc,  double *xc, int faceindex,
                                     complex<double> *K, complex<double> *L);
    void HelmDGMPMLEMatrices3d(double *xyz,
                               int ndir, complex<double> *dirs,
                               int *nldirs, complex<double> *ldirs,
                               double kappa, int *sflags, double *xsc,
                               double *xc, PMLProps *pml,
                               complex<double> *K, complex<double> *L);
    void HelmDGMENeumV(double *xyz,
                       int ndir, complex<double> *dirs,
                       double kappa, complex<double> *incdir, int faceindex,
                       double *xc,
                       complex<double>* v);
    double * getCubeDir(int n);
    int ndir;
    virtual void dir(int,complex<double>*) {};
    int o;
    void init(int _nnodes, int* nodenums);
    int defaultLMType() const override { return 51; }
    bool dgmFlag() const override { return true; }
    bool condensedFlag() const override  {
        for(int i=0;i<nFaces();i++) if (bc[i]==3) return false&&condensedF;
        return true&&condensedF;
    }
    bool storeMatricesFlag() override { return storeMatricesF; }
// bool storeMatricesFlag() override { return false; }
    int nPolynomialDofs() const  override { return 0; }
    int nEnrichmentDofs() const override { return ndir; }
    virtual int *faceCornerI(int fi)=0;
    int *faceCorners(int fi) override {
        int *ix = faceCornerI(fi);
        int nf = nFaceCorners(fi);
        for(int i=0;i<nf;i++) {
            int tmp = ix[i]; ix[i] = nn[tmp];
        }
        return ix;
    }

    virtual void externalNormal(double *xyz, int faceindex, double *n)=0;
    virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f)=0;
    virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f)=0;
    virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f)=0;
    virtual int isFlatAndStraight(double *xyz,int faceindex)=0;

#ifndef SANDIA
    int polyDofType() const override { return DofSet::Helm; }
#endif
    int polyDofsPerNode() const override { return 1; }

    virtual void getRef(double *xyz,double *xy);
    void createM(complex<double>*)  override;
    void interfMatrix(int fi, DEMElement*,complex<double>*) override;
    void createRHS(complex<double>*) override;
    void createSol(double *xyz, complex<double>*, complex<double>*) override;
};

class HexDGMElement3d: public DGMHelm3d {
public:
    HexDGMElement3d(int _n, int* nodenums);
    int nGeomNodes() const override { return o*o*o; }
    int nFaces() const override { return 6; }
    int nFaceCorners(int fi) const override { return 4; }
    int *faceCornerI(int fi) override;
    void externalNormal(double *xyz, int faceindex, double *n) override;
    void volumeInt3d(double *xyz, IntegFunctionV3d &f) override;
    void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f) override;
    void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f) override;
    int isFlatAndStraight(double *xyz,int faceindex) override;
};


class TetraDGMElement3d: public DGMHelm3d {
public:
    TetraDGMElement3d(int _n, int* nodenums);
    int nGeomNodes() const override { return (o*(o+1)*(o+2))/6; }
    int nFaces() const override { return 4; }
    int nFaceCorners(int fi) const override { return 3; }
    int *faceCornerI(int fi) override;
    void externalNormal(double *xyz, int faceindex, double *n) override;
    void volumeInt3d(double *xyz, IntegFunctionV3d &f) override;
    void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f) override;
    void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f) override;
    int isFlatAndStraight(double *xyz,int faceindex) override;
};


class PrismDGMElement3d: public DGMHelm3d {
public:
    PrismDGMElement3d(int _n, int* nodenums);
    int nGeomNodes() const override { return o*(o*(o+1))/2; }
    int nFaces() const override { return 5; }
    int nFaceCorners(int fi) const override { return (fi<3)?3:4; }
    int *faceCornerI(int fi) override;
    void externalNormal(double *xyz, int faceindex, double *n) override;
    void volumeInt3d(double *xyz, IntegFunctionV3d &f) override;
    void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f) override;
    void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f) override;
    int isFlatAndStraight(double *xyz,int faceindex) override;
};


class PyramidDGMElement3d: public DGMHelm3d {
public:
    PyramidDGMElement3d(int _n, int* nodenums);
    int nGeomNodes() const override { return 5; }
    int nFaces() const override { return 5; }
    int nFaceCorners(int fi) const override { return (fi<2)?4:3; }
    int *faceCornerI(int fi) override;
    void externalNormal(double *xyz, int faceindex, double *n) override;
    void volumeInt3d(double *xyz, IntegFunctionV3d &f) override;
    void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f) override;
    void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f) override;
    int isFlatAndStraight(double *xyz,int faceindex) override;
};


class DGMHelm3d_6: public HexDGMElement3d {
public:
    DGMHelm3d_6(int _o, int* nodenums);

	int getElementType() const override { return 1150; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 51; }
};


class DGMHelm3d_6t: public TetraDGMElement3d {
public:
    DGMHelm3d_6t(int _o, int* nodenums);

	int getElementType() const override { return 1160; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 51; }
};


class DGMHelm3d_6p: public PrismDGMElement3d {
public:
    DGMHelm3d_6p(int _o, int* nodenums);
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 51; }
};

class DGMHelm3d_6pd: public PyramidDGMElement3d {
public:
    DGMHelm3d_6pd(int _o, int* nodenums);
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 51; }
};

class DGMHelm3d_26: public HexDGMElement3d {
public:
    DGMHelm3d_26(int _o, int* nodenums);

	int getElementType() const override { return 1151; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_26t: public TetraDGMElement3d {
public:
    DGMHelm3d_26t(int _o, int* nodenums);

	int getElementType() const override { return 1161; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_26p: public PrismDGMElement3d {
public:
    DGMHelm3d_26p(int _o, int* nodenums);
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_26pd: public PyramidDGMElement3d {
public:
    DGMHelm3d_26pd(int _o, int* nodenums);
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_56: public HexDGMElement3d {
public:
    DGMHelm3d_56(int _o, int* nodenums);

	int getElementType() const override { return 1152; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_56t: public TetraDGMElement3d {
public:
    DGMHelm3d_56t(int _o, int* nodenums);

	int getElementType() const override { return 1162; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_56p: public PrismDGMElement3d {
public:
    DGMHelm3d_56p(int _o, int* nodenums);
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_56pd: public PyramidDGMElement3d {
public:
    DGMHelm3d_56pd(int _o, int* nodenums);
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_98: public HexDGMElement3d {
public:
    DGMHelm3d_98(int _o, int* nodenums);

	int getElementType() const override { return 1153; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 54; }
};


class DEMHelm3d: public HexDGMElement3d {
public:
    DEMHelm3d(int _o, int* nodenums);
    int defaultLMType() const override { return 51; }
    bool dgmFlag() const override { return false; }
    int nPolynomialDofs() const  override { return o*o*o; }

    void createM(complex<double>*)  override;
    void interfMatrix(int fi, DEMElement*,complex<double>*) override;
    void createRHS(complex<double>*) override;
};


class DEMHelm3d_6: public DEMHelm3d {
public:
    DEMHelm3d_6(int _o, int* nodenums);

	int getElementType() const override { return 1170; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 51; }
};

class DEMHelm3d_26: public DEMHelm3d {
public:
    DEMHelm3d_26(int _o, int* nodenums);

	int getElementType() const override { return 1171; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 52; }
};

class DEMHelm3d_56: public DEMHelm3d {
public:
    DEMHelm3d_56(int _o, int* nodenums);

	int getElementType() const override { return 1172; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 53; }
};

class DEMHelm3d_98: public DEMHelm3d {
public:
    DEMHelm3d_98(int _o, int* nodenums);

	int getElementType() const override { return 1173; }
    void dir(int,complex<double>*) override;
    int defaultLMType() const override { return 54; }
};

#endif
