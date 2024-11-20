#ifndef _GAUSSRULES_H_
#define _GAUSSRULES_H_

class GaussRuleLine {
public:
        int nsubdiv;
        int ng;
        int ngauss;
        double *pool;
        double *xigauss;
        double *wgauss;

        GaussRuleLine(int,int);
        void init();
};


class GaussLobattoRuleLine {
public:
        int ngauss;
        double *xigauss;
        double *wgauss;

        GaussLobattoRuleLine(int);
        ~GaussLobattoRuleLine();
};


class GaussRuleTriLine: public GaussRuleLine {
public:
        GaussRuleTriLine(int,int);
        void init(int);
private:
        void createFaceInteg(int faceindex, double *_xigauss);
};



class GaussRuleTetra {
public:
        int ngauss;
        double *pool;
        double *xigauss;
        double *wgauss;

        GaussRuleTetra(int);
        void init();
private:
        void genFrom1D();
};


class GaussRuleTriangle {
public:
        int ngauss;
        double *pool;
        double *xigauss;
        double *wgauss;

        GaussRuleTriangle(int);
        void init(int faceindex=0);
private:
        void createFaceInteg(int faceindex, double *_xigauss);
        void genFrom1D(int faceindex);
};
#endif
