#ifndef _GAUSSRULE_H_
#define _GAUSSRULE_H_

class GaussRule {

public:
        int order;
        int ngauss;
        double *pool;
        double *xigauss;
        double *wgauss;

        GaussRule(int);
        void init();
};
#endif
