#include <Element.d/Helm.d/GaussRules.h>

#define GSUBDIV 1

#include <Utils.d/dbg_alloca.h>

template<class Scalar> void IsoParamUtils::zeroOut(int n,Scalar *v) {
 int i;
 for(i=0;i<n;i++) v[i] = 0;
}


template<class Scalar> void IsoParamUtils::symmetrize(int n, Scalar *K) {
   int i,j;
   for(j=0;j<n;j++)
     for(i=0;i<j;i++)
       K[j*n+i] = K[i*n+j];
}




template<class Scalar> class GalStiffFunction : public IntegFunctionV3d {
 int n;
 Scalar *K;
public:
 GalStiffFunction(int _n, Scalar*_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {
   int i,j;
   double wdet = w*det;
   for(j=0;j<ncolumns();j++)
     for(i=j;i<nrows();i++)
       K[j*nrows()+i] += wdet*(
        dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1]+dNdx[i][2]*dNdx[j][2]);
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};


template<class Scalar> class GalMassFunction : public IntegFunctionV3d {
 int n;
 Scalar *K;
public:
 GalMassFunction(int _n, Scalar*_K) {
   n = _n; K = _K;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {
   int i,j;
   double wdet = w*det;
   for(j=0;j<ncolumns();j++)
     for(i=j;i<nrows();i++)
       K[j*nrows()+i] += wdet*N[i]*N[j];
 }
 int nrows() { return n; }
 int ncolumns() { return n; }
};
