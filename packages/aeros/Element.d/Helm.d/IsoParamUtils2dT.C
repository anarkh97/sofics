#include <Element.d/Helm.d/GaussRules.h>
#include <Element.d/Helm.d/IntegFunction.h>

#define GSUBDIV 2

#ifdef WINDOWS
#include <alloca.h>
#endif

template<class Scalar> void IsoParamUtils2d::zeroOut(int n,Scalar *v) {
 int i;
 for(i=0;i<n;i++) v[i] = 0;
}


template<class Scalar> void IsoParamUtils2d::symmetrize(int n, Scalar *K) {
   int i,j;
   for(j=0;j<n;j++)
     for(i=0;i<j;i++)
       K[j*n+i] = K[i*n+j];
}

