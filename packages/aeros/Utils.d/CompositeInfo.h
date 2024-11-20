#ifndef _COMPOSITE_INFO_H_
#define _COMPOSITE_INFO_H_

#include <Utils.d/resize_array.h>

// Coefficient Data class for composites
class CoefData {
    public:
      // the first 6 rows of c store the constitutive matrix
      // the 7th row of c stores the coefficients of thermal expansion
      double c[7][6];
      bool coefFlag; // true:  the the constitutive matrix relates stress and strain, material properties are assume constant through the thickness
                     // false: the constitutive matrix relates forces and moments to mid-surface strains and curvatures
    public:
      void zero();
      void setCoef(int,int,double);
      void setCoef(double d[36]);
      double *values() { return (double *)c; }
};


// composite layer information
class LayInfo {
   public:
     int numLayers;
     int maxLayer;
     int type;
   private:
     ResizeArray<int> matids;
   public:
     double (*data)[12];
     double (*grad)[12];

     LayInfo(int _type);
 
     void add(int k, double *d, int m = -1);
     void setAllLayers(int, double (*)[12]);
     void setGrad();
     void zeroGrad();

     int nLayers() { return numLayers; }
     int getType() { return type; }

     double *values()  { return (double *) data; }
     double *gradval() { return (double *) grad; }

     void setLayerMaterialProperties(int k, double *d);
     int getLayerMaterialId(int k);
};

// composite layer material information
class LayMat {
  public:
    int m;
    double data[10]; // data[] = { E1  E2  nu12  G12  mu1,12 mu2,12 rho, cte1, cte2, ta }
    LayMat(int _m, double *d) { m = _m; for(int i=0; i<10; ++i) data[i] = d[i]; }
};

#endif
