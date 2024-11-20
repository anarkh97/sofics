#include <cstdio>
#include <Utils.d/CompositeInfo.h>
#include <iostream>
using namespace std;

void CoefData::setCoef(double data[36])
{
 for(int i = 0; i < 6; i++)
   for(int j = 0; j < 6; j++)
     c[i][j] = data[6*i + j];
}

void
CoefData::zero()
{
 for(int i = 0; i < 6; i++)
   for(int j = 0; j < 6; j++)
     c[i][j] = 0.0;
}

void
CoefData::setCoef(int i, int j, double v)
{
 if(i < 0 || j < 0 || i >= 6 || j >= 6) {
   fprintf(stderr,"ERROR: in coefficient data:"
                  " out of bound indices %d %d\n",i+1,j+1);
 }
 c[i][j] = c[j][i] = v;
}

LayInfo::LayInfo(int _type) : matids(0)
{
 type      = _type;
 maxLayer  = 4;
 numLayers = 0;
 data      = new double[maxLayer][12];
 grad      = new double[maxLayer][12];
}

void
LayInfo::add(int ln, double *v, int m)
{ 
 int i,j;

 //allocate
 if(ln >= maxLayer) {
   int newMaxLayer = (ln+4);
   double (*newData)[12] = new double[newMaxLayer][12];
   for(i=0; i < numLayers; ++i) 
     for(j=0; j < 12; ++j)
       newData[i][j] = data[i][j];
   delete [] data;
   data     = newData;
   maxLayer = newMaxLayer;
 }
 if(ln >= numLayers) numLayers = ln+1;
 for(i=0; i < 12; ++i)
   data[ln][i] = v[i];

 matids[ln] = m;
}

//------------------------------------------------------------

void LayInfo::setAllLayers(int nLayers, double (*d)[12])  {

  delete [] data;
  maxLayer = nLayers;
  numLayers = nLayers;
  data = d;

}

//------------------------------------------------------------

void
LayInfo::setGrad()
{
  delete [] grad;

  grad = new double[maxLayer][12];
}

void
LayInfo::zeroGrad()
{
  int i,j;
  for(i=0; i < numLayers; ++i)
    for(j=0; j < 12; ++j)
      grad[i][j] = 0.0;
}

//------------------------------------------------------------
#include <cmath>
#define PI 3.14159265358979

void 
LayInfo::setLayerMaterialProperties(int k, double *d)
{
  double E1 = d[0], E2 = d[1], nu12 = d[2], G12 = d[3], mx = d[4], my = d[5], rho = d[6];
  double cte1 = d[7], cte2 = d[8], ta = d[9];
  data[k][0] = E1;
  data[k][1] = E2;
  data[k][2] = nu12;
  data[k][3] = G12;
  data[k][4] = mx; 
  data[k][5] = my; 
  data[k][6] = rho;
  data[k][9] = cte1;
  data[k][10] = cte2;
  data[k][11] = ta;
}

int 
LayInfo::getLayerMaterialId(int k)
{
  return matids[k];
}

