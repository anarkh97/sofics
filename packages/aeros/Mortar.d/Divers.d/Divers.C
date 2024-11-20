// --------------------------------------------------------
// HB - 08/01/03
// --------------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>

void Order2TriangleQuadRule(int, double&, double&, double&, double&);
void Order4TriangleQuadRule(int, double&, double&, double&, double&);
void Order5TriangleQuadRule(int, double&, double&, double&, double&);
void Order6TriangleQuadRule(int, double&, double&, double&, double&);
void Order8TriangleQuadRule(int, double&, double&, double&, double&);

// Compute 3D triangle area
double
Compute3DTriArea(double **NodesCoord, int DimFlag=0)
{
   double *M1 = NodesCoord[0];
   double *M2 = NodesCoord[1];
   double *M3 = NodesCoord[2];

   double V12[3], V13[3];
   V12[0] = M2[0]-M1[0]; V12[1] = M2[1]-M1[1]; V12[2] = M2[2]-M1[2];
   V13[0] = M3[0]-M1[0]; V13[1] = M3[1]-M1[1]; V13[2] = M3[2]-M1[2];
   
   double W[3];
   W[0] = V12[1]*V13[2] - V12[2]*V13[1];
   W[1] = V12[2]*V13[0] - V12[0]*V13[2];
   W[2] = V12[0]*V13[1] - V12[1]*V13[0];

   return(0.5*sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]));
}

// Gives integration pt location & weight on the reference triangle
void
getGaussPtOnTriangle(int ngp, int igp, double& r, double& s, double& t, double& weight)
{
  // Note: the sum of the weights is 0.5 for these rules. This implies that for example in the case
  // of a three-noded triangle value of "J" used in the integration should be equal to two times the
  // Area. This is not necessarily a standard practice, it is more conventional to have the sum of 
  // the weights equal 1.
  switch(ngp){
     case  1: r = 1./3.; s = 1./3.; t = 1.-r-s; weight = 0.5; break;
     case  3: Order2TriangleQuadRule(igp, r, s, t, weight)  ; break;
     case  6: Order4TriangleQuadRule(igp, r, s, t, weight)  ; break; 
     case  7: Order5TriangleQuadRule(igp, r, s, t, weight)  ; break; 
     case 12: Order6TriangleQuadRule(igp, r, s, t, weight)  ; break; 
     case 16: Order8TriangleQuadRule(igp, r, s, t, weight)  ; break; 
     default: {std::cerr <<"In getGaussPtOnTriangle: selected "<<ngp
                    <<" Gauss pts quadrature rule is NOT implemented !!!"<< std::endl; 
               break;
              }
   }
}

void
Order2TriangleQuadRule(int igp, double& r, double& s, double& t, double& weight)
{
  double c1 = 1./6.;
  double c2 = 2./3.;
  double w0 = 1./6.;

  switch(igp){
    case 1: r = c1; s = c1; t = 1.-r-s; weight = w0; break;
    case 2: r = c2; s = c1; t = 1.-r-s; weight = w0; break;
    case 3: r = c1; s = c2; t = 1.-r-s; weight = w0; break;
    default: break;
  }
}
void
Order4TriangleQuadRule(int igp, double& r, double& s, double& t, double& weight)
{
  double a1 = 0.445948490915965;
  double b1 = 0.091576213509771;
  double a2 = 1-2*a1;
  double b2 = 1-2*b1;
  double w1 = 0.111690794839005;
  double w2 = 0.054975871827661;
  
  switch(igp){
    case 1: r = a1; s = a1; t = 1.-r-s; weight = w1; break;
    case 2: r = a2; s = a1; t = 1.-r-s; weight = w1; break;
    case 3: r = a1; s = a2; t = 1.-r-s; weight = w1; break;
    case 4: r = b1; s = b1; t = 1.-r-s; weight = w2; break;
    case 5: r = b2; s = b1; t = 1.-r-s; weight = w2; break;
    case 6: r = b1; s = b2; t = 1.-r-s; weight = w2; break;
    default: break;
  }
}

void
Order5TriangleQuadRule(int igp, double& r, double& s, double& t, double& weight)
{
  double a1 = (6.+sqrt(15.))/21.;
  double a2 = 1-2*a1;
  double b1 = 4./7. - a1;
  double b2 = 1-2*b1;
  double w1 = (155.+sqrt(15.))/2400.;
  double w2 = 31./240. - w1;
  
  switch(igp){
    case 1: r = 1./3.; s = 1./3.; t = 1.-r-s; weight = 9./80.; break;
    case 2: r = a1; s = a1; t = 1.-r-s; weight = w1; break;
    case 3: r = a2; s = a1; t = 1.-r-s; weight = w1; break;
    case 4: r = a1; s = a2; t = 1.-r-s; weight = w1; break;
    case 5: r = b1; s = b1; t = 1.-r-s; weight = w2; break;
    case 6: r = b2; s = b1; t = 1.-r-s; weight = w2; break;
    case 7: r = b1; s = b2; t = 1.-r-s; weight = w2; break;
    default: break;
  }
}

void
Order6TriangleQuadRule(int igp, double& r, double& s, double& t, double& weight)
{

  double w1 = 0.116786275726379;
  double l1 = 0.501426509658179;
  double l2 = 0.249286745170910;
  double l3 = 1. - l1 - l2;

  double w2 = 0.050844906370207;
  double h1 = 0.873821971016996;
  double h2 = 0.063089014491502;
  double h3 = 1. - h1 - h2;

  double w3 = 0.082851075618374;
  double k1 = 0.053145049844817;
  double k2 = 0.310352451033784;
  double k3 = 1. - k1 - k2;

  // divide weight by 1/2 
  // (unit triangle area = 1/2 and NOT 1)
  w1 *= 0.5; w2 *= 0.5; w3 *= 0.5; 
  
  switch(igp){
    case  1: r = l1; s = l2; t = l3; weight = w1; break;
    case  2: r = l2; s = l3; t = l1; weight = w1; break;
    case  3: r = l3; s = l1; t = l2; weight = w1; break;
    case  4: r = h1; s = h2; t = h3; weight = w2; break;
    case  5: r = h2; s = h3; t = h1; weight = w2; break;
    case  6: r = h3; s = h1; t = h2; weight = w2; break;
    case  7: r = k1; s = k2; t = k3; weight = w3; break;
    case  8: r = k1; s = k3; t = k2; weight = w3; break;
    case  9: r = k2; s = k1; t = k3; weight = w3; break;
    case 10: r = k3; s = k1; t = k2; weight = w3; break;
    case 11: r = k2; s = k3; t = k1; weight = w3; break;
    case 12: r = k3; s = k2; t = k1; weight = w3; break;
    default: break;
  }
}

void
Order8TriangleQuadRule(int igp, double& r, double& s, double& t, double& weight)
{

  double w1 = 0.144315607677787;
  double l1 = 1./3.;
  double l2 = l1 ;
  double l3 = 1. - l1 - l2;

  double w2 = 0.095091634267285;
  double h1 = 0.081414823414554;
  double h2 = 0.459292588292723;
  double h3 = 1. - h1 - h2;

  double w3 = 0.103217370534718;
  double k1 = 0.658861384496480;
  double k2 = 0.170569307751760;
  double k3 = 1. - k1 - k2;

  double w4 = 0.032458497623198;
  double t1 = 0.898905543365938;
  double t2 = 0.050547228317031;
  double t3 = 1. - t1 - t2;

  double w5 = 0.027230314174435;
  double d1 = 0.008394777409958;
  double d2 = 0.263112829634638;
  double d3 = 1. - d1 - d2;

  // divide weight by 1/2 
  // (unit triangle area = 1/2 and NOT 1)
  w1 *= 0.5; w2 *= 0.5; w3 *= 0.5; 
  w4 *= 0.5; w5 *= 0.5;
  
  switch(igp) {
    case  1: r = l1; s = l2; t = l3; weight = w1; break;
    case  2: r = h1; s = h2; t = h3; weight = w2; break;
    case  3: r = h2; s = h3; t = h1; weight = w2; break;
    case  4: r = h3; s = h1; t = h2; weight = w2; break;
    case  5: r = k1; s = k2; t = k3; weight = w3; break;
    case  6: r = k2; s = k3; t = k1; weight = w3; break;
    case  7: r = k3; s = k1; t = k2; weight = w3; break;
    case  8: r = t1; s = t2; t = t3; weight = w4; break;
    case  9: r = t2; s = t3; t = t1; weight = w4; break;
    case 10: r = t3; s = t1; t = t2; weight = w4; break;
    case 11: r = d1; s = d2; t = d3; weight = w5; break;
    case 12: r = d1; s = d3; t = d2; weight = w5; break;
    case 13: r = d2; s = d1; t = d3; weight = w5; break;
    case 14: r = d3; s = d1; t = d2; weight = w5; break;
    case 15: r = d2; s = d3; t = d1; weight = w5; break;
    case 16: r = d3; s = d2; t = d1; weight = w5; break;
    default: break;
  }
}

