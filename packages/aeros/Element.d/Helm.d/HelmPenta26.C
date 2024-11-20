// ---------------------------------------------------------------------
// HB - 05-24-05
// ---------------------------------------------------------------------
// 26 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
// EXPERIMENTAL ...
// ---------------------------------------------------------------------
#include <cstdio>
#include <Element.d/Helm.d/HelmPenta26.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

#define CHECK_JACOBIAN //HB: force check nullity & constant sign of jacobian over el.
//#define HELMPENTA26_DEBUG
extern bool useFull;
extern "C"      {
void    _FORTRAN(lgauss)(int &, int &, double *, double *);
}

extern void   Penta26ShapeFct(double Shape[26], double dShape[26][3], double m[3]);
extern double Penta26ShapeFct(double Shape[26], double DShape[26][3], double m[3], double X[26], double Y[26], double Z[26]);
extern double computePenta26DShapeFct(double dShape[26][3], double X[26], double Y[26], double Z[26], double (*DShape)[3] = 0);
extern void addBtBtoK3DHelm(FullSquareMatrix &K, double (*DShape)[3], double alpha, int nnodes, int* ls);
extern void addNtNtoM3DHelm(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls);

#ifdef CHECK_JACOBIAN //HB: for checking zero/small & constant sign of jacobian over the el.
extern int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);
extern void printShapeFct3D(double* Shape, double (*DShape)[3], int nnodes, char* mssg= 0, FILE* file=stderr);
#endif
                                                                                                                                        

HelmPenta26::HelmPenta26(int* nodenums)
{
  for(int i=0; i<26; i++)
    nn[i] = nodenums[i];
}

Element *
HelmPenta26::clone()
{
  return(new HelmPenta26(*this));
}

void
HelmPenta26::renum(const int *table)
{
  for(int i=0; i<26; i++) { nn[i] = table[nn[i]]; }
}

void
HelmPenta26::renum(EleRenumMap& table)
{
  for(int i=0; i<26; i++) { nn[i] = table[nn[i]]; }
}

double
HelmPenta26::getMass(const CoordSet& cs) const
{
  fprintf(stderr," *** WARNING: HelmPenta26::getMass NOT implemented. Return zero mass.\n");
  return(0.0);
}

FullSquareMatrix
HelmPenta26::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
  int status = 0;
  const int nnodes  = 26;
  const int ndofs   = 26;
  //const int numgauss= 4;

  double X[26], Y[26], Z[26];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[26] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};

  FullSquareMatrix M(ndofs,d);
  M.zero();
 
  // hard coded order 8 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double w1 = 0.116786275726379 ;
  double l1 = 0.501426509658179 ;
  double l2 = 0.249286745170910 ;
  double l3 = 1. - l1 - l2;

  double w2 = 0.050844906370207 ;
  double h1 = 0.873821971016996 ;
  double h2 = 0.063089014491502 ;
  double h3 = 1. - h1 - h2;
 
  double w3 = 0.082851075618374 ;
  double k1 = 0.053145049847817 ;
  double k2 = 0.310352451033784 ;
  double k3 = 1. - k1 - k2;
    
  // divide weight by 1/2
  // (unit triangle area = 1/2 and NOT 1)
  w1 *= 0.5; w2 *= 0.5; w3 *= 0.5;

  double TriGPt12[12][4]= {{l1, l2, l3, w1},
                           {l2, l3, l1, w1},
                           {l3, l1, l2, w1},
                           {h1, h2, h3, w2},
                           {h2, h3, h1, w2}, 
                           {h3, h1, h2, w2},
                           {k1, k2, k3, w3},
                           {k1, k3, k2, w3},
                           {k2, k1, k3, w3},
                           {k3, k1, k2, w3},
                           {k2, k3, k1, w3},
                           {k3, k2, k1, w3}};

  // integration: loop over Gauss pts
  double wxy,wz,w;
  int ngpz =  4; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 12; // numbder of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[26], DShape[26][3];
  double dOmega;//det of jacobian
  int jSign = 0;
  for(int iz=1;iz<=ngpz;iz++){ // routine lgauss uses fortran indexing
    _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz); // get z position & weight of the Gauss pt
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt12[ixy][0]; m[1] = TriGPt12[ixy][1]; wxy = TriGPt12[ixy][3];
      dOmega = Penta26ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef HELMPENTA26_DEBUG
      fprintf(stderr," *** In HelmPenta26::massMatrix: iz = %d, ixy = %d, J = %e\n",iz,ixy+1,dOmega);
#endif
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "HelmPenta26::massMatrix");
#endif
      w = dOmega*wxy*wz;
      addNtNtoM3DHelm(M, Shape, w, nnodes, ls);
    }
  }

  if(status != 0) {
    fprintf(stderr, " *** FATAL ERROR in HelmPenta26::stiffness");
    fprintf(stderr, " -> Penta26's nodes (+1): \n");
    for(int i=0; i<26; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
    exit(-1);
  }
  
  M /= getProperty()->rho;
  return(M);
}

FullSquareMatrix
HelmPenta26::stiffness(const CoordSet &cs, double *d, int flg) const
{
  int status = 0;
  const int nnodes  = 26;
  const int ndofs   = 26;
  //const int numgauss= 4;

  double X[26], Y[26], Z[26];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[26] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
  
  FullSquareMatrix K(ndofs,d);
  K.zero();

  // hard coded order 4 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double a1 = 0.445978490915965;
  double b1 = 0.091576213509771;
  double a2 = 1-2.*a1;
  double b2 = 1-2.*b1;
  double w1 = 0.111690797839005;
  double w2 = 0.054975871827661;
  double TriGPt6[6][4] = {{a1, a1, 1.-a1-a1, w1},
                          {a2, a1, 1.-a2-a1, w1},
                          {a1, a2, 1.-a1-a2, w1},
                          {b1, b1, 1.-b1-b1, w2},
                          {b2, b1, 1.-b2-b1, w2},
                          {b1, b2, 1.-b1-b2, w2}};

  // integration: loop over Gauss pts
  double wxy,wz,w;
  int ngpz = 3; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 6; // numbder of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[26], DShape[26][3];
  double dOmega;//det of jacobian
  int jSign = 0;
#ifdef HELMPENTA26_DEBUG
  double Vol = 0.0;
#endif
  for(int iz=1;iz<=ngpz;iz++){ // routine lgauss uses fortran indexing
    _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz); // get z position & weight of the Gauss pt
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt6[ixy][0]; m[1] = TriGPt6[ixy][1]; wxy = TriGPt6[ixy][3];
      dOmega = Penta26ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef HELMPENTA26_DEBUG
      fprintf(stderr," *** In HelmPenta26::stiffness: iz = %d, ixy = %d, J = %e\n",iz,ixy+1,dOmega);
      Vol += dOmega*wxy*wz;
#endif
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "HelmPenta::stiffness");
#endif
      w = dOmega*wxy*wz;
      addBtBtoK3DHelm(K, DShape, w, nnodes, ls);
    }
  }
#ifdef HELMPENTA26_DEBUG
  fprintf(stderr," *** In HelmPenta26::stiffness: volume = %e\n",Vol);
#endif

  if(status != 0) {
    fprintf(stderr, " *** FATAL ERROR in HelmPenta26::stiffness");
    fprintf(stderr, " -> Pent26's nodes (+1): \n");
    for(int i = 0; i < 6; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
    exit(-1);
  }
  K /= getProperty()->rho;
  return(K);
}
 

FullSquareMatrix
HelmPenta26::acousticm(CoordSet &cs, double *d)
{
  double dK[676], dM[676];

  FullSquareMatrix K = stiffness(cs, dK);
  FullSquareMatrix M = massMatrix(cs, dM, 1);

  double kappa = prop ->kappaHelm;
  double kk = kappa*kappa;
  
  for(int i=0;i<676;i++)
    d[i] = dK[i] - kk*dM[i];

  /* take advantage of the symmetry to save some computation: NEED TESTING ...
  for(int i=0;i<26;i++)
    for(int j=i;j<26;j++){
      d[26*i+j] = dK[26*i+j] - kk*dM[26*i+j];
      d[i+26*j] = d[26*i+j];
    }
  */
  FullSquareMatrix Ka(26,d);

  Ka /= getProperty()->rho;
  return(Ka);
}

int
HelmPenta26::numNodes() const {
  if(useFull)
    return(26); 
  else
    return(8);
}
int HelmPenta26::numDofs() const { return(26); }

int*
HelmPenta26::nodes(int *p) const
{
  if(useFull)
    {
      if(!p) p = new int[26];
      for(int i=0; i<26; i+=2) {
	p[i] = nn[i]; p[i+1] = nn[i+1];
      }
      return(p);
    }
  else
    {
      if(!p) p = new int[6];
      for(int i=0; i<6; i+=1) {
	p[i] = nn[i]; 
      }
      return(p);
    }
}

int*
HelmPenta26::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[26];
  for(int i=0; i<26; i++)
    dsa.number(nn[i],DofSet::Helm, p+i);

  return(p);
}

void
HelmPenta26::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::Helm);
}

// Treat as Penta6 (MECH)
//int HelmPenta26::getTopNumber() const { return(124); }
//int HelmPenta26::numTopNodes() { return(6); }

int
HelmPenta26::getTopNumber() const { return(194); }
int
HelmPenta26::numTopNodes() const { return(26); }

void
HelmPenta26::addFaces(PolygonSet *pset)
{
  fprintf(stderr," *** ERROR: HelmPenta26::addFaces NOT implemented. Abort.\n");
  exit(-1);
}


