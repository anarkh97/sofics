#include <Element.d/Helm.d/HelmPenta.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

extern "C"      {
void    _FORTRAN(lgauss)(int &, int &, double *, double *);
}
extern double Penta6ShapeFct(double Shape[6], double DShape[6][3], double m[3], double X[6], double Y[6], double Z[6]);
extern void addBtBtoK3DHelm(FullSquareMatrix &K, double (*DShape)[3], double alpha, int nnodes, int* ls);
extern void addNtNtoM3DHelm(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls);

#define CHECK_JACOBIAN
#ifdef CHECK_JACOBIAN //HB: for checking zero/small & constant sign of jacobian over the el.
extern int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);
extern void printShapeFct3D(double* Shape, double (*DShape)[3], int nnodes, char* mssg= 0, FILE* file=stderr);
#endif
                                                                                                                                        

HelmPenta::HelmPenta(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];
  nn[4] = nodenums[4];
  nn[5] = nodenums[5];
}

Element *
HelmPenta::clone()
{
 return new HelmPenta(*this);
}

void
HelmPenta::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
}

void
HelmPenta::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
}

double
HelmPenta::getMass(const CoordSet& cs) const
{
 return 0.0;
}

FullSquareMatrix
HelmPenta::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
  int status = 0;
  const int nnodes= 6;	
  //const int ndofs = 6;

  double X[6], Y[6], Z[6];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[6] = {0,1,2,3,4,5};
  FullSquareMatrix M(6,d);
  M.zero();
 
 // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};

  // integration: loop over Gauss pts
  double wxy,wz,w;
  int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 3; // numbder of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[6], DShape[6][3];
  double dOmega;//det of jacobian
  int jSign = 0;                                                                                                                                         
  for(int iz=1;iz<=ngpz;iz++){
    // get z position & weight of the Gauss pt
    _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
      dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "HelmPenta::massMatrix");
#endif
      w = dOmega*wxy*wz;
      addNtNtoM3DHelm(M, Shape, w, nnodes, ls);
    }
  }

  if(status != 0) {
    fprintf(stderr, " *** FATAL ERROR in HelmPenta::stiffness");
    fprintf(stderr, " -> Pent6's nodes: \n");
    for(int i = 0; i < 6; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
    exit(-1);
  }

  M /= getProperty()->rho; 
  return(M);
}
                                                                                                                                        
FullSquareMatrix
HelmPenta::stiffness(const CoordSet &cs, double *d, int flg) const
{
  int status = 0;
  const int nnodes= 6;	
  //const int ndofs = 6;

  double X[6], Y[6], Z[6];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[6] = {0,1,2,3,4,5};
  FullSquareMatrix K(6,d);
  K.zero();

  // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};
  
  // integration: loop over Gauss pts
  double wxy,wz,w;
  int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 3; // numbder of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[6], DShape[6][3];
  double dOmega;//det of jacobian
  int jSign = 0;                                                                                                                                         
  for(int iz=1;iz<=ngpz;iz++){
    // get z position & weight of the Gauss pt
    _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
      dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "HelmPenta::stiffness");
#endif
      w = dOmega*wxy*wz;
      addBtBtoK3DHelm(K, DShape, w, nnodes, ls);
    }
  }

  if(status != 0) {
    fprintf(stderr, " *** FATAL ERROR in HelmPenta::stiffness");
    fprintf(stderr, " -> Pent6's nodes: \n");
    for(int i = 0; i < 6; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
    exit(-1);
  }
  K /= getProperty()->rho; 
  return(K);
}
 

FullSquareMatrix
HelmPenta::acousticm(CoordSet &cs, double *d)
{
  double dK[36], dM[36];

  FullSquareMatrix K = stiffness(cs, dK);
  FullSquareMatrix M = massMatrix(cs, dM, 1);

  double kappa = prop ->kappaHelm;
  double kk = kappa*kappa;

  for (int i=0;i<36;i++)
    d[i] = dK[i] - kk*dM[i];

  FullSquareMatrix Ka(6,d);

  Ka /= getProperty()->rho; 
  return(Ka);
}

int
HelmPenta::numNodes() const {return(6); }

int*
HelmPenta::nodes(int *p) const
{
  if(!p) { p = new int[6]; }
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  p[3] = nn[3];
  p[4] = nn[4];
  p[5] = nn[5];
  
  return(p);
}

int
HelmPenta::numDofs() const { return(6); }

int*
HelmPenta::dofs(DofSetArray &dsa, int *p) const
{
  p[0] = dsa.locate(nn[0],DofSet::Helm);
  p[1] = dsa.locate(nn[1],DofSet::Helm);
  p[2] = dsa.locate(nn[2],DofSet::Helm);
  p[3] = dsa.locate(nn[3],DofSet::Helm);
  p[4] = dsa.locate(nn[4],DofSet::Helm);
  p[5] = dsa.locate(nn[5],DofSet::Helm);

  return(p);
}

void
HelmPenta::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn[0],DofSet::Helm);
  dsa.mark(nn[1],DofSet::Helm);
  dsa.mark(nn[2],DofSet::Helm);
  dsa.mark(nn[3],DofSet::Helm);
  dsa.mark(nn[4],DofSet::Helm);
  dsa.mark(nn[5],DofSet::Helm);
}

int
HelmPenta::getTopNumber() const { return(190); }

void
HelmPenta::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmPenta::addFaces not implemented.\n");
/*
        pset->addQuad(this,nn[0], nn[1], nn[2], nn[3]);
        pset->addQuad(this,nn[4], nn[5], nn[6], nn[7]); // *
        pset->addQuad(this,nn[3], nn[0], nn[4], nn[7]); // *
        pset->addQuad(this,nn[0], nn[1], nn[5], nn[4]); // *
        pset->addQuad(this,nn[2], nn[1], nn[5], nn[6]);
        pset->addQuad(this,nn[3], nn[2], nn[6], nn[7]);
*/
}


