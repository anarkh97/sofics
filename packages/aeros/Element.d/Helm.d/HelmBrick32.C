// ---------------------------------------------------------------------
// HB - 05-24-05
// ---------------------------------------------------------------------
// 32 nodes brick element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
// EXPERIMENTAL ...
// ---------------------------------------------------------------------
#include <Element.d/Helm.d/HelmBrick32.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

#define CHECK_JACOBIAN //HB: force check nullity & constant sign of jacobian over el.
//#define HELMBRICK32_DEBUG 
extern bool useFull;
extern "C"      {
void _FORTRAN(lgauss)(const int &, int &, double *, double *);
}

// HB: for stiffness matrix with ansitropic constitutive matrix and/or consistent mass matrix
extern void   Hexa32ShapeFct(double Shape[32], double dShape[32][3], double m[3]);
extern double Hexa32ShapeFct(double Shape[32], double DShape[32][3], double m[3], double X[32], double Y[32], double Z[32]);
extern double computeHexa32DShapeFct(double dShape[32][3], double X[32], double Y[32], double Z[32], double (*DShape)[3] = 0);

extern void addBtBtoK3DHelm(FullSquareMatrix &K, double (*DShape)[3], double alpha, int nnodes, int* ls);
extern void addNtNtoM3DHelm(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls);
                                                                                                                             
#ifdef CHECK_JACOBIAN //HB: for checking zero/small & constant sign of jacobian over the el.
extern int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);
extern void printShapeFct3D(double* Shape, double (*DShape)[3], int nnodes, char* mssg= 0, FILE* file=stderr);
#endif
                                                                                                                             

HelmBrick32::HelmBrick32(int* nodenums)
{
  for(int i=0; i<32; i++)
    nn[i] = nodenums[i];
}

Element *
HelmBrick32::clone()
{
  return(new HelmBrick32(*this));
}

void
HelmBrick32::renum(const int *table)
{
  for(int i=0; i<32; i++) { nn[i] = table[nn[i]]; }
}

void
HelmBrick32::renum(EleRenumMap& table)
{
  for(int i=0; i<32; i++) { nn[i] = table[nn[i]]; }
}

double
HelmBrick32::getMass(const CoordSet& cs) const
{
  fprintf(stderr," *** WARNING: HelmBrick32::getMass NOT implemented. Return zero mass.\n");
  return(0.0);
}

FullSquareMatrix
HelmBrick32::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
  //int status = 0;
  const int nnodes  = 32;
  const int ndofs   = 32;
  const int numgauss= 4;

  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,d);

  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In HelmBrick32::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[32] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                  17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
 
    // integration: loop over Gauss pts
    double x,y,z,wx,wy,wz,w;
    double m[3], Shape[32], DShape[32][3];
    double dOmega;//det of jacobian
    int jSign = 0;                                                                                                                                        
    for(int i=1;i<=numgauss;i++){
      _FORTRAN(lgauss)(numgauss,i,&x,&wx);
      for(int j=1;j<=numgauss;j++){
        _FORTRAN(lgauss)(numgauss,j,&y,&wy);
        for(int k=1;k<=numgauss;k++){
          _FORTRAN(lgauss)(numgauss,k,&z,&wz);
          m[0] = x; m[1] = y; m[2] = z;
          dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef HELMBRICK32_DEBUG
          fprintf(stderr," *** In HelmBrick32::massMatrix: i = %d, j = %d, k = %d, J = %e\n",i,j,k,dOmega);
#endif
#ifdef CHECK_JACOBIAN
          checkJacobian(&dOmega, &jSign, getGlNum()+1, "HelmBrick32::massMatrix");
#endif
          w = dOmega*wx*wy*wz;
          addNtNtoM3DHelm(M, Shape, w, nnodes, ls);
        }
      }
    }
  } else { // Lumped mass matrix
    fprintf(stderr," *** In HelmBrick32::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  M /= getProperty()->rho;

  return(M);
}


FullSquareMatrix
HelmBrick32::stiffness(const CoordSet &cs, double *d, int flg) const
{
  //int status = 0;
  const int nnodes  = 32;
  const int ndofs   = 32;
  const int numgauss= 4;

  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

#ifdef HELMBRICK32_DEBUG
  fprintf(stderr," HelmBrick32 nodes = \n");
  for(int i=0; i<32; i++)
    fprintf(stderr,"  node[%2d] = %6d = %e  %e  %e \n",i+1, nn[i], X[i], Y[i], Z[i]);
  fprintf(stderr,"\n");
#endif

  FullSquareMatrix K(ndofs,d);
  K.zero();
  int ls[32] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
                                                                                                                             
  // integration: loop over Gauss pts
  double x,y,z,wx,wy,wz,w;
  double m[3], Shape[32], DShape[32][3];
  double dOmega;//det of jacobian
  int jSign = 0;                                                                                                                                          
  for(int i=1;i<=numgauss;i++){
    _FORTRAN(lgauss)(numgauss,i,&x,&wx);
    for(int j=1;j<=numgauss;j++){
      _FORTRAN(lgauss)(numgauss,j,&y,&wy);
      for(int k=1;k<=numgauss;k++){
        _FORTRAN(lgauss)(numgauss,k,&z,&wz);
        m[0] = x; m[1] = y; m[2] = z;
        dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef BRICK32_DEBUG
        fprintf(stderr," *** In HelmBrick32::stiffness: i = %d, j = %d, k = %d, J = %e\n",i,j,k,dOmega);
#endif
#ifdef CHECK_JACOBIAN
        checkJacobian(&dOmega, &jSign, getGlNum()+1, "HelmBrick32::stiffness");
#endif
        w = dOmega*wx*wy*wz;
        addBtBtoK3DHelm(K, DShape, w, nnodes, ls);
      }
    }
  }

  K /= getProperty()->rho;

  return(K);
}

FullSquareMatrix
HelmBrick32::acousticm(CoordSet &cs, double *d)
{
  double dK[1024], dM[1024];

  FullSquareMatrix K = stiffness(cs, dK, 1);
  FullSquareMatrix M = massMatrix(cs, dM, 1);

  double kappa = prop ->kappaHelm;
  double kk = kappa*kappa;

  for(int i=0;i<1024;i++)
    d[i] = dK[i] - kk*dM[i];

  /* take advantage of the symmetry to save some computation: NEED TESTING ...
  for(int i=0;i<32;i++)
    for(int j=i;j<32;j++){
      d[32*i+j] = dK[32*i+j] - kk*dM[32*i+j];
      d[i+32*j] = d[32*i+j];
    }
  */
  FullSquareMatrix Ka(32,d);

  return(Ka);
}

int
HelmBrick32::numNodes() const {
  if(useFull)
    return(32); 
  else
    return(8);
}

int
HelmBrick32::numDofs() const { return(32); }

int*
HelmBrick32::nodes(int *p) const
{
  if(useFull)
    {
      if(!p) p = new int[32];
      for(int i=0; i<32; i+=2) {
	p[i] = nn[i]; p[i+1] = nn[i+1];
      }
      return(p);
    }
  else
    {
      if(!p) p = new int[8];
      for(int i=0; i<8; i+=1) {
	p[i] = nn[i]; 
      }
      return(p);
    }
}

int*
HelmBrick32::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[32];
  for(int i=0; i<32; i++)
    dsa.number(nn[i],DofSet::Helm, p+i);

  return(p);
}

void
HelmBrick32::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::Helm);
}

// Treat as a Hexa8 (MECH)
//int HelmBrick32::getTopNumber() const { return(117); }
//int HelmBrick32::numTopNodes() { return(8); }

int
HelmBrick32::getTopNumber() const { return(193); }
int
HelmBrick32::numTopNodes() const { return(32); }

void
HelmBrick32::addFaces(PolygonSet *pset)
{
  fprintf(stderr," *** ERROR: HelmBrick32::addFaces NOT implemented. Abort.\n");
  exit(-1);
}


