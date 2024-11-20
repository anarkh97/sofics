// ----------------------------------------------------------------
// HB - 03/01/04
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
//          y
//          ^
//          |           Shape functions
//          7           ---------------
//   4+-----+-----+3      Phi[1] = (1/4).(1-x).(1-y).(-1-x-y)
//    |           |       Phi[2] = (1/4).(1+x).(1-y).(-1+x-y) 
//    |           |       Phi[3] = (1/4).(1+x).(1+y).(-1+x+y) 
//   8+     x     +6 -> x Phi[4] = (1/4).(1-x).(1+y).(-1-x+y) 
//    |           |       Phi[5] = (1/2).(1-x.x).(1-y) 
//    |           |       Phi[6] = (1/2).(1+x).(1-y.y) 
//   1+-----+-----+2      Phi[7] = (1/2).(1-x.x).(1+y) 
//          5             Phi[8] = (1/2).(1-x).(1-y.y) 
// ----------------------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// STL
#include <map>

// FEM headers 
#include <Element.d/Element.h>
#include <Math.d/matrix.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad8.d/FaceQuad8.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// Extern routine
extern "C" {
void _FORTRAN(qgauss)(int &, int &, int &, int &,
                      double &,  double &, double &);
}

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceQuad8::RefCoords[8][2] = {{-1.0,-1.0},
                                     { 1.0,-1.0},
                                     { 1.0, 1.0},
                                     {-1.0, 1.0},
                                     { 0.0,-1.0},
                                     { 1.0, 0.0},
                                     { 0.0, 1.0},
                                     {-1.0, 0.0}};

double*
FaceQuad8::ViewRefCoords() { return(FaceQuad8::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceQuad8::FaceQuad8(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
  Nodes[3] = nodenums[3];
  Nodes[4] = nodenums[4];
  Nodes[5] = nodenums[5];
  Nodes[6] = nodenums[6];
  Nodes[7] = nodenums[7];
}

FaceElement *
FaceQuad8::clone()
{
  return new FaceQuad8(Nodes);
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad8::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
  Nodes[3] = OldToNewNodeIds[Nodes[3]];
  Nodes[4] = OldToNewNodeIds[Nodes[4]];
  Nodes[5] = OldToNewNodeIds[Nodes[5]];
  Nodes[6] = OldToNewNodeIds[Nodes[6]];
  Nodes[7] = OldToNewNodeIds[Nodes[7]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF LOCAL METHODS
// -------------------------------
int
FaceQuad8::nQuad4Nodes() { return 4; }

int
FaceQuad8::GetQuad4Node(int i) { return Nodes[i]; }

void
FaceQuad8::GetQuad4Nodes(int *p, int* renumTable)
{
  if(renumTable) {
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
    p[3] = renumTable[Nodes[3]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
    p[3] = Nodes[3];
  }
}

void
FaceQuad8::GetQuad4Nodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceQuad8::nNodes() const { return 8; }

int
FaceQuad8::GetNode(int i) const { return Nodes[i]; }

void
FaceQuad8::GetNodes(int *p, int* renumTable) const
{
  if(renumTable) {
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
    p[3] = renumTable[Nodes[3]];
    p[4] = renumTable[Nodes[4]];
    p[5] = renumTable[Nodes[5]];
    p[6] = renumTable[Nodes[6]];
    p[7] = renumTable[Nodes[7]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
    p[3] = Nodes[3];
    p[4] = Nodes[4];
    p[5] = Nodes[5];
    p[6] = Nodes[6];
    p[7] = Nodes[7];
  }
}

void
FaceQuad8::GetNodes(int *p, std::map<int,int>& renumTable) const
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
  p[4] = renumTable[Nodes[4]];
  p[5] = renumTable[Nodes[5]];
  p[6] = renumTable[Nodes[6]];
  p[7] = renumTable[Nodes[7]];
}

int 
FaceQuad8::GetNodeIndex(int gNode) const
{
  int i;
  bool found = false; 
  for(i=0; i<8; i++)
    if(gNode==Nodes[i]) { found = true; break; }
  if(!found)
    filePrint(stderr," *** WARNING: FaceQuad8::GetNodeIndex(): node (%6d) does not belong to this element\n", gNode);
  return i; 
}

int
FaceQuad8::GetFaceElemType() { return FaceElement::QUADFACEQ8; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceQuad8::GetACMEFaceElemType() { return ContactSearch::QUADFACEQ8; }
#else
int
FaceQuad8::GetACMEFaceElemType() { return 2; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int
FaceQuad8::nVertices() { return nQuad4Nodes(); }

int
FaceQuad8::GetVertex(int i) { return GetQuad4Node(i); }

void
FaceQuad8::GetVertices(int* p, int* renumTable) { GetQuad4Nodes(p, renumTable); }

void
FaceQuad8::GetVertices(int* p, std::map<int,int>& renumTable) { GetQuad4Nodes(p, renumTable); }

// As ACME doesn't support the Quad8 face element for
// FaceFaceInteraction (FFI), we will pass to it the Quad4 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Quad8 face element has straight edges, but its is an
// APPROXIMATION in the general case (i.e. curved edges/face).
#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceQuad8::GetACMEFFIFaceElemType() { return ContactSearch::QUADFACEL4; }
#else
int
FaceQuad8::GetACMEFFIFaceElemType() { return 1; }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad8::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  LocalToGlobalCoord<double,CoordSet>(M, m, cs);
}

void
FaceQuad8::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFctVal<double>(Shape, m);
}

double
FaceQuad8::GetJacobian(double *m, CoordSet &cs)
{
  return GetJacobian<double,CoordSet>(m, cs);
}

double
FaceQuad8::GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs)
{
  return GetIsoParamMappingNormalAndJacobian<double,CoordSet>(Normal, m, cs);
}

void
FaceQuad8::GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs)
{
  GetIsoParamMappingNormalJacobianProduct<double,CoordSet>(JNormal, m, cs);
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
FaceQuad8::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
FaceQuad8::Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}

void
FaceQuad8::Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m)
{
  Getd3ShapeFct<double>(d3Shapex, d3Shapey, d2Shapex2y, d2Shapexy2, m);
}

void
FaceQuad8::ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs)
{
  ComputedMdxAnddMdy<double,CoordSet>(dMdx, dMdy, m, cs);
}

void
FaceQuad8::Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs)
{
  Computed2Mdx2d2Mdy2Andd2Mdxdy<double,CoordSet>(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);
}

void
FaceQuad8::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs)
{
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2<double,CoordSet>(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);
}

void
FaceQuad8::GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
{
  GetdJNormal<double,CoordSet>(dJNormal, m, cs);
}

void
FaceQuad8::Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs)
{
  Getd2JNormal<double,CoordSet>(d2JNormal, m, cs);
}

void
FaceQuad8::ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs)
{
  ComputedJNormaldxAnddJNormaldy<double,CoordSet>(dJNormaldx, dJNormaldy, m, cs);
}

void
FaceQuad8::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs)
{
  Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy<double,CoordSet>(d2JNormaldx2, d2JNormaldy2, d2JNormaldxdy, m, cs);
}

void
FaceQuad8::ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs)
{
  ComputeddJNormaldxAndddJNormaldy<double,CoordSet>(ddJNormaldx, ddJNormaldy, m, cs);
}

void
FaceQuad8::GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs)
{
  GetUnitNormal<double,CoordSet>(UnitNormal, m, cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceQuad8::ComputeArea(CoordSet &cs, const int ngp=2)
{
  double Area = 0.0;
  double dA, Shape[4];
  double xi, eta, weight, m[2];

  for(int i=0;i<ngp;i++) {
    for(int j=0;j<ngp;j++) {
      _FORTRAN(qgauss)(ngp,i,ngp,j,&xi,&eta,&weight);
      m[0] = xi; m[1] = eta;
      //dA = ComputeDiffSurfNormaleAndJacobian(Normal, m, cs);
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      Area += weight*dA;
    }
  }
  return Area;
} */

// -----------------------------------------------------------------------------------------------------
//                                            MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
FullM
FaceQuad8::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  double dA, Shape[8];
  double xi, eta, weight, m[2];

  FullM Mass(8);
  Mass.zero();

  for(int igp=0;igp<ngp;igp++) {
    for(int jgp=0;jgp<ngp;jgp++) {
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      // upper part
      for(int i=0;i<8;i++)
        for(int j=i;j<8;j++)
          Mass[i][j] += rho*weight*dA*Shape[i]*Shape[j];
    }
  }
  // lower part by symmetry 
  for(int i=0;i<8;i++)
    for(int j=0;j<i;j++)
      Mass[i][j] = Mass[j][i];

  return(Mass);
}

void 
FaceQuad8::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  for(int i=0; i<8; i++) ShapeIntg[i] = 0.0;

  double dA, Shape[8];
  double xi, eta, weight, m[2];

  for(int igp=0;igp<ngp;igp++) {
    for(int jgp=0;jgp<ngp;jgp++) {
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      for(int i=0;i<8;i++)
        ShapeIntg[i] += rho*weight*dA*Shape[i];
    }
  }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceQuad8::printNodes() const
{
  filePrint(stderr,"   # Quad8 face el., nodes = %6d %6d %6d %6d %6d %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2],Nodes[3],
                                                                                    Nodes[4],Nodes[5],Nodes[6],Nodes[7]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad8::print() const
{
  printNodes();
}
