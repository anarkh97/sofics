// ----------------------------------------------------------------
// HB - 05/10/05
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
//          y
//          ^
//          |           Shape functions
//       10   9         ---------------
//   4+---+---+---+3      see FaceQuad12::GetShapeFct
//    |           |
//  11+           +8
//    |           | -> x
//  12+           +7
//    |           |
//   1+---+---+---+2
//        5   6
//
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
#include <Mortar.d/FaceElement.d/FaceQuad12.d/FaceQuad12.h>

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
double FaceQuad12::RefCoords[12][2] = {{-1.0  ,-1.0  },
                                       { 1.0  ,-1.0  },
                                       { 1.0  , 1.0  },
                                       {-1.0  , 1.0  },
                                       {-1./3.,-1.0  },
                                       { 1./3.,-1.0  },
                                       { 1.0  ,-1./3.},
                                       { 1.0  , 1./3.},
                                       { 1./3., 1.0  },
                                       {-1./3., 1.0  },
                                       {-1.0  , 1./3.},
                                       {-1.0  ,-1./3.}};

double*
FaceQuad12::ViewRefCoords() { return(FaceQuad12::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceQuad12::FaceQuad12(int* nodenums)
{
  for(int i=0; i<12; i++) { Nodes[i] = nodenums[i]; }
}

FaceElement *
FaceQuad12::clone()
{
  return new FaceQuad12(Nodes);
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad12::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  for(int i=0; i<12; i++) { Nodes[i] = OldToNewNodeIds[Nodes[i]]; }
}

// -----------------------------------------------------------------------------------------------------
//                              GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF LOCAL METHODS
// -------------------------------
int
FaceQuad12::nQuad4Nodes() { return 4; }

int
FaceQuad12::GetQuad4Node(int i) { return Nodes[i]; }

void
FaceQuad12::GetQuad4Nodes(int *p, int* renumTable)
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
FaceQuad12::GetQuad4Nodes(int *p, std::map<int,int>& renumTable)
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
FaceQuad12::nNodes() const { return 12; }

int
FaceQuad12::GetNode(int i) const { return Nodes[i]; }

void
FaceQuad12::GetNodes(int *p, int* renumTable) const
{
  if(renumTable) {
    for(int i=0; i<12; i++)
      p[i] = renumTable[Nodes[i]];
  } else {
    for(int i=0; i<12; i++)
      p[i] = Nodes[i];
  }
}

void
FaceQuad12::GetNodes(int *p, std::map<int,int>& renumTable) const
{
  for(int i=0; i<12; i++)
    p[i] = renumTable[Nodes[i]];
}

int 
FaceQuad12::GetNodeIndex(int gNode) const
{
  int i;
  bool found = false; 
  for(i=0; i<12; i++)
    if(gNode==Nodes[i]) { found = true; break; }
  if(!found)
    filePrint(stderr," *** WARNING: FaceQuad12::GetNodeIndex(): node (%6d) does not belong to this element\n", gNode);
  return i; 
}

int
FaceQuad12::GetFaceElemType() { return FaceElement::QUADFACEC12; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceQuad12::GetACMEFaceElemType() { return ContactSearch::QUADFACEQ8; } //XXX
#else
int
FaceQuad12::GetACMEFaceElemType() { return 2; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int
FaceQuad12::nVertices() { return nQuad4Nodes(); }

int
FaceQuad12::GetVertex(int i) { return GetQuad4Node(i); }

void
FaceQuad12::GetVertices(int* p, int* renumTable) { GetQuad4Nodes(p, renumTable); }

void
FaceQuad12::GetVertices(int* p, std::map<int,int>& renumTable) { GetQuad4Nodes(p, renumTable); }

// As ACME doesn't support the Quad12 face element for
// FaceFaceInteraction (FFI), we will pass to it the Quad4 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Quad12 face element has straight edges, but its is an
// APPROXIMATION in the general case (i.e. curved edges/face).
#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceQuad12::GetACMEFFIFaceElemType() { return ContactSearch::QUADFACEL4; }
#else
int
FaceQuad12::GetACMEFFIFaceElemType() { return 1; }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad12::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  LocalToGlobalCoord<double,CoordSet>(M, m, cs);
}

void
FaceQuad12::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFctVal<double>(Shape, m);
}

double
FaceQuad12::GetJacobian(double *m, CoordSet &cs)
{
  return GetJacobian<double,CoordSet>(m, cs);
}

double
FaceQuad12::GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs)
{
  return GetIsoParamMappingNormalAndJacobian<double,CoordSet>(Normal, m, cs);
}

void
FaceQuad12::GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs)
{
  GetIsoParamMappingNormalJacobianProduct<double,CoordSet>(JNormal, m, cs);
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
FaceQuad12::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
FaceQuad12::Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}

void
FaceQuad12::Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m)
{
  Getd3ShapeFct<double>(d3Shapex, d3Shapey, d2Shapex2y, d2Shapexy2, m);
}

void
FaceQuad12::ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs)
{
  ComputedMdxAnddMdy<double,CoordSet>(dMdx, dMdy, m, cs);
}

void
FaceQuad12::Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs)
{
  Computed2Mdx2d2Mdy2Andd2Mdxdy<double,CoordSet>(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);
}

void
FaceQuad12::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs)
{
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2<double,CoordSet>(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);
}

void
FaceQuad12::GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
{
  GetdJNormal<double,CoordSet>(dJNormal, m, cs);
}

void
FaceQuad12::Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs)
{
  Getd2JNormal<double,CoordSet>(d2JNormal, m, cs);
}

void
FaceQuad12::ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs)
{
  ComputedJNormaldxAnddJNormaldy<double,CoordSet>(dJNormaldx, dJNormaldy, m, cs);
}

void
FaceQuad12::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs)
{
  Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy<double,CoordSet>(d2JNormaldx2, d2JNormaldy2, d2JNormaldxdy, m, cs);
}

void
FaceQuad12::ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs)
{
  ComputeddJNormaldxAndddJNormaldy<double,CoordSet>(ddJNormaldx, ddJNormaldy, m, cs);
}

void
FaceQuad12::GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs)
{
  GetUnitNormal<double,CoordSet>(UnitNormal, m, cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceQuad12::ComputeArea(CoordSet &cs,const int ngp=2)
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
FaceQuad12::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  double dA, Shape[12];
  double weight, m[2];

  FullM Mass(12);
  Mass.zero();

  for(int igp=0;igp<ngp;igp++) {
    for(int jgp=0;jgp<ngp;jgp++) {
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,m[0],m[1],weight);
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      // upper part
      for(int i=0;i<12;i++)
        for(int j=i;j<12;j++)
          Mass[i][j] += rho*weight*dA*Shape[i]*Shape[j];
    }
  }
  // lower part by symmetry 
  for(int i=0;i<12;i++)
    for(int j=0;j<i;j++)
      Mass[i][j] = Mass[j][i];

  return(Mass);
}

void 
FaceQuad12::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  for(int i=0; i<12; i++) ShapeIntg[i] = 0.0;

  double dA, Shape[12];
  double weight, m[2];

  for(int igp=0;igp<ngp;igp++) {
    for(int jgp=0;jgp<ngp;jgp++) {
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,m[0],m[1],weight);
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      for(int i=0;i<12;i++)
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
FaceQuad12::printNodes() const
{
  filePrint(stderr,"   # Quad12 face el., nodes = ");
  for(int i=0; i<12; i++)
    filePrint(stderr," %6d",Nodes[i]);
  filePrint(stderr,"\n");
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad12::print() const
{
  printNodes();
}
