// ----------------------------------------------------------------
// HB - 08/15/03
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
//    s
//    ^
//    |
//   2+
//    |\             Shape functions
//    | \            ---------------
//    |  \             Phi[1] = r*(2.r-1)
//   5+   +4           Phi[2] = s*(2.s-1)
//    |    \           Phi[3] = t*(2.t-1), t = 1-r-s
//    |     \          Phi[4] = 4.r.s
//    |      \         Phi[5] = 4.s.t
//   3+---+---+1 -> r  Phi[6] = 4.t.r
//        6
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
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/FaceElement.d/FaceTri6.d/FaceTri6.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceTri6::RefCoords[6][2] = {{ 1.0 , 0.0 },
                                    { 0.0 , 1.0 },
                                    { 0.0 , 0.0 },
                                    { 0.5 , 0.5 },
                                    { 0.0 , 0.5 },
                                    { 0.5 , 0.0 }};

double*
FaceTri6::ViewRefCoords() { return(FaceTri6::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceTri6::FaceTri6(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
  Nodes[3] = nodenums[3];
  Nodes[4] = nodenums[4];
  Nodes[5] = nodenums[5];
}

FaceElement *
FaceTri6::clone()
{
  return new FaceTri6(Nodes);
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri6::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
  Nodes[3] = OldToNewNodeIds[Nodes[3]];
  Nodes[4] = OldToNewNodeIds[Nodes[4]];
  Nodes[5] = OldToNewNodeIds[Nodes[5]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
int
FaceTri6::nTri3Nodes() { return 3; }

int
FaceTri6::GetTri3Node(int i) { return Nodes[i]; }

void
FaceTri6::GetTri3Nodes(int *p, int* renumTable)
{
  if(renumTable) {
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
  }
}

void
FaceTri6::GetTri3Nodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceTri6::nNodes() const { return 6; }

int
FaceTri6::GetNode(int i) const { return Nodes[i]; }

void
FaceTri6::GetNodes(int *p, int* renumTable) const
{
  if(renumTable) {
    p[0] = renumTable[Nodes[0]];
    p[1] = renumTable[Nodes[1]];
    p[2] = renumTable[Nodes[2]];
    p[3] = renumTable[Nodes[3]];
    p[4] = renumTable[Nodes[4]];
    p[5] = renumTable[Nodes[5]];
  } else {
    p[0] = Nodes[0];
    p[1] = Nodes[1];
    p[2] = Nodes[2];
    p[3] = Nodes[3];
    p[4] = Nodes[4];
    p[5] = Nodes[5];
  }
}

void
FaceTri6::GetNodes(int *p, std::map<int,int>& renumTable) const
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
  p[4] = renumTable[Nodes[4]];
  p[5] = renumTable[Nodes[5]];
}

int 
FaceTri6::GetNodeIndex(int gNode) const
{
  int i;
  bool found = false; 
  for(i=0; i<6; i++)
    if(gNode==Nodes[i]) { found = true; break; }
  if(!found)
    filePrint(stderr," *** WARNING: FaceTri6::GetNodeIndex(): node (%6d) does not belong to this element\n", gNode);
  return i; 
}

int
FaceTri6::GetFaceElemType() { return FaceElement::TRIFACEQ6; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri6::GetACMEFaceElemType() { return ContactSearch::TRIFACEQ6; }
#else
int
FaceTri6::GetACMEFaceElemType() { return 4; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int
FaceTri6::nVertices() { return nTri3Nodes(); }

int
FaceTri6::GetVertex(int i) { return GetTri3Node(i); }

void
FaceTri6::GetVertices(int* p, int* renumTable) { GetTri3Nodes(p, renumTable); }

void
FaceTri6::GetVertices(int* p, std::map<int,int>& renumTable) { GetTri3Nodes(p, renumTable); }

// As ACME doesn't support the Tri6 face element for
// FaceFaceInteraction (FFI), we will pass to it the Tri3 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Tri6 face element has straight edges, but its is an
// APPROXIMATION in the general case (i.e. curved edges/face).
#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri6::GetACMEFFIFaceElemType() { return ContactSearch::TRIFACEL3; }
#else
int
FaceTri6::GetACMEFFIFaceElemType() { return 3; }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri6::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  LocalToGlobalCoord<double,CoordSet>(M, m, cs);
}

void
FaceTri6::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFctVal<double>(Shape, m);
}

double
FaceTri6::GetJacobian(double *m, CoordSet &cs)
{
  return GetJacobian<double,CoordSet>(m, cs);
}

double
FaceTri6::GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs)
{
  return GetIsoParamMappingNormalAndJacobian<double,CoordSet>(Normal, m, cs);
}

void
FaceTri6::GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs)
{
  GetIsoParamMappingNormalJacobianProduct<double,CoordSet>(JNormal, m, cs);
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
FaceTri6::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
FaceTri6::Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}

void
FaceTri6::Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m)
{
  Getd3ShapeFct<double>(d3Shapex, d3Shapey, d2Shapex2y, d2Shapexy2, m);
}

void
FaceTri6::ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs)
{
  ComputedMdxAnddMdy<double,CoordSet>(dMdx, dMdy, m, cs);
}

void
FaceTri6::Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs)
{
  Computed2Mdx2d2Mdy2Andd2Mdxdy<double,CoordSet>(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);
}

void
FaceTri6::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs)
{
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2<double,CoordSet>(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);
}

void
FaceTri6::GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
{
  GetdJNormal<double,CoordSet>(dJNormal, m, cs);
}

void
FaceTri6::Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs)
{
  Getd2JNormal<double,CoordSet>(d2JNormal, m, cs);
}

void
FaceTri6::ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs)
{
  ComputedJNormaldxAnddJNormaldy<double,CoordSet>(dJNormaldx, dJNormaldy, m, cs);
}

void
FaceTri6::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs)
{
  Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy<double,CoordSet>(d2JNormaldx2, d2JNormaldy2, d2JNormaldxdy, m, cs);
}

void
FaceTri6::ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs)
{
  ComputeddJNormaldxAndddJNormaldy<double,CoordSet>(ddJNormaldx, ddJNormaldy, m, cs);
}

void
FaceTri6::GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs)
{
  GetUnitNormal<double,CoordSet>(UnitNormal, m, cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceTri6::ComputeArea(CoordSet &cs,const int ngp=0)
{
  return GetJacobian(cs);
} */

// -----------------------------------------------------------------------------------------------------
//                                            MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
FullM
FaceTri6::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  FullM Mass(6);
  Mass.zero();
  fprintf(stderr," *** WARNING:  FaceTri6::ScalarMass(...): NOT IMPLEMENTED !!!\n");
  return Mass;
}

void 
FaceTri6::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  fprintf(stderr," *** WARNING: FaceTri6::IntegrateShapeFcts(...): NOT IMPLEMENTED !!!\n");
  for(int i=0;i<6;i++) { ShapeIntg[i] = 0.0; }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceTri6::printNodes() const
{
  filePrint(stderr,"   # Tri6 face el., nodes = %6d %6d %6d %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2],
                                                                           Nodes[3],Nodes[4],Nodes[5]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri6::print() const
{
  printNodes();
}
