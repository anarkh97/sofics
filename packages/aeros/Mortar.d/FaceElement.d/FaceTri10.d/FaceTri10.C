// ----------------------------------------------------------------
// HB - 05/05/05
// ----------------------------------------------------------------
// WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING 
//    s
//    ^
//    |
//   2+
//    |\             Shape functions
//    | \            ---------------
//   6+  +5            see: FaceTri10::GetShapeFct
//    |   \
//    | +  \
//   7+ 10  +4
//    |      \
//   3+-+---+-+1 -> r
//      8   9
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
#include <Mortar.d/FaceElement.d/FaceTri10.d/FaceTri10.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceTri10::RefCoords[10][2] = {{ 1.0  , 0.0   },
                                      { 0.0  , 1.0   },
                                      { 0.0  , 0.0   },
                                      { 2./3., 1./3. },
                                      { 1./3., 2./3. },
                                      { 0.0  , 2./3. },
                                      { 0.0  , 1./3. }, 
                                      { 1./3., 0.0   },
                                      { 2./3., 0.0   },
                                      { 1./3., 1./3. }};

double*
FaceTri10::ViewRefCoords() { return(FaceTri10::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceTri10::FaceTri10(int* nodenums)
{
  for(int i=0; i<10; i++) { Nodes[i] = nodenums[i]; }
}

FaceElement *
FaceTri10::clone()
{
  return new FaceTri10(Nodes);
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri10::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  for(int i=0; i<10; i++) { Nodes[i] = OldToNewNodeIds[Nodes[i]]; }
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
int
FaceTri10::nTri3Nodes() { return 3; }

int
FaceTri10::GetTri3Node(int i) { return Nodes[i]; }

void
FaceTri10::GetTri3Nodes(int *p, int* renumTable)
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
FaceTri10::GetTri3Nodes(int *p, std::map<int,int>& renumTable)
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceTri10::nNodes() const { return 10; }

int
FaceTri10::GetNode(int i) const { return Nodes[i]; }

void
FaceTri10::GetNodes(int *p, int* renumTable) const
{
  if(renumTable) {
    for(int i=0; i<10; i++)
       p[i] = renumTable[Nodes[i]];
  } else {
    for(int i=0; i<10; i++)
       p[i] = Nodes[i];
  }
}

void
FaceTri10::GetNodes(int *p, std::map<int,int>& renumTable) const
{
  for(int i=0; i<10; i++)
    p[i] = renumTable[Nodes[i]];
}

int 
FaceTri10::GetNodeIndex(int gNode) const
{
  int i;
  bool found = false; 
  for(i=0; i<10; i++)
    if(gNode==Nodes[i]) { found = true; break; }
  if(!found)
    filePrint(stderr," *** WARNING: FaceTri10::GetNodeIndex(): node (%6d) does not belong to this element\n", gNode);
  return i; 
}

int
FaceTri10::GetFaceElemType() { return FaceElement::TRIFACEC10; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri10::GetACMEFaceElemType() { return ContactSearch::TRIFACEQ6; } //XXX
#else
int
FaceTri10::GetACMEFaceElemType() { return 4; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int
FaceTri10::nVertices() { return nTri3Nodes(); }

int
FaceTri10::GetVertex(int i) { return GetTri3Node(i); }

void
FaceTri10::GetVertices(int* p, int* renumTable) { GetTri3Nodes(p, renumTable); }

void
FaceTri10::GetVertices(int* p, std::map<int,int>& renumTable) { GetTri3Nodes(p, renumTable); }

// As ACME doesn't support the Tri10 face element for
// FaceFaceInteraction (), we will pass to it the Tri3 face element
// made of its vertices for the (geometric) contact search.
// This is OK if the Tri10 face element has straight edges, but its is an
// APPROXIMATION in the general case (i.e. curved edges/face).
#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri10::GetACMEFFIFaceElemType() { return ContactSearch::TRIFACEL3; }
#else
int
FaceTri10::GetACMEFFIFaceElemType() { return 3; }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri10::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  LocalToGlobalCoord<double,CoordSet>(M, m, cs);
}

void
FaceTri10::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFctVal<double>(Shape, m);
}

double
FaceTri10::GetJacobian(double *m, CoordSet &cs)
{
  return GetJacobian<double,CoordSet>(m, cs);
}

double
FaceTri10::GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs)
{
  return GetIsoParamMappingNormalAndJacobian<double,CoordSet>(Normal, m, cs);
}

void
FaceTri10::GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs)
{
  GetIsoParamMappingNormalJacobianProduct<double,CoordSet>(JNormal, m, cs);
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
FaceTri10::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
FaceTri10::Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}

void
FaceTri10::Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m)
{
  Getd3ShapeFct<double>(d3Shapex, d3Shapey, d2Shapex2y, d2Shapexy2, m);
}

void
FaceTri10::ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs)
{
  ComputedMdxAnddMdy<double,CoordSet>(dMdx, dMdy, m, cs);
}

void
FaceTri10::Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs)
{
  Computed2Mdx2d2Mdy2Andd2Mdxdy<double,CoordSet>(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);
}

void
FaceTri10::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs)
{
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2<double,CoordSet>(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);
}

void
FaceTri10::GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
{
  GetdJNormal<double,CoordSet>(dJNormal, m, cs);
}

void
FaceTri10::Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs)
{
  Getd2JNormal<double,CoordSet>(d2JNormal, m, cs);
}

void
FaceTri10::ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs)
{
  ComputedJNormaldxAnddJNormaldy<double,CoordSet>(dJNormaldx, dJNormaldy, m, cs);
}

void
FaceTri10::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs)
{
  Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy<double,CoordSet>(d2JNormaldx2, d2JNormaldy2, d2JNormaldxdy, m, cs);
}

void
FaceTri10::ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs)
{
  ComputeddJNormaldxAndddJNormaldy<double,CoordSet>(ddJNormaldx, ddJNormaldy, m, cs);
}

void
FaceTri10::GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs)
{
  GetUnitNormal<double,CoordSet>(UnitNormal, m, cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceTri10::ComputeArea(CoordSet &cs,const int ngp=0)
{
  return GetJacobian(cs);
} */

// -----------------------------------------------------------------------------------------------------
//                                            MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
FullM
FaceTri10::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  FullM Mass(10);
  Mass.zero();
  fprintf(stderr," *** WARNING: FaceTri10::ScalarMass(): NOT IMPLEMENTED. Return zero mass matrix.\n");
  return(Mass);
}

void 
FaceTri10::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  fprintf(stderr," *** WARNING: FaceTri10::IntegrateShapeFcts(): NOT IMPLEMENTED. Return zero shape functions'integral.\n");
  for(int i=0;i<10;i++) { ShapeIntg[i] = 0.0; }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceTri10::printNodes() const
{
  filePrint(stderr,"   # Tri10 face el., nodes = %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n",
            Nodes[0],Nodes[1],Nodes[2],Nodes[3],Nodes[4],Nodes[5],Nodes[6],Nodes[7],Nodes[8],Nodes[9]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri10::print() const
{
  printNodes();
}
