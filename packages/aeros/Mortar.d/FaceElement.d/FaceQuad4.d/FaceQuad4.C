// ----------------------------------------------------------------
// HB - 05/06/03
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
#include <Mortar.d/FaceElement.d/FaceQuad4.d/FaceQuad4.h>
#include <Utils.d/dofset.h>
#include <Hetero.d/FlExchange.h>
#include <Element.d/State.h>

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
double FaceQuad4::RefCoords[4][2] = {{-1.0,-1.0},
                                     { 1.0,-1.0},
                                     { 1.0, 1.0},
                                     {-1.0, 1.0}};

double*
FaceQuad4::ViewRefCoords() { return(FaceQuad4::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceQuad4::FaceQuad4(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
  Nodes[3] = nodenums[3];
}

FaceElement *
FaceQuad4::clone()
{
  return new FaceQuad4(Nodes);
}

// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad4::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
  Nodes[3] = OldToNewNodeIds[Nodes[3]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceQuad4::nNodes() const { return 4; }

int
FaceQuad4::GetNode(int i) const { return Nodes[i]; }

void
FaceQuad4::GetNodes(int *p, int* renumTable) const
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
FaceQuad4::GetNodes(int *p, std::map<int,int>& renumTable) const
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
  p[3] = renumTable[Nodes[3]];
}

int 
FaceQuad4::GetNodeIndex(int gNode) const
{
  int i;
  bool found = false; 
  for(i=0; i<4; i++)
    if(gNode==Nodes[i]) { found = true; break; }
  if(!found)
    filePrint(stderr," *** WARNING: FaceQuad4::GetNodeIndex(): node (%6d) does not belong to this element\n", gNode);
  return i; 
}

int
FaceQuad4::GetFaceElemType() { return FaceElement::QUADFACEL4; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceQuad4::GetACMEFaceElemType() { return ContactSearch::QUADFACEL4; }
#else
int
FaceQuad4::GetACMEFaceElemType() { return 1; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int
FaceQuad4::nVertices() { return nNodes(); }

int
FaceQuad4::GetVertex(int i) { return GetNode(i); }

void
FaceQuad4::GetVertices(int* p, int* renumTable) { GetNodes(p, renumTable); }

void
FaceQuad4::GetVertices(int* p, std::map<int,int>& renumTable) { GetNodes(p, renumTable); }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceQuad4::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#else
int
FaceQuad4::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
/*
void
FaceQuad4::IsoParamInterpolation(double* V, double* m, double* NdVals, int nComps=1, int* NdMap=0)
{
  double Shape[4];
  GetShapeFct(Shape,m);

  double (*ndVals)[nCmps] = reinterpret_cast<double (*)[nCmps]>(NdVals);  
  if(NdMap)
    for(int i=0; i<nComps; i++)
      V[i] = Shape[0]*ndVals[NdMap[0]][i]+Shape[1]*ndVals[NdMap[1]][i]
            +Shape[2]*ndVals[NdMap[2]][i]+Shape[3]*ndVals[NdMap[3]][i];
  else
    for(int i=0; i<nComps; i++)
      V[i] = Shape[0]*ndVals[0][i]+Shape[1]*ndVals[1][i]
            +Shape[2]*ndVals[2][i]+Shape[3]*ndVals[3][i];
}

void
IsoParamInterpolation(double* V, double* Shape, double* NdVals, int nNds, int nComps=1, int* NdMap=0)
{
  double (*ndVals)[nCmps] = reinterpret_cast<double (*)[nCmps]>(NdVals);  
  std::fill(V,V+nComps,0.0);
  if(NdMap)
    for(int j=0; j<nNds; ++j)
      for(int i=0; i<nComps; ++i) 
        V[i] += Shape[j]*ndVals[NdMap[j]][i];
  else
    for(int j=0; j<nNds; ++j)
      for(int i=0; i<nComps; ++i) 
        V[i] += Shape[j]*ndVals[j][i];
}
*/

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad4::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  LocalToGlobalCoord<double,CoordSet>(M, m, cs);
}

void
FaceQuad4::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFctVal<double>(Shape, m);
}

double
FaceQuad4::GetJacobian(double *m, CoordSet &cs)
{
  return GetJacobian<double,CoordSet>(m, cs);
}

double
FaceQuad4::GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs)
{
  return GetIsoParamMappingNormalAndJacobian<double,CoordSet>(Normal, m, cs);
}

void
FaceQuad4::GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs)
{
  GetIsoParamMappingNormalJacobianProduct<double,CoordSet>(JNormal, m, cs);
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
FaceQuad4::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
FaceQuad4::Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}

void
FaceQuad4::Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m)
{
  Getd3ShapeFct<double>(d3Shapex, d3Shapey, d2Shapex2y, d2Shapexy2, m);
}

void
FaceQuad4::ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs)
{
  ComputedMdxAnddMdy<double,CoordSet>(dMdx, dMdy, m, cs);
}

void
FaceQuad4::Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs)
{
  Computed2Mdx2d2Mdy2Andd2Mdxdy<double,CoordSet>(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);
}

void
FaceQuad4::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs)
{
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2<double,CoordSet>(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);
}

void
FaceQuad4::GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
{
  GetdJNormal<double,CoordSet>(dJNormal, m, cs);
}

void
FaceQuad4::Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs)
{
  Getd2JNormal<double,CoordSet>(d2JNormal, m, cs);
}

void
FaceQuad4::ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs)
{
  ComputedJNormaldxAnddJNormaldy<double,CoordSet>(dJNormaldx, dJNormaldy, m, cs);
}

void
FaceQuad4::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs)
{
  Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy<double,CoordSet>(d2JNormaldx2, d2JNormaldy2, d2JNormaldxdy, m, cs);
}

void
FaceQuad4::ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs)
{
  ComputeddJNormaldxAndddJNormaldy<double,CoordSet>(ddJNormaldx, ddJNormaldy, m, cs);
}

void
FaceQuad4::GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs)
{
  GetUnitNormal<double,CoordSet>(UnitNormal, m, cs);
}

void
FaceQuad4::GetdUnitNormal(double dUnitNormal[][3], double* m, CoordSet& cs)
{
  GetdUnitNormal<double,CoordSet>(dUnitNormal, m, cs);
}

void
FaceQuad4::Getd2UnitNormal(double d2UnitNormal[][3], double* m, CoordSet& cs)
{
  Getd2UnitNormal<double,CoordSet>(d2UnitNormal, m, cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceQuad4::ComputeArea(CoordSet &cs, const int ngp=2)
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
FaceQuad4::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  double dA, Shape[4];
  double xi, eta, weight, m[2];

  FullM Mass(4);
  Mass.zero();

  for(int igp=0;igp<ngp;igp++) {
    for(int jgp=0;jgp<ngp;jgp++) {
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      // upper part
      for(int i=0;i<4;i++)
        for(int j=i;j<4;j++)
          Mass[i][j] += rho*weight*dA*Shape[i]*Shape[j];
    }
  }
  // lower part by symmetry 
  for(int i=0;i<4;i++)
    for(int j=0;j<i;j++)
      Mass[i][j] = Mass[j][i];

  return(Mass);
}

void 
FaceQuad4::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  for(int i=0; i<4; i++) ShapeIntg[i] = 0.0;

  double dA, Shape[4];
  double xi, eta, weight, m[2];

  for(int igp=0;igp<ngp;igp++) {
    for(int jgp=0;jgp<ngp;jgp++) {
      _FORTRAN(qgauss)(ngp,igp,ngp,jgp,xi,eta,weight);
      m[0] = xi; m[1] = eta;
      dA = GetShapeFctAndJacobian(Shape, m, cs);
      for(int i=0;i<4;i++)
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
FaceQuad4::printNodes() const
{
  filePrint(stderr,"   # Quad4 face el., nodes = %6d %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2],Nodes[3]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceQuad4::print() const
{
  printNodes();
}

// -----------------------------------------------------------------------------------------------------
//                                            FS COMMUNICATION
// -----------------------------------------------------------------------------------------------------
int* FaceQuad4::dofs(DofSetArray &dsa, int *p, int *fnId) const
{
  if(p == 0) p = new int[12];
  dsa.number(fnId[Nodes[0]], DofSet::XYZdisp, p);
  dsa.number(fnId[Nodes[1]], DofSet::XYZdisp, p+3);
  dsa.number(fnId[Nodes[2]], DofSet::XYZdisp, p+6);
  dsa.number(fnId[Nodes[3]], DofSet::XYZdisp, p+9);
  return p;
}

void FaceQuad4::computeDisp(CoordSet&, State &state, const InterpPoint &ip, double *res, 
                            GeomState*, int *fnId) 
{
  double xyz[4][6];
  state.getDV(fnId[Nodes[0]], xyz[0], xyz[0]+3);
  state.getDV(fnId[Nodes[1]], xyz[1], xyz[1]+3);
  state.getDV(fnId[Nodes[2]], xyz[2], xyz[2]+3);
  state.getDV(fnId[Nodes[3]], xyz[3], xyz[3]+3);

  double Shape[4];
  GetShapeFctVal(Shape, const_cast<double*>(ip.xy));

  for(int j=0; j<6; ++j)
    res[j] = Shape[0]*xyz[0][j] +
             Shape[1]*xyz[1][j] +
             Shape[2]*xyz[2][j] +
             Shape[3]*xyz[3][j];
}

void FaceQuad4::getFlLoad(const InterpPoint &ip, double *flF, double *resF) 
{
  double Shape[4];
  GetShapeFctVal(Shape, const_cast<double*>(ip.xy));

  for(int i = 0; i < 3; ++i) {
    resF[i]   = Shape[0]*flF[i]; 
    resF[3+i] = Shape[1]*flF[i];
    resF[6+i] = Shape[2]*flF[i];
    resF[9+i] = Shape[3]*flF[i];
  }
}
