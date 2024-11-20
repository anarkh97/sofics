/* --------------------------------------------------------------
 HB - 08/15/03
 ----------------------------------------------------------------
 WARNINGS: NEED TO FOLLOW THE ACME NODES NUMBERING !!
    s 
    ^
    | 
   2+
    |\             Shape functions
    | \            ---------------
    |  \             Phi[1] = r
    |   \            Phi[2] = s
    |    \           Phi[3] = t = 1.-r-s
    |     \
   3+------+1 -> r

 ----------------------------------------------------------------*/
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
#include <Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.h>
#include <Utils.d/dofset.h>
#include <Hetero.d/FlExchange.h>
#include <Element.d/State.h>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// -----------------------------------------------------------------------------------------------------
//                                  STATIC MEMBERS
// -----------------------------------------------------------------------------------------------------
// coords of the nodes in the ref./parametric domain
double FaceTri3::RefCoords[3][2] = {{ 1.0, 0.0 },
                                    { 0.0, 1.0 },
                                    { 0.0, 0.0 }};

double*
FaceTri3::ViewRefCoords() { return(FaceTri3::RefCoords[0]); }

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------

FaceTri3::FaceTri3(int* nodenums)
{
  Nodes[0] = nodenums[0];
  Nodes[1] = nodenums[1];
  Nodes[2] = nodenums[2];
}

FaceElement *
FaceTri3::clone()
{
  return new FaceTri3(Nodes);
}


// -----------------------------------------------------------------------------------------------------
//                                       SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri3::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  Nodes[0] = OldToNewNodeIds[Nodes[0]];
  Nodes[1] = OldToNewNodeIds[Nodes[1]];
  Nodes[2] = OldToNewNodeIds[Nodes[2]];
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
int
FaceTri3::nNodes() const { return 3; }

int
FaceTri3::GetNode(int i) const { return Nodes[i]; }

void
FaceTri3::GetNodes(int *p, int* renumTable) const
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
FaceTri3::GetNodes(int *p, std::map<int,int>& renumTable) const
{
  p[0] = renumTable[Nodes[0]];
  p[1] = renumTable[Nodes[1]];
  p[2] = renumTable[Nodes[2]];
}

int 
FaceTri3::GetNodeIndex(int gNode) const
{
  int i;
  bool found = false; 
  for(i=0; i<3; i++)
    if(gNode==Nodes[i]) { found = true; break; }
  if(!found)
    filePrint(stderr," *** WARNING: FaceTri3::GetNodeIndex(): node (%6d) does not belong to this element\n", gNode);
  return i; 
}

int
FaceTri3::GetFaceElemType() { return FaceElement::TRIFACEL3; }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri3::GetACMEFaceElemType() { return ContactSearch::TRIFACEL3; }
#else
int
FaceTri3::GetACMEFaceElemType() { return 3; }
#endif

// -> for dealing with quadratic face element (see FaceElement.h for more details)
int
FaceTri3::nVertices() { return nNodes(); }

int
FaceTri3::GetVertex(int i) { return GetNode(i); }

void
FaceTri3::GetVertices(int* p, int* renumTable) { GetNodes(p, renumTable); }

void
FaceTri3::GetVertices(int* p, std::map<int,int>& renumTable) { GetNodes(p, renumTable); }

#ifdef USE_ACME
ContactSearch::ContactFace_Type
FaceTri3::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#else
int
FaceTri3::GetACMEFFIFaceElemType() { return GetACMEFaceElemType(); }
#endif

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri3::LocalToGlobalCoord(double *M, double *m, CoordSet &cs)
{
  LocalToGlobalCoord<double,CoordSet>(M, m, cs);
}

void
FaceTri3::GetShapeFctVal(double *Shape, double *m)
{
  GetShapeFctVal<double>(Shape, m);
}

double
FaceTri3::GetJacobian(double *m, CoordSet &cs)
{
  return GetJacobian<double,CoordSet>(m, cs);
}

double
FaceTri3::GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet& cs)
{
  return GetIsoParamMappingNormalAndJacobian<double,CoordSet>(Normal, m, cs);
}

void
FaceTri3::GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet& cs)
{
  GetIsoParamMappingNormalJacobianProduct<double,CoordSet>(JNormal, m, cs);
}

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
FaceTri3::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
FaceTri3::Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}

void
FaceTri3::Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m)
{
  Getd3ShapeFct<double>(d3Shapex, d3Shapey, d2Shapex2y, d2Shapexy2, m);
}

void
FaceTri3::ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs)
{
  ComputedMdxAnddMdy<double,CoordSet>(dMdx, dMdy, m, cs);
}

void
FaceTri3::Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs)
{
  Computed2Mdx2d2Mdy2Andd2Mdxdy<double,CoordSet>(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);
}

void
FaceTri3::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs)
{
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2<double,CoordSet>(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);
}

void
FaceTri3::GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
{
  GetdJNormal<double,CoordSet>(dJNormal, m, cs);
}

void
FaceTri3::Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs)
{
  Getd2JNormal<double,CoordSet>(d2JNormal, m, cs);
}

void
FaceTri3::ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs)
{
  ComputedJNormaldxAnddJNormaldy<double,CoordSet>(dJNormaldx, dJNormaldy, m, cs);
}

void
FaceTri3::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs)
{
  Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy<double,CoordSet>(d2JNormaldx2, d2JNormaldy2, d2JNormaldxdy, m, cs);
}

void
FaceTri3::ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs)
{
  ComputeddJNormaldxAndddJNormaldy<double,CoordSet>(ddJNormaldx, ddJNormaldy, m, cs);
}

void
FaceTri3::GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs)
{
  GetUnitNormal<double,CoordSet>(UnitNormal, m, cs); 
}

void
FaceTri3::GetdUnitNormal(double dUnitNormal[][3], double* m, CoordSet& cs)
{
  GetdUnitNormal<double,CoordSet>(dUnitNormal, m, cs);
}

void
FaceTri3::Getd2UnitNormal(double d2UnitNormal[][3], double* m, CoordSet& cs)
{
  Getd2UnitNormal<double,CoordSet>(d2UnitNormal, m, cs);
}

void
FaceTri3::GlobalToLocalCoord(double *m, double *M, CoordSet &cs)
{
  GlobalToLocalCoord<double,CoordSet>(m, M, cs);
}

void
FaceTri3::ComputedmdXdmdYAnddmdZ(double *dmdX, double *dmdY, double *dmdZ, double *M, CoordSet &cs)
{
  ComputedmdXdmdYAnddmdZ<double,CoordSet>(dmdX, dmdY, dmdZ, M, cs);
}

void
FaceTri3::Computed2mdX2d2mdY2Etc(double *d2mdX2, double *d2mdY2, double *d2mdZ2, double *d2mdXdY, double *d2mdYdZ, double *d2mdXdZ, double *M, CoordSet &cs)
{
  Computed2mdX2d2mdY2Etc<double,CoordSet>(d2mdX2, d2mdY2, d2mdZ2, d2mdXdY, d2mdYdZ, d2mdXdZ, M, cs);
}

void
FaceTri3::GetdLocalCoords(double dLocalCoords[][2], double *M, CoordSet &cs)
{
  GetdLocalCoords<double,CoordSet>(dLocalCoords, M, cs);
}

void
FaceTri3::Getd2LocalCoords(double d2LocalCoords[][2], double *M, CoordSet &cs)
{
  Getd2LocalCoords<double,CoordSet>(d2LocalCoords, M, cs);
}

void
FaceTri3::GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ(double ddmdX[][2], double ddmdY[][2], double ddmdZ[][2], double *M, CoordSet &cs)
{
  GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ<double,CoordSet>(ddmdX, ddmdY, ddmdZ, M, cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
FaceTri3::ComputeArea(CoordSet &cs, const int ngp=0)
{
  return 0.5*GetJacobian((double*)NULL, cs);
} */

// -----------------------------------------------------------------------------------------------------
//                                            MASS MATRIX METHODS
// -----------------------------------------------------------------------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
FullM
FaceTri3::ScalarMass(CoordSet &cs, double rho, int ngp)
{
  FullM Mass(3);
  Mass.zero();
  double Area = 0.5*GetJacobian<double,CoordSet>((double*)NULL, cs);
  Area *= rho/24;
  Mass[0][0] = 2.*Area; Mass[0][1] =    Area; Mass[0][1] =    Area;
  Mass[1][0] =    Area; Mass[1][1] = 2.*Area; Mass[1][1] =    Area;
  Mass[2][0] =    Area; Mass[2][1] =    Area; Mass[2][2] = 2.*Area;

  return(Mass);
}

void 
FaceTri3::IntegrateShapeFcts(double* ShapeIntg, CoordSet& cs, double rho, int ngp)
{
  double Area = 0.5*GetJacobian<double,CoordSet>((double*)NULL, cs);
  Area *= rho/6;
  ShapeIntg[0] = Area;
  ShapeIntg[1] = Area;
  ShapeIntg[2] = Area;
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
FaceTri3::printNodes() const
{
  filePrint(stderr,"   # Tri3  face el., nodes = %6d %6d %6d\n",Nodes[0],Nodes[1],Nodes[2]);
}

// --------------------------------------
// IMPLEMENTATION OF PURE VIRTUAL METHODS
// --------------------------------------
void
FaceTri3::print() const
{
  printNodes();
}

// -----------------------------------------------------------------------------------------------------
//                                            FS COMMUNICATION (KW) 
// -----------------------------------------------------------------------------------------------------
int* FaceTri3::dofs(DofSetArray &dsa, int *p, int *fnId) const
{
  if(p == 0) p = new int[9];
  dsa.number(fnId[Nodes[0]], DofSet::XYZdisp, p);
  dsa.number(fnId[Nodes[1]], DofSet::XYZdisp, p+3);
  dsa.number(fnId[Nodes[2]], DofSet::XYZdisp, p+6);
  return p;
}

void FaceTri3::computeDisp(CoordSet&, State &state, const InterpPoint &ip, double *res, 
                           GeomState*, int *fnId) 
{
  const double *gp = ip.xy;
  double xyz[3][6];
  state.getDV(fnId[Nodes[0]], xyz[0], xyz[0]+3);
  state.getDV(fnId[Nodes[1]], xyz[1], xyz[1]+3);
  state.getDV(fnId[Nodes[2]], xyz[2], xyz[2]+3);

  for(int j=0; j<6; ++j)
    res[j] = gp[0]*xyz[0][j] + gp[1]*xyz[1][j] + (1.0-gp[0]-gp[1])*xyz[2][j]; //using ACME convention
}

void FaceTri3::getFlLoad(const InterpPoint &ip, double *flF, double *resF) 
{
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]    = gp[0] * flF[i]; //using ACME convention
    resF[3+i]  = gp[1] * flF[i];
    resF[6+i] = (1.0-gp[0]-gp[1]) * flF[i];
  }
}
