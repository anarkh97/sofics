// ----------------------------------------------------------------
// HB - 08/15/03
// ----------------------------------------------------------------
#ifndef _FACETRI3_H_
#define _FACETRI3_H_

// STL
#include <map>

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;

class FaceTri3: public FaceElement {
  private:
        int Nodes[3];
        static double RefCoords[3][2]; // coords of the nodes in the ref./parametric domain

  public:
        enum { NumberOfNodes=3 };

        // Constructors
        // ~~~~~~~~~~~~
        FaceTri3(int *);
        FaceElement* clone() override;

        // Setup & update methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of pure virtual methods
        void Renumber(std::map<int,int>& OldToNewNodeIds) override;

        // Get methods
        // ~~~~~~~~~~~
        // -> implementation of pure virtual methods
        int  nNodes() const override;
        void GetNodes(int*, int* renumTable) const override;
        void GetNodes(int*, std::map<int,int>& renumTable) const override;
        int GetNode(int) const override;
        int GetNodeIndex(int) const override;

        int GetFaceElemType() override;
#ifdef USE_ACME
        ContactSearch::ContactFace_Type GetACMEFaceElemType() override;
#else
        int GetACMEFaceElemType() override;
#endif
        // -> pure virtual method for dealing with quadratic face element
        //    (see FaceElement.h for more details)
        int  nVertices() override;
        int  GetVertex(int) override;
        void GetVertices(int*, int* renumTable) override;
        void GetVertices(int*, std::map<int,int>& renumTable) override;

#ifdef USE_ACME
        ContactSearch::ContactFace_Type GetACMEFFIFaceElemType() override;
#else
        int GetACMEFFIFaceElemType() override;
#endif
        // Mapping & shape fct methods        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // -> local methods 
        template<typename Scalar>
          void   GetShapeFctVal(Scalar*, Scalar*);
        template<typename Scalar>
          void   GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m);
        template<typename Scalar>
          void   Getd2ShapeFct(Scalar *d2Shapex, Scalar *d2Shapey, Scalar *d2Shapexy, Scalar *m);
        template<typename Scalar>
          void   Getd3ShapeFct(Scalar *d3Shapex, Scalar *d3Shapey, Scalar *d2Shapex2y, Scalar *d2Shapexy2, Scalar *m);
        template<typename Scalar, typename CoordSetT>
          void   ComputedMdxAnddMdy(Scalar* dMdx, Scalar* dMdy, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   Computed2Mdx2d2Mdy2Andd2Mdxdy(Scalar *d2Mdx2, Scalar *d2Mdy2, Scalar *d2Mdxdy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void   Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(Scalar *d3Mdx3, Scalar *d3Mdy3, Scalar *d3Mdx2dy, Scalar *d3Mdxdy2, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          Scalar GetJacobian(Scalar*, CoordSetT&);
        template<typename Scalar, typename CoordSetT>
          Scalar GetShapeFctAndJacobian(Scalar* Shape, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   LocalToGlobalCoord(Scalar* M, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GetIsoParamMappingNormalJacobianProduct(Scalar* JNormal, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          Scalar GetIsoParamMappingNormalAndJacobian(Scalar* Normal, Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GetdJNormal(Scalar dJNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   Getd2JNormal(Scalar d2JNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   ComputedJNormaldxAnddJNormaldy(Scalar *dJNormaldx, Scalar *dJNormaldy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void   Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(Scalar *d2JNormaldx2, Scalar *d2JNormaldy2, Scalar *d2JNormaldxdy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void   ComputeddJNormaldxAndddJNormaldy(Scalar ddJNormaldx[][3], Scalar ddJNormaldy[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GetUnitNormal(Scalar UnitNormal[3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GetdUnitNormal(Scalar dUnitNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   Getd2UnitNormal(Scalar d2Normal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void   GlobalToLocalCoord(Scalar *m, Scalar *M, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void   ComputedmdXdmdYAnddmdZ(Scalar *dmdX, Scalar *dmdY, Scalar *dmdZ, Scalar *M, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void   Computed2mdX2d2mdY2Etc(Scalar *d2mdX2, Scalar *d2mdY2, Scalar *d2mdZ2, Scalar *d2mdXdY, Scalar *d2mdYdZ, Scalar *d2mdXdZ, Scalar *M, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void  GetdLocalCoords(Scalar dLocalCoords[][2], Scalar *M, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void  Getd2LocalCoords(Scalar d2LocalCoords[][2], Scalar *M, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void  GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ(Scalar ddmdX[][2], Scalar ddmdY[][2], Scalar ddmdZ[][2], Scalar *M, CoordSetT &cs);

        // -> implementation of pure virtual fcts
        void   LocalToGlobalCoord(double*, double*, CoordSet&) override;
        void   GetShapeFctVal(double*, double*) override;
        double GetJacobian(double*, CoordSet&) override;
        double GetIsoParamMappingNormalAndJacobian(double*, double*, CoordSet&) override;
        void   GetIsoParamMappingNormalJacobianProduct(double*, double*, CoordSet&) override;

        // -> implementation of virtual fcts
        double* ViewRefCoords() override;
        void GetdShapeFct(double* dShapex, double* dShapey, double* m) override;
        void Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m) override;
        void Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m) override;
        void ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs) override;
        void Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs) override;
        void Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs) override;
        void GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs) override;
        void Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs) override;
        void ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs) override;
        void Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs) override;
        void ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs) override;
        void GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs) override;
        void GetdUnitNormal(double dUnitNormal[][3], double* m, CoordSet& cs) override;
        void Getd2UnitNormal(double d2UnitNormal[][3], double* m, CoordSet& cs) override;
        void GlobalToLocalCoord(double *m, double *M, CoordSet &cs) override;
        void ComputedmdXdmdYAnddmdZ(double *dmdX, double *dmdY, double *dmdZ, double *M, CoordSet &cs) override;
        void Computed2mdX2d2mdY2Etc(double *d2mdX2, double *d2mdY2, double *d2mdZ2, double *d2mdXdY, double *d2mdYdZ, double *d2mdXdZ, double *M, CoordSet &cs) override;
        void GetdLocalCoords(double dLocalCoords[][2], double *M, CoordSet &cs) override;
        void Getd2LocalCoords(double d2LocalCoords[][2], double *M, CoordSet &cs) override;
        void GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ(double ddmdX[][2], double ddmdY[][2], double ddmdZ[][2], double *M, CoordSet &cs) override;

        // Miscelleaneous methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> implementation of virtual fcts
        //double ComputeArea(CoordSet&, const int ngp=2);

        // Mass matrix methods
        // ~~~~~~~~~~~~~~~~~~~
        // -> implementation of pure virtual fcts
        FullM ScalarMass(CoordSet&, double rho, int ngp) override;
        void  IntegrateShapeFcts(double*, CoordSet&, double rho, int ngp) override;

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        // -> local fcts
        void printNodes() const;

        // -> implementation of pure virtual fcts
        void print() const override;

        int numDofs() const override {return 9;}
        int* dofs(DofSetArray &dsa, int *p, int *fnId) const override;
        void computeDisp(CoordSet&, State &state, const InterpPoint &ip, double *res, GeomState*, int *fnId) override;
        void getFlLoad(const InterpPoint &ip, double *flF, double *resF) override;

};

// -----------------------------------------------------------------------------------------------------
//                          MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
template<typename Scalar>
void
FaceTri3::GetShapeFctVal(Scalar *Shape, Scalar *m)
{
  // Computes the values of the shape functions at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P.
  // Outputs: Shape = [N₀(ξ,η), N₁(ξ,η), N₂(ξ,η)], the values of the shape functions at P.

  Scalar& r = m[0];
  Scalar& s = m[1];
  Scalar t = 1.-r-s;

  // !! idem ACME !! 
  Shape[0] = r;
  Shape[1] = s;
  Shape[2] = t;
}

template<typename Scalar>
void
FaceTri3::GetdShapeFct(Scalar *dShapex, Scalar *dShapey, Scalar *m)
{
  // Computes the values of the shape functions' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m       = [ξ, η], the local coordinates of P.
  // Outputs: dShapex = [∂N₀/∂ξ, ∂N₁/∂ξ, ∂N₂/∂ξ]
  //          dShapey = [∂N₀/∂η, ∂N₁/∂η, ∂N₂/∂η] 
  //                    where N₀, N₁ and N₂ are the shape functions.

  dShapex[0] = 1;
  dShapex[1] = 0;
  dShapex[2] = -1;

  dShapey[0] = 0;
  dShapey[1] = 1;
  dShapey[2] = -1;
}

template<typename Scalar>
void
FaceTri3::Getd2ShapeFct(Scalar *d2Shapex, Scalar *d2Shapey, Scalar *d2Shapexy, Scalar *m)
{
  // Computes the values of the shape functions' second partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m         = [ξ, η], the local coordinates of P.
  // Outputs: d2Shapex  = [∂²N₀/∂ξ², ∂²N₁/∂ξ², ∂²N₂/∂ξ²] 
  //          d2Shapey  = [∂²N₀/∂η², ∂²N₁/∂η², ∂²N₂/∂η²]
  //          d2Shapexy = [∂²N₀/∂ξ∂η, ∂²N₁/∂ξ∂η, ∂²N₂/∂ξ∂η]
  //                      where N₀, N₁ and N₂ are the shape functions.

  d2Shapex[0] = 0;
  d2Shapex[1] = 0;
  d2Shapex[2] = 0;

  d2Shapey[0] = 0;
  d2Shapey[1] = 0;
  d2Shapey[2] = 0;

  d2Shapexy[0] = 0;
  d2Shapexy[1] = 0;
  d2Shapexy[2] = 0;
}

template<typename Scalar>
void
FaceTri3::Getd3ShapeFct(Scalar *d3Shapex, Scalar *d3Shapey, Scalar *d3Shapex2y, Scalar *d3Shapexy2, Scalar *m)
{
  // Computes the values of the shape functions' third partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m          = [ξ, η], the local coordinates of P.
  // Outputs: d3Shapex   = [∂³N₀/∂ξ³, ∂³N₁/∂ξ³, ∂³N₂/∂ξ³]
  //          d3Shapey   = [∂³N₀/∂η³, ∂³N₁/∂η³, ∂³N₂/∂η³]
  //          d3Shapex2y = [∂³N₀/∂ξ²∂η, ∂³N₁/∂ξ²∂η, ∂³N₂/∂ξ²∂η]
  //          d3Shapexy2 = [∂³N₀/∂ξ∂η², ∂³N₁/∂ξ∂η², ∂³N₂/∂ξ∂η²]
  //                       where N₀, N₁ and N₂ are the shape functions.

  d3Shapex[0] = 0;
  d3Shapex[1] = 0;
  d3Shapex[2] = 0;

  d3Shapey[0] = 0;
  d3Shapey[1] = 0;
  d3Shapey[2] = 0;

  d3Shapex2y[0] = 0;
  d3Shapex2y[1] = 0;
  d3Shapex2y[2] = 0;

  d3Shapexy2[0] = 0;
  d3Shapexy2[1] = 0;
  d3Shapexy2[2] = 0;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::LocalToGlobalCoord(Scalar *M, Scalar *m, CoordSetT &cs)
{
  // Computes the global X-, Y- and Z-coordinates of a point P.
  //
  // Inputs:  m  = [ξ, η], the local coordinates of P, and
  //          cs = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: M  = [M₀(ξ,η,X,Y,Z), M₁(ξ,η,X,Y,Z), M₂(ξ,η,X,Y,Z)], the global X-, Y- and Z-coordinates of P.

  Scalar& r = m[0];
  Scalar& s = m[1];
  Scalar t = 1-r-s;

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // !! idem ACME !!
  M[0] = r*X[0]+s*X[1]+t*X[2];
  M[1] = r*Y[0]+s*Y[1]+t*Y[2];
  M[2] = r*Z[0]+s*Z[1]+t*Z[2];
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::ComputedMdxAnddMdy(Scalar *dMdx, Scalar *dMdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P, and
  //          cs    = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dMdx  = [∂M₀/∂ξ, ∂M₁/∂ξ, ∂M₂/∂ξ]
  //          dMdy  = [∂M₀/∂η, ∂M₁/∂η, ∂M₂/∂η]
  //                  where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  dMdx[0] = X[0] - X[2];
  dMdx[1] = Y[0] - Y[2];
  dMdx[2] = Z[0] - Z[2];

  dMdy[0] = X[1] - X[2];
  dMdy[1] = Y[1] - Y[2];
  dMdy[2] = Z[1] - Z[2];
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::Computed2Mdx2d2Mdy2Andd2Mdxdy(Scalar *d2Mdx2, Scalar *d2Mdy2, Scalar *d2Mdxdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' second partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs: m        = [ξ, η], the local coordinates of P, and
  //         cs       = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d2Mdx2  = [∂²M₀/∂ξ², ∂²M₁/∂ξ², ∂²M₂/∂ξ²]
  //          d2Mdy2  = [∂²M₀/∂η², ∂²M₁/∂η², ∂²M₂/∂η²]
  //          d2Mdxdy = [∂²M₀/∂ξ∂η, ∂²M₁/∂ξ∂η, ∂²M₂/∂ξ∂η]
  //                    where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  d2Mdx2[0] = 0;
  d2Mdx2[1] = 0;
  d2Mdx2[2] = 0;

  d2Mdy2[0] = 0;
  d2Mdy2[1] = 0;
  d2Mdy2[2] = 0;

  d2Mdxdy[0] = 0;
  d2Mdxdy[1] = 0;
  d2Mdxdy[2] = 0;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(Scalar *d3Mdx3, Scalar *d3Mdy3, Scalar *d3Mdx2dy, Scalar *d3Mdxdy2, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' third partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs: m         = [ξ, η], the local coordinates of P, and
  //         cs        = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d3Mdx3   = [∂³M₀/∂ξ³, ∂³M₁/∂ξ³, ∂³M₂/∂ξ³]
  //          d3Mdy3   = [∂³M₀/∂η³, ∂³M₁/∂η³, ∂³M₂/∂η³]
  //          d3Mdx2dy = [∂³M₀/∂ξ²∂η, ∂³M₁/∂ξ²∂η, ∂³M₂/∂ξ²∂η]
  //          d3Mdxdy2 = [∂³M₀/∂ξ∂η², ∂³M₁/∂ξ∂η², ∂³M₂/∂ξ∂η²]
  //                     where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  d3Mdx3[0] = 0;
  d3Mdx3[1] = 0;
  d3Mdx3[2] = 0;

  d3Mdy3[0] = 0;
  d3Mdy3[1] = 0;
  d3Mdy3[2] = 0;

  d3Mdx2dy[0] = 0;
  d3Mdx2dy[1] = 0;
  d3Mdx2dy[2] = 0;

  d3Mdxdy2[0] = 0;
  d3Mdxdy2[1] = 0;
  d3Mdxdy2[2] = 0;
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri3::GetJacobian(Scalar *m, CoordSetT &cs)
{
  // J = 2*Area = ||12 x 13||
  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  Scalar V12[3], V13[3];
  V12[0] = X[1]-X[0]; V12[1] = Y[1]-Y[0]; V12[2] = Z[1]-Z[0];
  V13[0] = X[2]-X[0]; V13[1] = Y[2]-Y[0]; V13[2] = Z[2]-Z[0];

  Scalar N[3];
  N[0] = V12[1]*V13[2] - V12[2]*V13[1];
  N[1] = V12[2]*V13[0] - V12[0]*V13[2];
  N[2] = V12[0]*V13[1] - V12[1]*V13[0];

  using std::sqrt;
  return(sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri3::GetShapeFctAndJacobian(Scalar *Shape, Scalar *m, CoordSetT &cs)
{
  Scalar& r = m[0];
  Scalar& s = m[1];
  Scalar t = 1-r-s;

  // !! idem ACME !!
  Shape[0] = r;
  Shape[1] = s;
  Shape[2] = t;

  return(GetJacobian(m, cs));
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GetIsoParamMappingNormalJacobianProduct(Scalar *JNormal, Scalar *m, CoordSetT &cs)
{
  // Computes the vector normal to the surface at a point P whose magnitude is the Jacobian determinant
  // of the local-to-global mapping.

  // Inputs:  m       = [ξ, η], the local coordinates of P, and
  //          cs      = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: JNormal = [∂M₀/∂ξ, ∂M₁/∂ξ, ∂M₂/∂ξ] x [∂M₀/∂η, ∂M₁/∂η, ∂M₂/∂η], the vector normal
  //                    to the surface at point P whose magnitude is the Jacobian determinant of
  //                    the mapping [ξ, η] →- [M₀, M₁, M₂] where M₀, M₁ and M₂ are the
  //                    global X-, Y- and Z-coordinates of P.

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute the normal vector, n = ∂M/∂ξ x ∂M/∂η
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::ComputedJNormaldxAnddJNormaldy(Scalar *dJNormaldx, Scalar *dJNormaldy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the partial derivatives w.r.t the local coordinates of the components of the vector
  // normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m          = [ξ, η], the local coordinates of some point P, and
  //          cs         = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dJNormaldx = [∂n₀/∂ξ, ∂n₁/∂ξ, ∂n₂/∂ξ]
  //          dJNormaldy = [∂n₀/∂η, ∂n₁/∂η, ∂n₂/∂η]
  //                       where n₀, n₁, and n₂ are the global X-, Y- and Z-components of the vector
  //                       normal to the surface at the point P whose magnitude is the Jacobian determinant
  //                       of the mapping [ξ, η] →- [M₀, M₁, M₂] where M₀, M₁ and M₂ are the
  //                       global X-, Y- and Z-coordinates of P.

  dJNormaldx[0] = 0;
  dJNormaldx[1] = 0;
  dJNormaldx[2] = 0;

  dJNormaldy[0] = 0;
  dJNormaldy[1] = 0;
  dJNormaldy[2] = 0;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(Scalar *d2JNormaldx2, Scalar *d2JNormaldy2, Scalar *d2JNormaldxdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the second partial derivatives w.r.t the local coordinates of the components of the vector
  // normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m             = [ξ, η], the local coordinates of some point P, and
  //          cs            = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d2JNormaldx2  = [∂²n₀/∂ξ², ∂²n₁/∂ξ², ∂²n₂/∂ξ²]
  //          d2JNormaldy2  = [∂²n₀/∂η², ∂²n₁/∂η², ∂²n₂/∂η²]
  //          d2JNormaldxdy = [∂²n₀/∂ξ∂η, ∂²n₁/∂ξ∂η, ∂²n₂/∂ξ∂η]
  //                          where n₀, n₁, and n₂ are the global X-, Y- and Z-components of the vector
  //                          normal to the surface at the point P whose magnitude is the Jacobian determinant
  //                          of the mapping [ξ, η] →- [M₀, M₁, M₂] where M₀, M₁ and M₂ are the
  //                          global X-, Y- and Z-coordinates of P.

  d2JNormaldx2[0] = 0;
  d2JNormaldx2[1] = 0;
  d2JNormaldx2[2] = 0;

  d2JNormaldy2[0] = 0;
  d2JNormaldy2[1] = 0;
  d2JNormaldy2[2] = 0;

  d2JNormaldxdy[0] = 0;
  d2JNormaldxdy[1] = 0;
  d2JNormaldxdy[2] = 0;
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceTri3::GetIsoParamMappingNormalAndJacobian(Scalar *Normal, Scalar *m, CoordSetT &cs)
{
  GetIsoParamMappingNormalJacobianProduct(Normal, m, cs);

  using std::sqrt;
  Scalar J = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(J != 0.0) {
    Normal[0] /= J; Normal[1] /= J; Normal[2] /= J;
  }
  return J;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GetdJNormal(Scalar dJNormal[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.

  // Inputs:  m           = [ξ, η], the local coordinates of P, and
  //          cs          = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dJNormal[0] = [∂n₀/∂X₀, ∂n₀/∂Y₀, ∂n₀/∂Z₀, ∂n₀/∂X₁, ∂n₀/∂Y₁, ∂n₀/∂Z₁, ∂n₀/∂X₂, ∂n₀/∂Y₂, ∂n₀/∂Z₂]
  //          dJNormal[1] = [∂n₁/∂X₀, ∂n₁/∂Y₀, ∂n₁/∂Z₀, ∂n₁/∂X₁, ∂n₁/∂Y₁, ∂n₁/∂Z₁, ∂n₁/∂X₂, ∂n₁/∂Y₂, ∂n₁/∂Z₂]
  //          dJNormal[2] = [∂n₂/∂X₀, ∂n₂/∂Y₀, ∂n₂/∂Z₀, ∂n₂/∂X₁, ∂n₂/∂Y₁, ∂n₂/∂Z₁, ∂n₂/∂X₂, ∂n₂/∂Y₂, ∂n₂/∂Z₂]
  //                        where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                        at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                        where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                        and X = [X₀, X₁, X₂], Y = [Y₀, Y₁, Y₂] and  Z = [Z₀, Z₁, Z₂] are
  //                        the nodes' global X-, Y- and Z- coordinates, respectively.

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // Compute edge vectors
  Scalar V12[3] = { X[1] - X[0], Y[1] - Y[0], Z[1] - Z[0] };
  Scalar V13[3] = { X[2] - X[0], Y[2] - Y[0], Z[2] - Z[0] };
  Scalar V23[3] = { X[2] - X[1], Y[2] - Y[1], Z[2] - Z[1] };

  dJNormal[0][0] = 0;
  dJNormal[1][0] = -V23[2];
  dJNormal[2][0] =  V23[1];
  dJNormal[3][0] = 0;
  dJNormal[4][0] =  V13[2];
  dJNormal[5][0] = -V13[1];
  dJNormal[6][0] = 0;
  dJNormal[7][0] = -V12[2];
  dJNormal[8][0] =  V12[1];

  dJNormal[0][1] =  V23[2];
  dJNormal[1][1] = 0;
  dJNormal[2][1] = -V23[0];
  dJNormal[3][1] = -V13[2];
  dJNormal[4][1] = 0;
  dJNormal[5][1] =  V13[0];
  dJNormal[6][1] =  V12[2];
  dJNormal[7][1] = 0;
  dJNormal[8][1] = -V12[0];

  dJNormal[0][2] = -V23[1];
  dJNormal[1][2] =  V23[0];
  dJNormal[2][2] = 0;
  dJNormal[3][2] =  V13[1];
  dJNormal[4][2] = -V13[0];
  dJNormal[5][2] = 0;
  dJNormal[6][2] = -V12[1];
  dJNormal[7][2] =  V12[0];
  dJNormal[8][2] = 0;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::ComputeddJNormaldxAndddJNormaldy(Scalar ddJNormaldx[][3], Scalar ddJNormaldy[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the mixed second partial derivatives of the components of the vector normal to the surface at a point P whose
  // magnitude is the Jacobian determinant of the local-to-global mapping, w.r.t the nodes' global X-, Y- and Z-coordinates and the local
  // coordinates.
  //
  // Inputs:  m                = [ξ, η], the local coordinates of P, and
  //          cs               = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: ddJNormaldx[][0] = [∂²n₀/∂ξ∂X₀, ∂²n₀/∂ξ∂Y₀, ∂²n₀/∂ξ∂Z₀, ∂²n₀/∂ξ∂X₁, ∂²n₀/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ∂²n₀/∂ξ∂X₂, ∂²n₀/∂ξ∂Y₂, ∂²n₀/∂ξ∂Z₂]
  //          ddJNormaldx[][1] = [∂²n₁/∂ξ∂X₀, ∂²n₁/∂ξ∂Y₀, ∂²n₁/∂ξ∂Z₀, ∂²n₁/∂ξ∂X₁, ∂²n₁/∂ξ∂Y₁, ∂²n₁/∂ξ∂Z₁, ∂²n₁/∂ξ∂X₂, ∂²n₁/∂ξ∂Y₂, ∂²n₁/∂ξ∂Z₂]
  //          ddJNormaldx[][2] = [∂²n₂/∂ξ∂X₀, ∂²n₂/∂ξ∂Y₀, ∂²n₂/∂ξ∂Z₀, ∂²n₂/∂ξ∂X₁, ∂²n₂/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ∂²n₂/∂ξ∂X₂, ∂²n₂/∂ξ∂Y₂, ∂²n₂/∂ξ∂Z₂]
  //          ddJNormaldy[][0] = [∂²n₀/∂η∂X₀, ∂²n₀/∂η∂Y₀, ∂²n₀/∂η∂Z₀, ∂²n₀/∂η∂X₁, ∂²n₀/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ∂²n₀/∂η∂X₂, ∂²n₀/∂η∂Y₂, ∂²n₀/∂η∂Z₂]
  //          ddJNormaldy[][1] = [∂²n₁/∂η∂X₀, ∂²n₁/∂η∂Y₀, ∂²n₁/∂η∂Z₀, ∂²n₁/∂η∂X₁, ∂²n₁/∂η∂Y₁, ∂²n₁/∂η∂Z₁, ∂²n₁/∂η∂X₂, ∂²n₁/∂η∂Y₂, ∂²n₁/∂η∂Z₂]
  //          ddJNormaldy[][2] = [∂²n₂/∂η∂X₀, ∂²n₂/∂η∂Y₀, ∂²n₂/∂η∂Z₀, ∂²n₂/∂η∂X₁, ∂²n₂/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ∂²n₂/∂η∂X₂, ∂²n₂/∂η∂Y₂, ∂²n₂/∂η∂Z₂]
  //                             where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                             at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                             where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                             and X = [X₀, X₁, X₂], Y = [Y₀, Y₁, Y₂] and  Z = [Z₀, Z₁, Z₂] are
  //                             the nodes' global X-, Y- and Z- coordinates, respectively.

  for(int i = 0; i < 3; ++i) {
    ddJNormaldx[3*i  ][0] = 0;
    ddJNormaldx[3*i  ][1] = 0;
    ddJNormaldx[3*i  ][2] = 0;
    ddJNormaldx[3*i+1][0] = 0;
    ddJNormaldx[3*i+1][1] = 0;
    ddJNormaldx[3*i+1][2] = 0;
    ddJNormaldx[3*i+2][0] = 0;
    ddJNormaldx[3*i+2][1] = 0;
    ddJNormaldx[3*i+2][2] = 0;
  }

  for(int i = 0; i < 3; ++i) {
    ddJNormaldy[3*i  ][0] = 0;
    ddJNormaldy[3*i  ][1] = 0;
    ddJNormaldy[3*i  ][2] = 0;
    ddJNormaldy[3*i+1][0] = 0;
    ddJNormaldy[3*i+1][1] = 0;
    ddJNormaldy[3*i+1][2] = 0;
    ddJNormaldy[3*i+2][0] = 0;
    ddJNormaldy[3*i+2][1] = 0;
    ddJNormaldy[3*i+2][2] = 0;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::Getd2JNormal(Scalar H[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the second partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m    = [ξ, η], the local coordinates of P, and
  //          cs   = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: H[0] = [∂²n₀/∂X₀²,   ∂²n₀/∂X₀∂Y₀, ∂²n₀/∂X₀∂Z₀, ∂²n₀/∂X₀∂X₁, ∂²n₀/∂X₀∂Y₁, ∂²n₀/∂X₀∂Z₁, ∂²n₀/∂X₀∂X₂, ∂²n₀/∂X₀∂Y₂, ∂²n₀/∂X₀∂Z₂,
  //                  ∂²n₀/∂Y₀∂X₀, ∂²n₀/∂Y₀²,   ∂²n₀/∂Y₀∂Z₀, ∂²n₀/∂Y₀∂X₁, ∂²n₀/∂Y₀∂Y₁, ∂²n₀/∂Y₀∂Z₁, ∂²n₀/∂Y₀∂X₂, ∂²n₀/∂Y₀∂Y₂, ∂²n₀/∂Y₀∂Z₂,
  //                                                                         ⋮ 
  //                  ∂²n₀/∂Z₂∂X₀, ∂²n₀/∂Z₂∂Y₀, ∂²n₀/∂Z₂∂Z₀, ∂²n₀/∂Z₂∂X₁, ∂²n₀/∂Z₂∂Y₁, ∂²n₀/∂Z₂∂Z₁, ∂²n₀/∂Z₂∂X₂, ∂²n₀/∂Z₂∂Y₂, ∂²n₀/∂Z₂²]
  //          H[1] = [∂²n₁/∂X₀²,   ∂²n₁/∂X₀∂Y₀, ∂²n₁/∂X₀∂Z₀, ∂²n₁/∂X₀∂X₁, ∂²n₁/∂X₀∂Y₁, ∂²n₁/∂X₀∂Z₁, ∂²n₁/∂X₀∂X₂, ∂²n₁/∂X₀∂Y₂, ∂²n₁/∂X₀∂Z₂,
  //                  ∂²n₁/∂Y₀∂X₀, ∂²n₁/∂Y₀²,   ∂²n₁/∂Y₀∂Z₀, ∂²n₁/∂Y₀∂X₁, ∂²n₁/∂Y₀∂Y₁, ∂²n₁/∂Y₀∂Z₁, ∂²n₁/∂Y₀∂X₂, ∂²n₁/∂Y₀∂Y₂, ∂²n₁/∂Y₀∂Z₂,
  //                                                                         ⋮ 
  //                  ∂²n₁/∂Z₂∂X₀, ∂²n₁/∂Z₂∂Y₀, ∂²n₁/∂Z₂∂Z₀, ∂²n₁/∂Z₂∂X₁, ∂²n₁/∂Z₂∂Y₁, ∂²n₁/∂Z₂∂Z₁, ∂²n₁/∂Z₂∂X₂, ∂²n₁/∂Z₂∂Y₂, ∂²n₁/∂Z₂²]
  //          H[2] = [∂²n₂/∂X₀²,   ∂²n₂/∂X₀∂Y₀, ∂²n₂/∂X₀∂Z₀, ∂²n₂/∂X₀∂X₁, ∂²n₂/∂X₀∂Y₁, ∂²n₂/∂X₀∂Z₁, ∂²n₂/∂X₀∂X₂, ∂²n₂/∂X₀∂Y₂, ∂²n₂/∂X₀∂Z₂,
  //                  ∂²n₂/∂Y₀∂X₀, ∂²n₂/∂Y₀²,   ∂²n₂/∂Y₀∂Z₀, ∂²n₂/∂Y₀∂X₁, ∂²n₂/∂Y₀∂Y₁, ∂²n₂/∂Y₀∂Z₁, ∂²n₂/∂Y₀∂X₂, ∂²n₂/∂Y₀∂Y₂, ∂²n₂/∂Y₀∂Z₂,
  //                                                                         ⋮ 
  //                  ∂²n₂/∂Z₂∂X₀, ∂²n₂/∂Z₂∂Y₀, ∂²n₂/∂Z₂∂Z₀, ∂²n₂/∂Z₂∂X₁, ∂²n₂/∂Z₂∂Y₁, ∂²n₂/∂Z₂∂Z₁, ∂²n₂/∂Z₂∂X₂, ∂²n₂/∂Z₂∂Y₂, ∂²n₂/∂Z₂²]
  //                 where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                 at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                 where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                 and X = [X₀, X₁, X₂], Y = [Y₀, Y₁, Y₂] and  Z = [Z₀, Z₁, Z₂] are
  //                 the nodes' global X-, Y- and Z- coordinates, respectively.

  Scalar &x1 = cs[Nodes[0]]->x, &y1 = cs[Nodes[0]]->y, &z1 = cs[Nodes[0]]->z;
  Scalar &x2 = cs[Nodes[1]]->x, &y2 = cs[Nodes[1]]->y, &z2 = cs[Nodes[1]]->z;
  Scalar &x3 = cs[Nodes[2]]->x, &y3 = cs[Nodes[2]]->y, &z3 = cs[Nodes[2]]->z;

  // XXX: this is constant w.r.t local and nodal coordinates so it could be a static member variable
  H[0][0] = 0;  H[9 ][0] = 0;  H[18][0] = 0;  H[27][0] = 0;  H[36][0] = 0;  H[45][0] = 0;  H[54][0] = 0;  H[63][0] = 0;  H[72][0] = 0;
  H[1][0] = 0;  H[10][0] = 0;  H[19][0] = 0;  H[28][0] = 0;  H[37][0] = 0;  H[46][0] = 1;  H[55][0] = 0;  H[64][0] = 0;  H[73][0] = -1;
  H[2][0] = 0;  H[11][0] = 0;  H[20][0] = 0;  H[29][0] = 0;  H[38][0] = -1; H[47][0] = 0;  H[56][0] = 0;  H[65][0] = 1;  H[74][0] = 0;
  H[3][0] = 0;  H[12][0] = 0;  H[21][0] = 0;  H[30][0] = 0;  H[39][0] = 0;  H[48][0] = 0;  H[57][0] = 0;  H[66][0] = 0;  H[75][0] = 0;
  H[4][0] = 0;  H[13][0] = 0;  H[22][0] = -1; H[31][0] = 0;  H[40][0] = 0;  H[49][0] = 0;  H[58][0] = 0;  H[67][0] = 0;  H[76][0] = 1;
  H[5][0] = 0;  H[14][0] = 1;  H[23][0] = 0;  H[32][0] = 0;  H[41][0] = 0;  H[50][0] = 0;  H[59][0] = 0;  H[68][0] = -1; H[77][0] = 0;
  H[6][0] = 0;  H[15][0] = 0;  H[24][0] = 0;  H[33][0] = 0;  H[42][0] = 0;  H[51][0] = 0;  H[60][0] = 0;  H[69][0] = 0;  H[78][0] = 0;
  H[7][0] = 0;  H[16][0] = 0;  H[25][0] = 1;  H[34][0] = 0;  H[43][0] = 0;  H[52][0] = -1; H[61][0] = 0;  H[70][0] = 0;  H[79][0] = 0;
  H[8][0] = 0;  H[17][0] = -1; H[26][0] = 0;  H[35][0] = 0;  H[44][0] = 1;  H[53][0] = 0;  H[62][0] = 0;  H[71][0] = 0;  H[80][0] = 0;

  H[0][1] = 0;  H[9 ][1] = 0;  H[18][1] = 0;  H[27][1] = 0;  H[36][1] = 0;  H[45][1] = -1; H[54][1] = 0;  H[63][1] = 0;  H[72][1] = 1;
  H[1][1] = 0;  H[10][1] = 0;  H[19][1] = 0;  H[28][1] = 0;  H[37][1] = 0;  H[46][1] = 0;  H[55][1] = 0;  H[64][1] = 0;  H[73][1] = 0; 
  H[2][1] = 0;  H[11][1] = 0;  H[20][1] = 0;  H[29][1] = 1;  H[38][1] = 0;  H[47][1] = 0;  H[56][1] = -1; H[65][1] = 0;  H[74][1] = 0;
  H[3][1] = 0;  H[12][1] = 0;  H[21][1] = 1;  H[30][1] = 0;  H[39][1] = 0;  H[48][1] = 0;  H[57][1] = 0;  H[66][1] = 0;  H[75][1] = -1;
  H[4][1] = 0;  H[13][1] = 0;  H[22][1] = 0;  H[31][1] = 0;  H[40][1] = 0;  H[49][1] = 0;  H[58][1] = 0;  H[67][1] = 0;  H[76][1] = 0;
  H[5][1] = -1; H[14][1] = 0;  H[23][1] = 0;  H[32][1] = 0;  H[41][1] = 0;  H[50][1] = 0;  H[59][1] = 1;  H[68][1] = 0;  H[77][1] = 0;
  H[6][1] = 0;  H[15][1] = 0;  H[24][1] = -1; H[33][1] = 0;  H[42][1] = 0;  H[51][1] = 1;  H[60][1] = 0;  H[69][1] = 0;  H[78][1] = 0;
  H[7][1] = 0;  H[16][1] = 0;  H[25][1] = 0;  H[34][1] = 0;  H[43][1] = 0;  H[52][1] = 0;  H[61][1] = 0;  H[70][1] = 0;  H[79][1] = 0;
  H[8][1] = 1;  H[17][1] = 0;  H[26][1] = 0;  H[35][1] = -1; H[44][1] = 0;  H[53][1] = 0;  H[62][1] = 0;  H[71][1] = 0;  H[80][1] = 0;

  H[0][2] = 0;  H[9 ][2] = 0;  H[18][2] = 0;  H[27][2] = 0;  H[36][2] = 1;  H[45][2] = 0;  H[54][2] = 0;  H[63][2] = -1; H[72][2] = 0;
  H[1][2] = 0;  H[10][2] = 0;  H[19][2] = 0;  H[28][2] = -1; H[37][2] = 0;  H[46][2] = 0;  H[55][2] = 1;  H[64][2] = 0;  H[73][2] = 0; 
  H[2][2] = 0;  H[11][2] = 0;  H[20][2] = 0;  H[29][2] = 0;  H[38][2] = 0;  H[47][2] = 0;  H[56][2] = 0;  H[65][2] = 0;  H[74][2] = 0;
  H[3][2] = 0;  H[12][2] = -1; H[21][2] = 0;  H[30][2] = 0;  H[39][2] = 0;  H[48][2] = 0;  H[57][2] = 0;  H[66][2] = 1;  H[75][2] = 0;
  H[4][2] = 1;  H[13][2] = 0;  H[22][2] = 0;  H[31][2] = 0;  H[40][2] = 0;  H[49][2] = 0;  H[58][2] = -1; H[67][2] = 0;  H[76][2] = 0;
  H[5][2] = 0;  H[14][2] = 0;  H[23][2] = 0;  H[32][2] = 0;  H[41][2] = 0;  H[50][2] = 0;  H[59][2] = 0;  H[68][2] = 0;  H[77][2] = 0;
  H[6][2] = 0;  H[15][2] = 1;  H[24][2] = 0;  H[33][2] = 0;  H[42][2] = -1; H[51][2] = 0;  H[60][2] = 0;  H[69][2] = 0;  H[78][2] = 0;
  H[7][2] = -1; H[16][2] = 0;  H[25][2] = 0;  H[34][2] = 1;  H[43][2] = 0;  H[52][2] = 0;  H[61][2] = 0;  H[70][2] = 0;  H[79][2] = 0;
  H[8][2] = 0;  H[17][2] = 0;  H[26][2] = 0;  H[35][2] = 0;  H[44][2] = 0;  H[53][2] = 0;  H[62][2] = 0;  H[71][2] = 0;  H[80][2] = 0;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GetUnitNormal(Scalar Normal[3], Scalar *m, CoordSetT &cs)
{
  // Computes the unit vector normal to the surface at a point P.

  GetIsoParamMappingNormalJacobianProduct(Normal, m, cs);

  using std::sqrt;
  Scalar J = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(J != 0.0) {
    Normal[0] /= J; Normal[1] /= J; Normal[2] /= J;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GetdUnitNormal(Scalar dNormal[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the partial derivatives w.r.t. the nodal X-, Y- and Z-coordinates of the unit vector normal to the surface at a point P.
  //
  // ñ₀ = n₀/sqrt(n₀n₀+n₁n₁+n₂n₂) = n₀/J
  // ∂ñ₀/∂X₀ = (∂n₀/∂X₀)/sqrt(n₀n₀+n₁n₁+n₂n₂) + n₀*∂[1/sqrt(n₀n₀+n₁n₁+n₂n₂)]/∂X₀
  //         = (∂n₀/∂X₀)/J + n₀*(-1/2)*pow(n₀n₀+n₁n₁+n₂n₂,-3/2)*(2n₀*∂n₀/∂X₀+2n₁*∂n₁/∂X₀+2n₂*∂n₂/∂X₀)
  //         = (∂n₀/∂X₀)/J - n₀/J³*(n₀*∂n₀/∂X₀+n₁*∂n₁/∂X₀+n₂*∂n₂/∂X₀)
  //         = (∂n₀/∂X₀)/J - n₀/J³*(n·∂n/∂X₀)
  // Similarly,
  // ∂ñ₁/∂X₀ = (∂n₁/∂X₀)/J - n₁/J³*(n·∂n/∂X₀)
  // ∂ñ₂/∂X₀ = (∂n₂/∂X₀)/J - n₂/J³*(n·∂n/∂X₀)

  // Compute the normal, n
  Scalar JNormal[3];
  GetIsoParamMappingNormalJacobianProduct(JNormal, m, cs);

  // Compute the partial derivatives of n
  Scalar dJNormal[9][3];
  GetdJNormal(dJNormal, m, cs);

  // Compute the partial derivatives of the unit vector, ñ
  Scalar J2 = JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2];
  if(J2 != 0) {
    using std::sqrt;
    Scalar J = sqrt(J2);
    Scalar J3 = J*J2;
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        Scalar v = 1/J3*(JNormal[0]*dJNormal[3*i+j][0] + JNormal[1]*dJNormal[3*i+j][1] + JNormal[2]*dJNormal[3*i+j][2]);
        dNormal[3*i+j][0] = dJNormal[3*i+j][0]/J - JNormal[0]*v;
        dNormal[3*i+j][1] = dJNormal[3*i+j][1]/J - JNormal[1]*v;
        dNormal[3*i+j][2] = dJNormal[3*i+j][2]/J - JNormal[2]*v;
      }
    }
  }
  else {
    for(int i = 0; i < 9; ++i)
      for(int j = 0; j < 3; ++j) 
        dNormal[i][j] = dJNormal[i][j];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::Getd2UnitNormal(Scalar d2Normal[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the second partial derivatives w.r.t. the nodal X-, Y- and Z-coordinates of the unit vector normal to the surface at a point P.
  //
  // ñ₀ = n₀/sqrt(n₀n₀+n₁n₁+n₂n₂) = n₀/J
  // ∂ñ₀/∂X₀ = (∂n₀/∂X₀)/sqrt(n₀n₀+n₁n₁+n₂n₂) + n₀*∂[1/sqrt(n₀n₀+n₁n₁+n₂n₂)]/∂X₀
  //         = (∂n₀/∂X₀)/J + n₀*(-1/2)*pow(n₀n₀+n₁n₁+n₂n₂,-3/2)*(2n₀*∂n₀/∂X₀+2n₁*∂n₁/∂X₀+2n₂*∂n₂/∂X₀)
  //         = (∂n₀/∂X₀)/J - n₀/J³*(n₀*∂n₀/∂X₀+n₁*∂n₁/∂X₀+n₂*∂n₂/∂X₀)
  //         = (∂n₀/∂X₀)/J - n₀/J³*(n·∂n/∂X₀)
  // ∂²ñ₀/∂X₀² = (∂²n₀/∂X₀²)/J - (∂n₀/∂X₀)/J³*(n·∂n/∂X₀) - (∂n₀/∂X₀)/J³*(n·∂n/∂X₀) - n₀*(∂[1/J³*(n·∂n/∂X₀)]/∂X₀)
  //           = (∂²n₀/∂X₀²)/J - 2*(∂n₀/∂X₀)/J³*(n·∂n/∂X₀) - n₀*( (∂[1/J³]/∂X₀)*(n·∂n/∂X₀) + 1/J³*(∂[n·∂n/∂X₀]/∂X₀) )
  //           = (∂²n₀/∂X₀²)/J - 2*(∂n₀/∂X₀)/J³*(n·∂n/∂X₀) - n₀*( 3*(∂[1/J]/∂X₀)*1/J²*(n·∂n/∂X₀) + 1/J³*(∂n/∂X₀·∂n/∂X₀ + n·∂²n/∂X₀²) )
  //           = (∂²n₀/∂X₀²)/J - 2*(∂n₀/∂X₀)/J³*(n·∂n/∂X₀) - n₀*( 3*(-1/J³*(n·∂n/∂X₀))*1/J²*(n·∂n/∂X₀) + 1/J³*(∂n/∂X₀·∂n/∂X₀ + n·∂²n/∂X₀²) )
  //           = (∂²n₀/∂X₀²)/J - 2*(∂n₀/∂X₀)/J³*(n·∂n/∂X₀) + n₀*(3/J⁵*(n·∂n/∂X₀)² - 1/J³*(∂n/∂X₀·∂n/∂X₀ + n·∂²n/∂X₀²))

  // Similarly,
  // ∂²ñ₀/∂X₀∂Y₀ = (∂²n₀/∂X₀∂Y₀)/J - (∂n₀/∂X₀)/J³*(n·∂n/∂Y₀) - (∂n₀/∂Y₀)/J³*(n·∂n/∂X₀) + n₀*(3/J⁵*(n·∂n/∂X₀)(n·∂n/∂Y₀) - 1/J³*(∂n/∂X₀·∂n/∂Y₀ + n·∂²n/∂X₀∂Y₀))
  // Compute the normal, n
  Scalar JNormal[3];
  GetIsoParamMappingNormalJacobianProduct(JNormal, m, cs);

  // Compute the partial derivatives of n
  Scalar dJNormal[9][3];
  GetdJNormal(dJNormal, m, cs);

  // Compute the second partial derivatives of n
  Scalar d2JNormal[81][3];
  Getd2JNormal(d2JNormal, m, cs);

  // Compute the partial derivatives of the unit vector, ñ
  Scalar J2 = JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2];
  if(J2 != 0) {
    using std::sqrt;
    Scalar J = sqrt(J2);
    Scalar J3 = J*J2;
    Scalar threeJ = 3*J;

    Scalar v[9];
    for(int i = 0; i < 9; ++i) v[i] = 1/J3*(JNormal[0]*dJNormal[i][0] + JNormal[1]*dJNormal[i][1] + JNormal[2]*dJNormal[i][2]);

    for(int i = 0; i < 9; ++i) {
      for(int j = i; j < 9; ++j) {
        Scalar a = dJNormal[i][0]*dJNormal[j][0] + dJNormal[i][1]*dJNormal[j][1] + dJNormal[i][2]*dJNormal[j][2];
        Scalar b = JNormal[0]*d2JNormal[9*i+j][0]+JNormal[1]*d2JNormal[9*i+j][1]+JNormal[2]* d2JNormal[9*i+j][2];
        Scalar c = threeJ*v[i]*v[j] - (a+b)/J3;
        d2Normal[9*i+j][0] = d2Normal[9*j+i][0] = d2JNormal[9*i+j][0]/J - dJNormal[j][0]*v[i] - dJNormal[i][0]*v[j] + JNormal[0]*c;
        d2Normal[9*i+j][1] = d2Normal[9*j+i][1] = d2JNormal[9*i+j][1]/J - dJNormal[j][1]*v[i] - dJNormal[i][1]*v[j] + JNormal[1]*c;
        d2Normal[9*i+j][2] = d2Normal[9*j+i][2] = d2JNormal[9*i+j][2]/J - dJNormal[j][2]*v[i] - dJNormal[i][2]*v[j] + JNormal[2]*c;
      }
    }

  }
  else {
    for(int i = 0; i < 81; ++i)
      for(int j = 0; j < 3; ++j)
        d2Normal[i][j] = d2JNormal[i][j];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GlobalToLocalCoord(Scalar *m, Scalar *M, CoordSetT &cs)
{
  // Computes the local coordinates of a point P.
  //
  // Inputs:  M  = [M₀, M₁, M₂], the global X-, Y- and Z-coordinates of P, and
  //          cs = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: m  = [ξ(M₀,M₁,M₂,X,Y,Z), η(M₀,M₁,M₂,X,Y,Z)], the local coordinates of P.

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // Compute the edge vectors
  Scalar V12[3] = { X[1] - X[0], Y[1] - Y[0], Z[1] - Z[0] };
  Scalar V13[3] = { X[2] - X[0], Y[2] - Y[0], Z[2] - Z[0] };

  // Compute the normal (as V12 x V13)
  Scalar N[3];
  N[0] = V12[1]*V13[2] - V12[2]*V13[1];
  N[1] = V12[2]*V13[0] - V12[0]*V13[2];
  N[2] = V12[0]*V13[1] - V12[1]*V13[0];

  // Compute vectors from P to each node
  Scalar VP1[3] = { X[0] - M[0], Y[0] - M[1], Z[0] - M[2] };
  Scalar VP2[3] = { X[1] - M[0], Y[1] - M[1], Z[1] - M[2] };
  Scalar VP3[3] = { X[2] - M[0], Y[2] - M[1], Z[2] - M[2] };

  // Compute the areas as the scalar triple products: a0 = (VP2 x VP3)·N, a1 = (VP3 x VP1)·N, a2 = (VP1 x VP2)·N
  Scalar a0 = (VP2[1]*VP3[2] - VP2[2]*VP3[1])*N[0] + (VP2[2]*VP3[0] - VP2[0]*VP3[2])*N[1] + (VP2[0]*VP3[1] - VP2[1]*VP3[0])*N[2];
  Scalar a1 = (VP3[1]*VP1[2] - VP3[2]*VP1[1])*N[0] + (VP3[2]*VP1[0] - VP3[0]*VP1[2])*N[1] + (VP3[0]*VP1[1] - VP3[1]*VP1[0])*N[2];
  Scalar a2 = (VP1[1]*VP2[2] - VP1[2]*VP2[1])*N[0] + (VP1[2]*VP2[0] - VP1[0]*VP2[2])*N[1] + (VP1[0]*VP2[1] - VP1[1]*VP2[0])*N[2];

  // Compute the local coordinates
  Scalar A = a0 + a1 + a2;
  m[0] = a0/A;
  m[1] = a1/A;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::ComputedmdXdmdYAnddmdZ(Scalar *dmdX, Scalar *dmdY, Scalar *dmdZ, Scalar *M, CoordSetT &cs)
{
  // Computes the values of the local coordinates' partial derivatives w.r.t the global coordinates at a point P.
  //
  // Inputs:  M     = [M₀, M₁, M₂], the global X-, Y- and Z-coordinates of P, and
  //          cs    = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dmdX  = [∂ξ/∂M₀, ∂η/∂M₀]
  //          dmdY  = [∂ξ/∂M₁, ∂η/∂M₁]
  //          dmdZ  = [∂ξ/∂M₂, ∂η/∂M₂]
  //                  where ξ, η are the local coordinates of P.

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // Compute edge vectors
  Scalar V12[3] = { X[1] - X[0], Y[1] - Y[0], Z[1] - Z[0] };
  Scalar V13[3] = { X[2] - X[0], Y[2] - Y[0], Z[2] - Z[0] };
  Scalar V23[3] = { X[2] - X[1], Y[2] - Y[1], Z[2] - Z[1] };

  // Compute the normal (as V12 x V13)
  Scalar N[3];
  N[0] = V12[1]*V13[2] - V12[2]*V13[1];
  N[1] = V12[2]*V13[0] - V12[0]*V13[2];
  N[2] = V12[0]*V13[1] - V12[1]*V13[0];

  // Compute vectors from P to each node
  Scalar VP1[3] = { X[0] - M[0], Y[0] - M[1], Z[0] - M[2] };
  Scalar VP2[3] = { X[1] - M[0], Y[1] - M[1], Z[1] - M[2] };
  Scalar VP3[3] = { X[2] - M[0], Y[2] - M[1], Z[2] - M[2] };

  // Compute the areas as the scalar triple products: a0 = (VP2 x VP3)·N, a1 = (VP3 x VP1)·N, a2 = (VP1 x VP2)·N
  Scalar a0 = (VP2[1]*VP3[2] - VP2[2]*VP3[1])*N[0] + (VP2[2]*VP3[0] - VP2[0]*VP3[2])*N[1] + (VP2[0]*VP3[1] - VP2[1]*VP3[0])*N[2];
  Scalar a1 = (VP3[1]*VP1[2] - VP3[2]*VP1[1])*N[0] + (VP3[2]*VP1[0] - VP3[0]*VP1[2])*N[1] + (VP3[0]*VP1[1] - VP3[1]*VP1[0])*N[2];
  Scalar a2 = (VP1[1]*VP2[2] - VP1[2]*VP2[1])*N[0] + (VP1[2]*VP2[0] - VP1[0]*VP2[2])*N[1] + (VP1[0]*VP2[1] - VP1[1]*VP2[0])*N[2];
  Scalar A = a0 + a1 + a2;

  // Compute the partial derivatives of the local coordinates w.r.t. the global coordinates at P
  dmdX[0] = (V23[2]*N[1] - V23[1]*N[2])/A;
  dmdX[1] = (V13[1]*N[2] - V13[2]*N[1])/A;
  dmdY[0] = (V23[0]*N[2] - V23[2]*N[0])/A;
  dmdY[1] = (V13[2]*N[0] - V13[0]*N[2])/A;
  dmdZ[0] = (V23[1]*N[0] - V23[0]*N[1])/A;
  dmdZ[1] = (V13[0]*N[1] - V13[1]*N[0])/A;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::Computed2mdX2d2mdY2Etc(Scalar *d2mdX2, Scalar *d2mdY2, Scalar *d2mdZ2, Scalar *d2mdXdY, Scalar *d2mdYdZ, Scalar *d2mdXdZ, Scalar *M, CoordSetT &cs)
{
  // Computes the values of the local coordinates' second partial derivatives w.r.t the global coordinates at a point P.
  //
  // Inputs:  M       = [M₀, M₁, M₂], the global X-, Y- and Z-coordinates of P, and
  //          cs      = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d2mdX2  = [∂²ξ/∂M₀, ∂²η/∂M₀²]
  //          d2mdY2  = [∂²ξ/∂M₁, ∂²η/∂M₁²]
  //          d2mdZ2  = [∂²ξ/∂M₂, ∂²η/∂M₂²]
  //          d2mdXdY = [∂²ξ/∂M₀∂M₁, ∂²η/∂M₀∂M₁]
  //          d2mdYdZ = [∂²ξ/∂M₁∂M₂, ∂²η/∂M₁∂M₂]
  //          d2mdXdZ = [∂²ξ/∂M₀∂M₂, ∂²η/∂M₀∂M₂]
  //                    where ξ, η are the local coordinates of P.

  d2mdX2[0] = d2mdX2[1] = 0;
  d2mdY2[0] = d2mdY2[1] = 0;
  d2mdZ2[0] = d2mdZ2[1] = 0;
  d2mdXdY[0] = d2mdXdY[1] = 0;
  d2mdYdZ[0] = d2mdYdZ[1] = 0;
  d2mdXdZ[0] = d2mdXdZ[1] = 0;
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GetdLocalCoords(Scalar dLocalCoords[][2], Scalar *M, CoordSetT &cs)
{
  // Computes the partial derivatives w.r.t. the nodal X-, Y- and Z-coordinates of the local coordinates at a point P.
  //
  // Inputs:  M     = [M₀, M₁, M₂], the global X-, Y- and Z-coordinates of P, and
  //          cs    = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dLocalCoords[0] = [∂ξ/∂X₀, ∂ξ/∂Y₀, ∂ξ/∂Z₀, ∂ξ/∂X₁, ∂ξ/∂Y₁, ∂ξ/∂Z₁, ∂ξ/∂X₂, ∂ξ/∂Y₂, ∂ξ/∂Z₂]
  //          dLocalCoords[1] = [∂η/∂X₀, ∂η/∂Y₀, ∂η/∂Z₀, ∂η/∂X₁, ∂η/∂Y₁, ∂η/∂Z₁, ∂η/∂X₂, ∂η/∂Y₂, ∂η/∂Z₂]
  //                            where ξ, η are the local coordinates of P,
  //                            and X = [X₀, X₁, X₂], Y = [Y₀, Y₁, Y₂] and  Z = [Z₀, Z₁, Z₂] are
  //                            the nodes' global X-, Y- and Z- coordinates, respectively.

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // Compute edge vectors
  Scalar V12[3] = { X[1] - X[0], Y[1] - Y[0], Z[1] - Z[0] };
  Scalar V13[3] = { X[2] - X[0], Y[2] - Y[0], Z[2] - Z[0] };
  Scalar V23[3] = { X[2] - X[1], Y[2] - Y[1], Z[2] - Z[1] };

  // Compute the normal (as V12 x V13)
  Scalar N[3];
  N[0] = V12[1]*V13[2] - V12[2]*V13[1];
  N[1] = V12[2]*V13[0] - V12[0]*V13[2];
  N[2] = V12[0]*V13[1] - V12[1]*V13[0];

  // Compute vectors from P to each node
  Scalar VP1[3] = { X[0] - M[0], Y[0] - M[1], Z[0] - M[2] };
  Scalar VP2[3] = { X[1] - M[0], Y[1] - M[1], Z[1] - M[2] };
  Scalar VP3[3] = { X[2] - M[0], Y[2] - M[1], Z[2] - M[2] };

  // Compute the cross products VP2 x VP3, VP3 x VP1 and VP1 x VP2
  Scalar VP2xVP3[3] = { VP2[1]*VP3[2] - VP2[2]*VP3[1], VP2[2]*VP3[0] - VP2[0]*VP3[2], VP2[0]*VP3[1] - VP2[1]*VP3[0] };
  Scalar VP3xVP1[3] = { VP3[1]*VP1[2] - VP3[2]*VP1[1], VP3[2]*VP1[0] - VP3[0]*VP1[2], VP3[0]*VP1[1] - VP3[1]*VP1[0] };
  Scalar VP1xVP2[3] = { VP1[1]*VP2[2] - VP1[2]*VP2[1], VP1[2]*VP2[0] - VP1[0]*VP2[2], VP1[0]*VP2[1] - VP1[1]*VP2[0] };

  // Compute the areas as the scalar triple products: a0 = (VP2 x VP3)·N, a1 = (VP3 x VP1)·N, a2 = (VP1 x VP2)·N
  Scalar a0 = VP2xVP3[0]*N[0] + VP2xVP3[1]*N[1] + VP2xVP3[2]*N[2];
  Scalar a1 = VP3xVP1[0]*N[0] + VP3xVP1[1]*N[1] + VP3xVP1[2]*N[2];
  Scalar a2 = VP1xVP2[0]*N[0] + VP1xVP2[1]*N[1] + VP1xVP2[2]*N[2];
  Scalar A = a0+a1+a2;

  Scalar da0[9] = { VP2xVP3[1]*V23[2]               - VP2xVP3[2]*V23[1]              ,
                   -VP2xVP3[0]*V23[2]               + VP2xVP3[2]*V23[0]              ,
                    VP2xVP3[0]*V23[1]               - VP2xVP3[1]*V23[0]              ,
                   -VP2xVP3[1]*V13[2] - VP3[2]*N[1] + VP2xVP3[2]*V13[1] + VP3[1]*N[2],
                    VP2xVP3[0]*V13[2] + VP3[2]*N[0] - VP2xVP3[2]*V13[0] - VP3[0]*N[2],
                   -VP2xVP3[0]*V13[1] - VP3[1]*N[0] + VP2xVP3[1]*V13[0] + VP3[0]*N[1],
                    VP2xVP3[1]*V12[2] + VP2[2]*N[1] - VP2xVP3[2]*V12[1] - VP2[1]*N[2],
                   -VP2xVP3[0]*V12[2] - VP2[2]*N[0] + VP2xVP3[2]*V12[0] + VP2[0]*N[2],
                    VP2xVP3[0]*V12[1] + VP2[1]*N[0] - VP2xVP3[1]*V12[0] - VP2[0]*N[1] };

  Scalar da1[9] = { VP3xVP1[1]*V23[2] + VP3[2]*N[1] - VP3xVP1[2]*V23[1] - VP3[1]*N[2],
                   -VP3xVP1[0]*V23[2] - VP3[2]*N[0] + VP3xVP1[2]*V23[0] + VP3[0]*N[2],
                    VP3xVP1[0]*V23[1] + VP3[1]*N[0] - VP3xVP1[1]*V23[0] - VP3[0]*N[1],
                   -VP3xVP1[1]*V13[2]               + VP3xVP1[2]*V13[1]              ,
                    VP3xVP1[0]*V13[2]               - VP3xVP1[2]*V13[0]              ,
                   -VP3xVP1[0]*V13[1]               + VP3xVP1[1]*V13[0]              ,
                    VP3xVP1[1]*V12[2] - VP1[2]*N[1] - VP3xVP1[2]*V12[1] + VP1[1]*N[2],
                   -VP3xVP1[0]*V12[2] + VP1[2]*N[0] + VP3xVP1[2]*V12[0] - VP1[0]*N[2],
                    VP3xVP1[0]*V12[1] - VP1[1]*N[0] - VP3xVP1[1]*V12[0] + VP1[0]*N[1] };

  Scalar da2[9] = { VP1xVP2[1]*V23[2] - VP2[2]*N[1] - VP1xVP2[2]*V23[1] + VP2[1]*N[2],
                   -VP1xVP2[0]*V23[2] + VP2[2]*N[0] + VP1xVP2[2]*V23[0] - VP2[0]*N[2],
                    VP1xVP2[0]*V23[1] - VP2[1]*N[0] - VP1xVP2[1]*V23[0] + VP2[0]*N[1],
                   -VP1xVP2[1]*V13[2] + VP1[2]*N[1] + VP1xVP2[2]*V13[1] - VP1[1]*N[2],
                    VP1xVP2[0]*V13[2] - VP1[2]*N[0] - VP1xVP2[2]*V13[0] + VP1[0]*N[2],
                   -VP1xVP2[0]*V13[1] + VP1[1]*N[0] + VP1xVP2[1]*V13[0] - VP1[0]*N[1],
                    VP1xVP2[1]*V12[2]               - VP1xVP2[2]*V12[1]              ,
                   -VP1xVP2[0]*V12[2]               + VP1xVP2[2]*V12[0]              ,
                    VP1xVP2[0]*V12[1]               - VP1xVP2[1]*V12[0]               };

  Scalar A2 = A*A;
  Scalar a12 = a1+a2, a02 = a0+a2;
  for(int i=0; i<9; ++i) {
    dLocalCoords[i][0] = (da0[i]*a12 - a0*(da1[i]+da2[i]))/A2;
    dLocalCoords[i][1] = (da1[i]*a02 - a1*(da0[i]+da2[i]))/A2;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ(Scalar ddLocalCoordsdX[][2], Scalar ddLocalCoordsdY[][2], Scalar ddLocalCoordsdZ[][2],
                                                              Scalar *M, CoordSetT &cs)
{
  // Computes the values of the mixed second partial derivatives of the local coordinates at a point P w.r.t the nodes' global
  // X-, Y- and Z- coordinates and the global X-, Y- and Z- coordinates of P.
  //
  // Inputs:  M                    = [M₀, M₁, M₂], the global X-, Y- and Z-coordinates of P, and
  //          cs                   = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: ddLocalCoordsdX[][0] = [∂²ξ/∂M₀∂X₀, ∂²ξ/∂M₀∂Y₀, ∂²ξ/∂M₀∂Z₀, ∂²ξ/∂M₀∂X₁, ∂²ξ/∂M₀∂Y₁, ∂²ξ/∂M₀∂Z₁, ∂²ξ/∂M₀∂X₂, ∂²ξ/∂M₀∂Y₂, ∂²ξ/∂M₀∂Z₂]
  //          ddLocalCoordsdX[][1] = [∂²η/∂M₀∂X₀, ∂²η/∂M₀∂Y₀, ∂²η/∂M₀∂Z₀, ∂²η/∂M₀∂X₁, ∂²η/∂M₀∂Y₁, ∂²η/∂M₀∂Z₁, ∂²η/∂M₀∂X₂, ∂²η/∂M₀∂Y₂, ∂²η/∂M₀∂Z₂]
  //          ddLocalCoordsdY[][0] = [∂²ξ/∂M₁∂X₀, ∂²ξ/∂M₁∂Y₀, ∂²ξ/∂M₁∂Z₀, ∂²ξ/∂M₁∂X₁, ∂²ξ/∂M₁∂Y₁, ∂²ξ/∂M₁∂Z₁, ∂²ξ/∂M₁∂X₂, ∂²ξ/∂M₁∂Y₂, ∂²ξ/∂M₁∂Z₂]
  //          ddLocalCoordsdY[][1] = [∂²η/∂M₁∂X₀, ∂²η/∂M₁∂Y₀, ∂²η/∂M₁∂Z₀, ∂²η/∂M₁∂X₁, ∂²η/∂M₁∂Y₁, ∂²η/∂M₁∂Z₁, ∂²η/∂M₁∂X₂, ∂²η/∂M₁∂Y₂, ∂²η/∂M₁∂Z₂]
  //          ddLocalCoordsdZ[][0] = [∂²ξ/∂M₂∂X₀, ∂²ξ/∂M₂∂Y₀, ∂²ξ/∂M₂∂Z₀, ∂²ξ/∂M₂∂X₁, ∂²ξ/∂M₂∂Y₁, ∂²ξ/∂M₂∂Z₁, ∂²ξ/∂M₂∂X₂, ∂²ξ/∂M₂∂Y₂, ∂²ξ/∂M₂∂Z₂]
  //          ddLocalCoordsdZ[][1] = [∂²η/∂M₂∂X₀, ∂²η/∂M₂∂Y₀, ∂²η/∂M₂∂Z₀, ∂²η/∂M₂∂X₁, ∂²η/∂M₂∂Y₁, ∂²η/∂M₂∂Z₁, ∂²η/∂M₂∂X₂, ∂²η/∂M₂∂Y₂, ∂²η/∂M₂∂Z₂]
  //                                 where ξ, η are the local coordinates of P,
  //                                 and X = [X₀, X₁, X₂], Y = [Y₀, Y₁, Y₂] and  Z = [Z₀, Z₁, Z₂] are
  //                                 the nodes' global X-, Y- and Z- coordinates, respectively.

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // Compute edge vectors
  Scalar V12[3] = { X[1] - X[0], Y[1] - Y[0], Z[1] - Z[0] };
  Scalar V13[3] = { X[2] - X[0], Y[2] - Y[0], Z[2] - Z[0] };
  Scalar V23[3] = { X[2] - X[1], Y[2] - Y[1], Z[2] - Z[1] };

  // Compute the normal (as V12 x V13)
  Scalar N[3];
  N[0] = V12[1]*V13[2] - V12[2]*V13[1];
  N[1] = V12[2]*V13[0] - V12[0]*V13[2];
  N[2] = V12[0]*V13[1] - V12[1]*V13[0];

  // Compute vectors from P to each node
  Scalar VP1[3] = { X[0] - M[0], Y[0] - M[1], Z[0] - M[2] };
  Scalar VP2[3] = { X[1] - M[0], Y[1] - M[1], Z[1] - M[2] };
  Scalar VP3[3] = { X[2] - M[0], Y[2] - M[1], Z[2] - M[2] };

  // Compute the cross products VP2 x VP3, VP3 x VP1 and VP1 x VP2
  Scalar VP2xVP3[3] = { VP2[1]*VP3[2] - VP2[2]*VP3[1], VP2[2]*VP3[0] - VP2[0]*VP3[2], VP2[0]*VP3[1] - VP2[1]*VP3[0] };
  Scalar VP3xVP1[3] = { VP3[1]*VP1[2] - VP3[2]*VP1[1], VP3[2]*VP1[0] - VP3[0]*VP1[2], VP3[0]*VP1[1] - VP3[1]*VP1[0] };
  Scalar VP1xVP2[3] = { VP1[1]*VP2[2] - VP1[2]*VP2[1], VP1[2]*VP2[0] - VP1[0]*VP2[2], VP1[0]*VP2[1] - VP1[1]*VP2[0] };

  // Compute the partial derivatives of these cross products w.r.t. the global coordinates of P
  Scalar dVP2xVP3dX[3] = { 0, -VP2[2] + VP3[2], -VP3[1] + VP2[1] };
  Scalar dVP2xVP3dY[3] = { -VP3[2] + VP2[2], 0, -VP2[0] + VP3[0] };
  Scalar dVP2xVP3dZ[3] = { -VP2[1] + VP3[1], -VP3[0] + VP2[0], 0 };

  Scalar dVP3xVP1dX[3] = { 0, -VP3[2] + VP1[2], -VP1[1] + VP3[1] };
  Scalar dVP3xVP1dY[3] = { -VP1[2] + VP3[2], 0, -VP3[0] + VP1[0] };
  Scalar dVP3xVP1dZ[3] = { -VP3[1] + VP1[1], -VP1[0] + VP3[0], 0 };

  Scalar dVP1xVP2dX[3] = { 0, -VP1[2] + VP2[2], -VP2[1] + VP1[1] };
  Scalar dVP1xVP2dY[3] = { -VP2[2] + VP1[2], 0, -VP1[0] + VP2[0] };
  Scalar dVP1xVP2dZ[3] = { -VP1[1] + VP2[1], -VP2[0] + VP1[0], 0 };

  // Compute the areas as the scalar triple products: a0 = (VP2 x VP3)·N, a1 = (VP3 x VP1)·N, a2 = (VP1 x VP2)·N
  Scalar a0 = VP2xVP3[0]*N[0] + VP2xVP3[1]*N[1] + VP2xVP3[2]*N[2];
  Scalar a1 = VP3xVP1[0]*N[0] + VP3xVP1[1]*N[1] + VP3xVP1[2]*N[2];
  Scalar a2 = VP1xVP2[0]*N[0] + VP1xVP2[1]*N[1] + VP1xVP2[2]*N[2];
  Scalar A = a0+a1+a2;

  // Compute the partial derivatives of the areas w.r.t. the nodal coordinates
  Scalar da0[9] = { VP2xVP3[1]*V23[2]               - VP2xVP3[2]*V23[1]              ,
                   -VP2xVP3[0]*V23[2]               + VP2xVP3[2]*V23[0]              ,
                    VP2xVP3[0]*V23[1]               - VP2xVP3[1]*V23[0]              ,
                   -VP2xVP3[1]*V13[2] - VP3[2]*N[1] + VP2xVP3[2]*V13[1] + VP3[1]*N[2],
                    VP2xVP3[0]*V13[2] + VP3[2]*N[0] - VP2xVP3[2]*V13[0] - VP3[0]*N[2],
                   -VP2xVP3[0]*V13[1] - VP3[1]*N[0] + VP2xVP3[1]*V13[0] + VP3[0]*N[1],
                    VP2xVP3[1]*V12[2] + VP2[2]*N[1] - VP2xVP3[2]*V12[1] - VP2[1]*N[2],
                   -VP2xVP3[0]*V12[2] - VP2[2]*N[0] + VP2xVP3[2]*V12[0] + VP2[0]*N[2],
                    VP2xVP3[0]*V12[1] + VP2[1]*N[0] - VP2xVP3[1]*V12[0] - VP2[0]*N[1] };

  Scalar da1[9] = { VP3xVP1[1]*V23[2] + VP3[2]*N[1] - VP3xVP1[2]*V23[1] - VP3[1]*N[2],
                   -VP3xVP1[0]*V23[2] - VP3[2]*N[0] + VP3xVP1[2]*V23[0] + VP3[0]*N[2],
                    VP3xVP1[0]*V23[1] + VP3[1]*N[0] - VP3xVP1[1]*V23[0] - VP3[0]*N[1],
                   -VP3xVP1[1]*V13[2]               + VP3xVP1[2]*V13[1]              ,
                    VP3xVP1[0]*V13[2]               - VP3xVP1[2]*V13[0]              ,
                   -VP3xVP1[0]*V13[1]               + VP3xVP1[1]*V13[0]              ,
                    VP3xVP1[1]*V12[2] - VP1[2]*N[1] - VP3xVP1[2]*V12[1] + VP1[1]*N[2],
                   -VP3xVP1[0]*V12[2] + VP1[2]*N[0] + VP3xVP1[2]*V12[0] - VP1[0]*N[2],
                    VP3xVP1[0]*V12[1] - VP1[1]*N[0] - VP3xVP1[1]*V12[0] + VP1[0]*N[1] };

  Scalar da2[9] = { VP1xVP2[1]*V23[2] - VP2[2]*N[1] - VP1xVP2[2]*V23[1] + VP2[1]*N[2],
                   -VP1xVP2[0]*V23[2] + VP2[2]*N[0] + VP1xVP2[2]*V23[0] - VP2[0]*N[2],
                    VP1xVP2[0]*V23[1] - VP2[1]*N[0] - VP1xVP2[1]*V23[0] + VP2[0]*N[1],
                   -VP1xVP2[1]*V13[2] + VP1[2]*N[1] + VP1xVP2[2]*V13[1] - VP1[1]*N[2],
                    VP1xVP2[0]*V13[2] - VP1[2]*N[0] - VP1xVP2[2]*V13[0] + VP1[0]*N[2],
                   -VP1xVP2[0]*V13[1] + VP1[1]*N[0] + VP1xVP2[1]*V13[0] - VP1[0]*N[1],
                    VP1xVP2[1]*V12[2]               - VP1xVP2[2]*V12[1]              ,
                   -VP1xVP2[0]*V12[2]               + VP1xVP2[2]*V12[0]              ,
                    VP1xVP2[0]*V12[1]               - VP1xVP2[1]*V12[0]               };

  // Compute the partial derivatives of da0 w.r.t. the global coordinates of P 
  Scalar dda0dX[9] = { dVP2xVP3dX[1]*V23[2]        - dVP2xVP3dX[2]*V23[1]       ,
                      -dVP2xVP3dX[0]*V23[2]        + dVP2xVP3dX[2]*V23[0]       ,
                       dVP2xVP3dX[0]*V23[1]        - dVP2xVP3dX[1]*V23[0]       ,
                      -dVP2xVP3dX[1]*V13[2]        + dVP2xVP3dX[2]*V13[1]       ,
                       dVP2xVP3dX[0]*V13[2]        - dVP2xVP3dX[2]*V13[0] + N[2],
                      -dVP2xVP3dX[0]*V13[1]        + dVP2xVP3dX[1]*V13[0] - N[1],
                       dVP2xVP3dX[1]*V12[2]        - dVP2xVP3dX[2]*V12[1]       ,
                      -dVP2xVP3dX[0]*V12[2]        + dVP2xVP3dX[2]*V12[0] - N[2],
                       dVP2xVP3dX[0]*V12[1]        - dVP2xVP3dX[1]*V12[0] + N[1] };

  Scalar dda0dY[9] = { dVP2xVP3dY[1]*V23[2]        - dVP2xVP3dY[2]*V23[1]       ,
                      -dVP2xVP3dY[0]*V23[2]        + dVP2xVP3dY[2]*V23[0]       ,
                       dVP2xVP3dY[0]*V23[1]        - dVP2xVP3dY[1]*V23[0]       ,
                      -dVP2xVP3dY[1]*V13[2]        + dVP2xVP3dY[2]*V13[1] - N[2],
                       dVP2xVP3dY[0]*V13[2]        - dVP2xVP3dY[2]*V13[0]       ,
                      -dVP2xVP3dY[0]*V13[1] + N[0] + dVP2xVP3dY[1]*V13[0]       ,
                       dVP2xVP3dY[1]*V12[2]        - dVP2xVP3dY[2]*V12[1] + N[2],
                      -dVP2xVP3dY[0]*V12[2]        + dVP2xVP3dY[2]*V12[0]       ,
                       dVP2xVP3dY[0]*V12[1] - N[0] - dVP2xVP3dY[1]*V12[0]        };

  Scalar dda0dZ[9] = { dVP2xVP3dZ[1]*V23[2]        - dVP2xVP3dZ[2]*V23[1]       ,
                      -dVP2xVP3dZ[0]*V23[2]        + dVP2xVP3dZ[2]*V23[0]       ,
                       dVP2xVP3dZ[0]*V23[1]        - dVP2xVP3dZ[1]*V23[0]       ,
                      -dVP2xVP3dZ[1]*V13[2] + N[1] + dVP2xVP3dZ[2]*V13[1]       ,
                       dVP2xVP3dZ[0]*V13[2] - N[0] - dVP2xVP3dZ[2]*V13[0]       ,
                      -dVP2xVP3dZ[0]*V13[1]        + dVP2xVP3dZ[1]*V13[0]       ,
                       dVP2xVP3dZ[1]*V12[2] - N[1] - dVP2xVP3dZ[2]*V12[1]       ,
                      -dVP2xVP3dZ[0]*V12[2] + N[0] + dVP2xVP3dZ[2]*V12[0]       ,
                       dVP2xVP3dZ[0]*V12[1]        - dVP2xVP3dZ[1]*V12[0]        };

  // Compute the partial derivatives of da1 w.r.t. the global coordinates of P
  Scalar dda1dX[9] = { dVP3xVP1dX[1]*V23[2]        - dVP3xVP1dX[2]*V23[1]       ,
                      -dVP3xVP1dX[0]*V23[2]        + dVP3xVP1dX[2]*V23[0] - N[2],
                       dVP3xVP1dX[0]*V23[1]        - dVP3xVP1dX[1]*V23[0] + N[1],
                      -dVP3xVP1dX[1]*V13[2]        + dVP3xVP1dX[2]*V13[1]       ,
                       dVP3xVP1dX[0]*V13[2]        - dVP3xVP1dX[2]*V13[0]       ,
                      -dVP3xVP1dX[0]*V13[1]        + dVP3xVP1dX[1]*V13[0]       ,
                       dVP3xVP1dX[1]*V12[2]        - dVP3xVP1dX[2]*V12[1]       ,
                      -dVP3xVP1dX[0]*V12[2]        + dVP3xVP1dX[2]*V12[0] + N[2],
                       dVP3xVP1dX[0]*V12[1]        - dVP3xVP1dX[1]*V12[0] - N[1] };

  Scalar dda1dY[9] = { dVP3xVP1dY[1]*V23[2]        - dVP3xVP1dY[2]*V23[1] + N[2],
                      -dVP3xVP1dY[0]*V23[2]        + dVP3xVP1dY[2]*V23[0]       ,
                       dVP3xVP1dY[0]*V23[1] - N[0] - dVP3xVP1dY[1]*V23[0]       ,
                      -dVP3xVP1dY[1]*V13[2]        + dVP3xVP1dY[2]*V13[1]       ,
                       dVP3xVP1dY[0]*V13[2]        - dVP3xVP1dY[2]*V13[0]       ,
                      -dVP3xVP1dY[0]*V13[1]        + dVP3xVP1dY[1]*V13[0]       ,
                       dVP3xVP1dY[1]*V12[2]        - dVP3xVP1dY[2]*V12[1] - N[2],
                      -dVP3xVP1dY[0]*V12[2]        + dVP3xVP1dY[2]*V12[0]       ,
                       dVP3xVP1dY[0]*V12[1] + N[0] - dVP3xVP1dY[1]*V12[0]        };

  Scalar dda1dZ[9] = { dVP3xVP1dZ[1]*V23[2] - N[1] - dVP3xVP1dZ[2]*V23[1]       ,
                      -dVP3xVP1dZ[0]*V23[2] + N[0] + dVP3xVP1dZ[2]*V23[0]       ,
                       dVP3xVP1dZ[0]*V23[1]        - dVP3xVP1dZ[1]*V23[0]       ,
                      -dVP3xVP1dZ[1]*V13[2]        + dVP3xVP1dZ[2]*V13[1]       ,
                       dVP3xVP1dZ[0]*V13[2]        - dVP3xVP1dZ[2]*V13[0]       ,
                      -dVP3xVP1dZ[0]*V13[1]        + dVP3xVP1dZ[1]*V13[0]       ,
                       dVP3xVP1dZ[1]*V12[2] + N[1] - dVP3xVP1dZ[2]*V12[1]       ,
                      -dVP3xVP1dZ[0]*V12[2] - N[0] + dVP3xVP1dZ[2]*V12[0]       ,
                       dVP3xVP1dZ[0]*V12[1]        - dVP3xVP1dZ[1]*V12[0]        };

  // Compute the partial derivatives of da2 w.r.t. the global coordinates of P
  Scalar dda2dX[9] = { dVP1xVP2dX[1]*V23[2]        - dVP1xVP2dX[2]*V23[1]       ,
                      -dVP1xVP2dX[0]*V23[2]        + dVP1xVP2dX[2]*V23[0] + N[2],
                       dVP1xVP2dX[0]*V23[1]        - dVP1xVP2dX[1]*V23[0] - N[1],
                      -dVP1xVP2dX[1]*V13[2]        + dVP1xVP2dX[2]*V13[1]       ,
                       dVP1xVP2dX[0]*V13[2]        - dVP1xVP2dX[2]*V13[0] - N[2],
                      -dVP1xVP2dX[0]*V13[1]        + dVP1xVP2dX[1]*V13[0] + N[1],
                       dVP1xVP2dX[1]*V12[2]        - dVP1xVP2dX[2]*V12[1]       ,
                      -dVP1xVP2dX[0]*V12[2]        + dVP1xVP2dX[2]*V12[0]       ,
                       dVP1xVP2dX[0]*V12[1]        - dVP1xVP2dX[1]*V12[0]        };

  Scalar dda2dY[9] = { dVP1xVP2dY[1]*V23[2]        - dVP1xVP2dY[2]*V23[1] - N[2],
                      -dVP1xVP2dY[0]*V23[2]        + dVP1xVP2dY[2]*V23[0]       ,
                       dVP1xVP2dY[0]*V23[1] + N[0] - dVP1xVP2dY[1]*V23[0]       ,
                      -dVP1xVP2dY[1]*V13[2]        + dVP1xVP2dY[2]*V13[1] + N[2],
                       dVP1xVP2dY[0]*V13[2]        - dVP1xVP2dY[2]*V13[0]       ,
                      -dVP1xVP2dY[0]*V13[1] - N[0] + dVP1xVP2dY[1]*V13[0]       ,
                       dVP1xVP2dY[1]*V12[2]        - dVP1xVP2dY[2]*V12[1]       ,
                      -dVP1xVP2dY[0]*V12[2]        + dVP1xVP2dY[2]*V12[0]       ,
                       dVP1xVP2dY[0]*V12[1]        - dVP1xVP2dY[1]*V12[0]        };

  Scalar dda2dZ[9] = { dVP1xVP2dZ[1]*V23[2] + N[1] - dVP1xVP2dZ[2]*V23[1]       ,
                      -dVP1xVP2dZ[0]*V23[2] - N[0] + dVP1xVP2dZ[2]*V23[0]       ,
                       dVP1xVP2dZ[0]*V23[1]        - dVP1xVP2dZ[1]*V23[0]       ,
                      -dVP1xVP2dZ[1]*V13[2] - N[1] + dVP1xVP2dZ[2]*V13[1]       ,
                       dVP1xVP2dZ[0]*V13[2] + N[0] - dVP1xVP2dZ[2]*V13[0]       ,
                      -dVP1xVP2dZ[0]*V13[1]        + dVP1xVP2dZ[1]*V13[0]       ,
                       dVP1xVP2dZ[1]*V12[2]        - dVP1xVP2dZ[2]*V12[1]       ,
                      -dVP1xVP2dZ[0]*V12[2]        + dVP1xVP2dZ[2]*V12[0]       ,
                       dVP1xVP2dZ[0]*V12[1]        - dVP1xVP2dZ[1]*V12[0]        };

  // Compute the partial derivatives of the a0, a1 and a2 w.r.t. the global coordinates at P
  Scalar da0dX = V23[2]*N[1] - V23[1]*N[2];
  Scalar da1dX = V13[1]*N[2] - V13[2]*N[1];
  Scalar da2dX = V12[2]*N[1] - V12[1]*N[2];
  Scalar da0dY = V23[0]*N[2] - V23[2]*N[0];
  Scalar da1dY = V13[2]*N[0] - V13[0]*N[2];
  Scalar da2dY = V12[0]*N[2] - V12[2]*N[0];
  Scalar da0dZ = V23[1]*N[0] - V23[0]*N[1];
  Scalar da1dZ = V13[0]*N[1] - V13[1]*N[0];
  Scalar da2dZ = V12[1]*N[0] - V12[0]*N[1];

  Scalar A2 = A*A;
  Scalar a12 = a1+a2, a02 = a0+a2;
  Scalar da12dX = da1dX+da2dX, da02dX = da0dX+da2dX;
  Scalar da12dY = da1dY+da2dY, da02dY = da0dY+da2dY;
  Scalar da12dZ = da1dZ+da2dZ, da02dZ = da0dZ+da2dZ;
  for(int i=0; i<9; ++i) {
    ddLocalCoordsdX[i][0] = (dda0dX[i]*a12 + da0[i]*da12dX - da0dX*(da1[i]+da2[i]) - a0*(dda1dX[i]+dda2dX[i]))/A2;
    ddLocalCoordsdX[i][1] = (dda1dX[i]*a02 + da1[i]*da02dX - da1dX*(da0[i]+da2[i]) - a1*(dda0dX[i]+dda2dX[i]))/A2;

    ddLocalCoordsdY[i][0] = (dda0dY[i]*a12 + da0[i]*da12dY - da0dY*(da1[i]+da2[i]) - a0*(dda1dY[i]+dda2dY[i]))/A2;
    ddLocalCoordsdY[i][1] = (dda1dY[i]*a02 + da1[i]*da02dY - da1dY*(da0[i]+da2[i]) - a1*(dda0dY[i]+dda2dY[i]))/A2;

    ddLocalCoordsdZ[i][0] = (dda0dZ[i]*a12 + da0[i]*da12dZ - da0dZ*(da1[i]+da2[i]) - a0*(dda1dZ[i]+dda2dZ[i]))/A2;
    ddLocalCoordsdZ[i][1] = (dda1dZ[i]*a02 + da1[i]*da02dZ - da1dZ*(da0[i]+da2[i]) - a1*(dda0dZ[i]+dda2dZ[i]))/A2;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceTri3::Getd2LocalCoords(Scalar d2LocalCoords[][2], Scalar *M, CoordSetT &cs)
{
  // Computes the second partial derivatives w.r.t. the nodal X-, Y- and Z-coordinates of the local coordinates at a point P.
  //
  // Inputs:  M     = [M₀, M₁, M₂], the global X-, Y- and Z-coordinates of P, and
  //          cs    = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d2LocalCoords[0] = [∂²ξ/∂X₀²,   ∂²ξ/∂X₀∂Y₀, ∂²ξ/∂X₀∂Z₀, ∂²ξ/∂X₀∂X₁, ∂²ξ/∂X₀∂Y₁, ∂²ξ/∂X₀∂Z₁, ∂²ξ/∂X₀∂X₂, ∂²ξ/∂X₀∂Y₂, ∂²ξ/∂X₀∂Z₂,
  //                              ∂²ξ/∂Y₀∂X₀, ∂²ξ/∂Y₀²,   ∂²ξ/∂Y₀∂Z₀, ∂²ξ/∂Y₀∂X₁, ∂²ξ/∂Y₀∂Y₁, ∂²ξ/∂Y₀∂Z₁, ∂²ξ/∂Y₀∂X₂, ∂²ξ/∂Y₀∂Y₂, ∂²ξ/∂Y₀∂Z₂,
  //                                                                                  ⋮ 
  //                              ∂²ξ/∂Z₂∂X₀, ∂²ξ/∂Z₂∂Y₀, ∂²ξ/∂Z₂∂Z₀, ∂²ξ/∂Z₂∂X₁, ∂²ξ/∂Z₂∂Y₁, ∂²ξ/∂Z₂∂Z₁, ∂²ξ/∂Z₂∂X₂, ∂²ξ/∂Z₂∂Y₂, ∂²ξ/∂Z₂²]
  //          d2LocalCoords[1] = [∂²η/∂X₀²,   ∂²η/∂X₀∂Y₀, ∂²η/∂X₀∂Z₀, ∂²η/∂X₀∂X₁, ∂²η/∂X₀∂Y₁, ∂²η/∂X₀∂Z₁, ∂²η/∂X₀∂X₂, ∂²η/∂X₀∂Y₂, ∂²η/∂X₀∂Z₂,
  //                              ∂²η/∂Y₀∂X₀, ∂²η/∂Y₀²,   ∂²η/∂Y₀∂Z₀, ∂²η/∂Y₀∂X₁, ∂²η/∂Y₀∂Y₁, ∂²η/∂Y₀∂Z₁, ∂²η/∂Y₀∂X₂, ∂²η/∂Y₀∂Y₂, ∂²η/∂Y₀∂Z₂,
  //                                                                                  ⋮ 
  //                              ∂²η/∂Z₂∂X₀, ∂²η/∂Z₂∂Y₀, ∂²η/∂Z₂∂Z₀, ∂²η/∂Z₂∂X₁, ∂²η/∂Z₂∂Y₁, ∂²η/∂Z₂∂Z₁, ∂²η/∂Z₂∂X₂, ∂²η/∂Z₂∂Y₂, ∂²η/∂Z₂²]
  //                              where ξ, η are the local coordinates of P,
  //                              and X = [X₀, X₁, X₂], Y = [Y₀, Y₁, Y₂] and  Z = [Z₀, Z₁, Z₂] are
  //                              the nodes' global X-, Y- and Z- coordinates, respectively.

  Scalar X[3], Y[3], Z[3];
  for(int i=0; i<3; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  // Compute edge vectors
  Scalar V12[3] = { X[1] - X[0], Y[1] - Y[0], Z[1] - Z[0] };
  Scalar V13[3] = { X[2] - X[0], Y[2] - Y[0], Z[2] - Z[0] };
  Scalar V23[3] = { X[2] - X[1], Y[2] - Y[1], Z[2] - Z[1] };

  // Compute the normal (as V12 x V13)
  Scalar N[3];
  N[0] = V12[1]*V13[2] - V12[2]*V13[1];
  N[1] = V12[2]*V13[0] - V12[0]*V13[2];
  N[2] = V12[0]*V13[1] - V12[1]*V13[0];

  // Compute vectors from P to each node
  Scalar VP1[3] = { X[0] - M[0], Y[0] - M[1], Z[0] - M[2] };
  Scalar VP2[3] = { X[1] - M[0], Y[1] - M[1], Z[1] - M[2] };
  Scalar VP3[3] = { X[2] - M[0], Y[2] - M[1], Z[2] - M[2] };

  // Compute the cross products VP2 x VP3, VP3 x VP1 and VP1 x VP2
  Scalar VP2xVP3[3] = { VP2[1]*VP3[2] - VP2[2]*VP3[1], VP2[2]*VP3[0] - VP2[0]*VP3[2], VP2[0]*VP3[1] - VP2[1]*VP3[0] };
  Scalar VP3xVP1[3] = { VP3[1]*VP1[2] - VP3[2]*VP1[1], VP3[2]*VP1[0] - VP3[0]*VP1[2], VP3[0]*VP1[1] - VP3[1]*VP1[0] };
  Scalar VP1xVP2[3] = { VP1[1]*VP2[2] - VP1[2]*VP2[1], VP1[2]*VP2[0] - VP1[0]*VP2[2], VP1[0]*VP2[1] - VP1[1]*VP2[0] };

  // Compute the areas as the scalar triple products: a0 = (VP2 x VP3)·N, a1 = (VP3 x VP1)·N, a2 = (VP1 x VP2)·N
  Scalar a0 = VP2xVP3[0]*N[0] + VP2xVP3[1]*N[1] + VP2xVP3[2]*N[2];
  Scalar a1 = VP3xVP1[0]*N[0] + VP3xVP1[1]*N[1] + VP3xVP1[2]*N[2];
  Scalar a2 = VP1xVP2[0]*N[0] + VP1xVP2[1]*N[1] + VP1xVP2[2]*N[2];
  Scalar A = a0+a1+a2;

  // Compute outer products
  Scalar VP1V12[3][3], VP2V12[3][3], VP3V12[3][3];
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j)
    { VP1V12[i][j] = VP1[i]*V12[j]; VP2V12[i][j] = VP2[i]*V12[j]; VP3V12[i][j] = VP3[i]*V12[j]; }
  Scalar VP1V13[3][3], VP2V13[3][3], VP3V13[3][3];
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j)
    { VP1V13[i][j] = VP1[i]*V13[j]; VP2V13[i][j] = VP2[i]*V13[j]; VP3V13[i][j] = VP3[i]*V13[j]; }
  Scalar VP1V23[3][3], VP2V23[3][3], VP3V23[3][3];
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) 
    { VP1V23[i][j] = VP1[i]*V23[j]; VP2V23[i][j] = VP2[i]*V23[j]; VP3V23[i][j] = VP3[i]*V23[j]; }

  // Compute first partial derivatives of a0: ∂[a0]/∂X₀, ∂[a0]/∂Y₀, etc
  Scalar da0[9] = { VP2xVP3[1]*V23[2]               - VP2xVP3[2]*V23[1]              ,
                   -VP2xVP3[0]*V23[2]               + VP2xVP3[2]*V23[0]              ,
                    VP2xVP3[0]*V23[1]               - VP2xVP3[1]*V23[0]              ,
                   -VP2xVP3[1]*V13[2] - VP3[2]*N[1] + VP2xVP3[2]*V13[1] + VP3[1]*N[2],
                    VP2xVP3[0]*V13[2] + VP3[2]*N[0] - VP2xVP3[2]*V13[0] - VP3[0]*N[2],
                   -VP2xVP3[0]*V13[1] - VP3[1]*N[0] + VP2xVP3[1]*V13[0] + VP3[0]*N[1],
                    VP2xVP3[1]*V12[2] + VP2[2]*N[1] - VP2xVP3[2]*V12[1] - VP2[1]*N[2],
                   -VP2xVP3[0]*V12[2] - VP2[2]*N[0] + VP2xVP3[2]*V12[0] + VP2[0]*N[2],
                    VP2xVP3[0]*V12[1] + VP2[1]*N[0] - VP2xVP3[1]*V12[0] - VP2[0]*N[1] };

  // Compute second partial derivatives of a0: ∂²[a0]/∂X₀², ∂²[a0]/∂X₀∂Y₀, etc
  // Note: due to symmetry only the lower triangular part is computed. The strictly upper triangular part is shown commented out.
  Scalar d2a0[9][9] = {

    { 0,
      0,   // 0,
      0,   // 0,
      0,   //-VP3V23[2][2] - VP3V23[1][1],
      0,   // VP2xVP3[2] + VP3V23[0][1],
      0,   //-VP2xVP3[1] + VP3V23[0][2],
      0,   // VP2V23[2][2] + VP2V23[1][1],
      0,   //-VP2xVP3[2] - VP2V23[0][1],
      0 }, // VP2xVP3[1] - VP2V23[0][2] },

    { 0,
      0,
      0,   // 0, 
      0,   //-VP2xVP3[2] + VP3V23[1][0],
      0,   //-VP3V23[2][2] - VP3V23[0][0],
      0,   // VP2xVP3[0] + VP3V23[1][2],
      0,   // VP2xVP3[2] - VP2V23[1][0],
      0,   // VP2V23[2][2] + VP2V23[0][0],
      0 }, //-VP2xVP3[0] - VP2V23[1][2] },

    { 0,
      0,
      0,
      0,   // VP2xVP3[1] + VP3V23[2][0],
      0,   //-VP2xVP3[0] + VP3V23[2][1],
      0,   // -VP3V23[1][1] - VP3V23[0][0],
      0,   //-VP2xVP3[1] - VP2V23[2][0],
      0,   // VP2xVP3[0] - VP2V23[2][1],
      0 }, // VP2V23[1][1] + VP2V23[0][0] },

    {-VP3V23[2][2] - VP3V23[1][1],
      VP3V23[1][0] - VP2xVP3[2],
      VP3V23[2][0] + VP2xVP3[1],
      VP3V13[2][2] + VP3V13[1][1] + VP3V13[2][2] + VP3V13[1][1],
      0,   //-VP3V13[1][0] - VP3V13[0][1],
      0,   //-VP3V13[2][0] - VP3V13[0][2],
      0,   //-VP3[2]*V12[2] - VP3V12[1][1] - VP2V13[2][2] - VP2V13[1][1],
      0,   // VP3V12[1][0] + N[2] + VP2xVP3[2] + VP2V13[0][1],
      0 }, // VP3V12[2][0] - N[1] - VP2xVP3[1] + VP2V13[0][2] },

    { VP3V23[0][1] + VP2xVP3[2],
     -VP3V23[2][2] - VP3V23[0][0],
      VP3V23[2][1] - VP2xVP3[0],
     -VP3V13[0][1] - VP3V13[1][0],
      VP3V13[2][2] + VP3V13[0][0] + VP3V13[2][2] + VP3V13[0][0],
      0,   //-VP3V13[2][1] - VP3V13[1][2],
      0,   // VP3V12[0][1] - N[2] - VP2xVP3[2] + VP2V13[1][0],
      0,   //-VP3[2]*V12[2] - VP3V12[0][0] - VP2V13[2][2] - VP2V13[0][0],
      0 }, // VP3V12[2][1] + N[0] + VP2xVP3[0] + VP2V13[1][2] },

    { VP3V23[0][2] - VP2xVP3[1],
      VP3V23[1][2] + VP2xVP3[0],
     -VP3V23[1][1] - VP3V23[0][0],
     -VP3V13[0][2] - VP3V13[2][0],
     -VP3V13[1][2] - VP3V13[2][1],
      VP3V13[1][1] + VP3V13[0][0] + VP3V13[1][1] + VP3V13[0][0],
      0,   // VP3V12[0][2] + N[1] + VP2xVP3[1] + VP2V13[2][0],
      0,   // VP3V12[1][2] - N[0] - VP2xVP3[0] + VP2V13[2][1],
      0 }, //-VP3V12[1][1] - VP3V12[0][0] - VP2V13[1][1] - VP2V13[0][0] },

    { VP2V23[2][2] + VP2V23[1][1],
     -VP2V23[1][0] + VP2xVP3[2],
     -VP2V23[2][0] - VP2xVP3[1],
     -VP2V13[2][2] - VP2V13[1][1] - VP3V12[2][2] - VP3V12[1][1],
      VP2V13[1][0] - N[2] - VP2xVP3[2] + VP3V12[0][1],
      VP2V13[2][0] + N[1] + VP2xVP3[1] + VP3V12[0][2],
      VP2V12[2][2] + VP2V12[1][1] + VP2V12[2][2] + VP2V12[1][1],
      0,   //-VP2V12[1][0] - VP2V12[0][1],
      0 }, //-VP2V12[2][0] - VP2V12[0][2] },

    {-VP2V23[0][1] - VP2xVP3[2],
      VP2V23[2][2] + VP2V23[0][0],
     -VP2V23[2][1] + VP2xVP3[0],
      VP2V13[0][1] + N[2] + VP2xVP3[2] + VP3V12[1][0],
     -VP2V13[2][2] - VP2V13[0][0] - VP3V12[2][2] - VP3V12[0][0],
      VP2V13[2][1] - N[0] - VP2xVP3[0] + VP3V12[1][2],
     -VP2V12[0][1] - VP2V12[1][0],
      VP2V12[2][2] + VP2V12[0][0] + VP2V12[2][2] + VP2V12[0][0],
      0 }, //-VP2V12[2][1] - VP2V12[1][2] },

    {-VP2V23[0][2] + VP2xVP3[1],
     -VP2V23[1][2] - VP2xVP3[0],
      VP2V23[1][1] + VP2V23[0][0],
      VP2V13[0][2] - N[1] - VP2xVP3[1] + VP3V12[2][0],
      VP2V13[1][2] + N[0] + VP2xVP3[0] + VP3V12[2][1],
     -VP2V13[1][1] - VP2V13[0][0] - VP3V12[1][1] - VP3V12[0][0],
     -VP2V12[0][2] - VP2V12[2][0],
     -VP2V12[1][2] - VP2V12[2][1],
      VP2V12[1][1] + VP2V12[0][0] + VP2V12[1][1] + VP2V12[0][0] },

  };

  // Compute first partial derivatives of a1: ∂[a1]/∂X₀, ∂[a1]/∂Y₀, etc
  Scalar da1[9] = { VP3xVP1[1]*V23[2] + VP3[2]*N[1] - VP3xVP1[2]*V23[1] - VP3[1]*N[2],
                   -VP3xVP1[0]*V23[2] - VP3[2]*N[0] + VP3xVP1[2]*V23[0] + VP3[0]*N[2],
                    VP3xVP1[0]*V23[1] + VP3[1]*N[0] - VP3xVP1[1]*V23[0] - VP3[0]*N[1],
                   -VP3xVP1[1]*V13[2] + VP3xVP1[2]*V13[1],
                    VP3xVP1[0]*V13[2] - VP3xVP1[2]*V13[0],
                   -VP3xVP1[0]*V13[1] + VP3xVP1[1]*V13[0],
                    VP3xVP1[1]*V12[2] - VP1[2]*N[1] - VP3xVP1[2]*V12[1] + VP1[1]*N[2],
                   -VP3xVP1[0]*V12[2] + VP1[2]*N[0] + VP3xVP1[2]*V12[0] - VP1[0]*N[2],
                    VP3xVP1[0]*V12[1] - VP1[1]*N[0] - VP3xVP1[1]*V12[0] + VP1[0]*N[1] };

  // Compute second partial derivatives of a1: ∂²[a1]/∂X₀², ∂²[a1]/∂X₀∂Y₀, etc
  // Note: due to symmetry only the lower triangular part is computed. The strictly upper triangular part is shown commented out.
  Scalar d2a1[9][9] = {

    { VP3V23[2][2] + VP3V23[1][1] + VP3V23[2][2] + VP3V23[1][1],
      0,   //-VP3V23[1][0] - VP3V23[0][1],
      0,   //-VP3V23[2][0] - VP3V23[0][2],
      0,   //-VP3V13[2][2] - VP3V13[1][1],
      0,   // VP3V13[1][0] + VP3xVP1[2],
      0,   // VP3V13[2][0] - VP3xVP1[1],
      0,   // VP3V12[2][2] + VP3V12[1][1] - VP1V23[2][2] - VP1V23[1][1],
      0,   //-VP3V12[1][0] - N[2] - VP3xVP1[2] + VP1V23[0][1],
      0 }, //-VP3V12[2][0] + N[1] + VP3xVP1[1] + VP1V23[0][2] },

    {-VP3V23[0][1] - VP3V23[1][0],
      VP3V23[2][2] + VP3V23[0][0] + VP3V23[2][2] + VP3V23[0][0],
      0,   //-VP3V23[2][1] - VP3V23[1][2],
      0,   // VP3V13[0][1] - VP3xVP1[2],
      0,   //-VP3V13[2][2] - VP3V13[0][0],
      0,   // VP3V13[2][1] + VP3xVP1[0],
      0,   //-VP3V12[0][1] + N[2] + VP3xVP1[2] + VP1V23[1][0],
      0,   // VP3V12[2][2] + VP3V12[0][0] - VP1V23[2][2] - VP1V23[0][0],
      0 }, //-VP3[2]*V12[1] - N[0] - VP3xVP1[0] + VP1V23[1][2] },

    {-VP3V23[0][2] - VP3V23[2][0],
     -VP3V23[1][2] - VP3V23[2][1],
      VP3V23[1][1] + VP3V23[0][0] + VP3V23[1][1] + VP3V23[0][0],
      0,   // VP3V13[0][2] + VP3xVP1[1],
      0,   // VP3V13[1][2] - VP3xVP1[0],
      0,   //-VP3V13[1][1] - VP3V13[0][0],
      0,   //-VP3V12[0][2] - N[1] - VP3xVP1[1] + VP1V23[2][0],
      0,   //-VP3V12[1][2] + N[0] + VP3xVP1[0] + VP1V23[2][1],
      0 }, // VP3V12[1][1] + VP3V12[0][0] - VP1V23[1][1] - VP1V23[0][0] },

    {-VP3V13[2][2] - VP3V13[1][1],
     -VP3xVP1[2] + VP3V13[0][1],
      VP3xVP1[1] + VP3V13[0][2],
      0,
      0,   // 0,
      0,   // 0,
      0,   // VP1V13[2][2] + VP1V13[1][1],
      0,   // VP3xVP1[2] - VP1V13[0][1],
      0 }, //-VP3xVP1[1] - VP1V13[0][2] },

    { VP3xVP1[2] + VP3V13[1][0],
     -VP3V13[2][2] - VP3V13[0][0],
     -VP3xVP1[0] + VP3V13[1][2],
      0,
      0,
      0,   // 0,
      0,   //-VP3xVP1[2] - VP1V13[1][0],
      0,   // VP1V13[2][2] + VP1V13[0][0],
      0 }, // VP3xVP1[0] - VP1V13[1][2] },

    {-VP3xVP1[1] + VP3V13[2][0],
      VP3xVP1[0] + VP3V13[2][1],
     -VP3V13[1][1] - VP3V13[0][0],
      0,
      0,
      0,
      0,   // VP3xVP1[1] - VP1V13[2][0],
      0,   //-VP3xVP1[0] - VP1V13[2][1],
      0 }, // VP1V13[1][1] + VP1V13[0][0] },

    {-VP1V23[2][2] - VP1V23[1][1] + VP3V12[2][2] + VP3V12[1][1],
      VP1V23[1][0] + N[2] + VP3xVP1[2] - VP3V12[0][1],
      VP1V23[2][0] - N[1] - VP3xVP1[1] - VP3V12[0][2],
      VP1V13[2][2] + VP1V13[1][1],
     -VP1V13[1][0] - VP3xVP1[2],
     -VP1V13[2][0] + VP3xVP1[1],
     -VP1V12[2][2] - VP1V12[1][1] - VP1V12[2][2] - VP1V12[1][1],
      0,   // VP1V12[1][0] + VP1V12[0][1],
      0 }, // VP1V12[2][0] + VP1V12[0][2] },

    { VP1V23[0][1] - N[2] - VP3xVP1[2] - VP3V12[1][0],
     -VP1V23[2][2] - VP1V23[0][0] + VP3V12[2][2] + VP3V12[0][0],
      VP1V23[2][1] + N[0] + VP3xVP1[0] - VP3V12[1][2],
     -VP1V13[0][1] + VP3xVP1[2],
      VP1V13[2][2] + VP1V13[0][0],
     -VP1V13[2][1] - VP3xVP1[0],
      VP1V12[0][1] + VP1V12[1][0],
     -VP1V12[2][2] - VP1V12[0][0] - VP1V12[2][2] - VP1V12[0][0],
      0 }, // VP1V12[2][1] + VP1V12[1][2] },

    { VP1V23[0][2] + N[1] + VP3xVP1[1] - VP3V12[2][0],
      VP1V23[1][2] - N[0] - VP3xVP1[0] - VP3V12[2][1],
     -VP1V23[1][1] - VP1V23[0][0] + VP3V12[1][1] + VP3V12[0][0],
     -VP1V13[0][2] - VP3xVP1[1],
     -VP1V13[1][2] + VP3xVP1[0],
      VP1V13[1][1] + VP1V13[0][0],
      VP1V12[0][2] + VP1V12[2][0],
      VP1V12[1][2] + VP1V12[2][1],
     -VP1V12[1][1] - VP1V12[0][0] - VP1V12[1][1] - VP1V12[0][0] }

  };

  // Compute first partial derivatives of a2: ∂[a2]/∂X₀, ∂[a2]/∂Y₀, etc
  Scalar da2[9] = { VP1xVP2[1]*V23[2] - VP2[2]*N[1] - VP1xVP2[2]*V23[1] + VP2[1]*N[2],
                   -VP1xVP2[0]*V23[2] + VP2[2]*N[0] + VP1xVP2[2]*V23[0] - VP2[0]*N[2],
                    VP1xVP2[0]*V23[1] - VP2[1]*N[0] - VP1xVP2[1]*V23[0] + VP2[0]*N[1],
                   -VP1xVP2[1]*V13[2] + VP1[2]*N[1] + VP1xVP2[2]*V13[1] - VP1[1]*N[2],
                    VP1xVP2[0]*V13[2] - VP1[2]*N[0] - VP1xVP2[2]*V13[0] + VP1[0]*N[2],
                   -VP1xVP2[0]*V13[1] + VP1[1]*N[0] + VP1xVP2[1]*V13[0] - VP1[0]*N[1],
                    VP1xVP2[1]*V12[2] - VP1xVP2[2]*V12[1],
                   -VP1xVP2[0]*V12[2] + VP1xVP2[2]*V12[0],
                    VP1xVP2[0]*V12[1] - VP1xVP2[1]*V12[0] };

  // Compute second partial derivatives of a2: ∂²[a2]/∂X₀², ∂²[a2]/∂X₀∂Y₀, etc
  // Note: due to symmetry only the lower triangular part is computed. The strictly upper triangular part is shown commented out.
  Scalar d2a2[9][9] = {

    {-VP2V23[2][2] - VP2V23[1][1] - VP2V23[2][2] - VP2V23[1][1],
      0,   // VP2V23[1][0] + VP2V23[0][1],
      0,   // VP2V23[2][0] + VP2V23[0][2],
      0,   // VP2V13[2][2] + VP2V13[1][1] + VP1V23[2][2] + VP1V23[1][1],
      0,   //-VP2V13[1][0] + N[2] + VP1xVP2[2] - VP1V23[0][1],
      0,   //-VP2V13[2][0] - N[1] - VP1xVP2[1] - VP1V23[0][2],
      0,   //-VP2V12[2][2] - VP2V12[1][1],
      0,   // VP2V12[1][0] - VP1xVP2[2],
      0 }, // VP2V12[2][0] + VP1xVP2[1] },

    { VP2V23[0][1] + VP2V23[1][0],
     -VP2V23[2][2] - VP2V23[0][0] - VP2V23[2][2] - VP2V23[0][0],
      0,   // VP2V23[2][1] + VP2V23[1][2],
      0,   //-VP2V13[0][1] - N[2] - VP1xVP2[2] - VP1V23[1][0],
      0,   // VP2V13[2][2] + VP2V13[0][0] + VP1V23[2][2] + VP1V23[0][0],
      0,   //-VP2V13[2][1] + N[0] + VP1xVP2[0] - VP1V23[1][2],
      0,   // VP2V12[0][1] + VP1xVP2[2],
      0,   //-VP2V12[2][2] - VP2V12[0][0],
      0 }, // VP2V12[2][1] - VP1xVP2[0] },

    { VP2V23[0][2] + VP2V23[2][0],
      VP2V23[1][2] + VP2V23[2][1],
     -VP2V23[1][1] - VP2V23[0][0] - VP2V23[1][1] - VP2V23[0][0],
      0,   //-VP2V13[0][2] + N[1] + VP1xVP2[1] - VP1V23[2][0],
      0,   //-VP2V13[1][2] - N[0] - VP1xVP2[0] - VP1V23[2][1],
      0,   // VP2V13[1][1] + VP2V13[0][0] + VP1V23[1][1] + VP1V23[0][0],
      0,   // VP2V12[0][2] - VP1xVP2[1],
      0,   // VP2V12[1][2] + VP1xVP2[0],
      0 }, //-VP2V12[1][1] - VP2V12[0][0] },

    { VP1V23[2][2] + VP1V23[1][1] + VP2V13[2][2] + VP2V13[1][1],
     -VP1V23[1][0] - N[2] - VP1xVP2[2] - VP2V13[0][1],
     -VP1V23[2][0] + N[1] + VP1xVP2[1] - VP2V13[0][2],
     -VP1V13[2][2] - VP1V13[1][1] - VP1V13[2][2] - VP1V13[1][1],
      0,   // VP1V13[1][0] + VP1V13[0][1],
      0,   // VP1V13[2][0] + VP1V13[0][2],
      0,   // VP1V12[2][2] + VP1V12[1][1],
      0,   //-VP1V12[1][0] + VP1xVP2[2],
      0 }, //-VP1V12[2][0] - VP1xVP2[1] },

    {-VP1V23[0][1] + N[2] + VP1xVP2[2] - VP2V13[1][0],
      VP1V23[2][2] + VP1V23[0][0] + VP2V13[2][2] + VP2V13[0][0],
     -VP1V23[2][1] - N[0] - VP1xVP2[0] - VP2V13[1][2],
      VP1V13[0][1] + VP1V13[1][0],
     -VP1V13[2][2] - VP1V13[0][0] - VP1V13[2][2] - VP1V13[0][0],
      0,   // VP1V13[2][1] + VP1V13[1][2],
      0,   //-VP1V12[0][1] - VP1xVP2[2],
      0,   // VP1V12[2][2] + VP1V12[0][0],
      0 }, //-VP1V12[2][1] + VP1xVP2[0] },

    {-VP1V23[0][2] - N[1] - VP1xVP2[1] - VP2V13[2][0],
     -VP1V23[1][2] + N[0] + VP1xVP2[0] - VP2V13[2][1],
      VP1V23[1][1] + VP1V23[0][0] + VP2V13[1][1] + VP2V13[0][0],
      VP1V13[0][2] + VP1V13[2][0],
      VP1V13[1][2] + VP1V13[2][1],
     -VP1V13[1][1] - VP1V13[0][0] - VP1V13[1][1] - VP1V13[0][0],
      0,   //-VP1V12[0][2] + VP1xVP2[1],
      0,   //-VP1V12[1][2] - VP1xVP2[0],
      0 }, // VP1V12[1][1] + VP1V12[0][0] },

    {-VP2V12[2][2] - VP2V12[1][1],
      VP1xVP2[2] + VP2V12[0][1],
     -VP1xVP2[1] + VP2V12[0][2],
      VP1V12[2][2] + VP1V12[1][1],
     -VP1xVP2[2] - VP1V12[0][1],
      VP1xVP2[1] - VP1V12[0][2],
      0,
      0,   // 0,
      0 }, // 0 },

    {-VP1xVP2[2] + VP2V12[1][0],
     -VP2V12[2][2] - VP2V12[0][0],
      VP1xVP2[0] + VP2V12[1][2],
      VP1xVP2[2] - VP1V12[1][0],
      VP1V12[2][2] + VP1V12[0][0],
     -VP1xVP2[0] - VP1V12[1][2],
      0,
      0,
      0 }, // 0 },

    { VP1xVP2[1] + VP2V12[2][0],
     -VP1xVP2[0] + VP2V12[2][1],
     -VP2V12[1][1] - VP2V12[0][0],
     -VP1xVP2[1] - VP1V12[2][0],
      VP1xVP2[0] - VP1V12[2][1],
      VP1V12[1][1] + VP1V12[0][0],
      0,
      0,
      0 }

  };

  Scalar A2 = A*A;
  Scalar A4 = A2*A2;
  Scalar a12 = a1+a2, a02 = a0+a2;
  Scalar da12[9], da02[9], toto[9];
  for(int i=0; i<9; ++i) {
    da12[i] = da1[i]+da2[i];
    da02[i] = da0[i]+da2[i];
    toto[i] = 2*(da0[i]+da1[i]+da2[i])*A;
  }

  for(int i=0; i<9; ++i) {
    for(int j=0; j<=i; ++j) {
      d2LocalCoords[9*j+i][0] =
      d2LocalCoords[9*i+j][0] = ( (d2a0[i][j]*a12 + da0[j]*da12[i] - da0[i]*da12[j] - a0*(d2a1[i][j]+d2a2[i][j]))*A2
                                  -(da0[j]*a12 - a0*da12[j])*toto[i] ) / A4;
      d2LocalCoords[9*j+i][1] =
      d2LocalCoords[9*i+j][1] = ( (d2a1[i][j]*a02 + da1[j]*da02[i] - da1[i]*da02[j] - a1*(d2a0[i][j]+d2a2[i][j]))*A2
                                  -(da1[j]*a02 - a1*da02[j])*toto[i] ) / A4;
    }
  }
}

#endif
