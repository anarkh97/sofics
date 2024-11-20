// ----------------------------------------------------------------
// HB - 05/06/03
// ----------------------------------------------------------------
#ifndef _FACEQUAD4_H_
#define _FACEQUAD4_H_

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

class FaceQuad4: public FaceElement {
  private:
        int Nodes[4];
        static double RefCoords[4][2]; // coords of the nodes in the ref./parametric domain

  public:
        enum { NumberOfNodes=4 };

        // Constructors
        // ~~~~~~~~~~~~
        FaceQuad4(int *);
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
          void GetdJNormal(Scalar dJNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void Getd2JNormal(Scalar d2JNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void ComputedJNormaldxAnddJNormaldy(Scalar *dJNormaldx, Scalar *dJNormaldy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(Scalar *d2JNormaldx2, Scalar *d2JNormaldy2, Scalar *d2JNormaldxdy, Scalar *m, CoordSetT &cs);
        template<typename Scalar, typename CoordSetT>
          void ComputeddJNormaldxAndddJNormaldy(Scalar ddJNormaldx[][3], Scalar ddJNormaldy[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void GetUnitNormal(Scalar UnitNormal[3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void GetdUnitNormal(Scalar dUnitNormal[][3], Scalar* m, CoordSetT& cs);
        template<typename Scalar, typename CoordSetT>
          void Getd2UnitNormal(Scalar dUnitNormal[][3], Scalar* m, CoordSetT& cs);

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
        void Getd2UnitNormal(double dUnitNormal[][3], double* m, CoordSet& cs) override;

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

        int numDofs() const override {return 12;}
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
FaceQuad4::GetShapeFctVal(Scalar *Shape, Scalar *m)
{
  // Computes the values of the shape functions at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P.
  // Outputs: Shape = [N₀(ξ,η), N₁(ξ,η), N₂(ξ,η), N₃(ξ,η)], the values of the shape functions at P.

  Scalar& x = m[0];
  Scalar& y = m[1];

  Scalar d1 = 0.5*(1.0+x);
  Scalar d2 = 0.5*(1.0+y);
  Scalar d3 = 1.0-d1;
  Scalar d4 = 1.0-d2;

  Shape[0] = d3*d4;
  Shape[1] = d4*d1;
  Shape[2] = d1*d2;
  Shape[3] = d2*d3;
}

template<typename Scalar>
void
FaceQuad4::GetdShapeFct(Scalar *dShapex, Scalar *dShapey, Scalar *m)
{
  // Computes the values of the shape functions' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m       = [ξ, η], the local coordinates of P.
  // Outputs: dShapex = [∂N₀/∂ξ, ∂N₁/∂ξ, ∂N₂/∂ξ, ∂N₃/∂ξ]
  //          dShapey = [∂N₀/∂η, ∂N₁/∂η, ∂N₂/∂η, ∂N₃/∂η]
  //                    where N₀, N₁, N₂ and N₃ are the shape functions.

  Scalar& x = m[0];
  Scalar& y = m[1];
  Scalar onequart = 1/4.;

  Scalar xm = 1-x;
  Scalar xp = 1+x;
  Scalar ym = 1-y;
  Scalar yp = 1+y;

  dShapex[0] = -onequart*ym;
  dShapex[1] =  onequart*ym;
  dShapex[2] =  onequart*yp;
  dShapex[3] = -onequart*yp;

  dShapey[0] = -onequart*xm;
  dShapey[1] = -onequart*xp;
  dShapey[2] =  onequart*xp;
  dShapey[3] =  onequart*xm;
}

template<typename Scalar>
void
FaceQuad4::Getd2ShapeFct(Scalar *d2Shapex, Scalar *d2Shapey, Scalar *d2Shapexy, Scalar *m)
{
  // Computes the values of the shape functions' second partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m         = [ξ, η], the local coordinates of P.
  // Outputs: d2Shapex  = [∂²N₀/∂ξ², ∂²N₁/∂ξ², ∂²N₂/∂ξ², ∂²N₃/∂ξ²] 
  //          d2Shapey  = [∂²N₀/∂η², ∂²N₁/∂η², ∂²N₂/∂η², ∂²N₃/∂η²]
  //          d2Shapexy = [∂²N₀/∂ξ∂η, ∂²N₁/∂ξ∂η, ∂²N₂/∂ξ∂η, ∂²N₃/∂ξ∂η]
  //                      where N₀, N₁, N₂ and N₃ are the shape functions.

  Scalar onequart = 1/4.;

  d2Shapex[0] = 0;
  d2Shapex[1] = 0;
  d2Shapex[2] = 0;
  d2Shapex[3] = 0;

  d2Shapey[0] = 0;
  d2Shapey[1] = 0;
  d2Shapey[2] = 0;
  d2Shapey[3] = 0;

  d2Shapexy[0] = onequart;
  d2Shapexy[1] = -onequart;
  d2Shapexy[2] = onequart;
  d2Shapexy[3] = -onequart;
}

template<typename Scalar>
void
FaceQuad4::Getd3ShapeFct(Scalar *d3Shapex, Scalar *d3Shapey, Scalar *d3Shapex2y, Scalar *d3Shapexy2, Scalar *m)
{
  // Computes the values of the shape functions' third partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m          = [ξ, η], the local coordinates of P.
  // Outputs: d3Shapex   = [∂³N₀/∂ξ³, ∂³N₁/∂ξ³, ∂³N₂/∂ξ³, ∂³N₃/∂ξ³]
  //          d3Shapey   = [∂³N₀/∂η³, ∂³N₁/∂η³, ∂³N₂/∂η³, ∂³N₃/∂η³]
  //          d3Shapex2y = [∂³N₀/∂ξ²∂η, ∂³N₁/∂ξ²∂η, ∂³N₂/∂ξ²∂η, ∂³N₃/∂ξ²∂η]
  //          d3Shapexy2 = [∂³N₀/∂ξ∂η², ∂³N₁/∂ξ∂η², ∂³N₂/∂ξ∂η², ∂³N₃/∂ξ∂η²]
  //                       where N₀, N₁, N₂ and N₃ are the shape functions.

  d3Shapex[0] = 0;
  d3Shapex[1] = 0;
  d3Shapex[2] = 0;
  d3Shapex[3] = 0;

  d3Shapey[0] = 0;
  d3Shapey[1] = 0;
  d3Shapey[2] = 0;
  d3Shapey[3] = 0;

  d3Shapex2y[0] = 0;
  d3Shapex2y[1] = 0;
  d3Shapex2y[2] = 0;
  d3Shapex2y[3] = 0;

  d3Shapexy2[0] = 0;
  d3Shapexy2[1] = 0;
  d3Shapexy2[2] = 0;
  d3Shapexy2[3] = 0;
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::LocalToGlobalCoord(Scalar *M, Scalar *m, CoordSetT &cs)
{
  // Computes the global X-, Y- and Z-coordinates of a point P.
  //
  // Inputs:  m  = [ξ, η], the local coordinates of P, and
  //          cs = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: M  = [M₀(ξ,η,X,Y,Z), M₁(ξ,η,X,Y,Z), M₂(ξ,η,X,Y,Z)], the global X-, Y- and Z-coordinates of P.

  Scalar Shape[4];
  GetShapeFctVal(Shape,m);

  Scalar X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  M[0] = Shape[0]*X[0]+Shape[1]*X[1]+Shape[2]*X[2]+Shape[3]*X[3];
  M[1] = Shape[0]*Y[0]+Shape[1]*Y[1]+Shape[2]*Y[2]+Shape[3]*Y[3];
  M[2] = Shape[0]*Z[0]+Shape[1]*Z[1]+Shape[2]*Z[2]+Shape[3]*Z[3];
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::ComputedMdxAnddMdy(Scalar *dMdx, Scalar *dMdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs:  m     = [ξ, η], the local coordinates of P, and
  //          cs    = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dMdx  = [∂M₀/∂ξ, ∂M₁/∂ξ, ∂M₂/∂ξ]
  //          dMdy  = [∂M₀/∂η, ∂M₁/∂η, ∂M₂/∂η]
  //                  where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[4], dShapey[4];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  dMdx[0] = dShapex[0]*X[0] + dShapex[1]*X[1] + dShapex[2]*X[2] + dShapex[3]*X[3];
  dMdx[1] = dShapex[0]*Y[0] + dShapex[1]*Y[1] + dShapex[2]*Y[2] + dShapex[3]*Y[3];
  dMdx[2] = dShapex[0]*Z[0] + dShapex[1]*Z[1] + dShapex[2]*Z[2] + dShapex[3]*Z[3];

  dMdy[0] = dShapey[0]*X[0] + dShapey[1]*X[1] + dShapey[2]*X[2] + dShapey[3]*X[3];
  dMdy[1] = dShapey[0]*Y[0] + dShapey[1]*Y[1] + dShapey[2]*Y[2] + dShapey[3]*Y[3];
  dMdy[2] = dShapey[0]*Z[0] + dShapey[1]*Z[1] + dShapey[2]*Z[2] + dShapey[3]*Z[3];
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::Computed2Mdx2d2Mdy2Andd2Mdxdy(Scalar *d2Mdx2, Scalar *d2Mdy2, Scalar *d2Mdxdy, Scalar *m, CoordSetT &cs)
{
  // Computes the values of the global X, Y- and Z-coordinates' second partial derivatives w.r.t the local coordinates at a point P.
  //
  // Inputs: m        = [ξ, η], the local coordinates of P, and
  //         cs       = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: d2Mdx2  = [∂²M₀/∂ξ², ∂²M₁/∂ξ², ∂²M₂/∂ξ²]
  //          d2Mdy2  = [∂²M₀/∂η², ∂²M₁/∂η², ∂²M₂/∂η²]
  //          d2Mdxdy = [∂²M₀/∂ξ∂η, ∂²M₁/∂ξ∂η, ∂²M₂/∂ξ∂η]
  //                    where M₀, M₁, and M₂ are the global X-, Y- and Z-coordinates of P.

  // Special implementation for FaceQuad4
  Scalar onequart = 1/4.;

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  d2Mdx2[0] = 0; 
  d2Mdx2[1] = 0;
  d2Mdx2[2] = 0;

  d2Mdy2[0] = 0;
  d2Mdy2[1] = 0;
  d2Mdy2[2] = 0;

  d2Mdxdy[0] = onequart*(X[0] - X[1] + X[2] - X[3]);
  d2Mdxdy[1] = onequart*(Y[0] - Y[1] + Y[2] - Y[3]);
  d2Mdxdy[2] = onequart*(Z[0] - Z[1] + Z[2] - Z[3]);

/* General implementation
  // Compute shape functions' 2nd derivatives w.r.t. the local coordinates  
  Scalar d2Shapex[4], d2Shapey[4], d2Shapexy[4];
  Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  d2Mdx2[0] = d2Shapex[0]*X[0] + d2Shapex[1]*X[1] + d2Shapex[2]*X[2] + d2Shapex[3]*X[3];
  d2Mdx2[1] = d2Shapex[0]*Y[0] + d2Shapex[1]*Y[1] + d2Shapex[2]*Y[2] + d2Shapex[3]*Y[3];
  d2Mdx2[2] = d2Shapex[0]*Z[0] + d2Shapex[1]*Z[1] + d2Shapex[2]*Z[2] + d2Shapex[3]*Z[3];

  d2Mdy2[0] = d2Shapey[0]*X[0] + d2Shapey[1]*X[1] + d2Shapey[2]*X[2] + d2Shapey[3]*X[3];
  d2Mdy2[1] = d2Shapey[0]*Y[0] + d2Shapey[1]*Y[1] + d2Shapey[2]*Y[2] + d2Shapey[3]*Y[3];
  d2Mdy2[2] = d2Shapey[0]*Z[0] + d2Shapey[1]*Z[1] + d2Shapey[2]*Z[2] + d2Shapey[3]*Z[3];

  d2Mdxdy[0] = d2Shapexy[0]*X[0] + d2Shapexy[1]*X[1] + d2Shapexy[2]*X[2] + d2Shapexy[3]*X[3];
  d2Mdxdy[1] = d2Shapexy[0]*Y[0] + d2Shapexy[1]*Y[1] + d2Shapexy[2]*Y[2] + d2Shapexy[3]*Y[3];
  d2Mdxdy[2] = d2Shapexy[0]*Z[0] + d2Shapexy[1]*Z[1] + d2Shapexy[2]*Z[2] + d2Shapexy[3]*Z[3];
*/
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(Scalar *d3Mdx3, Scalar *d3Mdy3, Scalar *d3Mdx2dy, Scalar *d3Mdxdy2, Scalar *m, CoordSetT &cs)
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
FaceQuad4::GetJacobian(Scalar *m, CoordSetT &cs)
{
  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // JNormal = ∂M/∂ξ x ∂M/∂η
  Scalar JNormal[3];
  JNormal[0] = dMdx[1]*dMdy[2] - dMdx[2]*dMdy[1];
  JNormal[1] = dMdx[2]*dMdy[0] - dMdx[0]*dMdy[2];
  JNormal[2] = dMdx[0]*dMdy[1] - dMdx[1]*dMdy[0];

  using std::sqrt;
  return(sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]));
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceQuad4::GetShapeFctAndJacobian(Scalar *Shape, Scalar *m, CoordSetT &cs)
{
  Scalar x = m[0];
  Scalar y = m[1];

  Scalar d1 = 0.5*(1.0+x);
  Scalar d2 = 0.5*(1.0+y);
  Scalar d3 = 1.0-d1;
  Scalar d4 = 1.0-d2;

  Shape[0] = d3*d4;
  Shape[1] = d4*d1;
  Shape[2] = d1*d2;
  Shape[3] = d2*d3;

  double X[4], Y[4], Z[4];
  for(int i=0; i<4; ++i) {
    X[i] = cs[Nodes[i]]->x;
    Y[i] = cs[Nodes[i]]->y;
    Z[i] = cs[Nodes[i]]->z;
  }

  Scalar a[4], b[4], c[4];
  a[0] = (Y[1]-Y[0])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[1]-Z[0]);
  a[1] = (Y[1]-Y[0])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[1]-Z[0]);
  a[2] = (Y[2]-Y[3])*(Z[2]-Z[1]) - (Y[2]-Y[1])*(Z[2]-Z[3]);
  a[3] = (Y[2]-Y[3])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[2]-Z[3]);

  b[0] = (Z[1]-Z[0])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[1]-X[0]);
  b[1] = (Z[1]-Z[0])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[1]-X[0]);
  b[2] = (Z[2]-Z[3])*(X[2]-X[1]) - (Z[2]-Z[1])*(X[2]-X[3]);
  b[3] = (Z[2]-Z[3])*(X[3]-X[0]) - (Z[3]-Z[0])*(X[2]-X[3]);

  c[0] = (X[1]-X[0])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[1]-Y[0]);
  c[1] = (X[1]-X[0])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[1]-Y[0]);
  c[2] = (X[2]-X[3])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[2]-Y[3]);
  c[3] = (X[2]-X[3])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[2]-Y[3]);

  Scalar N[3];
  N[0] = Shape[0]*a[0]+Shape[1]*a[1]+Shape[2]*a[2]+Shape[3]*a[3];
  N[1] = Shape[0]*b[0]+Shape[1]*b[1]+Shape[2]*b[2]+Shape[3]*b[3];
  N[2] = Shape[0]*c[0]+Shape[1]*c[1]+Shape[2]*c[2]+Shape[3]*c[3];

  using std::sqrt;
  return(0.25*sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]));
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::GetIsoParamMappingNormalJacobianProduct(Scalar *JNormal, Scalar *m, CoordSetT &cs)
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
FaceQuad4::ComputedJNormaldxAnddJNormaldy(Scalar *dJNormaldx, Scalar *dJNormaldy, Scalar *m, CoordSetT &cs)
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

  // Special implementation for the FaceQuad4

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ∂n/∂ξ and ∂n/∂η
  dJNormaldx[0] = dMdx[1]*d2Mdxdy[2] - dMdx[2]*d2Mdxdy[1];
  dJNormaldx[1] = dMdx[2]*d2Mdxdy[0] - dMdx[0]*d2Mdxdy[2];
  dJNormaldx[2] = dMdx[0]*d2Mdxdy[1] - dMdx[1]*d2Mdxdy[0];

  dJNormaldy[0] = d2Mdxdy[1]*dMdy[2] - d2Mdxdy[2]*dMdy[1];
  dJNormaldy[1] = d2Mdxdy[2]*dMdy[0] - d2Mdxdy[0]*dMdy[2];
  dJNormaldy[2] = d2Mdxdy[0]*dMdy[1] - d2Mdxdy[1]*dMdy[0];

/* General implementation
  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ∂n/∂ξ and ∂n/∂η
  dJNormaldx[0] = d2Mdx2[1]*dMdy[2] + dMdx[1]*d2Mdxdy[2] - (d2Mdx2[2]*dMdy[1] + dMdx[2]*d2Mdxdy[1]);
  dJNormaldx[1] = d2Mdx2[2]*dMdy[0] + dMdx[2]*d2Mdxdy[0] - (d2Mdx2[0]*dMdy[2] + dMdx[0]*d2Mdxdy[2]);
  dJNormaldx[2] = d2Mdx2[0]*dMdy[1] + dMdx[0]*d2Mdxdy[1] - (d2Mdx2[1]*dMdy[0] + dMdx[1]*d2Mdxdy[0]);

  dJNormaldy[0] = d2Mdxdy[1]*dMdy[2] + dMdx[1]*d2Mdy2[2] - (d2Mdxdy[2]*dMdy[1] + dMdx[2]*d2Mdy2[1]);
  dJNormaldy[1] = d2Mdxdy[2]*dMdy[0] + dMdx[2]*d2Mdy2[0] - (d2Mdxdy[0]*dMdy[2] + dMdx[0]*d2Mdy2[2]);
  dJNormaldy[2] = d2Mdxdy[0]*dMdy[1] + dMdx[0]*d2Mdy2[1] - (d2Mdxdy[1]*dMdy[0] + dMdx[1]*d2Mdy2[0]);
*/
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(Scalar *d2JNormaldx2, Scalar *d2JNormaldy2, Scalar *d2JNormaldxdy, Scalar *m, CoordSetT &cs)
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

  // Special implementation for the FaceQuad4

  // Compute ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ∂²n/∂ξ², ∂²n/∂η² and ∂²n/∂ξ∂η
  d2JNormaldx2[0] = 0;
  d2JNormaldx2[1] = 0;
  d2JNormaldx2[2] = 0;

  d2JNormaldy2[0] = 0;
  d2JNormaldy2[1] = 0;
  d2JNormaldy2[2] = 0;

  d2JNormaldxdy[0] = d2Mdxdy[1]*d2Mdxdy[2] - d2Mdxdy[2]*d2Mdxdy[1];
  d2JNormaldxdy[1] = d2Mdxdy[2]*d2Mdxdy[0] - d2Mdxdy[0]*d2Mdxdy[2];
  d2JNormaldxdy[2] = d2Mdxdy[0]*d2Mdxdy[1] - d2Mdxdy[1]*d2Mdxdy[0];

/* General implementation
  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

  // Compute ∂³M/∂ξ³, ∂³M/∂η³, ∂³M/∂ξ²∂η and ∂³M/∂ξ∂η²
  Scalar d3Mdx3[3], d3Mdy3[3], d3Mdx2dy[3], d3Mdxdy2[3];
  Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(d3Mdx3, d3Mdy3, d3Mdx2dy, d3Mdxdy2, m, cs);

  // Compute ∂²n/∂ξ², ∂²n/∂η² and ∂²n/∂ξ∂η
  d2JNormaldx2[0] = d2Mdx2[1]*d2Mdxdy[2] + d2Mdx2[1]*d2Mdxdy[2] - (d2Mdx2[2]*d2Mdxdy[1] + d2Mdx2[2]*d2Mdxdy[1])
                  + d3Mdx3[1]*dMdy[2] + dMdx[1]*d3Mdx2dy[2] - (d3Mdx3[2]*dMdy[1] + dMdx[2]*d3Mdx2dy[1]);
  d2JNormaldx2[1] = d2Mdx2[2]*d2Mdxdy[0] + d2Mdx2[2]*d2Mdxdy[0] - (d2Mdx2[0]*d2Mdxdy[2] + d2Mdx2[0]*d2Mdxdy[2])
                  + d3Mdx3[2]*dMdy[0] + dMdx[2]*d3Mdx2dy[0] - (d3Mdx3[0]*dMdy[2] + dMdx[0]*d3Mdx2dy[2]);
  d2JNormaldx2[2] = d2Mdx2[0]*d2Mdxdy[1] + d2Mdx2[0]*d2Mdxdy[1] - (d2Mdx2[1]*d2Mdxdy[0] + d2Mdx2[1]*d2Mdxdy[0])
                  + d3Mdx3[0]*dMdy[1] + dMdx[0]*d3Mdx2dy[1] - (d3Mdx3[1]*dMdy[0] + dMdx[1]*d3Mdx2dy[0]);

  d2JNormaldy2[0] = d2Mdxdy[1]*d2Mdy2[2] + d2Mdxdy[1]*d2Mdy2[2] - (d2Mdxdy[2]*d2Mdy2[1] + d2Mdxdy[2]*d2Mdy2[1])
                  + d3Mdxdy2[1]*dMdy[2] + dMdx[1]*d3Mdy3[2] - (d3Mdxdy2[2]*dMdy[1] + dMdx[2]*d3Mdy3[1]);
  d2JNormaldy2[1] = d2Mdxdy[2]*d2Mdy2[0] + d2Mdxdy[2]*d2Mdy2[0] - (d2Mdxdy[0]*d2Mdy2[2] + d2Mdxdy[0]*d2Mdy2[2])
                  + d3Mdxdy2[2]*dMdy[0] + dMdx[2]*d3Mdy3[0] - (d3Mdxdy2[0]*dMdy[2] + dMdx[0]*d3Mdy3[2]);
  d2JNormaldy2[2] = d2Mdxdy[0]*d2Mdy2[1] + d2Mdxdy[0]*d2Mdy2[1] - (d2Mdxdy[1]*d2Mdy2[0] + d2Mdxdy[1]*d2Mdy2[0])
                  + d3Mdxdy2[0]*dMdy[1] + dMdx[0]*d3Mdy3[1] - (d3Mdxdy2[1]*dMdy[0] + dMdx[1]*d3Mdy3[0]);

  d2JNormaldxdy[0] = d2Mdx2[1]*d2Mdy2[2] + d2Mdxdy[1]*d2Mdxdy[2] - (d2Mdx2[2]*d2Mdy2[1] + d2Mdxdy[2]*d2Mdxdy[1])
                   + d3Mdx2dy[1]*dMdy[2] + dMdx[1]*d3Mdxdy2[2] - (d3Mdx2dy[2]*dMdy[1] + dMdx[2]*d3Mdxdy2[1]);
  d2JNormaldxdy[1] = d2Mdx2[2]*d2Mdy2[0] + d2Mdxdy[2]*d2Mdxdy[0] - (d2Mdx2[0]*d2Mdy2[2] + d2Mdxdy[0]*d2Mdxdy[2])
                   + d3Mdx2dy[2]*dMdy[0] + dMdx[2]*d3Mdxdy2[0] - (d3Mdx2dy[0]*dMdy[2] + dMdx[0]*d3Mdxdy2[2]);
  d2JNormaldxdy[2] = d2Mdx2[0]*d2Mdy2[1] + d2Mdxdy[0]*d2Mdxdy[1] - (d2Mdx2[1]*d2Mdy2[0] + d2Mdxdy[1]*d2Mdxdy[0])
                   + d3Mdx2dy[0]*dMdy[1] + dMdx[0]*d3Mdxdy2[1] - (d3Mdx2dy[1]*dMdy[0] + dMdx[1]*d3Mdxdy2[0]);
*/
}

template<typename Scalar, typename CoordSetT>
Scalar
FaceQuad4::GetIsoParamMappingNormalAndJacobian(Scalar *Normal, Scalar *m, CoordSetT &cs)
{
  GetIsoParamMappingNormalJacobianProduct(Normal, m, cs);

  using std::sqrt;
  Scalar NormN = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);

  if(NormN != 0.0) {
    Normal[0] /= NormN; Normal[1] /= NormN; Normal[2] /= NormN;
  }
  return(NormN);
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::GetdJNormal(Scalar dJNormal[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.

  // Inputs:  m           = [ξ, η], the local coordinates of P, and
  //          cs          = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: dJNormal[0] = [∂n₀/∂X₀, ∂n₀/∂Y₀, ∂n₀/∂Z₀, ∂n₀/∂X₁, ∂n₀/∂Y₁, ∂n₀/∂Z₁, ⋯ , ∂n₀/∂X₃, ∂n₀/∂Y₃, ∂n₀/∂Z₃]
  //          dJNormal[1] = [∂n₁/∂X₀, ∂n₁/∂Y₀, ∂n₁/∂Z₀, ∂n₁/∂X₁, ∂n₁/∂Y₁, ∂n₁/∂Z₁, ⋯ , ∂n₁/∂X₃, ∂n₁/∂Y₃, ∂n₁/∂Z₃]
  //          dJNormal[2] = [∂n₂/∂X₀, ∂n₂/∂Y₀, ∂n₂/∂Z₀, ∂n₂/∂X₁, ∂n₂/∂Y₁, ∂n₂/∂Z₁, ⋯ , ∂n₂/∂X₃, ∂n₂/∂Y₃, ∂n₂/∂Z₃]
  //                        where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                        at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                        where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                        and X = [X₀, X₁, ⋯ , X₃], Y = [Y₀, Y₁, ⋯ , Y₃] and  Z = [Z₀, Z₁, ⋯ , Z₃] are
  //                        the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[4], dShapey[4];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute dJNormal
  for(int i = 0; i < 4; ++i) {
    dJNormal[3*i  ][0] = 0;
    dJNormal[3*i  ][1] = dMdx[2]*dShapey[i] - dShapex[i]*dMdy[2];
    dJNormal[3*i  ][2] = dShapex[i]*dMdy[1] - dMdx[1]*dShapey[i];

    dJNormal[3*i+1][0] = dShapex[i]*dMdy[2] - dMdx[2]*dShapey[i];
    dJNormal[3*i+1][1] = 0;
    dJNormal[3*i+1][2] = dMdx[0]*dShapey[i] - dShapex[i]*dMdy[0];

    dJNormal[3*i+2][0] = dMdx[1]*dShapey[i] - dShapex[i]*dMdy[1];
    dJNormal[3*i+2][1] = dShapex[i]*dMdy[0] - dMdx[0]*dShapey[i];
    dJNormal[3*i+2][2] = 0;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::ComputeddJNormaldxAndddJNormaldy(Scalar ddJNormaldx[][3], Scalar ddJNormaldy[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the mixed second partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates and the local coordinate
  // of the components of the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m                = [ξ, η], the local coordinates of P, and
  //          cs               = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: ddJNormaldx[][0] = [∂²n₀/∂ξ∂X₀, ∂²n₀/∂ξ∂Y₀, ∂²n₀/∂ξ∂Z₀, ∂²n₀/∂ξ∂X₁, ∂²n₀/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ⋯ , ∂²n₀/∂ξ∂X₃, ∂²n₀/∂ξ∂Y₃, ∂²n₀/∂ξ∂Z₃]
  //          ddJNormaldx[][1] = [∂²n₁/∂ξ∂X₀, ∂²n₁/∂ξ∂Y₀, ∂²n₁/∂ξ∂Z₀, ∂²n₁/∂ξ∂X₁, ∂²n₁/∂ξ∂Y₁, ∂²n₁/∂ξ∂Z₁, ⋯ , ∂²n₁/∂ξ∂X₃, ∂²n₁/∂ξ∂Y₃, ∂²n₁/∂ξ∂Z₃]
  //          ddJNormaldx[][2] = [∂²n₂/∂ξ∂X₀, ∂²n₂/∂ξ∂Y₀, ∂²n₂/∂ξ∂Z₀, ∂²n₂/∂ξ∂X₁, ∂²n₂/∂ξ∂Y₁, ∂²n₀/∂ξ∂Z₁, ⋯ , ∂²n₂/∂ξ∂X₃, ∂²n₂/∂ξ∂Y₃, ∂²n₂/∂ξ∂Z₃]
  //          ddJNormaldy[][0] = [∂²n₀/∂η∂X₀, ∂²n₀/∂η∂Y₀, ∂²n₀/∂η∂Z₀, ∂²n₀/∂η∂X₁, ∂²n₀/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ⋯ , ∂²n₀/∂η∂X₃, ∂²n₀/∂η∂Y₃, ∂²n₀/∂η∂Z₃]
  //          ddJNormaldy[][1] = [∂²n₁/∂η∂X₀, ∂²n₁/∂η∂Y₀, ∂²n₁/∂η∂Z₀, ∂²n₁/∂η∂X₁, ∂²n₁/∂η∂Y₁, ∂²n₁/∂η∂Z₁, ⋯ , ∂²n₁/∂η∂X₃, ∂²n₁/∂η∂Y₃, ∂²n₁/∂η∂Z₃]
  //          ddJNormaldy[][2] = [∂²n₂/∂η∂X₀, ∂²n₂/∂η∂Y₀, ∂²n₂/∂η∂Z₀, ∂²n₂/∂η∂X₁, ∂²n₂/∂η∂Y₁, ∂²n₀/∂η∂Z₁, ⋯ , ∂²n₂/∂η∂X₃, ∂²n₂/∂η∂Y₃, ∂²n₂/∂η∂Z₃]
  //                             where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                             at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                             where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                             and X = [X₀, X₁, ⋯ , X₃], Y = [Y₀, Y₁, ⋯ , Y₃] and  Z = [Z₀, Z₁, ⋯ , Z₃] are
  //                             the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' derivatives w.r.t. the local coordinates
  Scalar dShapex[4], dShapey[4];
  GetdShapeFct(dShapex, dShapey, m);
  
  // Compute shape functions' 2nd derivatives w.r.t. the local coordinates (note: for the Quad4, d2Shapex = d2Shapey = 0) 
  Scalar d2Shapex[4], d2Shapey[4], d2Shapexy[4];
  Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);

  // Compute ∂M/∂ξ and ∂M/∂η
  Scalar dMdx[3], dMdy[3];
  ComputedMdxAnddMdy(dMdx, dMdy, m, cs);

  // Compute ∂²M/∂ξ², ∂²M/∂η² and ∂²M/∂ξ∂η (note: for the Quad4, ∂²M/∂ξ² = ∂²M/∂η² = 0)
  Scalar d2Mdx2[3], d2Mdy2[3], d2Mdxdy[3];
  Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2, d2Mdy2, d2Mdxdy, m, cs);

/*
  // Compute ddJNormaldx and ddJNormaldy (general implementation)
  for(int i = 0; i < 4; ++i) {
    ddJNormaldx[3*i  ][0] = 0;
    ddJNormaldx[3*i  ][1] = (d2Mdx2[2]*dShapey[i] + dMdx[2]*d2Shapexy[i]) - (d2Shapex[i]*dMdy[2] + dShapex[i]*d2Mdxdy[2]);
    ddJNormaldx[3*i  ][2] = (d2Shapex[i]*dMdy[1] + dShapex[i]*d2Mdxdy[1]) - (d2Mdx2[1]*dShapey[i] + dMdx[1]*d2Shapexy[i]);
    ddJNormaldx[3*i+1][0] = (d2Shapex[i]*dMdy[2] + dShapex[i]*d2Mdxdy[2]) - (d2Mdx2[2]*dShapey[i] + dMdx[2]*d2Shapexy[i]);
    ddJNormaldx[3*i+1][1] = 0;
    ddJNormaldx[3*i+1][2] = (d2Mdx2[0]*dShapey[i] + dMdx[0]*d2Shapexy[i]) - (d2Shapex[i]*dMdy[0] + dShapex[i]*d2Mdxdy[0]);
    ddJNormaldx[3*i+2][0] = (d2Mdx2[1]*dShapey[i] + dMdx[1]*d2Shapexy[i]) - (d2Shapex[i]*dMdy[1] + dShapex[i]*d2Mdxdy[1]);
    ddJNormaldx[3*i+2][1] = (d2Shapex[i]*dMdy[0] + dShapex[i]*d2Mdxdy[0]) - (d2Mdx2[0]*dShapey[i] + dMdx[0]*d2Shapexy[i]);
    ddJNormaldx[3*i+2][2] = 0;
  }

  for(int i = 0; i < 4; ++i) {
    ddJNormaldy[3*i  ][0] = 0;
    ddJNormaldy[3*i  ][1] = (d2Mdxdy[2]*dShapey[i] + dMdx[2]*d2Shapey[i]) - (d2Shapexy[i]*dMdy[2] + dShapex[i]*d2Mdy2[2]);
    ddJNormaldy[3*i  ][2] = (d2Shapexy[i]*dMdy[1] + dShapex[i]*d2Mdy2[1]) - (d2Mdxdy[1]*dShapey[i] + dMdx[1]*d2Shapey[i]);
    ddJNormaldy[3*i+1][0] = (d2Shapexy[i]*dMdy[2] + dShapex[i]*d2Mdy2[2]) - (d2Mdxdy[2]*dShapey[i] + dMdx[2]*d2Shapey[i]);
    ddJNormaldy[3*i+1][1] = 0;
    ddJNormaldy[3*i+1][2] = (d2Mdxdy[0]*dShapey[i] + dMdx[0]*d2Shapey[i]) - (d2Shapexy[i]*dMdy[0] + dShapex[i]*d2Mdy2[0]);
    ddJNormaldy[3*i+2][0] = (d2Mdxdy[1]*dShapey[i] + dMdx[1]*d2Shapey[i]) - (d2Shapexy[i]*dMdy[1] + dShapex[i]*d2Mdy2[1]);
    ddJNormaldy[3*i+2][1] = (d2Shapexy[i]*dMdy[0] + dShapex[i]*d2Mdy2[0]) - (d2Mdxdy[0]*dShapey[i] + dMdx[0]*d2Shapey[i]);
    ddJNormaldy[3*i+2][2] = 0;
  }
*/
  // Compute ddJNormaldx and ddJNormaldy (special implementation for the FaceQuad4)
  for(int i = 0; i < 4; ++i) {
    ddJNormaldx[3*i  ][0] = 0;
    ddJNormaldx[3*i  ][1] = dMdx[2]*d2Shapexy[i] - dShapex[i]*d2Mdxdy[2];
    ddJNormaldx[3*i  ][2] = dShapex[i]*d2Mdxdy[1] - dMdx[1]*d2Shapexy[i];
    ddJNormaldx[3*i+1][0] = dShapex[i]*d2Mdxdy[2] - dMdx[2]*d2Shapexy[i];
    ddJNormaldx[3*i+1][1] = 0;
    ddJNormaldx[3*i+1][2] = dMdx[0]*d2Shapexy[i] - dShapex[i]*d2Mdxdy[0];
    ddJNormaldx[3*i+2][0] = dMdx[1]*d2Shapexy[i] - dShapex[i]*d2Mdxdy[1];
    ddJNormaldx[3*i+2][1] = dShapex[i]*d2Mdxdy[0] - dMdx[0]*d2Shapexy[i];
    ddJNormaldx[3*i+2][2] = 0;
  }

  for(int i = 0; i < 4; ++i) {
    ddJNormaldy[3*i  ][0] = 0;
    ddJNormaldy[3*i  ][1] = (d2Mdxdy[2]*dShapey[i]) - (d2Shapexy[i]*dMdy[2]);
    ddJNormaldy[3*i  ][2] = (d2Shapexy[i]*dMdy[1]) - (d2Mdxdy[1]*dShapey[i]);
    ddJNormaldy[3*i+1][0] = (d2Shapexy[i]*dMdy[2]) - (d2Mdxdy[2]*dShapey[i]);
    ddJNormaldy[3*i+1][1] = 0;
    ddJNormaldy[3*i+1][2] = (d2Mdxdy[0]*dShapey[i]) - (d2Shapexy[i]*dMdy[0]);
    ddJNormaldy[3*i+2][0] = (d2Mdxdy[1]*dShapey[i]) - (d2Shapexy[i]*dMdy[1]);
    ddJNormaldy[3*i+2][1] = (d2Shapexy[i]*dMdy[0]) - (d2Mdxdy[0]*dShapey[i]);
    ddJNormaldy[3*i+2][2] = 0;
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::Getd2JNormal(Scalar H[][3], Scalar* m, CoordSetT& cs)
{
  // Computes the values of the second partial derivatives w.r.t the nodes' global X-, Y- and Z-coordinates of the components of
  // the vector normal to the surface at a point P whose magnitude is the Jacobian determinant of the local-to-global mapping.
  //
  // Inputs:  m    = [ξ, η], the local coordinates of P, and
  //          cs   = a container holding the nodes and their global X-, Y- and Z-coordinates.
  // Outputs: H[0] = [∂²n₀/∂X₀²,   ∂²n₀/∂X₀∂Y₀, ∂²n₀/∂X₀∂Z₀, ∂²n₀/∂X₀∂X₁, ∂²n₀/∂X₀∂Y₁, ∂²n₀/∂X₀∂Z₁, ⋯ , ∂²n₀/∂X₀∂X₃, ∂²n₀/∂X₀∂Y₃, ∂²n₀/∂X₀∂Z₃,
  //                  ∂²n₀/∂Y₀∂X₀, ∂²n₀/∂Y₀²,   ∂²n₀/∂Y₀∂Z₀, ∂²n₀/∂Y₀∂X₁, ∂²n₀/∂Y₀∂Y₁, ∂²n₀/∂Y₀∂Z₁, ⋯ , ∂²n₀/∂Y₀∂X₃, ∂²n₀/∂Y₀∂Y₃, ∂²n₀/∂Y₀∂Z₃,
  //                                                                         ⋮ 
  //                  ∂²n₀/∂Z₃∂X₀, ∂²n₀/∂Z₃∂Y₀, ∂²n₀/∂Z₃∂Z₀, ∂²n₀/∂Z₃∂X₁, ∂²n₀/∂Z₃∂Y₁, ∂²n₀/∂Z₃∂Z₁, ⋯ , ∂²n₀/∂Z₃∂X₃, ∂²n₀/∂Z₃∂Y₃, ∂²n₀/∂Z₃²]
  //          H[1] = [∂²n₁/∂X₀²,   ∂²n₁/∂X₀∂Y₀, ∂²n₁/∂X₀∂Z₀, ∂²n₁/∂X₀∂X₁, ∂²n₁/∂X₀∂Y₁, ∂²n₁/∂X₀∂Z₁, ⋯ , ∂²n₁/∂X₀∂X₃, ∂²n₁/∂X₀∂Y₃, ∂²n₁/∂X₀∂Z₃,
  //                  ∂²n₁/∂Y₀∂X₀, ∂²n₁/∂Y₀²,   ∂²n₁/∂Y₀∂Z₀, ∂²n₁/∂Y₀∂X₁, ∂²n₁/∂Y₀∂Y₁, ∂²n₁/∂Y₀∂Z₁, ⋯ , ∂²n₁/∂Y₀∂X₃, ∂²n₁/∂Y₀∂Y₃, ∂²n₁/∂Y₀∂Z₃,
  //                                                                         ⋮ 
  //                  ∂²n₁/∂Z₃∂X₀, ∂²n₁/∂Z₃∂Y₀, ∂²n₁/∂Z₃∂Z₀, ∂²n₁/∂Z₃∂X₁, ∂²n₁/∂Z₃∂Y₁, ∂²n₁/∂Z₃∂Z₁, ⋯ , ∂²n₁/∂Z₃∂X₃, ∂²n₁/∂Z₃∂Y₃, ∂²n₁/∂Z₃²]
  //          H[2] = [∂²n₂/∂X₀²,   ∂²n₂/∂X₀∂Y₀, ∂²n₂/∂X₀∂Z₀, ∂²n₂/∂X₀∂X₁, ∂²n₂/∂X₀∂Y₁, ∂²n₂/∂X₀∂Z₁, ⋯ , ∂²n₂/∂X₀∂X₃, ∂²n₂/∂X₀∂Y₃, ∂²n₂/∂X₀∂Z₃,
  //                  ∂²n₂/∂Y₀∂X₀, ∂²n₂/∂Y₀²,   ∂²n₂/∂Y₀∂Z₀, ∂²n₂/∂Y₀∂X₁, ∂²n₂/∂Y₀∂Y₁, ∂²n₂/∂Y₀∂Z₁, ⋯ , ∂²n₂/∂Y₀∂X₃, ∂²n₂/∂Y₀∂Y₃, ∂²n₂/∂Y₀∂Z₃,
  //                                                                         ⋮ 
  //                  ∂²n₂/∂Z₃∂X₀, ∂²n₂/∂Z₃∂Y₀, ∂²n₂/∂Z₃∂Z₀, ∂²n₂/∂Z₃∂X₁, ∂²n₂/∂Z₃∂Y₁, ∂²n₂/∂Z₃∂Z₁, ⋯ , ∂²n₂/∂Z₃∂X₃, ∂²n₂/∂Z₃∂Y₃, ∂²n₂/∂Z₃²]
  //                 where n₀, n₁ and n₂ are the X-, Y- and Z- components of the vector normal to the surface
  //                 at the point P whose magnitude is the Jacobian determinant of the mapping [ξ, η] →- [M₀, M₁, M₂]
  //                 where M₀, M₁ and M₂ are the global X-, Y- and Z-coordinates of P,
  //                 and X = [X₀, X₁, ⋯ , X₃], Y = [Y₀, Y₁, ⋯ , Y₃] and  Z = [Z₀, Z₁, ⋯ , Z₃] are
  //                 the nodes' global X-, Y- and Z- coordinates, respectively.

  // Compute shape functions' first partial derivatives w.r.t. the local coordinates
  Scalar dShapex[4], dShapey[4];
  GetdShapeFct(dShapex, dShapey, m);

  // Compute d2JNormal
  for(int i = 0; i < 4; ++ i) {
    for(int j = 0; j < 4; ++j) {
      H[12*(3*j  )+3*i  ][0] = 0;
      H[12*(3*j  )+3*i  ][1] = 0;
      H[12*(3*j  )+3*i  ][2] = 0;

      H[12*(3*j  )+3*i+1][0] = 0;
      H[12*(3*j  )+3*i+1][1] = 0;
      H[12*(3*j  )+3*i+1][2] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];

      H[12*(3*j  )+3*i+2][0] = 0;
      H[12*(3*j  )+3*i+2][1] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];
      H[12*(3*j  )+3*i+2][2] = 0;

      H[12*(3*j+1)+3*i  ][0] = 0;
      H[12*(3*j+1)+3*i  ][1] = 0;
      H[12*(3*j+1)+3*i  ][2] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];

      H[12*(3*j+1)+3*i+1][0] = 0;
      H[12*(3*j+1)+3*i+1][1] = 0;
      H[12*(3*j+1)+3*i+1][2] = 0;

      H[12*(3*j+1)+3*i+2][0] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];
      H[12*(3*j+1)+3*i+2][1] = 0;
      H[12*(3*j+1)+3*i+2][2] = 0;

      H[12*(3*j+2)+3*i  ][0] = 0;
      H[12*(3*j+2)+3*i  ][1] = dShapex[j]*dShapey[i] - dShapex[i]*dShapey[j];
      H[12*(3*j+2)+3*i  ][2] = 0;

      H[12*(3*j+2)+3*i+1][0] = dShapex[i]*dShapey[j] - dShapex[j]*dShapey[i];
      H[12*(3*j+2)+3*i+1][1] = 0;
      H[12*(3*j+2)+3*i+1][2] = 0;

      H[12*(3*j+2)+3*i+2][0] = 0;
      H[12*(3*j+2)+3*i+2][1] = 0;
      H[12*(3*j+2)+3*i+2][2] = 0;
    }
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::GetUnitNormal(Scalar Normal[3], Scalar *m, CoordSetT &cs)
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
FaceQuad4::GetdUnitNormal(Scalar dNormal[][3], Scalar* m, CoordSetT& cs)
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
  Scalar dJNormal[12][3];
  GetdJNormal(dJNormal, m, cs);

  // Compute the partial derivatives of the unit vector, ñ
  Scalar J2 = JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2];
  if(J2 != 0) {
    using std::sqrt;
    Scalar J = sqrt(J2);
    Scalar J3 = J*J2;
    for(int i = 0; i < 4; ++i) {
      for(int j = 0; j < 3; ++j) {
        Scalar v = 1/J3*(JNormal[0]*dJNormal[3*i+j][0] + JNormal[1]*dJNormal[3*i+j][1] + JNormal[2]*dJNormal[3*i+j][2]);
        dNormal[3*i+j][0] = dJNormal[3*i+j][0]/J - JNormal[0]*v;
        dNormal[3*i+j][1] = dJNormal[3*i+j][1]/J - JNormal[1]*v;
        dNormal[3*i+j][2] = dJNormal[3*i+j][2]/J - JNormal[2]*v;
      }
    }
  }
  else {
    for(int i = 0; i < 12; ++i)
      for(int j = 0; j < 3; ++j)
        dNormal[i][j] = dJNormal[i][j];
  }
}

template<typename Scalar, typename CoordSetT>
void
FaceQuad4::Getd2UnitNormal(Scalar d2Normal[][3], Scalar* m, CoordSetT& cs)
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
  Scalar dJNormal[12][3];
  GetdJNormal(dJNormal, m, cs);

  // Compute the second partial derivatives of n
  Scalar d2JNormal[144][3];
  Getd2JNormal(d2JNormal, m, cs);

  // Compute the partial derivatives of the unit vector, ñ
  Scalar J2 = JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2];
  if(J2 != 0) {
    using std::sqrt;
    Scalar J = sqrt(J2);
    Scalar J3 = J*J2;
    Scalar threeJ = 3*J;

    Scalar v[12];
    for(int i = 0; i < 12; ++i) v[i] = 1/J3*(JNormal[0]*dJNormal[i][0] + JNormal[1]*dJNormal[i][1] + JNormal[2]*dJNormal[i][2]);

    for(int i = 0; i < 12; ++i) {
      for(int j = i; j < 12; ++j) {
        Scalar a = dJNormal[i][0]*dJNormal[j][0] + dJNormal[i][1]*dJNormal[j][1] + dJNormal[i][2]*dJNormal[j][2];
        Scalar b = JNormal[0]*d2JNormal[12*i+j][0]+JNormal[1]*d2JNormal[12*i+j][1]+JNormal[2]* d2JNormal[12*i+j][2];
        Scalar c = threeJ*v[i]*v[j] - (a+b)/J3;
        d2Normal[12*i+j][0] = d2Normal[12*j+i][0] = d2JNormal[12*i+j][0]/J - dJNormal[j][0]*v[i] - dJNormal[i][0]*v[j] + JNormal[0]*c;
        d2Normal[12*i+j][1] = d2Normal[12*j+i][1] = d2JNormal[12*i+j][1]/J - dJNormal[j][1]*v[i] - dJNormal[i][1]*v[j] + JNormal[1]*c;
        d2Normal[12*i+j][2] = d2Normal[12*j+i][2] = d2JNormal[12*i+j][2]/J - dJNormal[j][2]*v[i] - dJNormal[i][2]*v[j] + JNormal[2]*c;
      }
    }
  }
  else {
    for(int i = 0; i < 144; ++i)
      for(int j = 0; j < 3; ++j)
        d2Normal[i][j] = d2JNormal[i][j];
  }
}

#endif
