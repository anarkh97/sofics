// ---------------------------------------------------------------- 
// HB - 05/06/03
// ---------------------------------------------------------------- 
#ifndef _FACEELEMENT_H_
#define _FACEELEMENT_H_

#include <iostream>

// STL 
#include <map>

// FEM headers
#include <Element.d/Element.h>
#include <Utils.d/resize_array.h>
#include <Math.d/matrix.h>
#include <Mortar.d/MortarDefines.h>

class DofSetArray;
struct InterpPoint;
class Connectivity;

#ifdef USE_EIGEN3
#include <Eigen/Sparse>
#endif

// ACME headers
#ifdef USE_ACME
#include "ContactSearch.h"
#endif

class FaceElement {

	double Area;

public:
	// public data
	// ~~~~~~~~~~~	
	enum {QUADFACEL4=1, QUADFACEQ8, TRIFACEL3, TRIFACEQ6, QUADFACEQ9, QUADFACEC12, TRIFACEC10};

	// Constructors 
	// ~~~~~~~~~~~~
	FaceElement():Area(0.0) { }

	// for future
	virtual FaceElement* clone()=0;

	// Setup & Update methods
	// ~~~~~~~~~~~~~~~~~~~~~~
	virtual void Renumber(std::map<int,int>& OldToNewNodeIds)=0;

	// Get and set methods
	// ~~~~~~~~~~~
	// -> pure interface methods
	virtual int nNodes() const = 0;
	virtual void GetNodes(int*, int* renumTable=0) const = 0;
	virtual void GetNodes(int*, std::map<int,int>& renumTable) const = 0;
	virtual int GetNode(int) const = 0;

	virtual int GetNodeIndex(int) const =0;
	virtual int GetFaceElemType()=0;

	// -> for the case of quadratic elements which are NOT
	//    yet supported by ACME for FFI
	//    basicaly, these methods will return the "linear" nodes
	//    of the element: for example, in the case of a FaceTri6 element
	//    these methods will return the 3 first nodes which are in fact
	//    the 3 vertices nodes of the underlaying FaceTri3 element
	// -> So we will perform the FFI search using the underlaying "linear"
	//    face element: FaceQuad8 -> FaceQuad4
	//                  FaceTri6  -> FaceTri3
	//    THIS IS AN APPROXIMATION !! IT IS OK FOR STRAIGHT EDGES FACE
	//    ELEMENT, BUT WE MAKE ERRORS IN THE CASE OF CURVED EDGES.
	//    NOTE THAT THIS IS ONLY VALID FOR THE "GEOMETRIC" FFI SEARCH.
	//    FOR THE MORTAR INTEGRATION, WE WILL USE THE "TRUE" (QUADRATIC)
	//    MAPPING.
	// -> in the case of "linear" face element, these methods are equivalent
	//    to nNodes(), GetNode(int) & GetNodes(int*).
	virtual int nVertices()=0;
	virtual int GetVertex(int)=0;
	virtual void GetVertices(int*, int* renumTable=0)=0;
	virtual void GetVertices(int*, std::map<int,int>& renumTable)=0;

	// -> for ACME
#ifdef USE_ACME
	virtual ContactSearch::ContactFace_Type GetACMEFaceElemType()=0;
#else
	virtual int GetACMEFaceElemType()=0;
#endif
	// -> see above, the discussion about the "quadratic" face elements.
	//    This method will return the ACME face element type used for
	//    the FFI search: FaceQuad8 -> FaceQuad4
	//                    FaceTri6  -> FaceTri3
	//    For "linear" face elements, this method is equivalent to
	//    GetACMEFaceElemType().
#ifdef USE_ACME
	virtual ContactSearch::ContactFace_Type GetACMEFFIFaceElemType()=0;
#else
	virtual int GetACMEFFIFaceElemType()=0;
#endif

	// Mapping & shape fct methods
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// -> pure interface methods
	virtual void GetShapeFctVal(double* Shape, double* m)=0;
	virtual double GetJacobian(double* m, CoordSet &cs)=0;
	virtual double GetIsoParamMappingNormalAndJacobian(double* Normal, double* m, CoordSet &cs)=0;
	virtual void GetIsoParamMappingNormalJacobianProduct(double* JNormal, double* m, CoordSet &cs)=0;
	virtual void LocalToGlobalCoord(double* M, double* m, CoordSet &cs)=0;

	virtual FullM ScalarMass(CoordSet& , double rho, int ngp)=0;
	virtual void IntegrateShapeFcts(double*, CoordSet&, double rho, int ngp)=0;

	virtual double* ViewRefCoords();

	virtual void GetdShapeFct(double* dShapex, double* dShapey, double* m)
	{ fprintf(stderr,"function GetdShapeFct(double*, double*, double*) undefined for this type of element!\n"); }
	virtual void Getd2ShapeFct(double *d2Shapex, double *d2Shapey, double *d2Shapexy, double *m)
	{ fprintf(stderr,"function Getd2ShapeFct(double*, double*, double*, double*) undefined for this type of element!\n"); }
	virtual void Getd3ShapeFct(double *d3Shapex, double *d3Shapey, double *d2Shapex2y, double *d2Shapexy2, double *m)
	{ fprintf(stderr,"function Getd3ShapeFct(double*, double*, double*, double*, double*) undefined for this type of element!\n"); }
	virtual void ComputedMdxAnddMdy(double* dMdx, double* dMdy, double* m, CoordSet& cs)
	{ fprintf(stderr,"function ComputedMdxAnddMdy(double*, double*, double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void Computed2Mdx2d2Mdy2Andd2Mdxdy(double *d2Mdx2, double *d2Mdy2, double *d2Mdxdy, double *m, CoordSet &cs)
	{ fprintf(stderr,"function Computed2Mdx2d2Mdy2Andd2Mdxdy(double*, double*, double*, double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double *d3Mdx3, double *d3Mdy3, double *d3Mdx2dy, double *d3Mdxdy2, double *m, CoordSet &cs)
	{ fprintf(stderr,"function Computed3Mdx3d3Mdy3d3Mdx2dyAndd3Mdxdy2(double*, double*, double*, double*, double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void GetdJNormal(double dJNormal[][3], double* m, CoordSet& cs)
	{ fprintf(stderr,"function GetdJNormal(double [][3], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void Getd2JNormal(double d2JNormal[][3], double* m, CoordSet& cs)
	{ fprintf(stderr,"function Getd2JNormal(double [][3], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void ComputedJNormaldxAnddJNormaldy(double *dJNormaldx, double *dJNormaldy, double *m, CoordSet &cs)
	{ fprintf(stderr,"function ComputedJNormaldxAnddJNormaldy(double*, double*, double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double *d2JNormaldx2, double *d2JNormaldy2, double *d2JNormaldxdy, double *m, CoordSet &cs)
	{ fprintf(stderr,"function Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(double*, double*, double*, double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void ComputeddJNormaldxAndddJNormaldy(double ddJNormaldx[][3], double ddJNormaldy[][3], double* m, CoordSet& cs)
	{ fprintf(stderr,"function ComputeddJNormaldxAndddJNormaldy(double[][3], double[][3], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void GetUnitNormal(double UnitNormal[3], double* m, CoordSet& cs)
	{ fprintf(stderr,"function  GetUnitNormal(double[3], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void GetdUnitNormal(double dUnitNormal[][3], double* m, CoordSet& cs)
	{ fprintf(stderr,"function GetdUnitNormal(double[][3], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void Getd2UnitNormal(double d2UnitNormal[][3], double* m, CoordSet& cs)
	{ fprintf(stderr,"function Getd2UnitNormal(double[][3], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void GlobalToLocalCoord(double *m, double *M, CoordSet &cs)
	{ fprintf(stderr,"function GlobalToLocalCoord(double*, double*, CoordSet) undefined for this type of element!\n"); }
	virtual void ComputedmdXdmdYAnddmdZ(double *dmdX, double *dmdY, double *dmdZ, double *M, CoordSet &cs)
	{ fprintf(stderr,"function ComputedmdXdmdYAnddmdZ(double*, double*, double*, double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void Computed2mdX2d2mdY2Etc(double *d2mdX2, double *d2mdY2, double *d2mdZ2, double *d2mdXdY, double *d2mdYdZ, double *d2mdXdZ, double *M, CoordSet &cs)
	{ fprintf(stderr,"function Computed2mdX2d2mdY2Etc(double*, double*, double*, double*, double*, double*, double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void GetdLocalCoords(double dLocalCoords[][2], double *M, CoordSet &cs)
	{ fprintf(stderr,"function GetdLocalCoords(double[][2], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void Getd2LocalCoords(double d2LocalCoords[][2], double *M, CoordSet &cs)
	{ fprintf(stderr,"function Getd2LocalCoords(double[][2], double*, CoordSet&) undefined for this type of element!\n"); }
	virtual void GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ(double ddmdX[][2], double ddmdY[][2], double ddmdZ[][2], double *M, CoordSet &cs)
	{ fprintf(stderr,"function GetddLocalCoordsdXddLocalCoordsdYAndddLocalCoordsdZ(double [][2], double [][2], double [][2], double*, CoordSet&) undefined for this type of element!\n"); }

	// Print, display ... methods
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~
	virtual void print() const =0;

	virtual int numDofs() const { fprintf(stderr,"function numDofs() undefined for this type of element!\n"); return 0; }
	virtual int* dofs(DofSetArray &, int*, int*) const
	{ fprintf(stderr,"function dofs(...) undefined for this type of element!\n"); return 0; }
	virtual void computeDisp(CoordSet&, State&, const InterpPoint&, double*, GeomState*, int*)
	{ fprintf(stderr,"function computeDisp(...) undefined for this type of element!\n"); }
	virtual void getFlLoad(const InterpPoint&, double*, double*)
	{ fprintf(stderr,"function computeDisp(...) undefined for this type of element!\n"); }

	int findEle(Connectivity *nodeToElem, int *eleTouch,
				int *eleCount, int myNum, int *fnId);

};

inline
int getNumNodes(const FaceElement &e) { return e.nNodes(); }
inline
void getNodes(const FaceElement &t, int *nd) { return t.GetNodes(nd); }
#endif

