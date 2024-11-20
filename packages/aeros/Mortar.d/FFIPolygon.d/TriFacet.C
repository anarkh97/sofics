#ifndef _TRIFACET_C_ 
#define _TRIFACET_C_

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

// FEM headers
#include <Mortar.d/FFIPolygon.d/TriFacet.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/Vector.h>

// External routines
extern void getGaussPtOnTriangle(int, int, double&, double&, double&, double&);

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar, class FaceType>
TriFacetTemplate<Scalar,FaceType>::TriFacetTemplate()
{
  Initialize();
}

template<class Scalar, class FaceType>
TriFacetTemplate<Scalar,FaceType>::TriFacetTemplate(FaceType* FaceElem, const Scalar* m1, const Scalar* m2, const Scalar* m3)
{
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];

  // set ptr to face el.
  FaceEl = FaceElem;
}

// -----------------------------------------------------------------------------------------------------
//                                             DESTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar, class FaceType>
TriFacetTemplate<Scalar,FaceType>::~TriFacetTemplate()
{
  // set ptr to face el. to NULL
  FaceEl = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                    INITIALIZATION & CLEAR/CLEAN METHODS
// -----------------------------------------------------------------------------------------------------
template<class Scalar, class FaceType>
void
TriFacetTemplate<Scalar,FaceType>::Initialize()
{
  // set local coord. of triangular vertices in face el. to 0.
  LocalCoordOnFaceEl[0][0] = 0.0; LocalCoordOnFaceEl[0][1] = 0.0;
  LocalCoordOnFaceEl[1][0] = 0.0; LocalCoordOnFaceEl[1][1] = 0.0;
  LocalCoordOnFaceEl[2][0] = 0.0; LocalCoordOnFaceEl[2][1] = 0.0;
  // set ptr to face el. to NULL
  FaceEl = 0;
}
// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar, class FaceType>
void
TriFacetTemplate<Scalar,FaceType>::SetTriFacet(FaceType* FaceElem, const Scalar* m1, const Scalar* m2, const Scalar* m3)
{
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];
  // set ptr to face el.
  FaceEl = FaceElem;
}

template<class Scalar, class FaceType>
void
TriFacetTemplate<Scalar,FaceType>::SetTriFacet(FaceType* FaceElem, const Scalar* m1, const Scalar* m2, const Scalar* m3,
                                               const Scalar* dm1, const Scalar* dm2, const Scalar* dm3, int nbDer)
{
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];

  // set ptr to face el.
  FaceEl = FaceElem;
#ifdef USE_EIGEN3
  // set derivatives of local coord.
  for(int j=0; j<2; ++j) {

    for(int i=0; i<3; ++i) {
      PartialLocalCoordOnFaceEl(i,j).resize(nbDer);
    }

    for(int k=0; k<nbDer; ++k) {
      PartialLocalCoordOnFaceEl(0,j)[k] = dm1[nbDer*j+k];
      PartialLocalCoordOnFaceEl(1,j)[k] = dm2[nbDer*j+k];
      PartialLocalCoordOnFaceEl(2,j)[k] = dm3[nbDer*j+k];
    }
  }
#endif
}

template<class Scalar, class FaceType>
void
TriFacetTemplate<Scalar,FaceType>::SetTriFacet(FaceType* FaceElem, const Scalar* m1, const Scalar* m2, const Scalar* m3,
                                               const Scalar* dm1, const Scalar* dm2, const Scalar* dm3,
                                               const Scalar* d2m1, const Scalar* d2m2, const Scalar* d2m3, int nbDer)
{
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];

  // set ptr to face el.
  FaceEl = FaceElem;
#ifdef USE_EIGEN3
  // set derivatives of local coord.
  int nbDer2 = nbDer*nbDer;
  for(int j=0; j<2; ++j) {

    for(int i=0; i<3; ++i) {
      PartialLocalCoordOnFaceEl(i,j).resize(nbDer);
      SecondPartialLocalCoordOnFaceEl(i,j).resize(nbDer,nbDer);
    }

    for(int k=0; k<nbDer; ++k) {
      PartialLocalCoordOnFaceEl(0,j)[k] = dm1[nbDer*j+k];
      PartialLocalCoordOnFaceEl(1,j)[k] = dm2[nbDer*j+k];
      PartialLocalCoordOnFaceEl(2,j)[k] = dm3[nbDer*j+k];
      for(int l=k; l<nbDer; ++l) { // note: only upper triangular part needs to be stored
        SecondPartialLocalCoordOnFaceEl(0,j)(k,l) = d2m1[nbDer2*j+nbDer*k+l];
        SecondPartialLocalCoordOnFaceEl(1,j)(k,l) = d2m2[nbDer2*j+nbDer*k+l];
        SecondPartialLocalCoordOnFaceEl(2,j)(k,l) = d2m3[nbDer2*j+nbDer*k+l];
      }
    }
  }
#endif
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
template<class Scalar, class FaceType>
void
TriFacetTemplate<Scalar,FaceType>::Print()
{
  std::cerr << " * TriFacet vertices coord. (on ref. face elem.):\n";
  for(int i=0; i<3; i++)
    std::cerr << "  -> vertex " << i+1 << ": x = " << LocalCoordOnFaceEl[i][0]
              << ", y = " << LocalCoordOnFaceEl[i][1] << std::endl;
  std::cerr << " * TriFacet ptr to face elem " << FaceEl << std::endl;
}

/*
void
TriFacet::PrintVerticesXYZ(FILE* file=stderr, CoordSet& cs, int& firstVertId)
{
  double mOnFaceEl[2];
  double XYZ[3];
  for(int i=0; i<3; i++){
     double m[2] = {LocalCoordOnFaceEl[i][0],LocalCoordOnFaceEl[i][0]};
     LocalToLocalCoordOnFaceEl(m,mOnFaceEl);    
     FaceEl->LocalToGlobalCoord(XYZ,mOnFaceEl,cs);
     fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
     firstVertId++;
  }
}

void
TriFacet::PrintTriFacetTopo(FILE* file=stderr, int& TriFacetId, int& firstVertId)
{
  fprintf(file," %6d  %6e  %6e  %6e\n",TriFacetId,firstVertId,firstVertId+1,firstVertId+2);
  TriFacetId++;
  firstVertId += 3;
}
*/

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
TriFacet::ComputeApproxArea(CoordSet &cs)
{
  double M1[3], M2[3], M3[3];
  double *m;
  m = &LocalCoordOnFaceEl[0][0];
  (*FaceEl).LocalToGlobalCoord(M1, m, cs);
  m = &LocalCoordOnFaceEl[1][0];
  (*FaceEl).LocalToGlobalCoord(M2, m, cs);
  m = &LocalCoordOnFaceEl[2][0];
  (*FaceEl).LocalToGlobalCoord(M3, m, cs);

  Vector V1(M1,3), V2(M2,3), V3(M3,3);
  Vector V12 = V2 - V1;  
  Vector V13 = V3 - V1; 
  Vector W = V12.cross(V13);

  double S = 0.5*W.magnitude();

  SetArea(S);
  return S;
}*/


// -----------------------------------------------------------------------------------------------------
//                         REF. TRIANGULAR FACET -> REF. FACE El. MAPPING METHODS 
// -----------------------------------------------------------------------------------------------------
/*void
TriFacet::MappingShapeFctAndDerivative(double* Shape, double* dShape, double* m)
{
  if(Shape==0)  Shape  = new double[3];
  if(dShape==0) dShape = new double[2][3];

  Shape[0] = m[0]; Shape[1] = m[1]; Shape[2] = 1.-m[0]-m[1];
  dShape[0][0] =  1.; dShape[0][1] =  0.; 
  dShape[1][0] =  0.; dShape[1][1] =  1.; 
  dShape[2][0] = -1.; dShape[2][1] = -1.; 
}*/

// Return jacobian of the mapping ref. TriFacet -> ref. face el.
template<class Scalar, class FaceType>
Scalar
TriFacetTemplate<Scalar,FaceType>::MappingJacobian()
{
  // J = 2.Area = ||12 x 13|| 
  Scalar X[3] = {LocalCoordOnFaceEl[0][0],LocalCoordOnFaceEl[1][0],LocalCoordOnFaceEl[2][0]};
  Scalar Y[3] = {LocalCoordOnFaceEl[0][1],LocalCoordOnFaceEl[1][1],LocalCoordOnFaceEl[2][1]};

  Scalar J = (X[1]-X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0]);

  using std::abs;
  //if(J < 0) std::cerr << "here in TriFacet::MappingJacobian, J < 0\n";
  return abs(J);// to avoid swapping the vertex to have a positive jacobian
                // this is OK because we just use it as the differential area
                // in the integration formula
}

// Map the given point m in the ref. Trifacet to the associated point in the ref. face el. 
template<class Scalar, class FaceType>
void
TriFacetTemplate<Scalar,FaceType>::LocalToLocalCoordOnFaceEl(const double* m, Scalar* mOnFaceEl)
{
  double r = m[0];
  double s = m[1];
  double t = 1 - r - s;

  mOnFaceEl[0] = t*LocalCoordOnFaceEl[0][0] + r*LocalCoordOnFaceEl[1][0] + s*LocalCoordOnFaceEl[2][0];
  mOnFaceEl[1] = t*LocalCoordOnFaceEl[0][1] + r*LocalCoordOnFaceEl[1][1] + s*LocalCoordOnFaceEl[2][1];
}

// Return the jacobian on the real face el. at the point associated to 
// the given point m in the ref. TriFace
template<class Scalar, class FaceType>
template<typename CoordSetType>
Scalar
TriFacetTemplate<Scalar,FaceType>::GetJacobianOnFaceEl(const Scalar* m, CoordSetType &cs)
{
  // J = J(mapping ref. face el. -> real face el.) * J(mapping TriFacet -> ref. face el.)
  Scalar JTriFacetMapping = MappingJacobian();
   
  Scalar mOnFaceEl[2];
  LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
  Scalar JFaceElMapping = FaceEl->GetJacobian(mOnFaceEl, cs);

  return JFaceElMapping*JTriFacetMapping;
}

// Return the isoparametric mapping normal and jacobian on the real face el. at the point 
// associated to the given point m in the ref. TriFacet
template<class Scalar, class FaceType>
template<typename CoordSetType>
Scalar
TriFacetTemplate<Scalar,FaceType>::GetIsoParamMappingNormalAndJacobianOnFaceEl(Scalar* Normal, const Scalar* m, CoordSetType &cs)
{
  // J = J(mapping ref. face el. -> real face el.) * J(mapping TriFacet -> ref. face el.)
  Scalar JTriFacetMapping = MappingJacobian();

  Scalar mOnFaceEl[2];
  LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
  Scalar JFaceElMapping = FaceEl->GetIsoParamMappingNormalAndJacobian(Normal, mOnFaceEl, cs);

  return JFaceElMapping*JTriFacetMapping;
}

// ------------------------------------------------------------------------------------------------------
//                         INTEGRATION OF SHAPE FUNCTIONS PRODUCT METHODS 
// ------------------------------------------------------------------------------------------------------
template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenFullM<Scalar>
TriFacetTemplate<Scalar,FaceType>::IntegrateShapeFctProduct(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                            CoordSetType &cs, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// NOTE: 
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facet 
  FriendFaceType* FaceElem = FriendFacet.GetPtrFaceEl();

  int nMortarShapeFct     = MortarEl->nNodes();
  int nShapeFctOnFaceElem = FaceElem->nNodes();

  GenFullM<Scalar> MatShapeFctProd(nMortarShapeFct,nShapeFctOnFaceElem);
  MatShapeFctProd.zero();

  Scalar* MortarShape      = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar)); 
  Scalar* ShapeOnFaceElem  = (Scalar*) dbg_alloca(nShapeFctOnFaceElem*sizeof(Scalar)); 
  Scalar mOnMortarEl[2], mOnFaceEl[2];
  double m[2], r, s, t, weight;

  int igp, i, j;
   
  for(igp=1; igp<=ngp; igp++){

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    // Jacobian on the face element supporting the CURRENT triangular facet
    Scalar dA = GetJacobianOnFaceEl(m, cs); 
 
    // Get local coord. on each face element 
    // -> for the mortar elem. (see the NOTE section) 
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
 
    // -> for the element supporting the given triangular facet 
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl); 
     
    // Compute shape fcts
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
 
    // -> for the element supporting the given triangular facet
    FaceElem->GetShapeFctVal(ShapeOnFaceElem, mOnFaceEl);
      
    // Shape fcts product & integration
    for(i=0;i<nMortarShapeFct;i++)
      for(j=0;j<nShapeFctOnFaceElem;j++)
        MatShapeFctProd[i][j] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j];
  }  
 
  return MatShapeFctProd; 
}

template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenFullM<Scalar>
TriFacetTemplate<Scalar,FaceType>::IntegrateNormalShapeFctProduct(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                                  CoordSetType& cs, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// times the normal (see Notes (2)).
// -> Mij = ∫_{current TriFacet}[Mortar.Shape(i)*FriendFacet.Shape(j/3)*(current TriFacet->FaceElem).normal(j%3)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH 
//                                       THE GIVEN (FRIEND) TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facet
  FriendFaceType* FaceElem = FriendFacet.GetPtrFaceEl();

  int nMortarShapeFct     = MortarEl->nNodes();
  int nShapeFctOnFaceElem = FaceElem->nNodes();

  GenFullM<Scalar> MatShapeFctProd(nMortarShapeFct,3*nShapeFctOnFaceElem);
  MatShapeFctProd.zero();

  Scalar* MortarShape      = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* ShapeOnFaceElem  = (Scalar*) dbg_alloca(nShapeFctOnFaceElem*sizeof(Scalar));
  Scalar mOnMortarEl[2], mOnFaceEl[2];
  double m[2], r, s, t, weight;
  Scalar Normal[3];
  int igp, i, j;

  for(igp=1; igp<=ngp; igp++){

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    // Jacobian on the face element supporting the CURRENT triangular facet
    Scalar dA = GetIsoParamMappingNormalAndJacobianOnFaceEl(Normal, m, cs);

    // Get local coord. on each face element
    // -> for the mortar elem. (see the NOTE section)
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl);

    // -> for the element supporting the given triangular facet
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl);

    // Compute shape fcts
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);

    // -> for the element supporting the given triangular facet
    FaceElem->GetShapeFctVal(ShapeOnFaceElem, mOnFaceEl);

    // Shape fcts product & integration
    for(i=0;i<nMortarShapeFct;i++)
      for(j=0;j<nShapeFctOnFaceElem;j++){
        MatShapeFctProd[i][3*j  ] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[0];
        MatShapeFctProd[i][3*j+1] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[1];
        MatShapeFctProd[i][3*j+2] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[2];
      } 
  }

  return MatShapeFctProd;
}

template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenVector<Scalar>
TriFacetTemplate<Scalar,FaceType>::IntegrateNormalGeoGap(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                         CoordSetType& cs, CoordSetType& cs2, int ngp, double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the gap vector between the element associated with the current triangular facet and
// the element associated to the given (Friend) triangular facet times the normal (see Notes (2)).
// -> Gi = ∫_{current TriFacet}[Mortar(i)*Gap·normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x 1
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facets
  FriendFaceType* FriendFaceEl = FriendFacet.GetPtrFaceEl();
  int nMortarShapeFct   = MortarEl->nNodes();

  GenVector<Scalar> NormalGeoGaps(nMortarShapeFct,Scalar(0.0));

  Scalar* MortarShape = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar mOnMortarEl[2], mOnFriendFaceEl[2];
  double m[2], r, s, t, weight;
  Scalar J, JNormal[3], MOnFaceEl[3], MOnFriendFaceEl[3], g;
  int igp, i;

  Scalar dA = MappingJacobian();
  if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
    FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, (Scalar*) NULL, cs);
    if(offset != 0) J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
  }

  for(igp=1; igp<=ngp; igp++) {

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    // Get local coord. on each face element
    // -> for the mortar elem. (see the NOTE section)
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl); 

    // -> for the element supporting the given triangular facet
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl);

    // Jacobian and J*Normal on the face element supporting the CURRENT triangular facet
    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
      if(offset != 0) {
        J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
      }
    }

    // Compute shape fcts
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);

    // Get global coord. on the element supporting the given triangular facets
    FaceEl->LocalToGlobalCoord(MOnFaceEl, mOnMortarEl, cs);
    FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs2);

    // Shape fcts product & integration
    g = (MOnFaceEl[0]-MOnFriendFaceEl[0])*JNormal[0]
       +(MOnFaceEl[1]-MOnFriendFaceEl[1])*JNormal[1]
       +(MOnFaceEl[2]-MOnFriendFaceEl[2])*JNormal[2];
    if(offset != 0) g += J*offset;

    for(i=0; i<nMortarShapeFct; i++) {
      NormalGeoGaps[i] += weight*MortarShape[i]*g;
    }
  }

  NormalGeoGaps *= dA;

  return NormalGeoGaps;
}

template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenFullM<Scalar>
TriFacetTemplate<Scalar, FaceType>::IntegrateGradNormalGeoGap(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                              CoordSetType& cs, CoordSetType& cs2, int ngp, double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the gradients of the products of each shape function defined
// by the given MortarElement and the gap vector between the element associated with the current
// triangular facet and the element associated to the given (Friend) triangular facet times the normal
// (see Notes (2)).
// -> C.row(i) = ∇ ∫_{current TriFacet}[Mortar(i)*Gap·normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH THE GIVEN
//                                       (FRIEND) TRIFACET + THE NUMBER OF DOFs OF THE FACE ELEMENT
//                                       ASSOCIATED WITH THE CURRENT TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE Mortar Element LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facet
  FriendFaceType* FriendFaceEl = FriendFacet.GetPtrFaceEl();

  int nMortarShapeFct         = MortarEl->nNodes();
  int nShapeFctOnFriendFaceEl = FriendFaceEl->nNodes();
  int nShapeFctOnFaceEl       = FaceEl->nNodes();

  GenFullM<Scalar> MatShapeFctProd(nMortarShapeFct, 3*nShapeFctOnFriendFaceEl + 3*nShapeFctOnFaceEl);
  MatShapeFctProd.zero();

  Scalar* MortarShape         = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* ShapeOnFaceEl       = (Scalar*) dbg_alloca(nShapeFctOnFaceEl*sizeof(Scalar));
  Scalar* ShapeOnFriendFaceEl = (Scalar*) dbg_alloca(nShapeFctOnFriendFaceEl*sizeof(Scalar));
  Scalar mOnMortarEl[2], mOnFriendFaceEl[2];
  double m[2], r, s, t, weight;
  Scalar J, JNormal[3];
  Scalar MOnFriendFaceEl[3], MOnFaceEl[3], (*dJNormal)[3];
  dJNormal = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*9*sizeof(Scalar));

  // temporary storage for pre-computation of partial derivatives
  Scalar (*PartialOnFriendFaceEl)[3], (*PartialOnFaceEl)[3];
  PartialOnFriendFaceEl = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFriendFaceEl*3*sizeof(Scalar));
  PartialOnFaceEl       = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*3*sizeof(Scalar));

  Scalar dA = MappingJacobian();
  if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
    // Get J*Normal on the face element supporting the CURRENT triangular facet
    FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, (Scalar*) NULL, cs);
    // Get the derivative of J*Normal
    FaceEl->GetdJNormal(dJNormal, (Scalar*) NULL, cs);
    // Get the norm of J*Normal
    if(offset != 0) {
      J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
    }
  }

  int igp, i, j;
  
  for(igp=1; igp<=ngp; igp++) {

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    // Get local coord. on each face element
    // -> for the mortar elem. (see the NOTE section)
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl);

    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      // Jacobian*Normal on the face element supporting the CURRENT triangular facet
      FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
      // Get the derivative of Jacobian*Normal
      FaceEl->GetdJNormal(dJNormal, mOnMortarEl, cs);
      // Get the norm of J*Normal
      if(offset != 0) {
        J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
      }
    }

    // Compute shape fcts
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
    // -> for the element supporting this triangular facet
    FaceEl->GetShapeFctVal(ShapeOnFaceEl, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFaceEl->GetShapeFctVal(ShapeOnFriendFaceEl, mOnFriendFaceEl);

    // Get global coord. on the element supporting this triangular facet
    FaceEl->LocalToGlobalCoord(MOnFaceEl, mOnMortarEl, cs);

    // Get global coord. on the element supporting the given triangular facet
    FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs2);

    // Shape fcts product & integration (partially optimized implementation)
    for(j=0; j<nShapeFctOnFriendFaceEl; j++) {
      PartialOnFriendFaceEl[j][0] = ShapeOnFriendFaceEl[j]*JNormal[0];
      PartialOnFriendFaceEl[j][1] = ShapeOnFriendFaceEl[j]*JNormal[1];
      PartialOnFriendFaceEl[j][2] = ShapeOnFriendFaceEl[j]*JNormal[2];
    }
    for(j=0; j<nShapeFctOnFaceEl; j++) {
      PartialOnFaceEl[j][0] = (ShapeOnFaceEl[j]*JNormal[0]+(MOnFaceEl[0]-MOnFriendFaceEl[0])*dJNormal[3*j  ][0]
                                                          +(MOnFaceEl[1]-MOnFriendFaceEl[1])*dJNormal[3*j  ][1]
                                                          +(MOnFaceEl[2]-MOnFriendFaceEl[2])*dJNormal[3*j  ][2]);
      PartialOnFaceEl[j][1] = (ShapeOnFaceEl[j]*JNormal[1]+(MOnFaceEl[0]-MOnFriendFaceEl[0])*dJNormal[3*j+1][0]
                                                          +(MOnFaceEl[1]-MOnFriendFaceEl[1])*dJNormal[3*j+1][1]
                                                          +(MOnFaceEl[2]-MOnFriendFaceEl[2])*dJNormal[3*j+1][2]);
      PartialOnFaceEl[j][2] = (ShapeOnFaceEl[j]*JNormal[2]+(MOnFaceEl[0]-MOnFriendFaceEl[0])*dJNormal[3*j+2][0]
                                                          +(MOnFaceEl[1]-MOnFriendFaceEl[1])*dJNormal[3*j+2][1]
                                                          +(MOnFaceEl[2]-MOnFriendFaceEl[2])*dJNormal[3*j+2][2]);
      if(offset != 0) {
        // note: for TRIFACEL3 this term can be pre-computed
        PartialOnFaceEl[j][0] += 1/J*(JNormal[0]*dJNormal[3*j  ][0]+JNormal[1]*dJNormal[3*j  ][1]+JNormal[2]*dJNormal[3*j  ][2])*offset;
        PartialOnFaceEl[j][1] += 1/J*(JNormal[0]*dJNormal[3*j+1][0]+JNormal[1]*dJNormal[3*j+1][1]+JNormal[2]*dJNormal[3*j+1][2])*offset;
        PartialOnFaceEl[j][2] += 1/J*(JNormal[0]*dJNormal[3*j+2][0]+JNormal[1]*dJNormal[3*j+2][1]+JNormal[2]*dJNormal[3*j+2][2])*offset;
      }
    }
    for(i=0; i<nMortarShapeFct; i++) {
      // -N
      for(j=0; j<nShapeFctOnFriendFaceEl; j++) {
        MatShapeFctProd[i][3*j  ] -= weight*MortarShape[i]*PartialOnFriendFaceEl[j][0];
        MatShapeFctProd[i][3*j+1] -= weight*MortarShape[i]*PartialOnFriendFaceEl[j][1];
        MatShapeFctProd[i][3*j+2] -= weight*MortarShape[i]*PartialOnFriendFaceEl[j][2];
      }
      // M
      for(j=0; j<nShapeFctOnFaceEl; j++) {
        int k = nShapeFctOnFriendFaceEl+j;
        MatShapeFctProd[i][3*k  ] += weight*MortarShape[i]*PartialOnFaceEl[j][0];
        MatShapeFctProd[i][3*k+1] += weight*MortarShape[i]*PartialOnFaceEl[j][1];
        MatShapeFctProd[i][3*k+2] += weight*MortarShape[i]*PartialOnFaceEl[j][2];
      }
    }
  }

  MatShapeFctProd *= dA;

#ifdef USE_EIGEN3
  if(PartialLocalCoordOnFaceEl(0,0).size() > 0) {
    GenFullM<Scalar> LocalGradNormalGeoGap(nMortarShapeFct, 12);
    LocalGradNormalGeoGap = IntegrateLocalGradNormalGeoGap(MortarEl, FriendFacet, cs, cs2, ngp, offset);

    for(int i=0; i<3; ++i)
      for(int j=0; j<2; ++j) {
        for(int k=0; k<nMortarShapeFct; ++k)
          for(int l=0; l<3*nShapeFctOnFriendFaceEl+3*nShapeFctOnFaceEl; ++l)
            MatShapeFctProd[k][l] += LocalGradNormalGeoGap[k][2*i+j]*FriendFacet.GetPartialLocalCoordOnFaceEl(i,j)[l]
                                    +LocalGradNormalGeoGap[k][6+2*i+j]*PartialLocalCoordOnFaceEl(i,j)[l];
      }
  }
#endif

  return MatShapeFctProd;
}

template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenFullM<Scalar>
TriFacetTemplate<Scalar, FaceType>::IntegrateHessNormalGeoGap(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                              CoordSetType& cs, CoordSetType& cs2, double* mu, int ngp, double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the Hessians of the products of each shape function defined
// by the given MortarElement and the gap vector between the element associated with the current
// triangular facet and the element associated to the given (Friend) triangular facet times the normal
// (see Notes (2)), contracted with the Lagrange multipliers
// -> H.row(j) = ∑_{i=0}^{i=N} mu(i) * ∇ [∇ ∫_{current TriFacet}[Mortar(i)*Gap·normal(current TriFacet->FaceElem)]](j)
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH THE GIVEN (FRIEND) TRIFACET + THE NUMBER
//       OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH THE CURRENT TRIFACET) x (THE NUMBER OF DOFs OF THE
//       FACE ELEMENT ASSOCIATED WITH THE GIVEN (FRIEND) TRIFACET + THE NUMBER OF DOFs OF THE FACE ELEMENT
//       ASSOCIATED WITH THE CURRENT TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facet
  FriendFaceType* FriendFaceEl = FriendFacet.GetPtrFaceEl();

  int nMortarShapeFct         = MortarEl->nNodes();
  int nShapeFctOnFriendFaceEl = FriendFaceEl->nNodes();
  int nShapeFctOnFaceEl       = FaceEl->nNodes();
  int nn3 = 3*(nShapeFctOnFriendFaceEl+nShapeFctOnFaceEl);

  GenFullM<Scalar> MatShapeFctProd(nn3, nn3);
  MatShapeFctProd.zero();

  // Check if all of the Lagrange multipliers are zero. 
  bool isNonZero = false;
  for(int i=0; i<MortarEl->nNodes(); ++i) {
    if(mu[i] != 0) { isNonZero = true; break; } 
  }
  if(!isNonZero) return MatShapeFctProd;

  Scalar* MortarShape         = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* ShapeOnFaceEl       = (Scalar*) dbg_alloca(nShapeFctOnFaceEl*sizeof(Scalar));
  Scalar* ShapeOnFriendFaceEl = (Scalar*) dbg_alloca(nShapeFctOnFriendFaceEl*sizeof(Scalar));
  Scalar mOnMortarEl[2], mOnFriendFaceEl[2];
  double m[2], r, s, t, weight;
  Scalar JNormal[3], J, J3;
  Scalar MOnFriendFaceEl[3], MOnFaceEl[3], (*dJNormal)[3], (*d2JNormal)[3];
  dJNormal = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*9*sizeof(Scalar));
  d2JNormal = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*nShapeFctOnFaceEl*27*sizeof(Scalar));

  // temporary storage for pre-computation of partial derivatives
  GenFullM<Scalar> SecondPartialOnFriendFaceEl(3*nShapeFctOnFriendFaceEl, 3*nShapeFctOnFaceEl);
  GenFullM<Scalar> SecondPartialOnFaceEl(3*nShapeFctOnFaceEl, 3*nShapeFctOnFaceEl);

  Scalar dA = MappingJacobian();
  if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
    // Get the derivative of J*Normal on the face element supporting the CURRENT triangular facet
    FaceEl->GetdJNormal(dJNormal, mOnMortarEl, cs);
    // Get the second derivative of J*Normal
    FaceEl->Getd2JNormal(d2JNormal, mOnMortarEl, cs);
    // Get the norm of J*Normal
    if(offset != 0) {
      FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
      J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
      J3 = J*J*J;
    }
  }

  int igp, i, j, k;
  
  for(igp=1; igp<=ngp; igp++) {

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    // Get local coord. on each face element
    // -> for the mortar elem. (see the NOTE section)
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl);

    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      // Get the derivative of Jacobian*Normal on the face element supporting the CURRENT triangular facet
      FaceEl->GetdJNormal(dJNormal, mOnMortarEl, cs);
      // Get the second derivative of Jacobian*Normal
      FaceEl->Getd2JNormal(d2JNormal, mOnMortarEl, cs);
      // Get the norm of J*Normal
      if(offset != 0) {
        FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
        J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
        J3 = J*J*J;
      }
    }

    // Compute shape fcts
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
    // -> for the element supporting this triangular facet
    FaceEl->GetShapeFctVal(ShapeOnFaceEl, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFaceEl->GetShapeFctVal(ShapeOnFriendFaceEl, mOnFriendFaceEl);

    // Get global coord. on the element supporting this triangular facet
    FaceEl->LocalToGlobalCoord(MOnFaceEl, mOnMortarEl, cs);

    // Get global coord. on the element supporting the given triangular facet
    FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs2);

    // Shape fcts product & integration (partially optimized implementation)
    for(j=0; j<nShapeFctOnFriendFaceEl; j++) {
      for(k=0; k<nShapeFctOnFaceEl; k++) {
        for(int p=0; p<3; ++p)
          for(int q=0; q<3; ++q)
            SecondPartialOnFriendFaceEl[3*j+p][3*k+q] = ShapeOnFriendFaceEl[j]*dJNormal[3*k+q][p];
      }
    }
    for(j=0; j<nShapeFctOnFaceEl; j++) {
      for(k=j; k<nShapeFctOnFaceEl; k++) { // symmetric
        for(int p=0; p<3; ++p)
          for(int q=0; q<3; ++q) {
            SecondPartialOnFaceEl[3*j+p][3*k+q] = (ShapeOnFaceEl[j]*dJNormal[3*k+q][p]+ShapeOnFaceEl[k]*dJNormal[3*j+p][q]
                                                   +(MOnFaceEl[0]-MOnFriendFaceEl[0])*d2JNormal[3*nShapeFctOnFaceEl*(3*j+p)+3*k+q][0]
                                                   +(MOnFaceEl[1]-MOnFriendFaceEl[1])*d2JNormal[3*nShapeFctOnFaceEl*(3*j+p)+3*k+q][1]
                                                   +(MOnFaceEl[2]-MOnFriendFaceEl[2])*d2JNormal[3*nShapeFctOnFaceEl*(3*j+p)+3*k+q][2]);
            if(offset != 0) {
              // note: for TRIFACEL3 this term is constant
              SecondPartialOnFaceEl[3*j+p][3*k+q] += (-1/J3*(JNormal[0]*dJNormal[3*k+q][0]+JNormal[1]*dJNormal[3*k+q][1]+JNormal[2]*dJNormal[3*k+q][2])
                                                           *(JNormal[0]*dJNormal[3*j+p][0]+JNormal[1]*dJNormal[3*j+p][1]+JNormal[2]*dJNormal[3*j+p][2])
                                                      +1/J*(dJNormal[3*k+q][0]*dJNormal[3*j+p][0]+JNormal[0]*d2JNormal[3*nShapeFctOnFaceEl*(3*j+p)+3*k+q][0]
                                                           +dJNormal[3*k+q][1]*dJNormal[3*j+p][1]+JNormal[1]*d2JNormal[3*nShapeFctOnFaceEl*(3*j+p)+3*k+q][1]
                                                           +dJNormal[3*k+q][2]*dJNormal[3*j+p][2]+JNormal[2]*d2JNormal[3*nShapeFctOnFaceEl*(3*j+p)+3*k+q][2]))*offset;
            }
          }
      }
    }

    for(i=0; i<nMortarShapeFct; i++) {
      if(mu[i] == 0) continue;
      // mu[i] * (Jacobian of -N[i])
      for(j=0; j<nShapeFctOnFriendFaceEl; j++) {
        // note: friend-friend block is zero
        // here is the friend-face block:
        for(k=0; k<nShapeFctOnFaceEl; k++) {
          int K = nShapeFctOnFriendFaceEl+k;
          MatShapeFctProd[3*j  ][3*K  ] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j  ][3*k  ];
          MatShapeFctProd[3*j  ][3*K+1] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j  ][3*k+1];
          MatShapeFctProd[3*j  ][3*K+2] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j  ][3*k+2];
          MatShapeFctProd[3*j+1][3*K  ] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j+1][3*k  ];
          MatShapeFctProd[3*j+1][3*K+1] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j+1][3*k+1];
          MatShapeFctProd[3*j+1][3*K+2] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j+1][3*k+2];
          MatShapeFctProd[3*j+2][3*K  ] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j+2][3*k  ];
          MatShapeFctProd[3*j+2][3*K+1] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j+2][3*k+1];
          MatShapeFctProd[3*j+2][3*K+2] -= mu[i]*weight*MortarShape[i]*SecondPartialOnFriendFaceEl[3*j+2][3*k+2];
        }
      }
      // mu[i] * (Jacobian of M[i])
      for(j=0; j<nShapeFctOnFaceEl; j++) {
        int J = nShapeFctOnFriendFaceEl+j;
        // note: face-friend block is transpose of friend-face block and will be added at the end
        // here is the face-face block:
        for(k=j; k<nShapeFctOnFaceEl; k++) { // compute the upper-triangular part only, due to symmetry
          int K = nShapeFctOnFriendFaceEl+k;
          MatShapeFctProd[3*J  ][3*K  ] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j  ][3*k  ];
          MatShapeFctProd[3*J  ][3*K+1] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j  ][3*k+1];
          MatShapeFctProd[3*J  ][3*K+2] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j  ][3*k+2];
          MatShapeFctProd[3*J+1][3*K  ] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j+1][3*k  ];
          MatShapeFctProd[3*J+1][3*K+1] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j+1][3*k+1];
          MatShapeFctProd[3*J+1][3*K+2] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j+1][3*k+2];
          MatShapeFctProd[3*J+2][3*K  ] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j+2][3*k  ];
          MatShapeFctProd[3*J+2][3*K+1] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j+2][3*k+1];
          MatShapeFctProd[3*J+2][3*K+2] += mu[i]*weight*MortarShape[i]*SecondPartialOnFaceEl[3*j+2][3*k+2];
        }
      }
    }
  }

  MatShapeFctProd *= dA;

#ifdef USE_EIGEN3
  if(PartialLocalCoordOnFaceEl(0,0).size() > 0 && SecondPartialLocalCoordOnFaceEl(0,0).size() > 0) {
    GenFullM<Scalar> LocalGradNormalGeoGap(nMortarShapeFct, 12);
    LocalGradNormalGeoGap = IntegrateLocalGradNormalGeoGap(MortarEl, FriendFacet, cs, cs2, ngp, offset);
    GenFullM<Scalar> LocalHessNormalGeoGap(12, 12);
    LocalHessNormalGeoGap = IntegrateLocalHessNormalGeoGap(MortarEl, FriendFacet, cs, cs2, mu, ngp, offset);
    GenFullM<Scalar> LocalGradGradNormalGeoGap(nn3, 12);
    LocalGradGradNormalGeoGap = IntegrateLocalGradGradNormalGeoGap(MortarEl, FriendFacet, cs, cs2, mu, ngp, offset);

    // consider a function of the form: f(x(t),t)
    // df/dt = ∂f/∂t*dt/dt + ∂f/∂x*dx/dt = ∂f/∂t + ∂f/∂x*dx/dt = ∂f/∂t + ∂f/∂x*∂x/∂t
    // d²f/dt² = ∂[df/dt]/∂t + ∂[df/dt]/∂x*∂x/∂t
    //         = ∂²f/∂t² + ∂[∂f/∂x*∂x/∂t]/∂t + ∂²f/∂x∂t*∂x/∂t + ∂[∂f/∂x*∂x/∂t]∂x*∂x/∂t
    //         = ∂²f/∂t² + ∂f/∂x*∂²x/∂t² + 2(∂²f/∂x∂t*∂x/∂t) + ∂²f/∂x²*(∂x/∂t)²
    // where: ∂²f/∂t² is MatShapeFctProd (already computed in this function) 
    //        ∂²x/∂t² is Facet.SecondPartialLocalCoordOnFaceEl
    //        ∂²f/∂x² is LocalHessNormalGeoGap
    //        ∂f/∂x is LocalGradNormalGeoGap
    //        ∂x/∂t is Facet.PartialLocalCoordOnFaceEl
    //        ∂²f/∂x∂t is LocalGradGradNormalGeoGap.

    // f = ∑ mu[i]*G[i]
    // x(t) are the LocalCoordOnFaceEls
    // t are the global nodal coordinates in CoordSets

    // ∂f/∂x*∂²x/∂t² 
    for(int k=0; k<nMortarShapeFct; ++k)
      if(mu[k] != 0)
        for(int l=0; l<nn3; ++l)
          for(int m=l; m<nn3; ++m) 
            for(int i=0,r1=0,r2=6; i<3; ++i)
              for(int j=0; j<2; ++j, ++r1, ++r2) // r1 = 2*i+j, r2 = 6+2*i+j
                MatShapeFctProd[l][m] += mu[k]*(LocalGradNormalGeoGap[k][r1]*FriendFacet.GetSecondPartialLocalCoordOnFaceEl(i,j)(l,m)
                                               +LocalGradNormalGeoGap[k][r2]*SecondPartialLocalCoordOnFaceEl(i,j)(l,m));

    // ∂²f/∂x²*(∂x/∂t)²
    for(int l=0; l<nn3; ++l)
      for(int p=0,s1=0,s2=6; p<3; ++p)
        for(int q=0; q<2; ++q, ++s1, ++s2) { // s1 = 2*p+q, s2 = 6+2*p+q
          Scalar toto1 = 0, toto2 = 0;
          for(int i=0,r1=0,r2=6; i<3; ++i)
            for(int j=0; j<2; ++j, ++r1, ++r2) { // r1 = 2*i+j, r2 = 6+2*i+j
              toto1 += LocalHessNormalGeoGap[r2][s2]*PartialLocalCoordOnFaceEl(i,j)[l]
                      +LocalHessNormalGeoGap[r1][s2]*FriendFacet.GetPartialLocalCoordOnFaceEl(i,j)[l];
              toto2 += LocalHessNormalGeoGap[r2][s1]*PartialLocalCoordOnFaceEl(i,j)[l]
                      +LocalHessNormalGeoGap[r1][s1]*FriendFacet.GetPartialLocalCoordOnFaceEl(i,j)[l];
            }
          for(int m=l; m<nn3; ++m) 
            MatShapeFctProd[l][m] += toto1*PartialLocalCoordOnFaceEl(p,q)[m] + toto2*FriendFacet.GetPartialLocalCoordOnFaceEl(p,q)[m];
        }

    // 2(∂²f/∂x∂t*∂x/∂t)
    for(int l=0; l<nn3; ++l)
      for(int m=l; m<nn3; ++m)
        for(int i=0,r1=0,r2=6; i<3; ++i) // r1 = 2*i+j, r2 = 6+2*i+j
          for(int j=0; j<2; ++j, ++r1, ++r2) 
            MatShapeFctProd[l][m] += LocalGradGradNormalGeoGap[l][r2]*PartialLocalCoordOnFaceEl(i,j)[m]
                                    +LocalGradGradNormalGeoGap[l][r1]*FriendFacet.GetPartialLocalCoordOnFaceEl(i,j)[m]
                                    +LocalGradGradNormalGeoGap[m][r2]*PartialLocalCoordOnFaceEl(i,j)[l]
                                    +LocalGradGradNormalGeoGap[m][r1]*FriendFacet.GetPartialLocalCoordOnFaceEl(i,j)[l];
  }
#endif

  // now, fill in the lower-triangular part of the symmetric matrix
  for(j=0; j<nn3; ++j)
    for(k=0; k<j; k++) MatShapeFctProd[j][k] = MatShapeFctProd[k][j];

  return MatShapeFctProd;
}

template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenFullM<Scalar>
TriFacetTemplate<Scalar, FaceType>::IntegrateLocalGradNormalGeoGap(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                                   CoordSetType& cs, CoordSetType& cs2, int ngp, double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the local gradients of the products of each shape function 
// defined by the given MortarElement and the gap vector between the element associated with the current
// triangular facet and the element associated to the given (Friend) triangular facet times the normal
// (see Notes (2)).
// -> C.row(i) = ∇ ∫_{current TriFacet}[Mortar(i)*Gap·normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH THE GIVEN
//                                       (FRIEND) TRIFACET + THE NUMBER OF DOFs OF THE FACE ELEMENT
//                                       ASSOCIATED WITH THE CURRENT TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE Mortar Element LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
//     (3) CURRENTLY ASSUMING cs and cs2 ARE CONSTANT
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facet
  FriendFaceType* FriendFaceEl = FriendFacet.GetPtrFaceEl();

  int nMortarShapeFct         = MortarEl->nNodes();
  int nShapeFctOnFriendFaceEl = FriendFaceEl->nNodes();
  int nShapeFctOnFaceEl       = FaceEl->nNodes();

  GenVector<Scalar> NormalGeoGaps(nMortarShapeFct,Scalar(0.0));
  GenFullM<Scalar> MatShapeFctProd(nMortarShapeFct, 12);
  MatShapeFctProd.zero();

  Scalar* MortarShape         = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* dShapexOnMortarEl   = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* dShapeyOnMortarEl   = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar mOnMortarEl[2], mOnFriendFaceEl[2];
  double m[2], r, s, t, weight;
  Scalar J, JNormal[3];
  Scalar MOnFriendFaceEl[3], MOnFaceEl[3];
  Scalar dMdxOnFriendFaceEl[3], dMdyOnFriendFaceEl[3], dMdxOnFaceEl[3], dMdyOnFaceEl[3];

  // temporary storage for pre-computation of partial derivatives
  Scalar (*PartialMortarShapePartialLocalCoordsOnFaceEl)[6];
  PartialMortarShapePartialLocalCoordsOnFaceEl = (Scalar (*)[6]) dbg_alloca(nMortarShapeFct*6*sizeof(Scalar));
  Scalar PartialgPartialLocalCoordsOnFriendFaceEl[6], PartialgPartialLocalCoordsOnFaceEl[6];
  Scalar dJNormaldx[3], dJNormaldy[3];

  if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
    // Get J*Normal on the face element supporting the CURRENT triangular facet
    FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, (Scalar*) NULL, cs);
    // Get the norm of J*Normal
    if(offset != 0) {
      J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
    }
  }

  int igp, i, j;
  
  for(igp=1; igp<=ngp; igp++) {

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    // Get local coord. on each face element
    // -> for the mortar elem. (see the NOTE section)
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl);
/*
    --> dmOnMortarEl/dLocalCoordsOnFaceEl = [ t 0 r 0 s 0 ]
                                            [ 0 t 0 r 0 s ]
        dmOnFriendFaceEl/dFriendFacet.LocalCoordsOnFaceEl = [ t 0 r 0 s 0 ]
                                                            [ 0 t 0 r 0 s ]
*/
    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      // Jacobian*Normal on the face element supporting the CURRENT triangular facet
      FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
      // Get the norm of J*Normal
      if(offset != 0) {
        J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
      }
    }

    // Compute shape fcts & derivatives w.r.t m
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
    MortarEl->GetdShapeFct(dShapexOnMortarEl, dShapeyOnMortarEl, mOnMortarEl);

    for(i=0; i<nMortarShapeFct; i++) {
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][0] = dShapexOnMortarEl[i]*t;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][1] = dShapeyOnMortarEl[i]*t;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][2] = dShapexOnMortarEl[i]*r;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][3] = dShapeyOnMortarEl[i]*r;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][4] = dShapexOnMortarEl[i]*s;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][5] = dShapeyOnMortarEl[i]*s;
    }

    // Get global coord. on the element supporting this triangular facet & derivatives w.r.t. m
    FaceEl->LocalToGlobalCoord(MOnFaceEl, mOnMortarEl, cs);

    // Get global coord. on the element supporting the given triangular facet & derivatives w.r.t m
    FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs2);
    FriendFaceEl->ComputedMdxAnddMdy(dMdxOnFriendFaceEl, dMdyOnFriendFaceEl, mOnFriendFaceEl, cs2);

    //  ∂[g]/∂[FriendFacet.LocalCoordsOnFaceEl]
    Scalar a = (dMdxOnFriendFaceEl[0]*JNormal[0]+dMdxOnFriendFaceEl[1]*JNormal[1]+dMdxOnFriendFaceEl[2]*JNormal[2]);
    Scalar b = (dMdyOnFriendFaceEl[0]*JNormal[0]+dMdyOnFriendFaceEl[1]*JNormal[1]+dMdyOnFriendFaceEl[2]*JNormal[2]);
    PartialgPartialLocalCoordsOnFriendFaceEl[0] = -a*t;
    PartialgPartialLocalCoordsOnFriendFaceEl[1] = -b*t;
    PartialgPartialLocalCoordsOnFriendFaceEl[2] = -a*r;
    PartialgPartialLocalCoordsOnFriendFaceEl[3] = -b*r;
    PartialgPartialLocalCoordsOnFriendFaceEl[4] = -a*s;
    PartialgPartialLocalCoordsOnFriendFaceEl[5] = -b*s;

    // ∂[g]/∂[LocalCoordsOnFaceEl]  (note that the ∂[MOnFaceEl]/∂[LocalCoordsOnFaceEl] term is always zero)
    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      Scalar gvec[3] = { MOnFaceEl[0]-MOnFriendFaceEl[0], MOnFaceEl[1]-MOnFriendFaceEl[1], MOnFaceEl[2]-MOnFriendFaceEl[2] };
      if(offset != 0) {
        gvec[0] += offset/J*JNormal[0];
        gvec[1] += offset/J*JNormal[1];
        gvec[2] += offset/J*JNormal[2];
      }
      // ∂[JNormal]/∂[LocalCoordsOnFaceEl] = ∂[JNormal]/∂[mOnMortarEl]*∂[mOnMortarEl]/∂[LocalCoordsOnFaceEl]
      FaceEl->ComputedJNormaldxAnddJNormaldy(dJNormaldx, dJNormaldy, mOnMortarEl, cs);
      Scalar c = gvec[0]*dJNormaldx[0] + gvec[1]*dJNormaldx[1] + gvec[2]*dJNormaldx[2];
      Scalar d = gvec[0]*dJNormaldy[0] + gvec[1]*dJNormaldy[1] + gvec[2]*dJNormaldy[2];
      PartialgPartialLocalCoordsOnFaceEl[0] = c*t;
      PartialgPartialLocalCoordsOnFaceEl[1] = d*t;
      PartialgPartialLocalCoordsOnFaceEl[2] = c*r;
      PartialgPartialLocalCoordsOnFaceEl[3] = d*r;
      PartialgPartialLocalCoordsOnFaceEl[4] = c*s;
      PartialgPartialLocalCoordsOnFaceEl[5] = d*s;
    }
    else {
      for(j=0; j<6; j++) PartialgPartialLocalCoordsOnFaceEl[j] = 0;
    }

    // Shape fcts product & integration
    Scalar g = (MOnFaceEl[0]-MOnFriendFaceEl[0])*JNormal[0] +(MOnFaceEl[1]-MOnFriendFaceEl[1])*JNormal[1] +(MOnFaceEl[2]-MOnFriendFaceEl[2])*JNormal[2];
    if(offset != 0) g += J*offset;

    for(i=0; i<nMortarShapeFct; i++) {
      NormalGeoGaps[i] += weight*MortarShape[i]*g;
      for(j=0; j<6; j++) { // ∂[NormalGeoGaps]/∂[FriendFacet.LocalCoordsOnFaceEl]
        MatShapeFctProd[i][j] += weight*MortarShape[i]*PartialgPartialLocalCoordsOnFriendFaceEl[j];
      }
      for(j=0; j<6; ++j) { // ∂[NormalGeoGaps]/∂[LocalCoordsOnFaceEl] 
        MatShapeFctProd[i][6+j] += weight*(MortarShape[i]*PartialgPartialLocalCoordsOnFaceEl[j] + PartialMortarShapePartialLocalCoordsOnFaceEl[i][j]*g);
      }
    }
  }

  Scalar dA = (LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0])*(LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1])
             -(LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0])*(LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]);

  Scalar PartialdAPartialLocalCoordsOnFaceEl[6] = {
     -(LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1]) + (LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]),
     -(LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0]) + (LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0]),
      (LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1]),
     -(LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0]),
     -(LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]),
      (LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0]) };

  if(dA < 0) {
    dA = -dA;
    for(i=0; i<6; ++i) PartialdAPartialLocalCoordsOnFaceEl[i] = -PartialdAPartialLocalCoordsOnFaceEl[i];
  }

  MatShapeFctProd *= dA;
  for(i=0; i<nMortarShapeFct; ++i)
    for(j=0; j<6; ++j) MatShapeFctProd[i][6+j] += NormalGeoGaps[i]*PartialdAPartialLocalCoordsOnFaceEl[j];

  return MatShapeFctProd;
}

template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenFullM<Scalar>
TriFacetTemplate<Scalar, FaceType>::IntegrateLocalHessNormalGeoGap(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                                   CoordSetType& cs, CoordSetType& cs2, double* mu, int ngp, double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the local Hessians of the products of each shape function
// defined by the given MortarElement and the gap vector between the element associated with the current
// triangular facet and the element associated to the given (Friend) triangular facet times the normal
// (see Notes (2)), contracted with the Lagrange multipliers
// -> H.row(j) = ∑_{i=0}^{i=N} mu(i) * ∇ [∇ ∫_{current TriFacet}[Mortar(i)*Gap·normal(current TriFacet->FaceElem)]](j)
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH THE GIVEN (FRIEND) TRIFACET + THE NUMBER
//       OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH THE CURRENT TRIFACET) x (THE NUMBER OF DOFs OF THE
//       FACE ELEMENT ASSOCIATED WITH THE GIVEN (FRIEND) TRIFACET + THE NUMBER OF DOFs OF THE FACE ELEMENT
//       ASSOCIATED WITH THE CURRENT TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
//     (3) CURRENTLY ASSUMING cs and cs2 ARE CONSTANT
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facet
  FriendFaceType* FriendFaceEl = FriendFacet.GetPtrFaceEl();

  int nMortarShapeFct         = MortarEl->nNodes();
  int nShapeFctOnFriendFaceEl = FriendFaceEl->nNodes();
  int nShapeFctOnFaceEl       = FaceEl->nNodes();

  GenFullM<Scalar> MatShapeFctProd(12, 12);
  MatShapeFctProd.zero();

  // Check if all of the Lagrange multipliers are zero. 
  bool isNonZero = false;
  for(int i=0; i<MortarEl->nNodes(); ++i) {
    if(mu[i] != 0) { isNonZero = true; break; }
  }
  if(!isNonZero) return MatShapeFctProd;

  GenVector<Scalar> NormalGeoGaps(nMortarShapeFct,Scalar(0.0));
  GenFullM<Scalar> PartialNormalGeoGapsPartialLocalCoordsOnFriendFaceEl(nMortarShapeFct,6);
  PartialNormalGeoGapsPartialLocalCoordsOnFriendFaceEl.zero();
  GenFullM<Scalar> PartialNormalGeoGapsPartialLocalCoordsOnFaceEl(nMortarShapeFct,6);
  PartialNormalGeoGapsPartialLocalCoordsOnFaceEl.zero();

  Scalar* MortarShape         = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* dShapexOnMortarEl   = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* dShapeyOnMortarEl   = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* d2ShapexOnMortarEl  = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* d2ShapeyOnMortarEl  = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* d2ShapexyOnMortarEl = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar mOnMortarEl[2], mOnFriendFaceEl[2];
  double m[2], r, s, t, weight, rr, ss, tt, rs, sr, rt, tr, st, ts;
  Scalar J, JNormal[3];
  Scalar MOnFriendFaceEl[3], MOnFaceEl[3];
  Scalar dMdxOnFriendFaceEl[3], dMdyOnFriendFaceEl[3], dMdxOnFaceEl[3], dMdyOnFaceEl[3];
  Scalar d2Mdx2OnFriendFaceEl[3], d2Mdy2OnFriendFaceEl[3], d2MdxdyOnFriendFaceEl[3],
         d2Mdx2OnFaceEl[3], d2Mdy2OnFaceEl[3], d2MdxdyOnFaceEl[3];

  // temporary storage for pre-computation of partial derivatives
  Scalar (*PartialMortarShapePartialLocalCoordsOnFaceEl)[6];
  PartialMortarShapePartialLocalCoordsOnFaceEl = (Scalar (*)[6]) dbg_alloca(nMortarShapeFct*6*sizeof(Scalar));
  Scalar PartialgPartialLocalCoordsOnFriendFaceEl[6], PartialgPartialLocalCoordsOnFaceEl[6];
  Scalar dJNormaldx[3], dJNormaldy[3], d2JNormaldx2[3], d2JNormaldy2[3], d2JNormaldxdy[3];
  Scalar SecondPartialgPartialLocalCoordsOnFriendFaceEl2[6][6], SecondPartialgPartialLocalCoordsOnFaceEl2[6][6]; 
  Scalar SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[6][6];
  Scalar (*SecondPartialMortarShapePartialLocalCoordsOnFaceEl2)[6][6];
  SecondPartialMortarShapePartialLocalCoordsOnFaceEl2 = (Scalar (*)[6][6]) dbg_alloca(nMortarShapeFct*36*sizeof(Scalar));

  if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
    // Get J*Normal on the face element supporting the CURRENT triangular facet
    FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, (Scalar*) NULL, cs);
    // Get the norm of J*Normal
    if(offset != 0) {
      J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
    }
  }

  int igp, i, j, k;
  
  for(igp=1; igp<=ngp; igp++) {

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    rr = r*r; ss = s*s; tt = t*t; rs = sr = r*s; rt = tr = r*t; st = ts = s*t;

    // Get local coord. on each face element
    // -> for the mortar elem. (see the NOTE section)
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl);

    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      // Jacobian*Normal on the face element supporting the CURRENT triangular facet
      FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
      FaceEl->ComputedJNormaldxAnddJNormaldy(dJNormaldx, dJNormaldy, mOnMortarEl, cs);
      FaceEl->Computed2JNormaldx2d2JNormaldy2Andd2JNormaldxdy(d2JNormaldx2, d2JNormaldy2, d2JNormaldxdy, mOnMortarEl, cs);
      // Get the norm of J*Normal
      if(offset != 0) {
        J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
      }
    }

    // Compute shape fcts & derivatives w.r.t m
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
    MortarEl->GetdShapeFct(dShapexOnMortarEl, dShapeyOnMortarEl, mOnMortarEl);

    for(i=0; i<nMortarShapeFct; i++) {
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][0] = dShapexOnMortarEl[i]*t;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][1] = dShapeyOnMortarEl[i]*t;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][2] = dShapexOnMortarEl[i]*r;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][3] = dShapeyOnMortarEl[i]*r;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][4] = dShapexOnMortarEl[i]*s;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][5] = dShapeyOnMortarEl[i]*s;
    }
    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      MortarEl->Getd2ShapeFct(d2ShapexOnMortarEl, d2ShapeyOnMortarEl, d2ShapexyOnMortarEl, mOnMortarEl);

      for(i=0; i<nMortarShapeFct; i++) {
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][0][0] =  d2ShapexOnMortarEl[i]*tt;

        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][1][0] = d2ShapexyOnMortarEl[i]*tt;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][1][1] =  d2ShapeyOnMortarEl[i]*tt;

        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][2][0] =  d2ShapexOnMortarEl[i]*rt;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][2][1] = d2ShapexyOnMortarEl[i]*rt;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][2][2] =  d2ShapexOnMortarEl[i]*rr;

        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][3][0] = d2ShapexyOnMortarEl[i]*rt;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][3][1] =  d2ShapeyOnMortarEl[i]*rt;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][3][2] = d2ShapexyOnMortarEl[i]*rr;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][3][3] =  d2ShapeyOnMortarEl[i]*rr;

        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][4][0] =  d2ShapexOnMortarEl[i]*st;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][4][1] = d2ShapexyOnMortarEl[i]*st;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][4][2] =  d2ShapexOnMortarEl[i]*sr;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][4][3] = d2ShapexyOnMortarEl[i]*sr;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][4][4] =  d2ShapexOnMortarEl[i]*ss;

        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][5][0] = d2ShapexyOnMortarEl[i]*st;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][5][1] =  d2ShapeyOnMortarEl[i]*st;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][5][2] = d2ShapexyOnMortarEl[i]*sr;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][5][3] =  d2ShapeyOnMortarEl[i]*sr;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][5][4] = d2ShapexyOnMortarEl[i]*ss;
        SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][5][5] =  d2ShapeyOnMortarEl[i]*ss;
      }
    }

    // Get global coord. on the element supporting this triangular facet & derivatives w.r.t. m
    FaceEl->LocalToGlobalCoord(MOnFaceEl, mOnMortarEl, cs);

    // Get global coord. on the element supporting the given triangular facet & derivatives w.r.t m
    FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs2);
    FriendFaceEl->ComputedMdxAnddMdy(dMdxOnFriendFaceEl, dMdyOnFriendFaceEl, mOnFriendFaceEl, cs2);

    //  ∂[g]/∂[FriendFacet.LocalCoordsOnFaceEl]
    Scalar a = (dMdxOnFriendFaceEl[0]*JNormal[0]+dMdxOnFriendFaceEl[1]*JNormal[1]+dMdxOnFriendFaceEl[2]*JNormal[2]);
    Scalar b = (dMdyOnFriendFaceEl[0]*JNormal[0]+dMdyOnFriendFaceEl[1]*JNormal[1]+dMdyOnFriendFaceEl[2]*JNormal[2]);
    PartialgPartialLocalCoordsOnFriendFaceEl[0] = -a*t;
    PartialgPartialLocalCoordsOnFriendFaceEl[1] = -b*t;
    PartialgPartialLocalCoordsOnFriendFaceEl[2] = -a*r;
    PartialgPartialLocalCoordsOnFriendFaceEl[3] = -b*r;
    PartialgPartialLocalCoordsOnFriendFaceEl[4] = -a*s;
    PartialgPartialLocalCoordsOnFriendFaceEl[5] = -b*s;

    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      Scalar a2 = (dMdxOnFriendFaceEl[0]*dJNormaldx[0]+dMdxOnFriendFaceEl[1]*dJNormaldx[1]+dMdxOnFriendFaceEl[2]*dJNormaldx[2]);
      Scalar b2 = (dMdyOnFriendFaceEl[0]*dJNormaldy[0]+dMdyOnFriendFaceEl[1]*dJNormaldy[1]+dMdyOnFriendFaceEl[2]*dJNormaldy[2]);
      Scalar ab = (dMdxOnFriendFaceEl[0]*dJNormaldy[0]+dMdxOnFriendFaceEl[1]*dJNormaldy[1]+dMdxOnFriendFaceEl[2]*dJNormaldy[2]);
      Scalar ba = (dMdyOnFriendFaceEl[0]*dJNormaldx[0]+dMdyOnFriendFaceEl[1]*dJNormaldx[1]+dMdyOnFriendFaceEl[2]*dJNormaldx[2]);
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[0][0] = -a2*tt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[0][1] = -ab*tt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[0][2] = -a2*tr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[0][3] = -ab*tr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[0][4] = -a2*ts;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[0][5] = -ab*ts;

      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[1][0] = -ba*tt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[1][1] = -b2*tt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[1][2] = -ba*tr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[1][3] = -b2*tr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[1][4] = -ba*ts;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[1][5] = -b2*ts;

      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[2][0] = -a2*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[2][1] = -ab*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[2][2] = -a2*rr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[2][3] = -ab*rr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[2][4] = -a2*rs;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[2][5] = -ab*rs;

      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[3][0] = -ba*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[3][1] = -b2*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[3][2] = -ba*rr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[3][3] = -b2*rr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[3][4] = -ba*rs;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[3][5] = -b2*rs;

      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[4][0] = -a2*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[4][1] = -ab*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[4][2] = -a2*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[4][3] = -ab*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[4][4] = -a2*ss;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[4][5] = -ab*ss;

      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[5][0] = -ba*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[5][1] = -b2*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[5][2] = -ba*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[5][3] = -b2*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[5][4] = -ba*ss;
      SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[5][5] = -b2*ss;
    }
    else {
      for(j=0; j<6; ++j) for(k=0; k<6; ++k) SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k] = 0;
    }

    if(FriendFaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      FriendFaceEl->Computed2Mdx2d2Mdy2Andd2Mdxdy(d2Mdx2OnFriendFaceEl, d2Mdy2OnFriendFaceEl, d2MdxdyOnFriendFaceEl, mOnFriendFaceEl, cs2);
      Scalar a2 = (d2Mdx2OnFriendFaceEl[0]*JNormal[0]+d2Mdx2OnFriendFaceEl[1]*JNormal[1]+d2Mdx2OnFriendFaceEl[2]*JNormal[2]);
      Scalar b2 = (d2Mdy2OnFriendFaceEl[0]*JNormal[0]+d2Mdy2OnFriendFaceEl[1]*JNormal[1]+d2Mdy2OnFriendFaceEl[2]*JNormal[2]);
      Scalar ab = (d2MdxdyOnFriendFaceEl[0]*JNormal[0]+d2MdxdyOnFriendFaceEl[1]*JNormal[1]+d2MdxdyOnFriendFaceEl[2]*JNormal[2]);
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[0][0] = -a2*tt;
    
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[1][0] = -ab*tt;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[1][1] = -b2*tt;

      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[2][0] = -a2*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[2][1] = -ab*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[2][2] = -a2*rr;

      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[3][0] = -ab*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[3][1] = -b2*rt;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[3][2] = -ab*rr;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[3][3] = -b2*rr;

      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[4][0] = -a2*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[4][1] = -ab*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[4][2] = -a2*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[4][3] = -ab*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[4][4] = -a2*ss;

      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[5][0] = -ab*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[5][1] = -b2*st;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[5][2] = -ab*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[5][3] = -b2*sr;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[5][4] = -ab*ss;
      SecondPartialgPartialLocalCoordsOnFriendFaceEl2[5][5] = -b2*ss;
    }

    // ∂[g]/∂[LocalCoordsOnFaceEl]  (note that the ∂[MOnFaceEl]/∂[LocalCoordsOnFaceEl]*JNormal term is always zero)
    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      Scalar gvec[3] = { MOnFaceEl[0]-MOnFriendFaceEl[0], MOnFaceEl[1]-MOnFriendFaceEl[1], MOnFaceEl[2]-MOnFriendFaceEl[2] };
      if(offset != 0) {
        gvec[0] += offset/J*JNormal[0];
        gvec[1] += offset/J*JNormal[1];
        gvec[2] += offset/J*JNormal[2];
      }
      // ∂[JNormal]/∂[LocalCoordsOnFaceEl] = ∂[JNormal]/∂[mOnMortarEl]*∂[mOnMortarEl]/∂[LocalCoordsOnFaceEl]
      Scalar c = gvec[0]*dJNormaldx[0] + gvec[1]*dJNormaldx[1] + gvec[2]*dJNormaldx[2];
      Scalar d = gvec[0]*dJNormaldy[0] + gvec[1]*dJNormaldy[1] + gvec[2]*dJNormaldy[2];
      PartialgPartialLocalCoordsOnFaceEl[0] = c*t;
      PartialgPartialLocalCoordsOnFaceEl[1] = d*t;
      PartialgPartialLocalCoordsOnFaceEl[2] = c*r;
      PartialgPartialLocalCoordsOnFaceEl[3] = d*r;
      PartialgPartialLocalCoordsOnFaceEl[4] = c*s;
      PartialgPartialLocalCoordsOnFaceEl[5] = d*s;
      // ∂^2[g]/∂[LocalCoordsOnFaceEl]^2 = 
      FaceEl->ComputedMdxAnddMdy(dMdxOnFaceEl, dMdyOnFaceEl, mOnMortarEl, cs);
      Scalar cc = gvec[0]*d2JNormaldx2[0] + gvec[1]*d2JNormaldx2[1] + gvec[2]*d2JNormaldx2[2] +
                  dMdxOnFaceEl[0]*dJNormaldx[0] + dMdxOnFaceEl[1]*dJNormaldx[1] + dMdxOnFaceEl[2]*dJNormaldx[2];
      Scalar dd = gvec[0]*d2JNormaldy2[0] + gvec[1]*d2JNormaldy2[1] + gvec[2]*d2JNormaldy2[2] +
                  dMdyOnFaceEl[0]*dJNormaldy[0] + dMdyOnFaceEl[1]*dJNormaldy[1] + dMdyOnFaceEl[2]*dJNormaldy[2];
      Scalar cd = gvec[0]*d2JNormaldxdy[0] + gvec[1]*d2JNormaldxdy[1] + gvec[2]*d2JNormaldxdy[2] +
                  dMdyOnFaceEl[0]*dJNormaldx[0] + dMdyOnFaceEl[1]*dJNormaldx[1] + dMdyOnFaceEl[2]*dJNormaldx[2];
      Scalar dc = gvec[0]*d2JNormaldxdy[0] + gvec[1]*d2JNormaldxdy[1] + gvec[2]*d2JNormaldxdy[2] +
                  dMdxOnFaceEl[0]*dJNormaldy[0] + dMdxOnFaceEl[1]*dJNormaldy[1] + dMdxOnFaceEl[2]*dJNormaldy[2];
      if(offset != 0) {
        // Normal = (1/J)*JNormal
        // dNormal/dx = d[1/J]/dx*JNormal + (1/J)*dJNormal/dx
        // 1/J = (J2)^-1/2
        // --> d[1/J]/dx = -0.5*(J2)^{-3/2}*d[J2]/dx = -/J3*[JNormal[0]*dJNormaldx[0]+JNormal[1]*dJNormaldx[1]+JNormal[2]*dJNormaldx[2])
        Scalar InvJ = 1/J;
        Scalar InvJ3 = 1/(J*J*J);
        Scalar dInvJdx = -InvJ3*(JNormal[0]*dJNormaldx[0]+JNormal[1]*dJNormaldx[1]+JNormal[2]*dJNormaldx[2]);
        Scalar dInvJdy = -InvJ3*(JNormal[0]*dJNormaldy[0]+JNormal[1]*dJNormaldy[1]+JNormal[2]*dJNormaldy[2]);
        Scalar dNormaldx[3] = { InvJ*dJNormaldx[0] + dInvJdx*JNormal[0],
                                InvJ*dJNormaldx[1] + dInvJdx*JNormal[1],
                                InvJ*dJNormaldx[2] + dInvJdx*JNormal[2] };
        Scalar dNormaldy[3] = { InvJ*dJNormaldy[0] + dInvJdy*JNormal[0],
                                InvJ*dJNormaldy[1] + dInvJdy*JNormal[1],
                                InvJ*dJNormaldy[2] + dInvJdy*JNormal[2] };
        cc += offset*(dNormaldx[0]*dJNormaldx[0] + dNormaldx[1]*dJNormaldx[1] + dNormaldx[2]*dJNormaldx[2]);
        dd += offset*(dNormaldy[0]*dJNormaldy[0] + dNormaldy[1]*dJNormaldy[1] + dNormaldy[2]*dJNormaldy[2]);
        cd += offset*(dNormaldy[0]*dJNormaldx[0] + dNormaldy[1]*dJNormaldx[1] + dNormaldy[2]*dJNormaldx[2]);
        dc += offset*(dNormaldx[0]*dJNormaldy[0] + dNormaldx[1]*dJNormaldy[1] + dNormaldx[2]*dJNormaldy[2]);
      }
      SecondPartialgPartialLocalCoordsOnFaceEl2[0][0] = cc*tt;

      SecondPartialgPartialLocalCoordsOnFaceEl2[1][0] = dc*tt;
      SecondPartialgPartialLocalCoordsOnFaceEl2[1][1] = dd*tt;

      SecondPartialgPartialLocalCoordsOnFaceEl2[2][0] = cc*rt;
      SecondPartialgPartialLocalCoordsOnFaceEl2[2][1] = cd*rt;
      SecondPartialgPartialLocalCoordsOnFaceEl2[2][2] = cc*rr;

      SecondPartialgPartialLocalCoordsOnFaceEl2[3][0] = dc*rt;
      SecondPartialgPartialLocalCoordsOnFaceEl2[3][1] = dd*rt;
      SecondPartialgPartialLocalCoordsOnFaceEl2[3][2] = dc*rr;
      SecondPartialgPartialLocalCoordsOnFaceEl2[3][3] = dd*rr;

      SecondPartialgPartialLocalCoordsOnFaceEl2[4][0] = cc*st;
      SecondPartialgPartialLocalCoordsOnFaceEl2[4][1] = cd*st;
      SecondPartialgPartialLocalCoordsOnFaceEl2[4][2] = cc*sr;
      SecondPartialgPartialLocalCoordsOnFaceEl2[4][3] = cd*sr;
      SecondPartialgPartialLocalCoordsOnFaceEl2[4][4] = cc*ss;

      SecondPartialgPartialLocalCoordsOnFaceEl2[5][0] = dc*st;
      SecondPartialgPartialLocalCoordsOnFaceEl2[5][1] = dd*st;
      SecondPartialgPartialLocalCoordsOnFaceEl2[5][2] = dc*sr;
      SecondPartialgPartialLocalCoordsOnFaceEl2[5][3] = dd*sr;
      SecondPartialgPartialLocalCoordsOnFaceEl2[5][4] = dc*ss;
      SecondPartialgPartialLocalCoordsOnFaceEl2[5][5] = dd*ss;
    }
    else {
      for(j=0; j<6; j++) PartialgPartialLocalCoordsOnFaceEl[j] = 0;
    }

    // Shape fcts product & integration
    Scalar g = (MOnFaceEl[0]-MOnFriendFaceEl[0])*JNormal[0] +(MOnFaceEl[1]-MOnFriendFaceEl[1])*JNormal[1] +(MOnFaceEl[2]-MOnFriendFaceEl[2])*JNormal[2];
    if(offset != 0) g += J*offset;

    for(i=0; i<nMortarShapeFct; i++) {
      if(mu[i] == 0) continue;
      NormalGeoGaps[i] += weight*MortarShape[i]*g;
      for(j=0; j<6; j++) { // ∂[NormalGeoGaps]/∂[FriendFacet.LocalCoordsOnFaceEl]
        PartialNormalGeoGapsPartialLocalCoordsOnFriendFaceEl[i][j] += weight*MortarShape[i]*PartialgPartialLocalCoordsOnFriendFaceEl[j];
      }
      for(j=0; j<6; ++j) { // ∂[NormalGeoGaps]/∂[LocalCoordsOnFaceEl] 
        PartialNormalGeoGapsPartialLocalCoordsOnFaceEl[i][j] += weight*(MortarShape[i]*PartialgPartialLocalCoordsOnFaceEl[j]+
                                                                        PartialMortarShapePartialLocalCoordsOnFaceEl[i][j]*g);
      }
      // working with lower triangular part only, due to symmetry
      if(FriendFaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
        for(j=0; j<6; j++) { // ∂^2[NormalGeoGaps]/∂[FriendFacet.LocalCoordsOnFaceEl]^2
          for(k=0; k<=j; ++k) {
            MatShapeFctProd[j][k] += weight*mu[i]*MortarShape[i]*SecondPartialgPartialLocalCoordsOnFriendFaceEl2[j][k];
          }
        }
      }
      for(j=0; j<6; j++) { // ∂^2[NormalGeoGaps]/∂[FriendFacet.LocalCoordsOnFaceEl]∂[LocalCoordsOnFaceEl]
        for(k=0; k<6; ++k) {
          (MatShapeFctProd[6+k][j] += weight*mu[i]*(PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialgPartialLocalCoordsOnFriendFaceEl[j]+
                                                   MortarShape[i]*SecondPartialgPartialLocalCoordsOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k]));
        }
      }
      if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
        for(j=0; j<6; ++j) { // ∂^2[NormalGeoGaps]/∂[LocalCoordsOnFaceEl]^2
          for(k=0; k<=j; ++k) {
            MatShapeFctProd[6+j][6+k] += weight*mu[i]*(MortarShape[i]*SecondPartialgPartialLocalCoordsOnFaceEl2[j][k] 
                                                     + SecondPartialMortarShapePartialLocalCoordsOnFaceEl2[i][j][k]*g
                                                     + PartialMortarShapePartialLocalCoordsOnFaceEl[i][j]*PartialgPartialLocalCoordsOnFaceEl[k]
                                                     + PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialgPartialLocalCoordsOnFaceEl[j]);
          }
        }
      }
    }
  }

  Scalar dA = (LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0])*(LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1])
             -(LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0])*(LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]);

  Scalar PartialdAPartialLocalCoordsOnFaceEl[6] = {
     -(LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1]) + (LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]),
     -(LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0]) + (LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0]),
      (LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1]),
     -(LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0]),
     -(LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]),
      (LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0]) };

  Scalar SecondPartialdAPartialLocalCoordsOnFaceEl2[6][6] = {
    {  0,  0,  0,  1,  0, -1 },
    {  0,  0, -1,  0,  1,  0 }, 
    {  0, -1,  0,  0,  0,  1 },
    {  1,  0,  0,  0, -1,  0 },
    {  0,  1,  0, -1,  0,  0 },
    { -1,  0,  1,  0,  0,  0 }
  }; 

  if(dA < 0) {
    dA = -dA;
    for(i=0; i<6; ++i) {
      PartialdAPartialLocalCoordsOnFaceEl[i] = -PartialdAPartialLocalCoordsOnFaceEl[i];
      for(j=0; j<6; ++j) SecondPartialdAPartialLocalCoordsOnFaceEl2[i][j] = -SecondPartialdAPartialLocalCoordsOnFaceEl2[i][j];
    }
  }

  MatShapeFctProd *= dA;
  for(i=0; i<nMortarShapeFct; ++i) { // work with lower-triangular part only, due to symmetry
    for(j=0; j<6; ++j) {
      for(k=0; k<6; ++k) {
        (MatShapeFctProd[6+k][j] += mu[i]*PartialNormalGeoGapsPartialLocalCoordsOnFriendFaceEl[i][j]*PartialdAPartialLocalCoordsOnFaceEl[k]);
      }
      for(k=0; k<=j; ++k) {
        MatShapeFctProd[6+j][6+k] += mu[i]*(NormalGeoGaps[i]*SecondPartialdAPartialLocalCoordsOnFaceEl2[j][k]+
                                            PartialNormalGeoGapsPartialLocalCoordsOnFaceEl[i][j]*PartialdAPartialLocalCoordsOnFaceEl[k] +
                                            PartialNormalGeoGapsPartialLocalCoordsOnFaceEl[i][k]*PartialdAPartialLocalCoordsOnFaceEl[j]);
      }
    }
  }

  // now, fill in the upper-triangular part of the symmetric matrix
  for(j=0; j<12; ++j)
    for(k=j; k<12; k++) MatShapeFctProd[j][k] = MatShapeFctProd[k][j];

  return MatShapeFctProd;
}

template<class Scalar, class FaceType>
template<typename CoordSetType, typename MortarType, typename FriendFaceType>
GenFullM<Scalar>
TriFacetTemplate<Scalar, FaceType>::IntegrateLocalGradGradNormalGeoGap(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                                       CoordSetType& cs, CoordSetType& cs2, double *mu, int ngp, double offset)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the gradients of the products of each shape function defined
// by the given MortarElement and the gap vector between the element associated with the current
// triangular facet and the element associated to the given (Friend) triangular facet times the normal
// (see Notes (2)).
// -> C.row(i) = ∇ ∫_{current TriFacet}[Mortar(i)*Gap·normal(current TriFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED WITH THE GIVEN
//                                       (FRIEND) TRIFACET + THE NUMBER OF DOFs OF THE FACE ELEMENT
//                                       ASSOCIATED WITH THE CURRENT TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE Mortar Element LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT TRIANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT TRIANGULAR FACET
//     (3) CURRENTLY ASSUMING LocalCoordOnFaceEl IS CONSTANT
// ******************************************************************************************************
{
  // Get ptr to the face element supporting the given triangular facet
  FriendFaceType* FriendFaceEl = FriendFacet.GetPtrFaceEl();

  int nMortarShapeFct         = MortarEl->nNodes();
  int nShapeFctOnFriendFaceEl = FriendFaceEl->nNodes();
  int nShapeFctOnFaceEl       = FaceEl->nNodes();

  GenFullM<Scalar> MatShapeFctProd(3*nShapeFctOnFriendFaceEl + 3*nShapeFctOnFaceEl, 12);
  MatShapeFctProd.zero();
  GenVector<Scalar> GradNormalGeoGaps(3*nShapeFctOnFriendFaceEl + 3*nShapeFctOnFaceEl, 0.0);

  Scalar* MortarShape           = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* dShapexOnMortarEl     = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* dShapeyOnMortarEl     = (Scalar*) dbg_alloca(nMortarShapeFct*sizeof(Scalar));
  Scalar* ShapeOnFaceEl         = (Scalar*) dbg_alloca(nShapeFctOnFaceEl*sizeof(Scalar));
  Scalar* dShapexOnFaceEl       = (Scalar*) dbg_alloca(nShapeFctOnFaceEl*sizeof(Scalar));
  Scalar* dShapeyOnFaceEl       = (Scalar*) dbg_alloca(nShapeFctOnFaceEl*sizeof(Scalar));
  Scalar* ShapeOnFriendFaceEl   = (Scalar*) dbg_alloca(nShapeFctOnFriendFaceEl*sizeof(Scalar));
  Scalar* dShapexOnFriendFaceEl = (Scalar*) dbg_alloca(nShapeFctOnFriendFaceEl*sizeof(Scalar));
  Scalar* dShapeyOnFriendFaceEl = (Scalar*) dbg_alloca(nShapeFctOnFriendFaceEl*sizeof(Scalar));
  Scalar mOnMortarEl[2], mOnFriendFaceEl[2];
  double m[2], r, s, t, weight;
  Scalar J, JNormal[3], dJNormaldx[3], dJNormaldy[3];
  Scalar MOnFriendFaceEl[3], MOnFaceEl[3], (*dJNormal)[3], (*ddJNormaldx)[3], (*ddJNormaldy)[3];
  dJNormal = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*9*sizeof(Scalar));
  ddJNormaldx = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*9*sizeof(Scalar));
  ddJNormaldy = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*9*sizeof(Scalar));
  Scalar dMdxOnFriendFaceEl[3], dMdyOnFriendFaceEl[3], dMdxOnFaceEl[3], dMdyOnFaceEl[3];

  // temporary storage for pre-computation of partial derivatives
  Scalar (*PartialOnFriendFaceEl)[3], (*PartialOnFaceEl)[3];
  PartialOnFriendFaceEl = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFriendFaceEl*3*sizeof(Scalar));
  PartialOnFaceEl       = (Scalar (*)[3]) dbg_alloca(nShapeFctOnFaceEl*3*sizeof(Scalar));
  Scalar (*PartialMortarShapePartialLocalCoordsOnFaceEl)[6];
  PartialMortarShapePartialLocalCoordsOnFaceEl = (Scalar (*)[6]) dbg_alloca(nMortarShapeFct*6*sizeof(Scalar));

  Scalar (*PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl)[3][6];
  Scalar (*PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl)[3][6];
  Scalar (*PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl)[3][6];
  Scalar (*PartialPartialOnFaceElPartialLocalCoordsOnFaceEl)[3][6];
  PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl = (Scalar (*)[3][6]) dbg_alloca(nShapeFctOnFriendFaceEl*18*sizeof(Scalar));
  PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl       = (Scalar (*)[3][6]) dbg_alloca(nShapeFctOnFriendFaceEl*18*sizeof(Scalar));
  PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl       = (Scalar (*)[3][6]) dbg_alloca(nShapeFctOnFaceEl*18*sizeof(Scalar));
  PartialPartialOnFaceElPartialLocalCoordsOnFaceEl             = (Scalar (*)[3][6]) dbg_alloca(nShapeFctOnFaceEl*18*sizeof(Scalar));

  //Scalar dA = MappingJacobian();
  if(FaceEl->GetFaceElemType() == FaceElement::TRIFACEL3) {
    // Get J*Normal on the face element supporting the CURRENT triangular facet
    FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, (Scalar*) NULL, cs);
    // Get the derivative of J*Normal
    FaceEl->GetdJNormal(dJNormal, (Scalar*) NULL, cs);
    // Get the norm of J*Normal
    if(offset != 0) {
      J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
    }
  }

  int igp, i, j, k, l;
  
  for(igp=1; igp<=ngp; igp++) {

    getGaussPtOnTriangle(ngp,igp,r,s,t,weight);
    m[0] = r; m[1] = s;

    // Get local coord. on each face element
    // -> for the mortar elem. (see the NOTE section)
    LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFriendFaceEl);

    if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
      // Jacobian*Normal on the face element supporting the CURRENT triangular facet
      FaceEl->GetIsoParamMappingNormalJacobianProduct(JNormal, mOnMortarEl, cs);
      // Get the derivative of Jacobian*Normal
      FaceEl->GetdJNormal(dJNormal, mOnMortarEl, cs);
      // Get the norm of J*Normal
      if(offset != 0) {
        J = sqrt(JNormal[0]*JNormal[0]+JNormal[1]*JNormal[1]+JNormal[2]*JNormal[2]);
      }
    }

    // Compute shape fcts and their derivatives
    // -> mortar elem.
    MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
    MortarEl->GetdShapeFct(dShapexOnMortarEl, dShapeyOnMortarEl, mOnMortarEl);

    for(i=0; i<nMortarShapeFct; i++) {
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][0] = dShapexOnMortarEl[i]*t;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][1] = dShapeyOnMortarEl[i]*t;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][2] = dShapexOnMortarEl[i]*r;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][3] = dShapeyOnMortarEl[i]*r;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][4] = dShapexOnMortarEl[i]*s;
      PartialMortarShapePartialLocalCoordsOnFaceEl[i][5] = dShapeyOnMortarEl[i]*s;
    }

    // -> for the element supporting this triangular facet
    FaceEl->GetShapeFctVal(ShapeOnFaceEl, mOnMortarEl);
    FaceEl->GetdShapeFct(dShapexOnFaceEl, dShapeyOnFaceEl, mOnMortarEl);
    // -> for the element supporting the given triangular facet
    FriendFaceEl->GetShapeFctVal(ShapeOnFriendFaceEl, mOnFriendFaceEl);
    FriendFaceEl->GetdShapeFct(dShapexOnFriendFaceEl, dShapeyOnFriendFaceEl, mOnFriendFaceEl);

    // Get global coord. on the element supporting this triangular facet
    FaceEl->LocalToGlobalCoord(MOnFaceEl, mOnMortarEl, cs);
    FaceEl->ComputedMdxAnddMdy(dMdxOnFaceEl, dMdyOnFaceEl, mOnMortarEl, cs);

    // Get global coord. on the element supporting the given triangular facet
    FriendFaceEl->LocalToGlobalCoord(MOnFriendFaceEl, mOnFriendFaceEl, cs2);
    FriendFaceEl->ComputedMdxAnddMdy(dMdxOnFriendFaceEl, dMdyOnFriendFaceEl, mOnFriendFaceEl, cs2);

    // Shape fcts product & integration (partially optimized implementation)
    for(j=0; j<nShapeFctOnFriendFaceEl; j++) {
      for(k=0; k<3; ++k) {
        PartialOnFriendFaceEl[j][k] = ShapeOnFriendFaceEl[j]*JNormal[k];

        PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][k][0] = dShapexOnFriendFaceEl[j]*JNormal[k]*t;
        PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][k][1] = dShapeyOnFriendFaceEl[j]*JNormal[k]*t;
        PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][k][2] = dShapexOnFriendFaceEl[j]*JNormal[k]*r;
        PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][k][3] = dShapeyOnFriendFaceEl[j]*JNormal[k]*r;
        PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][k][4] = dShapexOnFriendFaceEl[j]*JNormal[k]*s;
        PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][k][5] = dShapeyOnFriendFaceEl[j]*JNormal[k]*s;

        if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
          FaceEl->ComputedJNormaldxAnddJNormaldy(dJNormaldx, dJNormaldy, mOnMortarEl, cs);
          PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k][0] = ShapeOnFriendFaceEl[j]*dJNormaldx[k]*t;
          PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k][1] = ShapeOnFriendFaceEl[j]*dJNormaldy[k]*t;
          PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k][2] = ShapeOnFriendFaceEl[j]*dJNormaldx[k]*r;
          PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k][3] = ShapeOnFriendFaceEl[j]*dJNormaldy[k]*r;
          PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k][4] = ShapeOnFriendFaceEl[j]*dJNormaldx[k]*s;
          PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k][5] = ShapeOnFriendFaceEl[j]*dJNormaldy[k]*s;
        }
        else {
          for(l=0; l<6; ++l) PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][k][l] = 0;
        }
      }
    }

    // note that the ∂[MOnFaceEl]/∂[LocalCoordsOnFaceEl] term is always zero
    for(j=0; j<nShapeFctOnFaceEl; j++) {
      for(k=0; k<3; ++k) {
        PartialOnFaceEl[j][k] = (ShapeOnFaceEl[j]*JNormal[k]+(MOnFaceEl[0]-MOnFriendFaceEl[0])*dJNormal[3*j+k][0]
                                                            +(MOnFaceEl[1]-MOnFriendFaceEl[1])*dJNormal[3*j+k][1]
                                                            +(MOnFaceEl[2]-MOnFriendFaceEl[2])*dJNormal[3*j+k][2]);
        if(offset != 0) {
          // note: for TRIFACEL3 this term can be pre-computed
          PartialOnFaceEl[j][k] += 1/J*(JNormal[0]*dJNormal[3*j+k][0]+JNormal[1]*dJNormal[3*j+k][1]+JNormal[2]*dJNormal[3*j+k][2])*offset;
        }

        Scalar a = dMdxOnFriendFaceEl[0]*dJNormal[3*j+k][0]+dMdxOnFriendFaceEl[1]*dJNormal[3*j+k][1]+dMdxOnFriendFaceEl[2]*dJNormal[3*j+k][2];
        Scalar b = dMdyOnFriendFaceEl[0]*dJNormal[3*j+k][0]+dMdyOnFriendFaceEl[1]*dJNormal[3*j+k][1]+dMdyOnFriendFaceEl[2]*dJNormal[3*j+k][2];
        PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][k][0] = -a*t;
        PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][k][1] = -b*t;
        PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][k][2] = -a*r;
        PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][k][3] = -b*r;
        PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][k][4] = -a*s;
        PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][k][5] = -b*s;

        Scalar c = dShapexOnFaceEl[j]*JNormal[k] + dMdxOnFaceEl[0]*dJNormal[3*j+k][0]+dMdxOnFaceEl[1]*dJNormal[3*j+k][1]+dMdxOnFaceEl[2]*dJNormal[3*j+k][2];
        Scalar d = dShapeyOnFaceEl[j]*JNormal[k] + dMdyOnFaceEl[0]*dJNormal[3*j+k][0]+dMdyOnFaceEl[1]*dJNormal[3*j+k][1]+dMdyOnFaceEl[2]*dJNormal[3*j+k][2];
        PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][0] = c*t;
        PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][1] = d*t;
        PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][2] = c*r;
        PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][3] = d*r;
        PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][4] = c*s;
        PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][5] = d*s;
      }

      if(FaceEl->GetFaceElemType() != FaceElement::TRIFACEL3) {
        FaceEl->ComputeddJNormaldxAndddJNormaldy(ddJNormaldx, ddJNormaldy, mOnMortarEl, cs);
        for(k=0; k<3; ++k) {
          Scalar e = ShapeOnFaceEl[j]*dJNormaldx[k] + (MOnFaceEl[0]-MOnFriendFaceEl[0])*ddJNormaldx[3*j+k][0] 
                   + (MOnFaceEl[1]-MOnFriendFaceEl[1])*ddJNormaldx[3*j+k][1] + (MOnFaceEl[2]-MOnFriendFaceEl[2])*ddJNormaldx[3*j+k][2];
          Scalar f = ShapeOnFaceEl[j]*dJNormaldy[k] + (MOnFaceEl[0]-MOnFriendFaceEl[0])*ddJNormaldy[3*j+k][0]
                   + (MOnFaceEl[1]-MOnFriendFaceEl[1])*ddJNormaldy[3*j+k][1] + (MOnFaceEl[2]-MOnFriendFaceEl[2])*ddJNormaldy[3*j+k][2];
          PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][0] += e*t;
          PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][1] += f*t;
          PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][2] += e*r;
          PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][3] += f*r;
          PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][4] += e*s;
          PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][5] += f*s;
        }
        if(offset != 0) {
          Scalar InvJ = 1/J;
          Scalar InvJ3 = 1/(J*J*J);
          Scalar dInvJdx = -InvJ3*(JNormal[0]*dJNormaldx[0]+JNormal[1]*dJNormaldx[1]+JNormal[2]*dJNormaldx[2]);
          Scalar dInvJdy = -InvJ3*(JNormal[0]*dJNormaldy[0]+JNormal[1]*dJNormaldy[1]+JNormal[2]*dJNormaldy[2]);
          Scalar Normal[3] = { InvJ*JNormal[0], InvJ*JNormal[1], InvJ*JNormal[2] };
          Scalar dNormaldx[3] = { InvJ*dJNormaldx[0] + dInvJdx*JNormal[0],
                                  InvJ*dJNormaldx[1] + dInvJdx*JNormal[1],
                                  InvJ*dJNormaldx[2] + dInvJdx*JNormal[2] };
          Scalar dNormaldy[3] = { InvJ*dJNormaldy[0] + dInvJdy*JNormal[0],
                                  InvJ*dJNormaldy[1] + dInvJdy*JNormal[1],
                                  InvJ*dJNormaldy[2] + dInvJdy*JNormal[2] };
          for(k=0; k<3; ++k) {
            Scalar ee = offset*(dNormaldx[0]*dJNormal[3*j+k][0] + Normal[0]*ddJNormaldx[3*j+k][0] 
                               +dNormaldx[1]*dJNormal[3*j+k][1] + Normal[1]*ddJNormaldx[3*j+k][1]
                               +dNormaldx[2]*dJNormal[3*j+k][2] + Normal[2]*ddJNormaldx[3*j+k][2]);
            Scalar ff = offset*(dNormaldy[0]*dJNormal[3*j+k][0] + Normal[0]*ddJNormaldy[3*j+k][0] 
                               +dNormaldy[1]*dJNormal[3*j+k][1] + Normal[1]*ddJNormaldy[3*j+k][1]
                               +dNormaldy[2]*dJNormal[3*j+k][2] + Normal[2]*ddJNormaldy[3*j+k][2]);
            PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][0] += ee*t;
            PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][1] += ff*t;
            PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][2] += ee*r;
            PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][3] += ff*r;
            PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][4] += ee*s;
            PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][k][5] += ff*s;
          }
        }
      }
    }

    for(i=0; i<nMortarShapeFct; i++) {

      for(j=0; j<nShapeFctOnFriendFaceEl; j++) {
        GradNormalGeoGaps[3*j  ] -= weight*mu[i]*MortarShape[i]*PartialOnFriendFaceEl[j][0];
        GradNormalGeoGaps[3*j+1] -= weight*mu[i]*MortarShape[i]*PartialOnFriendFaceEl[j][1];
        GradNormalGeoGaps[3*j+2] -= weight*mu[i]*MortarShape[i]*PartialOnFriendFaceEl[j][2];
      }
      for(j=0; j<nShapeFctOnFaceEl; j++) {
        int J = nShapeFctOnFriendFaceEl+j;
        GradNormalGeoGaps[3*J  ] += weight*mu[i]*MortarShape[i]*PartialOnFaceEl[j][0];
        GradNormalGeoGaps[3*J+1] += weight*mu[i]*MortarShape[i]*PartialOnFaceEl[j][1];
        GradNormalGeoGaps[3*J+2] += weight*mu[i]*MortarShape[i]*PartialOnFaceEl[j][2];
      }

      for(j=0; j<nShapeFctOnFriendFaceEl; j++) {
        for(k=0; k<6; ++k) {
          MatShapeFctProd[3*j  ][k] -= weight*mu[i]*( MortarShape[i]*PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][0][k] );
          MatShapeFctProd[3*j+1][k] -= weight*mu[i]*( MortarShape[i]*PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][1][k] );
          MatShapeFctProd[3*j+2][k] -= weight*mu[i]*( MortarShape[i]*PartialPartialOnFriendFaceElPartialLocalCoordsOnFriendFaceEl[j][2][k] );
        }
        for(k=0; k<6; ++k) {
          MatShapeFctProd[3*j  ][6+k] -= weight*mu[i]*( PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialOnFriendFaceEl[j][0] 
                                                       +MortarShape[i]*PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][0][k] );
          MatShapeFctProd[3*j+1][6+k] -= weight*mu[i]*( PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialOnFriendFaceEl[j][1]
                                                       +MortarShape[i]*PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][1][k] );
          MatShapeFctProd[3*j+2][6+k] -= weight*mu[i]*( PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialOnFriendFaceEl[j][2]
                                                       +MortarShape[i]*PartialPartialOnFriendFaceElPartialLocalCoordsOnFaceEl[j][2][k] );
        }
      }

      for(j=0; j<nShapeFctOnFaceEl; j++) {
        int J = nShapeFctOnFriendFaceEl+j;
        for(k=0; k<6; ++k) {
          MatShapeFctProd[3*J  ][k] += weight*mu[i]*( MortarShape[i]*PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][0][k] );
          MatShapeFctProd[3*J+1][k] += weight*mu[i]*( MortarShape[i]*PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][1][k] );
          MatShapeFctProd[3*J+2][k] += weight*mu[i]*( MortarShape[i]*PartialPartialOnFaceElPartialLocalCoordsOnFriendFaceEl[j][2][k] );
        }
        for(k=0; k<6; ++k) {
          MatShapeFctProd[3*J  ][6+k] += weight*mu[i]*( PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialOnFaceEl[j][0] 
                                                       +MortarShape[i]*PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][0][k] );
          MatShapeFctProd[3*J+1][6+k] += weight*mu[i]*( PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialOnFaceEl[j][1]      
                                                       +MortarShape[i]*PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][1][k] );
          MatShapeFctProd[3*J+2][6+k] += weight*mu[i]*( PartialMortarShapePartialLocalCoordsOnFaceEl[i][k]*PartialOnFaceEl[j][2]      
                                                       +MortarShape[i]*PartialPartialOnFaceElPartialLocalCoordsOnFaceEl[j][2][k] );
        }
      }    
    }
  }

  Scalar dA = (LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0])*(LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1])
             -(LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0])*(LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]);

  Scalar PartialdAPartialLocalCoordsOnFaceEl[6] = {
     -(LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1]) + (LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]),
     -(LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0]) + (LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0]),
      (LocalCoordOnFaceEl[2][1]-LocalCoordOnFaceEl[0][1]),
     -(LocalCoordOnFaceEl[2][0]-LocalCoordOnFaceEl[0][0]),
     -(LocalCoordOnFaceEl[1][1]-LocalCoordOnFaceEl[0][1]),
      (LocalCoordOnFaceEl[1][0]-LocalCoordOnFaceEl[0][0]) };

  if(dA < 0) {
    dA = -dA;
    for(i=0; i<6; ++i) PartialdAPartialLocalCoordsOnFaceEl[i] = -PartialdAPartialLocalCoordsOnFaceEl[i];
  }

  MatShapeFctProd *= dA;
  for(i=0; i<3*nShapeFctOnFriendFaceEl + 3*nShapeFctOnFaceEl; ++i)
    for(j=0; j<6; ++j) MatShapeFctProd[i][6+j] += GradNormalGeoGaps[i]*PartialdAPartialLocalCoordsOnFaceEl[j];

  return MatShapeFctProd;
}
#endif
