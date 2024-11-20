// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <alloca.h>
#include <cmath>

// FEM headers
#include <Element.d/Element.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FFIPolygon.d/QuadFacet.h>

// External routines
extern void getGaussPtOnQuadrangle(int, int, double&, double&, double&, double&);

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------
QuadFacet::QuadFacet()
{
  Initialize();
}

QuadFacet::QuadFacet(FaceElement* FaceElem, const double* m1, const double* m2, const double* m3, const double* m4)
{
  // set area to zero
  Area   = 0.0;
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];
  LocalCoordOnFaceEl[3][0] = m4[0]; LocalCoordOnFaceEl[3][1] = m4[1];
  // set ptr to face el.
  FaceEl = FaceElem;
}

// -----------------------------------------------------------------------------------------------------
//                                             DESTRUCTORS 
// -----------------------------------------------------------------------------------------------------
QuadFacet::~QuadFacet()
{
  // set area to zero
  Area   = 0.0;
  // set ptr to face el. to NULL
  FaceEl = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                    INITIALIZATION & CLEAR/CLEAN METHODS
// -----------------------------------------------------------------------------------------------------
void
QuadFacet::Initialize()
{
  // set area to zero
  Area   = 0.0;
  // set local coord. of triangular vertices in face el. to 0.
  LocalCoordOnFaceEl[0][0] = 0.0; LocalCoordOnFaceEl[0][1] = 0.0;
  LocalCoordOnFaceEl[1][0] = 0.0; LocalCoordOnFaceEl[1][1] = 0.0;
  LocalCoordOnFaceEl[2][0] = 0.0; LocalCoordOnFaceEl[2][1] = 0.0;
  LocalCoordOnFaceEl[3][0] = 0.0; LocalCoordOnFaceEl[3][1] = 0.0;
  // set ptr to face el. to NULL
  FaceEl = 0;
}
// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS 
// -----------------------------------------------------------------------------------------------------
void
QuadFacet::SetQuadFacet(FaceElement* FaceElem, const double* m1, const double* m2, const double* m3, const double* m4)
{
  // set area to zero
  Area = 0.0;
  // set local coord. of triangular vertices in face el.
  LocalCoordOnFaceEl[0][0] = m1[0]; LocalCoordOnFaceEl[0][1] = m1[1];
  LocalCoordOnFaceEl[1][0] = m2[0]; LocalCoordOnFaceEl[1][1] = m2[1];
  LocalCoordOnFaceEl[2][0] = m3[0]; LocalCoordOnFaceEl[2][1] = m3[1];
  LocalCoordOnFaceEl[3][0] = m4[0]; LocalCoordOnFaceEl[3][1] = m4[1];
  // set ptr to face el.
  FaceEl = FaceElem;
}

void
QuadFacet::SetArea(double S) { Area = S; };

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
void
QuadFacet::Print()
{
  fprintf(stderr," * QuadFacet area = %e\n", Area);
  fprintf(stderr," * QuadFacet vertices coord. (on ref. face elem.):\n");
  for(int i=0; i<4; i++)
    fprintf(stderr,"  -> vertex %d: x = %e, y = %e\n",i+1,LocalCoordOnFaceEl[i][0],LocalCoordOnFaceEl[i][1]);
  fprintf(stderr," * QuadFacet ptr to face elem %p\n",FaceEl);
}

/*
void
QuadFacet::PrintVerticesXYZ(FILE* file=stderr, CoordSet& cs, int& firstVertId)
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
QuadFacet::PrintQuadFacetTopo(FILE* file=stderr, int& QuadFacetId, int& firstVertId)
{
  fprintf(file," %6d  %6e  %6e  %6e\n",QuadFacetId,firstVertId,firstVertId+1,firstVertId+2);
  QuadFacetId++;
  firstVertId += 3;
}
*/

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*double
QuadFacet::ComputeApproxArea(CoordSet &cs)
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
//                         REF. QUADRANGULAR FACET -> REF. FACE El. MAPPING METHODS 
// -----------------------------------------------------------------------------------------------------
/*void
QuadFacet::MappingShapeFctAndDerivative(double* Shape, double* dShape, double* m)
{
  if(Shape==0)  Shape  = new double[3];
  if(dShape==0) dShape = new double[2][3];

  Shape[0] = m[0]; Shape[1] = m[1]; Shape[2] = 1.-m[0]-m[1];
  dShape[0][0] =  1.; dShape[0][1] =  0.; 
  dShape[1][0] =  0.; dShape[1][1] =  1.; 
  dShape[2][0] = -1.; dShape[2][1] = -1.; 
}*/

// Return jacobian of the mapping ref. QuadFacet -> ref. face el.
double
QuadFacet::MappingJacobian(const double* m)
{
  double d1 = 0.5*(1.0+m[0]);
  double d2 = 0.5*(1.0+m[1]);
  double d3 = 1.0-d1;
  double d4 = 1.0-d2;

  double Shape0 = d3*d4;
  double Shape1 = d4*d1;
  double Shape2 = d1*d2;
  double Shape3 = d2*d3;

  double X[4] = {LocalCoordOnFaceEl[0][0],LocalCoordOnFaceEl[1][0],LocalCoordOnFaceEl[2][0],LocalCoordOnFaceEl[3][0]};
  double Y[4] = {LocalCoordOnFaceEl[0][1],LocalCoordOnFaceEl[1][1],LocalCoordOnFaceEl[2][1],LocalCoordOnFaceEl[3][1]};

  double c0 = (X[1]-X[0])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[1]-Y[0]);
  double c1 = (X[1]-X[0])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[1]-Y[0]);
  double c2 = (X[2]-X[3])*(Y[2]-Y[1]) - (X[2]-X[1])*(Y[2]-Y[3]);
  double c3 = (X[2]-X[3])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[2]-Y[3]);

  double J = 0.25*(Shape0*c0+Shape1*c1+Shape2*c2+Shape3*c3);

  return fabs(J);// to avoid swapping the vertex to have a positive jacobian
                 // this is OK because we just use it as the differential area
                 // in the integration formula
}

// Map the given point m in the ref. Trifacet to the associated point in the ref. face el. 
void
QuadFacet::LocalToLocalCoordOnFaceEl(const double* m, double* mOnFaceEl)
{
  double d1 = 0.5*(1.0+m[0]);
  double d2 = 0.5*(1.0+m[1]);
  double d3 = 1.0-d1;
  double d4 = 1.0-d2;

  double Shape0 = d3*d4;
  double Shape1 = d4*d1;
  double Shape2 = d1*d2;
  double Shape3 = d2*d3;

  mOnFaceEl[0] = Shape0*LocalCoordOnFaceEl[0][0] + Shape1*LocalCoordOnFaceEl[1][0] 
               + Shape2*LocalCoordOnFaceEl[2][0] + Shape3*LocalCoordOnFaceEl[3][0];
  mOnFaceEl[1] = Shape0*LocalCoordOnFaceEl[0][1] + Shape1*LocalCoordOnFaceEl[1][1] 
               + Shape2*LocalCoordOnFaceEl[2][1] + Shape3*LocalCoordOnFaceEl[3][1];
}

// Return the jacobian on the real face el. at the point associated to 
// the given point m in the ref. QuadFacet 
double
QuadFacet::GetJacobianOnFaceEl(const double* m, CoordSet &cs)
{
  // J = J(mapping ref. face el. -> real face el.) * J(mapping QuadFacet -> ref. face el.)
  double JQuadFacetMapping = MappingJacobian(m);
   
  double mOnFaceEl[2];
  LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
  double JFaceElMapping = FaceEl->GetJacobian(mOnFaceEl, cs);

  return JFaceElMapping*JQuadFacetMapping;
}

// Return the isoparametric mapping normal and jacobian on the real face el. at the point 
// associated to the given point m in the ref. QuadFacet
double
QuadFacet::GetIsoParamMappingNormalAndJacobianOnFaceEl(double* Normal, const double* m, CoordSet &cs)
{
  // J = J(mapping ref. face el. -> real face el.) * J(mapping QuadFacet -> ref. face el.)
  double JQuadFacetMapping = MappingJacobian(m);

  double mOnFaceEl[2];
  LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
  double JFaceElMapping = FaceEl->GetIsoParamMappingNormalAndJacobian(Normal, mOnFaceEl, cs);

  return JFaceElMapping*JQuadFacetMapping;
}

// ------------------------------------------------------------------------------------------------------
//                         INTEGRATION OF SHAPE FUNCTIONS PRODUCT METHODS 
// ------------------------------------------------------------------------------------------------------
FullM
QuadFacet::IntegrateShapeFctProduct(MortarElement* MortarEl, QuadFacet& FriendFacet, CoordSet &cs, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// NOTE: 
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT QUADRANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT QUADRANGULAR FACET)
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet 
   FaceElement* FaceElem = FriendFacet.GetPtrFaceEl();

   int nMortarShapeFct     = MortarEl->nNodes();
   int nShapeFctOnFaceElem = FaceElem->nNodes();

   //cerr << "In QuadFacet::IntegrateShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct     = " << nMortarShapeFct << endl;
   //cerr << " -> nShapeFctOnFaceElem = " << nShapeFctOnFaceElem << endl;
   //cerr << " -> ngp at input        = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,nShapeFctOnFaceElem);
   MatShapeFctProd.zero();

   double* MortarShape      = (double*) dbg_alloca(nMortarShapeFct*sizeof(double)); 
   double* ShapeOnFaceElem  = (double*) dbg_alloca(nShapeFctOnFaceElem*sizeof(double)); 
   double m[2], mOnMortarEl[2], mOnFaceEl[2];
   double r, s, t, weight;

   int igp, i, j;
   // FOR TEST
   //ngp = 1; r = 1/3.; s = 1/3.; t = 1-r-s; weight = 0.5;
   
   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnQuadrangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Jacobian on the face element supporting the CURRENT triangular facet
     double dA = GetJacobianOnFaceEl(m, cs); 
     //cerr << " # dA = " << dA << endl;
 
     // Get local coord. on each face element 
     // -> for the mortar elem. (see the NOTE section) 
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;
 
     // -> for the element supporting the given triangular facet 
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl); 
     //cerr << " # mOnFaceEl:   x = " << mOnFaceEl[0] << ", y = " << mOnFaceEl[1] << endl;
     
     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //for(i=0;i<nMortarShapeFct;i++)
     //  cerr << " MortarShape[" << i << "]     = " << MortarShape[i] << endl;
 
     // -> for the element supporting the given triangular facet
     FaceElem->GetShapeFctVal(ShapeOnFaceElem, mOnFaceEl);
     //for(j=0;j<nShapeFctOnFaceElem;j++)
     //  cerr << " ShapeOnFaceElem[" << j << "] = " << ShapeOnFaceElem[j] << endl;
      
     // Shape fcts product & integration
     for(i=0;i<nMortarShapeFct;i++)
       for(j=0;j<nShapeFctOnFaceElem;j++)
         MatShapeFctProd[i][j] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j];
   }  
 
  //MatShapeFctProd.print("M[QuadFacet]="); 
  return MatShapeFctProd; 
}

FullM
QuadFacet::IntegrateNormalShapeFctProduct(MortarElement* MortarEl, QuadFacet& FriendFacet, CoordSet &cs, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// times the normal (see Notes (2)).
// -> Mij = Intg[current QuadFacet][Mortar(i).FriendFacet.Shape(j).normal(current QuadFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED TO 
//                                       THE GIVEN (FRIEND) TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT QUADRANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT QUADRANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT QUADRANGULAR FACET
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet
   FaceElement* FaceElem = FriendFacet.GetPtrFaceEl();

   int nMortarShapeFct     = MortarEl->nNodes();
   int nShapeFctOnFaceElem = FaceElem->nNodes();
   //cerr << "In QuadFacet::IntegrateShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct     = " << nMortarShapeFct << endl;
   //cerr << " -> nShapeFctOnFaceElem = " << nShapeFctOnFaceElem << endl;
   //cerr << " -> ngp at input        = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,3*nShapeFctOnFaceElem);
   MatShapeFctProd.zero();

   double* MortarShape      = (double*) dbg_alloca(nMortarShapeFct*sizeof(double));
   double* ShapeOnFaceElem  = (double*) dbg_alloca(nShapeFctOnFaceElem*sizeof(double));
   double m[2], mOnMortarEl[2], mOnFaceEl[2];
   double r, s, t, weight;
   double Normal[3];
   int igp, i, j;
   // FOR TEST
   //ngp = 1; r = 1/3.; s = 1/3.; t = 1-r-s; weight = 0.5;
  
   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnQuadrangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Jacobian on the face element supporting the CURRENT triangular facet
     double dA = GetIsoParamMappingNormalAndJacobianOnFaceEl(Normal, m, cs);
     //cerr << " # dA = " << dA << endl;
     //cerr << " # normal = " << Normal[0] <<" "<< Normal[1] <<" "<< Normal[2] << endl;

     // Get local coord. on each face element
     // -> for the mortar elem. (see the NOTE section)
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;

     // -> for the element supporting the given triangular facet
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
     //cerr << " # mOnFaceEl:   x = " << mOnFaceEl[0] << ", y = " << mOnFaceEl[1] << endl;

     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //for(i=0;i<nMortarShapeFct;i++)
     //  cerr << " MortarShape[" << i << "]     = " << MortarShape[i] << endl;

     // -> for the element supporting the given triangular facet
     FaceElem->GetShapeFctVal(ShapeOnFaceElem, mOnFaceEl);
     //for(j=0;j<nShapeFctOnFaceElem;j++)
     //  cerr << " ShapeOnFaceElem[" << j << "] = " << ShapeOnFaceElem[j] << endl;

     // Shape fcts product & integration
     for(i=0;i<nMortarShapeFct;i++)
       for(j=0;j<nShapeFctOnFaceElem;j++){
         MatShapeFctProd[i][3*j  ] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[0];
         MatShapeFctProd[i][3*j+1] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[1];
         MatShapeFctProd[i][3*j+2] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[2];
       } 
  }

  //MatShapeFctProd.print("M[QuadFacet]=");
  return MatShapeFctProd;
}

FullM
QuadFacet::IntegrateGradNormalShapeFctProduct(MortarElement* MortarEl, QuadFacet& FriendFacet, CoordSet& cs,
                                             CoordSet& cs1, QuadFacet& FriendFacet2, CoordSet& cs2, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// times the normal (see Notes (2)).
// -> Mij = Intg[current QuadFacet][Mortar(i).FriendFacet.Shape(j).normal(current QuadFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED TO 
//                                       THE GIVEN (FRIEND) TRIFACET) 
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT QUADRANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT QUADRANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT QUADRANGULAR FACET
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet
   FaceElement* FaceElem = FriendFacet.GetPtrFaceEl();

   int nMortarShapeFct     = MortarEl->nNodes();
   int nShapeFctOnFaceElem = FaceElem->nNodes();
   //cerr << "In QuadFacet::IntegrateShapeFctProduct" << endl;
   //cerr << " -> nMortarShapeFct     = " << nMortarShapeFct << endl;
   //cerr << " -> nShapeFctOnFaceElem = " << nShapeFctOnFaceElem << endl;
   //cerr << " -> ngp at input        = " << ngp << endl;

   FullM MatShapeFctProd(nMortarShapeFct,3*nShapeFctOnFaceElem);
   MatShapeFctProd.zero();

   double* MortarShape      = (double*) dbg_alloca(nMortarShapeFct*sizeof(double));
   double* ShapeOnFaceElem  = (double*) dbg_alloca(nShapeFctOnFaceElem*sizeof(double));
   double m[2], mOnMortarEl[2], mOnFaceEl[2];
   double r, s, t, weight;
   double Normal[3];
   int igp, i, j;
   // FOR TEST
   //ngp = 1; r = 1/3.; s = 1/3.; t = 1-r-s; weight = 0.5;
  
   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnQuadrangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Jacobian on the face element supporting the CURRENT triangular facet
     double dA = GetIsoParamMappingNormalAndJacobianOnFaceEl(Normal, m, cs);
     //cerr << " # dA = " << dA << endl;
     //cerr << " # normal = " << setprecision(10) << Normal[0] <<" "<< Normal[1] <<" "<< Normal[2] << endl;

     // Get local coord. on each face element
     // -> for the mortar elem. (see the NOTE section)
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;

     // -> for the element supporting the given triangular facet
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
     //cerr << " # mOnFaceEl:   x = " << mOnFaceEl[0] << ", y = " << mOnFaceEl[1] << endl;

     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);
     //for(i=0;i<nMortarShapeFct;i++)
     //  cerr << " MortarShape[" << i << "]     = " << MortarShape[i] << endl;

     // -> for the element supporting the given triangular facet
     FaceElem->GetShapeFctVal(ShapeOnFaceElem, mOnFaceEl);
     //for(j=0;j<nShapeFctOnFaceElem;j++)
     //  cerr << " ShapeOnFaceElem[" << j << "] = " << ShapeOnFaceElem[j] << endl;

     // Shape fcts product & integration
     for(i=0;i<nMortarShapeFct;i++)
       for(j=0;j<nShapeFctOnFaceElem;j++){
         MatShapeFctProd[i][3*j  ] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[0];
         MatShapeFctProd[i][3*j+1] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[1];
         MatShapeFctProd[i][3*j+2] += weight*dA*MortarShape[i]*ShapeOnFaceElem[j]*Normal[2];
       } 
  }

  //MatShapeFctProd.print("M[QuadFacet]=");
  return MatShapeFctProd;
}

// Compute normal "geometrical" gap 
Vector
QuadFacet::IntegrateNormalGeoGagsProduct(MortarElement* MortarEl, QuadFacet& FriendFacet, CoordSet& cs, CoordSet& cs1, int ngp)
// ******************************************************************************************************
// Integrate on the CURRENT triangular facet the product of the shape functions defined by the given
// MortarElement and the shape functions of the element associated to the given (Friend) triangular facet
// times the normal (see Notes (2)).
// -> Mij = Intg[current QuadFacet][Mortar(i).FriendFacet.Shape(j).normal(current QuadFacet->FaceElem)]
// -> THE SIZE OF THE OUTPUT MATRIX IS:
//      (NUMBER OF MORTAR SHAPE FCTS) x (THE NUMBER OF DOFs OF THE FACE ELEMENT ASSOCIATED TO
//                                       THE GIVEN (FRIEND) TRIFACET)
// NOTES:
//     (1) WE ASSUME THAT THE MortarElement LIVES ON THE GEOMETRY ASSOCIATED WITH THE CURRENT QUADRANGULAR
//         FACET (I.E. THE ELEMENT SUPPORTING THE CURRENT QUADRANGULAR FACET)
//     (2) THE NORMAL IS THE NORMAL OF THE ELEMENT SUPPORTING THE CURRENT QUADRANGULAR FACET
// ******************************************************************************************************
{
   // Get ptr to the face element supporting the given triangular facet
   FaceElement* FaceElem = FriendFacet.GetPtrFaceEl();
   int nMortarShapeFct   = MortarEl->nNodes();

   Vector NormalGeoGaps(nMortarShapeFct,0.0);

   double* MortarShape = (double*) dbg_alloca(nMortarShapeFct*sizeof(double));
   double m[2], mOnMortarEl[2], mOnFaceEl[2];
   double r, s, t, weight;
   double Normal[3], MOnFaceEl[3];
   int igp, i;
   // FOR TEST
   //ngp = 1; r = 1/3.; s = 1/3.; t = 1-r-s; weight = 0.5;

   for(igp=1; igp<=ngp; igp++){
     //cerr << " # Gauss point " << igp << endl;

     getGaussPtOnQuadrangle(ngp,igp,r,s,t,weight);
     //cerr << " # r = " << r << ", s = " << s << ", w = " << weight << endl;
     m[0] = r; m[1] = s;

     // Jacobian on the face element supporting the CURRENT triangular facet
     double dA = GetIsoParamMappingNormalAndJacobianOnFaceEl(Normal, m, cs);
     //cerr << " # dA = " << dA << endl;
     //cerr << " # normal = " << Normal[0] <<" "<< Normal[1] <<" "<< Normal[2] << endl;

     // Get local coord. on each face element
     // -> for the mortar elem. (see the NOTE section)
     LocalToLocalCoordOnFaceEl(m, mOnMortarEl);
     //cerr << " # mOnMortarEl: x = " << mOnMortarEl[0] << ", y = " << mOnMortarEl[1] << endl;

     // -> for the element supporting the given triangular facet
     FriendFacet.LocalToLocalCoordOnFaceEl(m, mOnFaceEl);
     //cerr << " # mOnFaceEl:   x = " << mOnFaceEl[0] << ", y = " << mOnFaceEl[1] << endl;

     // Compute shape fcts
     // -> mortar elem.
     MortarEl->GetShapeFctVal(MortarShape, mOnMortarEl);

     // Get global coord. on the element supporting the given triangular facet
     FaceElem->LocalToGlobalCoord(MOnFaceEl, mOnFaceEl, cs1);

     // Shape fcts product & integration
     double sdot = weight*dA*(MOnFaceEl[0]*Normal[0]+MOnFaceEl[1]*Normal[1]+MOnFaceEl[2]*Normal[2]);
     for(i=0;i<nMortarShapeFct;i++)
       NormalGeoGaps[i] += sdot*MortarShape[i];
  }

  //NormalGeoGaps.print("NormalGeoGaps[QuadFacet]=");
  return NormalGeoGaps;
}
