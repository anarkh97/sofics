#ifndef _FFIPOLYGON_C_ 
#define _FFIPOLYGON_C_

// Std C/C++ headers
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FFIPolygon.d/FFIPolygon.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarDefines.h>

#include <Acme.d/search/Contact_Defines.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class MasterFaceType, class SlaveFaceType, class MortarType>
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::FFIPolygonTemplate()
: Facets()
{
  Initialize();
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::FFIPolygonTemplate(MasterFaceType* MasterFaceEl, SlaveFaceType* SlaveFaceEl)
: Facets()
{
  Initialize();
  MasterFace = MasterFaceEl;
  SlaveFace = SlaveFaceEl;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::FFIPolygonTemplate(MasterFaceType* MasterFaceEl, SlaveFaceType* SlaveFaceEl, int nVert, double* ACME_FFI_Data,
                     int ACME_FFI_Derivatives_Order)
: Facets()
{
  Initialize();
  SetFFIPolygon(MasterFaceEl, SlaveFaceEl, nVert, ACME_FFI_Data, ACME_FFI_Derivatives_Order);
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTORS 
// -----------------------------------------------------------------------------------------------------
template<class MasterFaceType, class SlaveFaceType, class MortarType>
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::~FFIPolygonTemplate()
{
#ifdef HB_ACME_FFI_DEBUG
  if(VertexLlCoordOnSFaceEl) { delete VertexLlCoordOnSFaceEl; }
  if(VertexLlCoordOnMFaceEl) { delete VertexLlCoordOnMFaceEl; }
#endif
  if(dM) delete [] dM;
  if(dN) delete [] dN;
}

// -----------------------------------------------------------------------------------------------------
//                                    INITIALIZATION & CLEAR/CLEAN METHODS
// -----------------------------------------------------------------------------------------------------
template<class MasterFaceType, class SlaveFaceType, class MortarType>
void 
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::Initialize()
{
  Area      = 0.0;
  nVertices = 0;
  MasterFace= 0;
  SlaveFace = 0;

#ifdef HB_ACME_FFI_DEBUG
  VertexLlCoordOnSFaceEl = 0;
  VertexLlCoordOnMFaceEl = 0;
#endif

  dM = 0;
  dN = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS 
// -----------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS 
// -----------------------------------------------------------------------------------------------------
template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::SetPtrMasterFace(MasterFaceType* PtrMasterFace) { MasterFace = PtrMasterFace; }

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::SetPtrSlaveFace(SlaveFaceType* PtrSlaveFace) { SlaveFace = PtrSlaveFace; }


template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::SetFFIPolygon(MasterFaceType* MasterFaceEl, SlaveFaceType* SlaveFaceEl, int nVert, double* ACME_FFI_Data,
                int ACME_FFI_Derivatives_Order)
{
  nVertices  = nVert;
  MasterFace = MasterFaceEl;
  SlaveFace  = SlaveFaceEl;

  // Create & fill polygon triangularization
  int LocalCoordOffset = 2*nVertices+1;
  ACME_FFI_LocalCoordData = &ACME_FFI_Data[LocalCoordOffset];
  switch(ACME_FFI_Derivatives_Order) {
    case 0 :
      CreateTriangularization(ACME_FFI_LocalCoordData);
      break;
    case 1 :
      CreateTriangularizationExp(ACME_FFI_LocalCoordData);
      break;
    case 2 :
      CreateTriangularizationExp2(ACME_FFI_LocalCoordData);
      break;
  }
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS 
// -----------------------------------------------------------------------------------------------------
template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::Print()
{
  fprintf(stderr," * -------------------------------- \n");
  fprintf(stderr," * FFIPolygon data:\n");
  fprintf(stderr," + number of vertices: %d\n", GetnVertices()); 
  fprintf(stderr," + number of TriFacet: %d\n", GetnFacets()); 
  fprintf(stderr," + ptr to master elem: %p\n", MasterFace);
  if(MasterFace) MasterFace->print();
  fprintf(stderr," + ptr to slave elem : %p\n", SlaveFace);
  if(SlaveFace) SlaveFace->print();
  fprintf(stderr," + TriFacets data on master face:\n");  
  for(int iTri=0; iTri<GetnFacets(); iTri++){
    fprintf(stderr," -> TriFacet %d (on master) \n",iTri+1);
    MasterFacet(Facets,iTri).Print();    
  } 
  fprintf(stderr," + TriFacets data on slave face:\n");  
  for(int iTri=0; iTri<GetnFacets(); iTri++){
    fprintf(stderr," -> TriFacet %d (on slave) \n",iTri+1);
    SlaveFacet(Facets,iTri).Print();    
  } 
  fprintf(stderr," * -------------------------------- \n");
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::PrintM()
{
  GetPtrSlaveFace()->print(); 
  GetPtrSlaveFace()->print(); 
  M.print("M[FFI] = ");
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::PrintN()
{
  GetPtrSlaveFace()->print(); 
  GetPtrMasterFace()->print(); 
  N.print("N[FFI] = ");
}

#ifdef HB_ACME_FFI_DEBUG
template<class MasterFaceType, class SlaveFaceType, class MortarType>
void 
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::PrintSlaveVertices(FILE* file, CoordSet& cs, int& firstVertId)
{
  double XYZ[3];
  for(int iVert=0; iVert<nVertices; iVert++){
    double m[2] = {VertexLlCoordOnSFaceEl[iVert][0],VertexLlCoordOnSFaceEl[iVert][1]};
    SlaveFace->LocalToGlobalCoord(XYZ,m,cs);
    fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
    firstVertId++;
  }
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::PrintMasterVertices(FILE* file, CoordSet& cs, int& firstVertId)
{
  double XYZ[3];
  for(int iVert=0; iVert<nVertices; iVert++){
    double m[2] = {VertexLlCoordOnMFaceEl[iVert][0],VertexLlCoordOnMFaceEl[iVert][1]};
    MasterFace->LocalToGlobalCoord(XYZ,m,cs);
    fprintf(file," %6d  %6e  %6e  %6e\n",firstVertId,XYZ[0],XYZ[1],XYZ[2]);
    firstVertId++;
  }
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::PrintFFIPolygonTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode)
{
   int firstVertId = VertOffset;
   switch(nVertices) {
     case 4 : { 
       fprintf(file," %6d  %3d  %6d  %6d  %6d  %6d\n",EdgeOffset,188,VertOffset,VertOffset+1,VertOffset+2,VertOffset+3);
       VertOffset+=4; EdgeOffset++;
     } break;
     case 3 : {
       fprintf(file," %6d  %3d  %6d  %6d  %6d\n",EdgeOffset,108,VertOffset,VertOffset+1,VertOffset+2);
       VertOffset+=3; EdgeOffset++;
     } break;
     default : {
       for(int iVert=0; iVert<nVertices-1; iVert++){
         fprintf(file," %6d  %3d  %6d  %6d\n",EdgeOffset,elCode,VertOffset,VertOffset+1);
         VertOffset++; EdgeOffset++;
       }
       fprintf(file," %6d  %3d  %6d  %6d\n",EdgeOffset,elCode,VertOffset,firstVertId);
       VertOffset++; EdgeOffset++;
     } break;
   }
}
#endif

// -----------------------------------------------------------------------------------------------------
//                                       TRIANGULARIZATION METHODS 
// -----------------------------------------------------------------------------------------------------
template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::CreateTriangularization(double* ACME_FFI_LocalCoordData, bool use_acme_centroid)
{
  // Allocate space for Facets
  // By default assume as many triangular facet as polygon vertices
  // -> can be overwritten in for optimization XXX
  Facets.assign(nVertices, Facet_pair_t(TriFacetTemplate<double,MasterFaceType>(), TriFacetTemplate<double,SlaveFaceType>()));

  // Create master & slave polygon centroids
  double MasterCentroid[2] = {0.0, 0.0};
  double SlaveCentroid[2]  = {0.0, 0.0};

#ifdef HB_ACME_FFI_DEBUG
  VertexLlCoordOnSFaceEl = new double[nVertices][2];
  VertexLlCoordOnMFaceEl = new double[nVertices][2];
#endif
  for(int iVert=0; iVert<nVertices; ++iVert) {
     int offset = 4*iVert;
#ifdef HB_ACME_FFI_DEBUG
     VertexLlCoordOnMFaceEl[iVert][0] = ACME_FFI_LocalCoordData[offset  ];
     VertexLlCoordOnMFaceEl[iVert][1] = ACME_FFI_LocalCoordData[offset+1];
     VertexLlCoordOnSFaceEl[iVert][0] = ACME_FFI_LocalCoordData[offset+2];
     VertexLlCoordOnSFaceEl[iVert][1] = ACME_FFI_LocalCoordData[offset+3];
#endif 
     MasterCentroid[0] += ACME_FFI_LocalCoordData[offset  ]; 
     MasterCentroid[1] += ACME_FFI_LocalCoordData[offset+1]; 
     SlaveCentroid[0]  += ACME_FFI_LocalCoordData[offset+2]; 
     SlaveCentroid[1]  += ACME_FFI_LocalCoordData[offset+3]; 
  }
  SlaveCentroid[0]  /= nVertices; SlaveCentroid[1]  /= nVertices;
  MasterCentroid[0] /= nVertices; MasterCentroid[1] /= nVertices;

  // Create master & slave triangularization
  for(size_t i=0, ii=Facets.size(), offsetm1=0; i<ii; i++, offsetm1+=4) {
     size_t offsetm2 = (i==ii-1) ? 0 : offsetm1+4;
     //HB: note that ACME Master <=> our Slave, ... 
     if(use_acme_centroid) {
       size_t offsetmc = 4*Facets.size();
       MasterFacet(Facets, i).SetTriFacet(MasterFace, &ACME_FFI_LocalCoordData[offsetmc],   &ACME_FFI_LocalCoordData[offsetm1],
                                                                                            &ACME_FFI_LocalCoordData[offsetm2]);
       SlaveFacet(Facets , i).SetTriFacet(SlaveFace,  &ACME_FFI_LocalCoordData[offsetmc+2], &ACME_FFI_LocalCoordData[offsetm1+2],
                                                                                            &ACME_FFI_LocalCoordData[offsetm2+2]);
     }
     else {
       MasterFacet(Facets, i).SetTriFacet(MasterFace, MasterCentroid, &ACME_FFI_LocalCoordData[offsetm1],
                                                                      &ACME_FFI_LocalCoordData[offsetm2]);
       SlaveFacet(Facets , i).SetTriFacet(SlaveFace,  SlaveCentroid,  &ACME_FFI_LocalCoordData[offsetm1+2],
                                                                      &ACME_FFI_LocalCoordData[offsetm2+2]);
    }
  } 
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::CreateTriangularizationExp(double* ACME_FFI_LocalCoordData, bool no_ffi_derivatives)
{
  int nbDer = 3*(MasterFace->nNodes() + SlaveFace->nNodes());
  // Allocate space for Facets
  Facets.assign(nVertices, Facet_pair_t(TriFacetTemplate<double,MasterFaceType>(), TriFacetTemplate<double,SlaveFaceType>()));
  // Create master & slave polygon centroids
  double MasterCentroid[2] = { 0.0, 0.0 };
  double SlaveCentroid[2]  = { 0.0, 0.0 };
  double *MasterCentroidDerivatives, *SlaveCentroidDerivatives;
  MasterCentroidDerivatives = new double[2*nbDer];
  SlaveCentroidDerivatives = new double[2*nbDer];
  for(int i=0; i<2*nbDer; ++i) {
    MasterCentroidDerivatives[i] = 0.0;
    SlaveCentroidDerivatives[i]  = 0.0;
  }

  for(int iVert=0; iVert<nVertices; ++iVert) {
    int offset = 4*iVert;
    MasterCentroid[0] += ACME_FFI_LocalCoordData[offset  ];
    MasterCentroid[1] += ACME_FFI_LocalCoordData[offset+1];
    SlaveCentroid[0]  += ACME_FFI_LocalCoordData[offset+2];
    SlaveCentroid[1]  += ACME_FFI_LocalCoordData[offset+3];
    for(int j=0; j<nbDer; ++j) {
      MasterCentroidDerivatives[j]       += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset  )+j];
      MasterCentroidDerivatives[nbDer+j] += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset+1)+j];
      SlaveCentroidDerivatives[j]        += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset+2)+j];
      SlaveCentroidDerivatives[nbDer+j]  += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset+3)+j];
    }
  }
  MasterCentroid[0] /= nVertices; MasterCentroid[1] /= nVertices;
  SlaveCentroid[0]  /= nVertices; SlaveCentroid[1]  /= nVertices;
  for(int i=0; i<2*nbDer; ++i) {
    MasterCentroidDerivatives[i] /= nVertices;
    SlaveCentroidDerivatives[i]  /= nVertices;
  }
  // Create master & slave triangularization
  for(size_t i=0, ii=Facets.size(), offsetm1=0; i<ii; i++, offsetm1+=4) {
    size_t offsetm2 = (i==ii-1) ? 0 : offsetm1+4;
    //HB: note that ACME Master <=> our Slave, ... 
    MasterFacet(Facets, i).SetTriFacet(MasterFace, MasterCentroid, &ACME_FFI_LocalCoordData[offsetm1], &ACME_FFI_LocalCoordData[offsetm2],
                                       MasterCentroidDerivatives, &ACME_FFI_LocalCoordData[4*nVertices+nbDer*offsetm1],
                                       &ACME_FFI_LocalCoordData[4*nVertices+nbDer*offsetm2], nbDer);
    SlaveFacet(Facets , i).SetTriFacet(SlaveFace,  SlaveCentroid,  &ACME_FFI_LocalCoordData[offsetm1+2], &ACME_FFI_LocalCoordData[offsetm2+2],
                                       SlaveCentroidDerivatives, &ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offsetm1+2)],
                                       &ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offsetm2+2)], nbDer);
  }
  delete [] MasterCentroidDerivatives;
  delete [] SlaveCentroidDerivatives;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::CreateTriangularizationExp2(double* ACME_FFI_LocalCoordData, bool no_ffi_derivatives)
{
  int nbDer = 3*(MasterFace->nNodes() + SlaveFace->nNodes());
  int nbDer2 = nbDer*nbDer;
  // Allocate space for Facets
  Facets.assign(nVertices, Facet_pair_t(TriFacetTemplate<double,MasterFaceType>(), TriFacetTemplate<double,SlaveFaceType>()));
  // Create master & slave polygon centroids
  double MasterCentroid[2] = { 0.0, 0.0 };
  double SlaveCentroid[2]  = { 0.0, 0.0 };
  double *MasterCentroidDerivatives, *SlaveCentroidDerivatives;
  double *MasterCentroid2ndDerivatives, *SlaveCentroid2ndDerivatives;
  MasterCentroidDerivatives = new double[2*nbDer];
  SlaveCentroidDerivatives = new double[2*nbDer];
  for(int i=0; i<2*nbDer; ++i) {
    MasterCentroidDerivatives[i] = 0.0;
    SlaveCentroidDerivatives[i]  = 0.0;
  }
  MasterCentroid2ndDerivatives = new double[2*nbDer2];
  SlaveCentroid2ndDerivatives = new double[2*nbDer2];
  for(int i=0; i<2*nbDer2; ++i) {
    MasterCentroid2ndDerivatives[i] = 0.0;
    SlaveCentroid2ndDerivatives[i]  = 0.0;
  }

  for(int iVert=0; iVert<nVertices; ++iVert) {
    int offset = 4*iVert;
    MasterCentroid[0] += ACME_FFI_LocalCoordData[offset  ];
    MasterCentroid[1] += ACME_FFI_LocalCoordData[offset+1];
    SlaveCentroid[0]  += ACME_FFI_LocalCoordData[offset+2];
    SlaveCentroid[1]  += ACME_FFI_LocalCoordData[offset+3];
    for(int j=0; j<nbDer; ++j) {
      MasterCentroidDerivatives[j]       += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset  )+j];
      MasterCentroidDerivatives[nbDer+j] += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset+1)+j];
      SlaveCentroidDerivatives[j]        += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset+2)+j];
      SlaveCentroidDerivatives[nbDer+j]  += ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offset+3)+j];
    }
    for(int j=0; j<nbDer2; ++j) {
      MasterCentroid2ndDerivatives[j]        += ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*(offset  )+j];
      MasterCentroid2ndDerivatives[nbDer2+j] += ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*(offset+1)+j];
      SlaveCentroid2ndDerivatives[j]         += ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*(offset+2)+j];
      SlaveCentroid2ndDerivatives[nbDer2+j]  += ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*(offset+3)+j];
    }
  }
  MasterCentroid[0] /= nVertices; MasterCentroid[1] /= nVertices;
  SlaveCentroid[0]  /= nVertices; SlaveCentroid[1]  /= nVertices;
  for(int i=0; i<2*nbDer; ++i) {
    MasterCentroidDerivatives[i] /= nVertices;
    SlaveCentroidDerivatives[i]  /= nVertices;
  }
  for(int i=0; i<2*nbDer2; ++i) {
    MasterCentroid2ndDerivatives[i] /= nVertices;
    SlaveCentroid2ndDerivatives[i]  /= nVertices;
  }
  // Create master & slave triangularization
  for(size_t i=0, ii=Facets.size(), offsetm1=0; i<ii; i++, offsetm1+=4) {
    size_t offsetm2 = (i==ii-1) ? 0 : offsetm1+4;
    //HB: note that ACME Master <=> our Slave, ... 
    MasterFacet(Facets, i).SetTriFacet(MasterFace, MasterCentroid, &ACME_FFI_LocalCoordData[offsetm1], &ACME_FFI_LocalCoordData[offsetm2],
                                       MasterCentroidDerivatives, &ACME_FFI_LocalCoordData[4*nVertices+nbDer*offsetm1],
                                       &ACME_FFI_LocalCoordData[4*nVertices+nbDer*offsetm2], MasterCentroid2ndDerivatives,
                                       &ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*offsetm1],
                                       &ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*offsetm2], nbDer);
    SlaveFacet(Facets , i).SetTriFacet(SlaveFace,  SlaveCentroid,  &ACME_FFI_LocalCoordData[offsetm1+2], &ACME_FFI_LocalCoordData[offsetm2+2],
                                       SlaveCentroidDerivatives, &ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offsetm1+2)],
                                       &ACME_FFI_LocalCoordData[4*nVertices+nbDer*(offsetm2+2)], SlaveCentroid2ndDerivatives,
                                       &ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*(offsetm1+2)],
                                       &ACME_FFI_LocalCoordData[4*nVertices*(1+nbDer)+nbDer2*(offsetm2+2)], nbDer);
  }
  delete [] MasterCentroidDerivatives;
  delete [] SlaveCentroidDerivatives;
  delete [] MasterCentroid2ndDerivatives;
  delete [] SlaveCentroid2ndDerivatives;
}

// -----------------------------------------------------------------------------------------------------
//                                            MISCELLEANEOUS METHODS 
// -----------------------------------------------------------------------------------------------------
/*
template<class MasterFaceType, class SlaveFaceType, class MortarType>
double
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>::ComputeArea()
{
   double Area = 0.0;
   for(size_t i=0, ii=GetnFacets(); i<ii; i++) {
      Area += TriFacet[i]->ComputeArea();
   }
   return Area;
}*/

// -----------------------------------------------------------------------------------------------------
//                         INTEGRATION OF SHAPE FUNCTIONS PRODUCT METHODS 
// -----------------------------------------------------------------------------------------------------
template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnMaster_MasterShapeFctProduct(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the MASTER SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the MASTER face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
  int nMortarShapeFct = MortarEl->nNodes();
  int nMasterShapeFct = MasterFace->nNodes();

  FullM MatShapeFctProd(nMortarShapeFct,nMasterShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += MasterFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, MasterFacet(Facets, i),
                                                                       cs, ngp);
  }

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnMaster_SlaveShapeFctProduct(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the MASTER SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the SLAVE face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & SLAVE face element node ordering
// *****************************************************************************************************
{
  int nMortarShapeFct = MortarEl->nNodes();
  int nSlaveShapeFct  = SlaveFace->nNodes();

  FullM MatShapeFctProd(nMortarShapeFct,nSlaveShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += MasterFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, SlaveFacet(Facets, i),
                                                                       cs, ngp);
  } 

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnSlave_MasterShapeFctProduct(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the MASTER face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
  int nMortarShapeFct = MortarEl->nNodes();
  int nMasterShapeFct = MasterFace->nNodes();

  FullM MatShapeFctProd(nMortarShapeFct,nMasterShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += SlaveFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, MasterFacet(Facets, i),
                                                                      cs, ngp);
  }

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnSlave_SlaveShapeFctProduct(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the shape functions defined by the given
// MortarElement and the shape functions of the SLAVE face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face 
//         element of THIS FFIPolygon object  
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering 
//         of the MORTAR & SLAVE face element node ordering
// *****************************************************************************************************
{
  int nMortarShapeFct = MortarEl->nNodes();
  int nSlaveShapeFct  = SlaveFace->nNodes();
   
  FullM MatShapeFctProd(nMortarShapeFct,nSlaveShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += SlaveFacet(Facets, i).IntegrateShapeFctProduct(MortarEl, SlaveFacet(Facets, i),
                                                                      cs, ngp);
  }

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void 
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::ComputeM(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Slave shape fcts product mass-like matrix M
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygonTemplate::IntegrateOnSlave_SlaveShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering 
//         of the MORTAR & SLAVE face element node ordering
//         -> see FFIPolygonTemplate::IntegrateOnSlave_SlaveShapeFctProduct for more details
// *****************************************************************************************************
{
  M = IntegrateOnSlave_SlaveShapeFctProduct(MortarEl, cs, ngp);
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::ComputeN(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Master shape fcts product mass-like matrix N
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygonTemplate::IntegrateOnSlave_MasterShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering 
//         of the MORTAR & MASTER face element node ordering
//         -> see FFIPolygonTemplate::IntegrateOnSlave_MasterShapeFctProduct for more details
// *****************************************************************************************************
{
  N = IntegrateOnSlave_MasterShapeFctProduct(MortarEl, cs, ngp);
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnSlave_MasterNormalShapeFctProduct(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the normal and the shape functions defined by the given
// MortarElement and the shape functions of the MASTER face element 
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
   int nMortarShapeFct = MortarEl->nNodes();
   int nMasterShapeFct = MasterFace->nNodes();

  FullM MatShapeFctProd(nMortarShapeFct,3*nMasterShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += SlaveFacet(Facets, i).IntegrateNormalShapeFctProduct(MortarEl, MasterFacet(Facets, i),
                                                                            cs, ngp);
  }

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the product of the normal and the shape functions defined by the given
// MortarElement and the shape functions of the SLAVE face element
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
  int nMortarShapeFct = MortarEl->nNodes();
  int nSlaveShapeFct  = SlaveFace->nNodes();                                                                                                                                          
  FullM MatShapeFctProd(nMortarShapeFct,3*nSlaveShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += SlaveFacet(Facets, i).IntegrateNormalShapeFctProduct(MortarEl, SlaveFacet(Facets, i),
                                                                            cs, ngp);
  }

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::ComputeNormalM(MortarType* MortarEl, CoordSet &cs, int ngp)
{
  M = IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarEl, cs, ngp);
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::ComputeNormalN(MortarType* MortarEl, CoordSet &cs, int ngp)
// *****************************************************************************************************
// Compute the FFI contribution to the Mortar-Master shape fcts product mass-like matrix N
// NOTE:
//     (1) by DEFAULT integrate on the SLAVE SIDE
//         -> see FFIPolygonTemplate::IntegrateOnSlave_MasterShapeFctProduct for more details
//     (2) the rows & colums ordering of the matrix M is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
//         -> see FFIPolygonTemplate::IntegrateOnSlave_MasterShapeFctProduct for more details
// *****************************************************************************************************
{
  N = IntegrateOnSlave_MasterNormalShapeFctProduct(MortarEl, cs, ngp);
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnSlave_GradNormalGeoGap(MortarType* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp, double offset)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the contribution to the gradient of the gap function term due to the
// products of the normal and the shape functions defined by the given MortarElement and the shape
// functions of both the MASTER and the SLAVE face elements i.e. C = [-N M]
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
  int nMortarShapeFct = MortarEl->nNodes();
  int nMasterShapeFct = MasterFace->nNodes();
  int nSlaveShapeFct  = SlaveFace->nNodes();

  FullM MatShapeFctProd(nMortarShapeFct,3*nMasterShapeFct+3*nSlaveShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += SlaveFacet(Facets, i).IntegrateGradNormalGeoGap(MortarEl, MasterFacet(Facets, i),
                                                                       SlaveCoords, MasterCoords, ngp, offset);
  }

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
FullM
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::IntegrateOnSlave_HessNormalGeoGap(MortarType* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, double *mu, int ngp, double offset)
// *****************************************************************************************************
// Integrate on the SLAVE SIDE the contribution to the hessian of the gap function contracted with mu due to the
// products of the normal and the shape functions defined by the given MortarElement and the shape
// functions of both the MASTER and the SLAVE face elements i.e. H = sum_i mu[i]*(Jacobian of [-N[i] M[i] ])
// NOTE:
//     (1) you HAVE to ensure that the MortarElement is ASSOCIATED with the SLAVE (or MASTER) face
//         element of THIS FFIPolygon object
//     (2) the rows & colums ordering of the shape fcts product matrix is ASSOCIATED with the ordering
//         of the MORTAR & MASTER face element node ordering
// *****************************************************************************************************
{
  int nMortarShapeFct = MortarEl->nNodes();
  int nMasterShapeFct = MasterFace->nNodes();
  int nSlaveShapeFct  = SlaveFace->nNodes();

  FullM MatShapeFctProd(3*nMasterShapeFct+3*nSlaveShapeFct, 3*nMasterShapeFct+3*nSlaveShapeFct);
  MatShapeFctProd.zero();

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    MatShapeFctProd += SlaveFacet(Facets, i).IntegrateHessNormalGeoGap(MortarEl, MasterFacet(Facets, i),
                                                                       SlaveCoords, MasterCoords, mu, ngp, offset);
  }

  return MatShapeFctProd;
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::ComputeNormalGeoGap(MortarType* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp, double offset)
{
  int nMortarShapeFct = MortarEl->nNodes();

  NormalGeoGaps.reset(nMortarShapeFct,0.0);

  // Loop over triangularization
  for(size_t i=0, ii=GetnFacets(); i<ii; ++i) {
    NormalGeoGaps += SlaveFacet(Facets, i).IntegrateNormalGeoGap(MortarEl, MasterFacet(Facets, i),
                                                                 SlaveCoords, MasterCoords, ngp, offset);
  }
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::ComputeGradNormalGeoGap(MortarType* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, int ngp, double offset)
{
  FullM C = IntegrateOnSlave_GradNormalGeoGap(MortarEl, SlaveCoords, MasterCoords, ngp, offset);

  int nMortarShapeFct = MortarEl->nNodes();
  M.setNewSize(nMortarShapeFct,3*SlaveFace->nNodes());
  N.setNewSize(nMortarShapeFct,3*MasterFace->nNodes());
  for(int i=0; i<nMortarShapeFct; ++i) {
    for(int j=0; j<3*SlaveFace->nNodes(); ++j) M[i][j] = C(i,3*MasterFace->nNodes()+j);
    for(int j=0; j<3*MasterFace->nNodes(); ++j) N[i][j] = -C(i,j);
  }
}

template<class MasterFaceType, class SlaveFaceType, class MortarType>
void
FFIPolygonTemplate<MasterFaceType,SlaveFaceType,MortarType>
::ComputeHessNormalGeoGap(MortarType* MortarEl, CoordSet &SlaveCoords, CoordSet &MasterCoords, double *mu, int ngp, double offset)
{
  H = IntegrateOnSlave_HessNormalGeoGap(MortarEl, SlaveCoords, MasterCoords, mu, ngp, offset);
}
#endif
