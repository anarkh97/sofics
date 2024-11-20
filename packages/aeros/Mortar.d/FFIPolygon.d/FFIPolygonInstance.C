#include <Mortar.d/FFIPolygon.d/FFIPolygon.C>

template
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::FFIPolygonTemplate();

template
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::FFIPolygonTemplate(FaceElement*, FaceElement*);

template
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::FFIPolygonTemplate(FaceElement*, FaceElement*, int, double*, int);

template
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::~FFIPolygonTemplate();

template
void 
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::Initialize();

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::SetPtrMasterFace(FaceElement*);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::SetPtrSlaveFace(FaceElement*);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::SetFFIPolygon(FaceElement*, FaceElement*, int, double*, int);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::Print();

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::PrintM();

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::PrintN();

#ifdef HB_ACME_FFI_DEBUG
template
void 
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::PrintSlaveVertices(FILE*, CoordSet&, int&);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::PrintMasterVertices(FILE*, CoordSet&, int&);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::PrintFFIPolygonTopo(FILE*, int&, int&, int);
#endif

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::CreateTriangularization(double*, bool);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::CreateTriangularizationExp(double*, bool);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::CreateTriangularizationExp2(double*, bool);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnMaster_MasterShapeFctProduct(MortarElement*, CoordSet&, int);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnMaster_SlaveShapeFctProduct(MortarElement*, CoordSet&, int);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnSlave_MasterShapeFctProduct(MortarElement*, CoordSet&, int);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnSlave_SlaveShapeFctProduct(MortarElement*, CoordSet&, int);

template
void 
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::ComputeM(MortarElement*, CoordSet&, int);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::ComputeN(MortarElement*, CoordSet&, int);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnSlave_MasterNormalShapeFctProduct(MortarElement*, CoordSet&, int);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarElement*, CoordSet&, int);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::ComputeNormalM(MortarElement*, CoordSet&, int);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::ComputeNormalN(MortarElement*, CoordSet&, int);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::ComputeNormalGeoGap(MortarElement*, CoordSet&, CoordSet&, int, double);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnSlave_GradNormalGeoGap(MortarElement*, CoordSet&, CoordSet&, int, double);

template
FullM
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::IntegrateOnSlave_HessNormalGeoGap(MortarElement*, CoordSet&, CoordSet&, double*, int, double);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::ComputeGradNormalGeoGap(MortarElement*, CoordSet&, CoordSet&, int, double);

template
void
FFIPolygonTemplate<FaceElement,FaceElement,MortarElement>
::ComputeHessNormalGeoGap(MortarElement*, CoordSet&, CoordSet&, double*, int, double);

