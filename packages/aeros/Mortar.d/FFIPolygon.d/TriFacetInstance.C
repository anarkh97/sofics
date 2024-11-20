#include <Mortar.d/FFIPolygon.d/TriFacet.C>
#include <Mortar.d/MortarElement.d/MortarElement.h>

template
TriFacetTemplate<double,FaceElement>
::TriFacetTemplate();

template
TriFacetTemplate<double,FaceElement>
::TriFacetTemplate(FaceElement*, const double*, const double*, const double*);

template
TriFacetTemplate<double,FaceElement>
::~TriFacetTemplate();

template
void
TriFacetTemplate<double,FaceElement>
::Initialize();

template
void
TriFacetTemplate<double,FaceElement>
::SetTriFacet(FaceElement*, const double*, const double*, const double*);

template
void
TriFacetTemplate<double,FaceElement>
::SetTriFacet(FaceElement*, const double* m1, const double* m2, const double* m3,
              const double* dm1, const double* dm2, const double* dm3, int nbDer);

template
void
TriFacetTemplate<double,FaceElement>
::SetTriFacet(FaceElement*, const double* m1, const double* m2, const double* m3,
              const double* dm1, const double* dm2, const double* dm3,
              const double* d2m1, const double* d2m2, const double* d2m3, int nbDer);

template
void
TriFacetTemplate<double,FaceElement>
::Print();

template
double
TriFacetTemplate<double,FaceElement>
::MappingJacobian();

template
void
TriFacetTemplate<double,FaceElement>
::LocalToLocalCoordOnFaceEl(const double*, double*);

template
double
TriFacetTemplate<double,FaceElement>
::GetJacobianOnFaceEl<CoordSet>(const double*, CoordSet&);

template
double
TriFacetTemplate<double,FaceElement>
::GetIsoParamMappingNormalAndJacobianOnFaceEl<CoordSet>(double*, const double*, CoordSet&);

template
FullM
TriFacetTemplate<double,FaceElement>
::IntegrateShapeFctProduct<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, int);

template
FullM
TriFacetTemplate<double,FaceElement>
::IntegrateNormalShapeFctProduct<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, int);

template
Vector
TriFacetTemplate<double,FaceElement>
::IntegrateNormalGeoGap<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, CoordSet&, int, double);

template
FullM 
TriFacetTemplate<double,FaceElement>
::IntegrateGradNormalGeoGap<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, CoordSet&, int, double);

template
FullM 
TriFacetTemplate<double,FaceElement>
::IntegrateHessNormalGeoGap<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, CoordSet&, double*, int, double);

template
FullM
TriFacetTemplate<double,FaceElement>
::IntegrateLocalGradNormalGeoGap<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, CoordSet&, int, double);

template
FullM
TriFacetTemplate<double,FaceElement>
::IntegrateLocalHessNormalGeoGap<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, CoordSet&, double*, int, double);

template
FullM
TriFacetTemplate<double,FaceElement>
::IntegrateLocalGradGradNormalGeoGap<CoordSet,MortarElement,FaceElement>(MortarElement*, TriFacet&, CoordSet&, CoordSet&, double*, int, double);

