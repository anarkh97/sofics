// ---------------------------------------------------------------- 
// HB - 05/29/03
// HB - Modification 03/02/04
// ---------------------------------------------------------------- 
#ifndef _FFIPOLYGON_H_ 
#define _FFIPOLYGON_H_

// STL
#include <vector>
#include <utility> // std::pair

// FEM
#include <Math.d/matrix.h>
#include <Math.d/Vector.h>
#include <Mortar.d/FFIPolygon.d/TriFacet.h>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

class CoordSet;
class SurfaceEntity;

template<class MasterFaceType, class SlaveFaceType, class MortarType>
class FFIPolygonTemplate {
  private:
        double Area;                 // (approximate) area 
        int nVertices;               // number of polygon vertices

        MasterFaceType* MasterFace;     // ptr to the associated master face el.
        SlaveFaceType* SlaveFace;       // ptr to the associated slave face el.

        typedef std::pair<TriFacetTemplate<double,MasterFaceType>,               // first  -> master
                          TriFacetTemplate<double,SlaveFaceType> > Facet_pair_t; // second -> slave
        std::vector<Facet_pair_t> Facets;

        FullM M;                     // store the (FFI contribution to) M matrix (Mortar-Slave)
        FullM N;                     // store the (FFI contribution to) N matrix (Mortar-Master)

        Vector NormalGeoGaps;        // store the (FFI contribution to) the "geometrical" normal gaps (Slave-Master)

        FullM *dM;                   // store the (FFI contribution to) derivative of the M matrix (Mortar-Slave)
        FullM *dN;                   // store the (FFI contribution to) derivative of the N matrix (Mortar-Master)

        FullM H;                     // store the (FFI contribution to) derivative of [-N,M] contracted with mu (vector of Lagrange multipliers)
#ifdef HB_ACME_FFI_DEBUG
        double (*VertexLlCoordOnSFaceEl)[2];
        double (*VertexLlCoordOnMFaceEl)[2];
#endif
        double* ACME_FFI_LocalCoordData;

        // Helper methods
        // ~~~~~~~~~~~~~~
        static TriFacetTemplate<double,MasterFaceType>& MasterFacet(std::vector<Facet_pair_t>& FacetSet, size_t i)
        {
          return FacetSet[i].first;
        }
        static TriFacetTemplate<double,SlaveFaceType>& SlaveFacet(std::vector<Facet_pair_t>& FacetSet, size_t i)
        {
          return FacetSet[i].second;
        }

  public:
        // Constructors
        // ~~~~~~~~~~~~
        FFIPolygonTemplate();
        FFIPolygonTemplate(MasterFaceType* MasterFaceEl, SlaveFaceType* SlaveFaceEl);
        FFIPolygonTemplate(MasterFaceType* MasterFaceEl, SlaveFaceType* SlaveFaceEl, int nVert, double* ACME_FFI_Data,
                           int ACME_FFI_Derivatives_Order = 0);
        
        // Destructor 
        // ~~~~~~~~~~
        ~FFIPolygonTemplate();

        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Set methods
        // ~~~~~~~~~~~ 
        void SetFFIPolygon(MasterFaceType* MasterFaceEl, SlaveFaceType* SlaveFaceEl, int nVert, double* ACME_FFI_Data,
                           int ACME_FFI_Derivatives_Order);
        void SetPtrMasterFace(MasterFaceType*); 
        void SetPtrSlaveFace(SlaveFaceType*);
        
        // Get methods  
        // ~~~~~~~~~~~
        int GetnVertices() { return nVertices; }

        int GetnFacets() { return Facets.size(); }

        MasterFaceType* GetPtrMasterFace() { return MasterFace; }

        SlaveFaceType* GetPtrSlaveFace() { return SlaveFace; }

        FullM* GetPtrM() { return &M; }

        FullM* GetPtrN() { return &N; }

        Vector* GetPtrNormalGeoGaps() { return &NormalGeoGaps; }

        FullM* GetdM() { return dM; }

        FullM* GetdN() { return dN; }

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void Print(); 
        void PrintM();
        void PrintN();
#ifdef HB_ACME_FFI_DEBUG
        void PrintSlaveVertices(FILE* file, CoordSet& cs, int& firstVertId);
        void PrintMasterVertices(FILE* file, CoordSet& cs, int& firstVertId);
        void PrintFFIPolygonTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode);
#endif
        // Triangularization methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        void CreateTriangularization(double*, bool use_acme_centroid = false);
        void CreateTriangularizationExp(double*, bool use_acme_centroid = false);
        void CreateTriangularizationExp2(double*, bool use_acme_centroid = false);

        // Integration methods
        // ~~~~~~~~~~~~~~~~~~~ 
        FullM IntegrateOnMaster_MasterShapeFctProduct(MortarType*, CoordSet&, int ngp=6);
        FullM IntegrateOnMaster_SlaveShapeFctProduct(MortarType*, CoordSet&, int ngp=6);

        FullM IntegrateOnSlave_MasterShapeFctProduct(MortarType*, CoordSet&, int ngp=6);
        FullM IntegrateOnSlave_SlaveShapeFctProduct(MortarType*, CoordSet&, int ngp=6);

        FullM IntegrateOnSlave_MasterNormalShapeFctProduct(MortarType*, CoordSet &, int ngp=6);
        FullM IntegrateOnSlave_SlaveNormalShapeFctProduct(MortarType*, CoordSet &, int ngp=6);

        FullM IntegrateOnSlave_GradNormalGeoGap(MortarType*, CoordSet &, CoordSet &, int ngp=6, double offset=0.);
        FullM IntegrateOnSlave_HessNormalGeoGap(MortarType*, CoordSet &, CoordSet &, double*, int ngp=6, double offset=0.);

        void ComputeM(MortarType*, CoordSet&, int ngp=6);
        void ComputeN(MortarType*, CoordSet&, int ngp=6);

        void ComputeNormalM(MortarType*, CoordSet&, int ngp=6);
        void ComputeNormalN(MortarType*, CoordSet&, int ngp=6);

        void ComputeNormalGeoGap(MortarType*, CoordSet&, CoordSet&, int ngp=6, double offset=0.);
        void ComputeGradNormalGeoGap(MortarType*, CoordSet&, CoordSet&, int ngp=6, double offset=0.);
        void ComputeHessNormalGeoGap(MortarType*, CoordSet&, CoordSet&, double*, int ngp=6, double offset=0.);
};

class FaceElement;
class MortarElement;
typedef FFIPolygonTemplate<FaceElement,FaceElement,MortarElement> FFIPolygon;

#endif
