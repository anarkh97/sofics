// --------------------------------------------------------------
// HB - 05/29/03
// --------------------------------------------------------------
#ifndef _TRIFACET_H_ 
#define _TRIFACET_H_

#include <Mortar.d/MortarDefines.h>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

template <class Scalar> class GenFullM;
template <class Scalar> class GenVector;

template <class Scalar, class FaceType>
class TriFacetTemplate {
  private:
        Scalar LocalCoordOnFaceEl[3][2]; // coord. of the vertices in the ref. face el.
#ifdef USE_EIGEN3
        Eigen::Array<Eigen::Matrix<Scalar,Eigen::Dynamic,1,Eigen::ColMajor,Eigen::Dynamic>,3,2> PartialLocalCoordOnFaceEl;
        Eigen::Array<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor,Eigen::Dynamic,Eigen::Dynamic>,3,2> SecondPartialLocalCoordOnFaceEl;
#endif
        FaceType* FaceEl;                // ptr to the associated face el. 

  public:
        // Constructors
        // ~~~~~~~~~~~~
        TriFacetTemplate();
        TriFacetTemplate(FaceType*, const Scalar* m1, const Scalar* m2, const Scalar* m3);

        // Destructor
        // ~~~~~~~~~~
        ~TriFacetTemplate();

        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Get methods  
        // ~~~~~~~~~~~~
        int nVertices() { return 3; }
        int nNodes() { return nVertices(); }
        FaceType* GetPtrFaceEl() { return FaceEl; }
        Scalar GetLocalCoordOnFaceElem(int i, int j) { return LocalCoordOnFaceEl[i][j]; }
#ifdef USE_EIGEN3
        Eigen::Matrix<Scalar,Eigen::Dynamic,1,Eigen::ColMajor,Eigen::Dynamic>&
          GetPartialLocalCoordOnFaceEl(int i, int j) { return PartialLocalCoordOnFaceEl(i,j); }
        Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor,Eigen::Dynamic,Eigen::Dynamic>&
          GetSecondPartialLocalCoordOnFaceEl(int i, int j) { return SecondPartialLocalCoordOnFaceEl(i,j); }
#endif
        // Set methods
        // ~~~~~~~~~~~ 
        void SetTriFacet(FaceType*, const Scalar* m1, const Scalar* m2, const Scalar* m3);
        void SetTriFacet(FaceType*, const Scalar* m1, const Scalar* m2, const Scalar* m3, 
                         const Scalar* dm1, const Scalar* dm2, const Scalar* dm3, int nbDer);
        void SetTriFacet(FaceType*, const Scalar* m1, const Scalar* m2, const Scalar* m3,
                         const Scalar* dm1, const Scalar* dm2, const Scalar* dm3,
                         const Scalar* d2m1, const Scalar* d2m2, const Scalar* d2m3, int nbDer);

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void Print();
      
        // Mapping methods
        // ~~~~~~~~~~~~~~~ 
        Scalar MappingJacobian();
        Scalar MappingJacobian(const Scalar* m) { return MappingJacobian(); }
        void LocalToLocalCoordOnFaceEl(const double* m, Scalar* mOnFaceEl);
        template<typename CoordSetType>
          Scalar GetJacobianOnFaceEl(const Scalar* m, CoordSetType& cs);
        template<typename CoordSetType>
          Scalar GetIsoParamMappingNormalAndJacobianOnFaceEl(Scalar* Normal, const Scalar* m, CoordSetType& cs);

        // Integration methods
        // ~~~~~~~~~~~~~~~~~~~ 
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenFullM<Scalar> IntegrateShapeFctProduct(MortarType*, TriFacetTemplate<Scalar,FriendFaceType>&, CoordSetType&, int ngp=6);
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenFullM<Scalar> IntegrateNormalShapeFctProduct(MortarType*, TriFacetTemplate<Scalar,FriendFaceType>&, CoordSetType&, int ngp=6);
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenVector<Scalar> IntegrateNormalGeoGap(MortarType*, TriFacetTemplate<Scalar,FriendFaceType>&, CoordSetType&, CoordSetType&,
                                                  int ngp=6, double offset=0.);
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenFullM<Scalar> IntegrateGradNormalGeoGap(MortarType*, TriFacetTemplate<Scalar,FriendFaceType>&, CoordSetType&, CoordSetType&,
                                                     int ngp=6, double offset=0.);
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenFullM<Scalar> IntegrateHessNormalGeoGap(MortarType*, TriFacetTemplate<Scalar,FriendFaceType>&, CoordSetType&, CoordSetType&,
                                                     double*, int ngp=6, double offset=0.);
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenFullM<Scalar> IntegrateLocalGradNormalGeoGap(MortarType*, TriFacetTemplate<Scalar,FriendFaceType>&, CoordSetType&, CoordSetType&,
                                                          int ngp=6, double offset=0.);
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenFullM<Scalar> IntegrateLocalHessNormalGeoGap(MortarType*, TriFacetTemplate<Scalar,FriendFaceType>&, CoordSetType&, CoordSetType&,
                                                          double*, int ngp=6, double offset=0.);
        template<typename CoordSetType, typename MortarType, typename FriendFaceType>
          GenFullM<Scalar> IntegrateLocalGradGradNormalGeoGap(MortarType* MortarEl, TriFacetTemplate<Scalar,FriendFaceType>& FriendFacet,
                                                              CoordSetType& cs, CoordSetType& cs2, double *mu, int ngp=6, double offset=0.);
};

class FaceElement;
typedef TriFacetTemplate<double, FaceElement> TriFacet;

#endif
