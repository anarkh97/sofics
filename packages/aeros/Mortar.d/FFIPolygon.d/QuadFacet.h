// --------------------------------------------------------------	
// HB - 02/18/07
// --------------------------------------------------------------	
#ifndef _QUADFACET_H_ 
#define _QUADFACET_H_

class FaceElement;
class MortarElement;
class CoordSet;

template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

class QuadFacet {
  private:
	double Area;                     // (approximate) area 
        double LocalCoordOnFaceEl[4][2]; // coord. of the vertices in the ref. face el.
        FaceElement* FaceEl;             // ptr to the associated face el. 

  public:
        // Constructors
        // ~~~~~~~~~~~~
        QuadFacet();
        QuadFacet(FaceElement*, const double* m1, const double* m2, const double* m3, const double* m4);

        // Destructor
        // ~~~~~~~~~~
        ~QuadFacet();

        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Get methods  
        // ~~~~~~~~~~~~
        int nVertices() { return 4; }

        int nNodes() { return nVertices(); }

        double GetArea() { return Area; };

        FaceElement* GetPtrFaceEl() { return FaceEl; }

        // Set methods
        // ~~~~~~~~~~~ 
        void SetQuadFacet(FaceElement*, const double* m1, const double* m2, const double* m3);
        void SetArea(double);

        // Print, display methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void Print();
      
        // Mapping methods
        // ~~~~~~~~~~~~~~~ 
        double MappingJacobian(const double* m);
        void LocalToLocalCoordOnFaceEl(const double* m, double*  mOnFaceEl);
        double GetJacobianOnFaceEl(const double* m, CoordSet& cs);
        double GetIsoParamMappingNormalAndJacobianOnFaceEl(double* Normal, const double* m, CoordSet& cs);

        // Integration methods
        // ~~~~~~~~~~~~~~~~~~~ 
        FullM IntegrateShapeFctProduct(MortarElement*, QuadFacet&, CoordSet&, int ngp=6);

        FullM IntegrateNormalShapeFctProduct(MortarElement*, QuadFacet&, CoordSet&, int ngp=6);
        FullM IntegrateGradNormalShapeFctProduct(MortarElement*, QuadFacet&, CoordSet&, CoordSet&, QuadFacet&, CoordSet&, int ngp=6)

        // Compute normal "geometrical" gap
        Vector IntegrateNormalGeoGagsProduct(MortarElement* MortarEl, QuadFacet& FriendQuadFacet, CoordSet& cs, CoordSet& cs1, int ngp=6);
};
#endif
