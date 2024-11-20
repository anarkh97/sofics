// ---------------------------------------------------------------- 
// HB - 08/26/03
// ---------------------------------------------------------------- 
#ifndef _DUALMORTARTRI3_H_
#define _DUALMORTARTRI3_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.h>

class CoordSet;
//template <class Scalar> class GenFullM;
//typedef GenFullM<double> FullM;
class FaceElement;

class DualMortarTri3: public MortarElement {
  private:
        // !!! NO NEED TO COMPUTE THE Alpha COEFFS:
        // -> THEY ARE A PRIORI INDEPENDANT OF THE FACE ELEMENT GEOMETRY
        //double Alpha[3]; // coeffs for computing dual shape fcts 
                           // from the std ones: 
                           //  Phidual[i] = Sum{alpha[i,j].PhiStd[j]} 
  public:
        // Constructors
        // ~~~~~~~~~~~~
        DualMortarTri3();
        DualMortarTri3(FaceElement*);
        DualMortarTri3(FaceElement*, CoordSet&);
        
        DualMortarTri3(double, FaceElement*);  
        DualMortarTri3(double, FaceElement*, CoordSet&);  
        
        // Destructor 
        // ~~~~~~~~~~
        virtual ~DualMortarTri3();
        
        // Get methods
        // ~~~~~~~~~~~
        // -> implementation of virtual methods
        int nNodes() override;
        int nMortarShapeFct() override;
        bool GetDualFlag() override { return true; }

        // Shape fct methods
        // ~~~~~~~~~~~~~~~~~
        // -> local methods
        //void ComputeDualCoeffs(CoordSet&);
        template<typename Scalar>
          void GetShapeFctVal(Scalar* Shape, Scalar* m);

        // -> implementation of virtual methods
        void GetShapeFctVal(double* Shape, double* m) override;
};

template<typename Scalar>
void
DualMortarTri3::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  Scalar StdShape[3];
  static_cast<FaceTri3*>(GetPtrMasterFace())->GetShapeFctVal(StdShape, m);

  Shape[0] = 3.*StdShape[0] -    StdShape[1] -    StdShape[2];
  Shape[1] =  - StdShape[0] + 3.*StdShape[1] -    StdShape[2];
  Shape[2] =  - StdShape[0] -    StdShape[1] + 3.*StdShape[2];
}

#endif
