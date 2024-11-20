// ---------------------------------------------------------------- 
// HB - 08/13/03
// ---------------------------------------------------------------- 
#ifndef _DUALMORTARTRI6_H_
#define _DUALMORTARTRI6_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceTri6.d/FaceTri6.h>

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class FaceElement;

class DualMortarTri6: public MortarElement {
  private:
        FullM* Alpha; // coeffs for computing dual shape fcts 
                      // from the std ones: Phidual[i] = Sum{alpha[i,j].PhiStd[j]} 
  public:
        // Constructors
        // ~~~~~~~~~~~~
        DualMortarTri6();
        DualMortarTri6(FaceElement*);
        DualMortarTri6(FaceElement*, CoordSet&);
        
        DualMortarTri6(double, FaceElement*);  
        DualMortarTri6(double, FaceElement*, CoordSet&);  
        
        // Destructor 
        // ~~~~~~~~~~
        virtual ~DualMortarTri6();

        // Initialization & clean/clear methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();
 
        // Get methods
        // ~~~~~~~~~~~
        // -> implementation of virtual methods
        int nNodes() override;
        int nMortarShapeFct() override;
        bool GetDualFlag() override { return true; }

        // Shape fct methods
        // ~~~~~~~~~~~~~~~~~
        // -> local methods
        void ComputeDualCoeffs(CoordSet&);
        template<typename Scalar>
          void GetShapeFctVal(Scalar* Shape, Scalar* m);

        // -> implementation of virtual methods
        void GetShapeFctVal(double* Shape, double* m) override;
};

// !!! ASSUME THAT DUAL MORTAR COEFFS Alpha HAVE ALREADY BEEN COMPUTED !!!
template<typename Scalar>
void
DualMortarTri6::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  Scalar StdShape[6];
  static_cast<FaceTri6*>(GetPtrMasterFace())->GetShapeFctVal(StdShape, m);
  for(int i=0; i<6; ++i) {
    Shape[i] = 0;
    for(int j=0; j<6; j++)
      Shape[i] += (*Alpha)[i][j]*StdShape[j]; 
  }
}

#endif
