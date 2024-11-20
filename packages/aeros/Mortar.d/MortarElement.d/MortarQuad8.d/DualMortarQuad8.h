// ---------------------------------------------------------------- 
// HB - 08/13/03
// ---------------------------------------------------------------- 
#ifndef _DUALMORTARQUAD8_H_
#define _DUALMORTARQUAD8_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad8.d/FaceQuad8.h>

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class FaceElement;

class DualMortarQuad8: public MortarElement {
  private:
        FullM* Alpha; // coeffs for computing dual shape fcts 
                      // from the std ones: Phidual[i] = Sum{alpha[i,j].PhiStd[j]} 
  public:
        // Constructors
        // ~~~~~~~~~~~~
        DualMortarQuad8();
        DualMortarQuad8(FaceElement*);
        DualMortarQuad8(FaceElement*, CoordSet&);
        
        DualMortarQuad8(double, FaceElement*);  
        DualMortarQuad8(double, FaceElement*, CoordSet&);  
        
        // Destructor 
        // ~~~~~~~~~~
        virtual ~DualMortarQuad8();

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
DualMortarQuad8::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  Scalar StdShape[8];
  static_cast<FaceQuad8*>(GetPtrMasterFace())->GetShapeFctVal(StdShape, m);
  for(int i=0; i<8; ++i) {
    Shape[i] = 0;
    for(int j=0; j<8; j++)
      Shape[i] += (*Alpha)[i][j]*StdShape[j]; 
  }
}

#endif
