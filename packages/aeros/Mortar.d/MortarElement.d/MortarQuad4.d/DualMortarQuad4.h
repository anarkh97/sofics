// ---------------------------------------------------------------- 
// HB - 08/13/03
// ---------------------------------------------------------------- 
#ifndef _DUALMORTARQUAD4_H_
#define _DUALMORTARQUAD4_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad4.d/FaceQuad4.h>

class CoordSet;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class FaceElement;

class DualMortarQuad4: public MortarElement {
  private:
        FullM* Alpha; // coeffs for computing dual shape fcts 
                      // from the std ones: Phidual[i] = Sum{alpha[i,j].PhiStd[j]} 
  public:
        // Constructors
        // ~~~~~~~~~~~~
        DualMortarQuad4();
        DualMortarQuad4(FaceElement*);
        DualMortarQuad4(FaceElement*, CoordSet&);
        
        DualMortarQuad4(double, FaceElement*);  
        DualMortarQuad4(double, FaceElement*, CoordSet&);  
        
        // Destructor 
        // ~~~~~~~~~~
        virtual ~DualMortarQuad4();

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
DualMortarQuad4::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  Scalar StdShape[4];
  static_cast<FaceQuad4*>(GetPtrMasterFace())->GetShapeFctVal(StdShape, m);
  for(int i=0; i<4; ++i) {
    Shape[i] = 0;
    for(int j=0; j<4; j++)
      Shape[i] += (*Alpha)[i][j]*StdShape[j]; 
  }
}

#endif
