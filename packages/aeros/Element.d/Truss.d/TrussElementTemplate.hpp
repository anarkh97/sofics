#ifndef _TRUSSELEMENTSEMITEMPLATE_HPP_
#define _TRUSSELEMENTSEMITEMPLATE_HPP_

#include <string>

template<typename doublereal>
class TrussElementTemplate 
{
#ifdef USE_EIGEN3
  public:
    void stiffnessMatrix(doublereal *_stiffness, doublereal *_x,doublereal *_y, doublereal *_z,
                         doublereal E, doublereal A, doublereal preload); 
    void vmst(doublereal *_stress, doublereal *_x,doublereal *_y, doublereal *_z, doublereal *_elDisp,
              doublereal E, doublereal A, doublereal W, doublereal Ta, doublereal preload, 
              int strInd, int surface=0, doublereal *_ndTemps=0, doublereal ylayer=0, doublereal zlayer=0, int avgnum=1);

  protected:
#endif
};

#endif
