#ifndef _SHEARPANELTEMPLATE_HPP_
#define _SHEARPANELTEMPLATE_HPP_

#include <Element.d/Quad4.d/QuadElementTemplate.hpp>

template<typename doublereal>
class ShearPanelTemplate: public QuadElementTemplate<doublereal> 
{
  public: 
    ShearPanelTemplate() {}
    ~ShearPanelTemplate() {}

    void spstress(doublereal *_xg, doublereal *_yg, doublereal *_zg, doublereal *_v,
                  doublereal G, doublereal E, doublereal F1, doublereal F2,
                  doublereal *_stress, doublereal *_strain,
                  doublereal &vmssig, doublereal &vmseps);
    void vmssWRTdisp(doublereal *_xg, doublereal *_yg, doublereal *_zg, doublereal *_v,
                     doublereal G, doublereal E, doublereal *_vmsWRTdisp, doublereal &vmssig);
};
#endif
