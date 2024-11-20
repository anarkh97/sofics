#ifndef _TRIANGLE3ELEMENTTEMPLATE_HPP_
#define _TRIANGLE3ELEMENTTEMPLATE_HPP_

template<typename doublereal>
class Triangle3ElementTemplate
{
  public:
    Triangle3ElementTemplate() {}
    ~Triangle3ElementTemplate() {}

    void sands4(doublereal *_x, doublereal *_y, doublereal *_v,
                doublereal *_stress, doublereal *_strain,
                doublereal E, doublereal nu);

    void vms4WRTdisp(doublereal *_x, doublereal *_y, doublereal *_v,
                     doublereal *_vmsWRTdisp,
                     doublereal E, doublereal nu);

};

#endif

