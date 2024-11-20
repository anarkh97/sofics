#ifndef _TETRAELEMENTTEMPLATE_HPP_
#define _TETRAELEMENTTEMPLATE_HPP_

template<typename doublereal>
class TetraElementTemplate 
{
  public:
    TetraElementTemplate() {}

    ~TetraElementTemplate() {
    }

    void sands23(int elm, doublereal *_x, doublereal *_y, doublereal *_z, doublereal e, doublereal xnu, doublereal *_u, 
                 doublereal *_stress, doublereal *_strain, int maxgus, int maxstr, int msize, int outerr, bool vmflg, bool strainFlg);  

    void vmelmv(doublereal *_stress, int maxgus, int maxstr, int msize, int elm, int nno);

    void strainvm(doublereal *_strain, int maxgus, int maxstr, int msize, int numnod);      

};

#endif

