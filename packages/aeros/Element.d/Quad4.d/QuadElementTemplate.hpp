#ifndef _QUADELEMENTTEMPLATE_HPP_
#define _QUADELEMENTTEMPLATE_HPP_

template<typename doublereal>
class QuadElementTemplate
{
  public:
    QuadElementTemplate() {}
    ~QuadElementTemplate() {}

    void sands2(char *escm, doublereal *_x, doublereal *_y, doublereal *_c, doublereal *_v,
                doublereal *_stress, doublereal *_strain,
                int maxgus, int maxstr, int elm, int msize, bool vmflg,
                bool strainFlg, doublereal tc=0, doublereal tref=0, doublereal *_ndtemps=0);

    void getcmt(doublereal rip, doublereal e, doublereal nu, doublereal *_c);

    void vms2WRTdisp(char *escm, doublereal *_x, doublereal *_y, doublereal *_c, doublereal *_v,
                     doublereal *_vmsWRTdisp, doublereal *_stressWRTdisp,
                     int maxgus, int maxstr, int elm, int msize, bool vmflg,
                     bool strainFlg, doublereal tc, doublereal tref, doublereal *_ndtemps);

  protected:
    void q4shpe(doublereal xi, doublereal eta, doublereal* x, doublereal* y,
                doublereal* s, doublereal* sx, doublereal* sy, doublereal& det);

    void vmelmv(doublereal *_stress, int maxgus, int maxstr, int msize, int elm, int nno);

    void strainvm(doublereal *_strain, int maxgus, int maxstr, int msize, int numnod);
};

#endif
