#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/ShearPanelTemplate.cpp>

template
void
ShearPanelTemplate<double>
::spstress(double *_xg, double *_yg, double *_zg, double *_v,
           double G, double E, double F1, double F2,
           double *_stress, double *_strain,
           double &vmssig, double &vmseps);

template
void
ShearPanelTemplate<double>
::vmssWRTdisp(double *_xg, double *_yg, double *_zg, double *_v,
              double G, double E,
              double *_vmsWRTdisp, double &vmssig);

#endif

