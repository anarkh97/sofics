#ifdef USE_EIGEN3

#include <Element.d/Beam.d/BeamElementTemplate.cpp>


template
void
BeamElementTemplate<double>
::sands6(double area, double e, int elm,
         double *_stress, int maxsze, int maxgus, int maxstr,
         double *_eframe,
         double ix, double iy, double iz,
         double nu, double *_x, double *_y, double *_z,
         double *_ug, double alpha, double tref, double *_temp);

template
void
BeamElementTemplate<double>
::sands7(int elm, double A, double E,
         double *eframe, double Ix, double Iy, double Iz,
         double alphay, double alphaz, double C1, double nu,
         double *x, double *y, double *z, double *ug,
         double *stress, int numel, int maxgus, int maxstr,
         int msize, double alpha, double tref, double *temp);

/*
template
void
BeamElementTemplate<double>
::vmsWRTdisp(int elm, double A, double E,
             double *eframe, double Ix, double Iy, double Iz,
             double alphay, double alphaz, double C1, double nu,
             double *x, double *y, double *z, double *ug,
             double *vmsWRTdisp,
             double alpha, double tref, double *temp);
*/

template
void
BeamElementTemplate<double>
::transform(double *_l, double *_g, double *_str);

#endif
