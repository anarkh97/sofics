#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/NoBendingTriangle.hpp>

#define SHELLELEMENTTEMPLATE_INSTANTIATION_HELPER(Membrane,Bending) \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andescrd(int elm, double *x, double *y, double *z, double *rot, \
           double *xlp, double *ylp, double *zlp, double &area); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesfrm(int elm, double *x, double *y, double *z, double *aframe, \
           double *cframe); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesgf(int elm, double *x, double *y, double *z, double *_gravityForce, \
          double *gamma, int gravflg, double rhoh); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesgfWRTcoord(int elm, double *x, double *y, double *z, \
                  double *J, double *gamma, int gravflg, double rhoh); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesmm(int elm, double *x, double *y, double *z, double *emass, \
          double rhoh, int mflg); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesms(int elm, double *x, double *y, double *z, double *gamma, \
          double *grvfor, double &totmas, double rhoh); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesmsWRTcoord(int elm, double *x, double *y, double *z, double *J, \
                  double rhoh); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesstf(int elm, double *estiff, double *fint, double nu, \
           double *x, double *y, double *z, double *v, int ctyp, \
           ShellMaterial<double> *gpmat, int flag, int tflg, double *ndtemps, \
           double dt, double *staten, double *statenp); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesstfWRTthic(int elm, double *destiffdh, double *dfinitdh, double nu, \
                  double *x, double *y, double *z, double *u, int ctyp, \
                  ShellMaterial<double> *gpmat, int flag, int tflg, \
                  double *ndtemps); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesstfWRTcoord(int elm, double *destiffdx[9], double E, double nu, \
                   double rho, double eh, double Ta, double W, double *cFrame, \
                   double *x, double *y, double *z, int ctyp, \
                   double *coefs, int flag, int tflg, double *ndtemps); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesvms(int elm, int maxstr, double nu, double *x, double *y, double *z, \
           double *v, double *stress, int ctyp, ShellMaterial<double> *nmat, \
           int strainflg, int surface, int sflg, double *ndtemps, \
           int flag, double *staten, double *statenp); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesvmsWRTdisp(int elm, double nu, double *x, double *y, double *z, \
                  double *v, double *dvmsdu, int ctyp, \
                  ShellMaterial<double> *nmat, int surface, int sflg, \
                  double *ndtemps); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesvmsWRTthic(int elm, double nu, double *x, double *y, double *z, \
                  double *v, double *dvmsdh, int ctyp, \
                  ShellMaterial<double> *nmat, int surface, int sflg, \
                  double *ndtemps); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesvmsWRTcoord(int elm, double E, double nu, double rho, double eh, \
                   double Ta, double W, double *cFrame, double *x, \
                   double *y, double *z, double *v, double *dvmsdx, \
                   int ctyp, double *coefs, int surface, int sflg, \
                   double *ndtemps); \
template \
double \
ShellElementTemplate<double,Membrane,Bending> \
::equivstr(double sxx, double syy, double szz, double sxy); \
template \
Eigen::Matrix<double,1,18> \
ShellElementTemplate<double,Membrane,Bending> \
::equivstrSensitivityWRTdisp(double vms, double sxx, double syy, double szz, \
                             double sxy, Eigen::Matrix<double,3,18> &dsigmadu); \
template \
double \
ShellElementTemplate<double,Membrane,Bending> \
::equivstrSensitivityWRTthic(double vms, double sxx, double syy, double szz, \
                             double sxy, Eigen::Matrix<double,3,1> &dsigmadh); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::transform(double *lframe, double *gframe, double *str); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesups(int elm, double *staten, double *statenp, double *X, double *Y, double *Z, double *_v, \
           ShellMaterial<double> *gpmat, ShellMaterial<double> *nmat, int sflg, int tflg, \
           double *ndtemps, double dt); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesare(int elm, double *X, double *Y, double *Z, double &area); \
template \
void \
ShellElementTemplate<double,Membrane,Bending> \
::andesden(int elm, double *X, double *Y, double *Z, \
           ShellMaterial<double> *gpmat, double *statenp, double &D);

SHELLELEMENTTEMPLATE_INSTANTIATION_HELPER(EffMembraneTriangle,AndesBendingTriangle);
SHELLELEMENTTEMPLATE_INSTANTIATION_HELPER(EffMembraneTriangle,NoBendingTriangle);

#endif
