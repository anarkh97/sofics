#ifndef _SHELLELEMENTTEMPLATE_HPP_
#define _SHELLELEMENTTEMPLATE_HPP_

#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

// Note: For full compatability with the legacy implementation of elements 8 and 20, use:
//       (a) mflg = 0 for andesmm. In this case the mass augmentation is not used. The recommended setting
//           is mflg = 1 which typically gives a larger critical time-step for explicit dynamics.
//       (b) tflg = 0 for andesstf. In this case the temperature field is assumed to be constant (i.e. the
//           average of the nodal temperatures is used at each gauss point. The recommended setting is tflg
//           = 1 in which case the temperature at the gauss point is interpolated from the nodes using the
//           shape functions of a 3-node triangle.
//       (c) sflg = 0 for andesvms. In this case the higher-order contribution to the B matrix is not used
//           for stress/strain recovery. The recommended setting is sflg = 1 which results in more accurate
//           stress/strain recovery but is also more expensive.

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
class ShellElementTemplate : public Membrane<doublereal>, public Bending<doublereal>
{
  public:
    static void
    andesgf(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *gravityForce,
            doublereal *gamma, int gravflg, doublereal rhoh);

    static void
    andesgfWRTcoord(int elm, doublereal *x, doublereal *y, doublereal *z,
                    doublereal *J, doublereal *gamma, int gravflg, doublereal rhoh);

    static void
    andesmm(int elm, doublereal *x, doublereal *y, doublereal *z,
            doublereal *emass, doublereal rhoh, int mflg);

    static void
    andesms(int elm, doublereal *x, doublereal *y, doublereal *z,
            doublereal *gamma, doublereal *grvfor, doublereal &totmas,
            doublereal rhoh);

    static void
    andesmsWRTcoord(int elm, doublereal *x, doublereal *y, doublereal *z,
                    doublereal *J, doublereal rhoh);

    static void
    andesstf(int elm, doublereal *estiff, doublereal *fint, doublereal nu,
             doublereal *x, doublereal *y, doublereal *z, doublereal *u,
             int ctyp, ShellMaterial<doublereal> *gpmat, int flag,
             int tflg = 1, doublereal *ndtemps = 0, doublereal dt = 0.,
             doublereal *staten = 0, doublereal *statenp = 0);

    static void 
    andesstfWRTthic(int elm, doublereal *destiffdh, doublereal *dfintdh, doublereal nu,
                    doublereal *x, doublereal *y, doublereal *z, doublereal *u,
                    int ctyp, ShellMaterial<doublereal> *gpmat, int flag,
                    int tflg = 1, doublereal *ndtemps = 0);

    static void
    andesstfWRTcoord(int elm, doublereal *destiffdx[9], doublereal E,
                     doublereal nu, doublereal rho, doublereal eh,
                     doublereal Ta, doublereal W, doublereal *cFrame,
                     doublereal *x, doublereal *y, doublereal *z,
                     int ctyp, doublereal *coefs, int flag, int tflg = 1,
                     doublereal *ndtemps = 0);

    static void
    andesvms(int elm, int maxstr, doublereal nu, doublereal *x, doublereal *y,
             doublereal *z, doublereal *u, doublereal *stress, int ctyp,
             ShellMaterial<doublereal> *nmat, int strainflg, int surface,
             int sflg, doublereal *ndtemps = 0, int flag = 1,
             doublereal *staten = 0, doublereal *statenp = 0);

    static void
    andesvmsWRTdisp(int elm, doublereal nu, doublereal *x, doublereal *y,
                    doublereal *z, doublereal *u, doublereal *dvmsdu,
                    int ctyp, ShellMaterial<doublereal> *nmat, int surface,
                    int sflg, doublereal *ndtemps = 0);

    static void
    andesvmsWRTthic(int elm, doublereal nu, doublereal *x, doublereal *y,
                    doublereal *z, doublereal *u, doublereal *dvmsdh,
                    int ctyp, ShellMaterial<doublereal> *nmat, int surface,
                    int sflg, doublereal *ndtemps = 0);

    static void
    andesvmsWRTcoord(int elm, doublereal E, doublereal nu, doublereal rho,
                     doublereal eh, doublereal Ta, doublereal W, doublereal *cFrame,
                     doublereal *x, doublereal *y, doublereal *z, doublereal *u,
                     doublereal *dvmsdx, int ctyp, doublereal *coefs,
                     int surface, int sflg, doublereal *ndtemps = 0);

    static void
    andesups(int elm, doublereal *staten, doublereal *statenp, doublereal *X, doublereal *Y,
             doublereal *Z, doublereal *v, ShellMaterial<doublereal> *gpmat,
             ShellMaterial<doublereal> *nmat, int sflg, int tflg, doublereal *ndtemps,
             doublereal dt = 0);

    // (AN) function updated ... now recieves the t_n+1 state
    static void
    andesden(int elm, doublereal *X, doublereal *Y, doublereal *Z,
             ShellMaterial<doublereal> *gpmat, doublereal *statenp, doublereal &D);

    static void
    andesfrm(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *aframe,
             doublereal *cframe);

    static void
    andescrd(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *rot,
             doublereal *xlp, doublereal *ylp, doublereal *zlp, doublereal &area);

  private:
    static doublereal
    equivstr(doublereal sxx, doublereal syy, doublereal szz, doublereal sxy);

    static Eigen::Matrix<doublereal,1,18>
    equivstrSensitivityWRTdisp(doublereal vms, doublereal sxx, doublereal syy,
                               doublereal szz, doublereal sxy,
                               Eigen::Matrix<doublereal,3,18> &dsigmadu);

    static doublereal
    equivstrSensitivityWRTthic(doublereal vms, doublereal sxx, doublereal syy,
                               doublereal szz, doublereal sxy,
                               Eigen::Matrix<doublereal,3,1> &dsigmadh);

    static void
    transform(doublereal *lframe, doublereal *gframe, doublereal *str);

    static void
    andesare(int elm, doublereal *X, doublereal *Y, doublereal *Z, doublereal &area);

};

#endif
