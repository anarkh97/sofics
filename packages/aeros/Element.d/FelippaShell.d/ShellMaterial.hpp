#ifndef _SHELLMATERIAL_HPP_
#define _SHELLMATERIAL_HPP_

#ifdef USE_EIGEN3
#include <iostream>
#include <vector>
#include <cstdio>
#include <stdexcept>
#include <Eigen/Core>

template<typename doublereal>
class ShellMaterial
{
  public:
    virtual ~ShellMaterial() {}
    virtual void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                         doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                         doublereal *staten = 0, doublereal *statenp = 0) = 0; 
                                                                          // Upsilon is the generalized "strains" {e,chi}
                                                                          // Sigma is the generalized "stress" {N,M}
                                                                          // D is the tangent constitutive matrix { Dm, Dmb; Dbm, Db }
                                                                          // temp is the temperature
    virtual doublereal* GetCoefOfConstitutiveLaw() { return NULL; }
    virtual doublereal GetShellThickness() = 0;
    virtual doublereal GetAreaDensity() = 0; // mass per unit area
    virtual doublereal GetAmbientTemperature() = 0;
    virtual void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                              doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                              doublereal *staten = 0, doublereal *statenp = 0) = 0;
    virtual int GetNumStates() { return 0; }
    virtual int GetNumLocalStates() { return 0; }
    virtual void SetState(doublereal *state) {}
    virtual void GetState(doublereal *state) {}
    virtual void UpdateState(doublereal *Upsilon, doublereal *staten, doublereal *statenp, int gp, doublereal temp, doublereal dt = 0) {}
    virtual std::vector<doublereal> GetLocalPlasticStrain(int nd, doublereal z, doublereal *statenp = 0) { return std::vector<doublereal>(); }
    virtual std::vector<doublereal> GetLocalBackStress(int nd, doublereal z, doublereal *statenp = 0) { return std::vector<doublereal>(); }
    virtual doublereal GetLocalEquivalentPlasticStrain(int nd, doublereal z, doublereal *statenp = 0) { return 0; }
    virtual doublereal GetLocalDamage(int nd, doublereal z, doublereal *statenp = 0) { return 0; }
    virtual doublereal GetDissipatedEnergy(int gp, doublereal *statenp = 0) { return 0; }
    virtual bool CheckFailure(doublereal *statenp = 0) { return false; } // used to initiate element deletion

    virtual void GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh,
                                                           doublereal *eframe, int gp, doublereal temp) = 0;
    virtual void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                           doublereal *eframe, int gp) = 0;
    virtual void GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh,
                                                                doublereal *, int) = 0;
    virtual void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                                doublereal *eframe, int gp) = 0;

    static Eigen::Matrix<doublereal,3,3>
    andesinvt(doublereal *_eframe, doublereal *_aframe, doublereal thetaf);
};

//     ------------------------------------------------ 
//     ISOTROPIC (ESSENTIALLY PLANE) MATERIAL           
//     NO COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS 
//     ------------------------------------------------ 
template<typename doublereal>
class ShellMaterialType0 : public ShellMaterial<doublereal>
{
    doublereal E;   // Young's modulus
    doublereal h;   // shell thickness
    doublereal nu;  // Poisson's ratio
    doublereal rho; // volume density
    doublereal Ta;  // ambient temperature
    doublereal w;   // coefficient of thermal expansion
  public:
    ShellMaterialType0(doublereal _E, doublereal _h, doublereal _nu, doublereal _rho,
                       doublereal _Ta = 0., doublereal _w = 0.) 
      : E(_E), h(_h), nu(_nu), rho(_rho), Ta(_Ta), w(_w) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                 doublereal *staten = 0, doublereal *statenp = 0);
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                      doublereal *staten = 0, doublereal *statenp = 0);
    doublereal GetShellThickness() { return h; }
    doublereal GetAreaDensity() { return rho*h; }
    doublereal GetAmbientTemperature() { return Ta; }

    void GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh,
                                                   doublereal *eframe, int gp, doublereal temp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh,
                                                        doublereal *, int);
    void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                        doublereal *eframe, int gp);
};

//     ------------------------------------------------- 
//     COMPOSITE MATERIAL WITH KNOWN CONSTITUTIVE MATRIX 
//     ------------------------------------------------- 
template<typename doublereal>
class ShellMaterialType1 : public ShellMaterial<doublereal>
{
    Eigen::Map<Eigen::Matrix<doublereal,6,6,Eigen::RowMajor> > coef; // coefficients of the constitutive law
    doublereal *aframe; // arbitrary 3x3 frame of the constitutive law
    doublereal rhoh;    // area density (mass per unit surface area)
    doublereal h;       // shell thickness
    doublereal Ta;      // ambient temperature
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Alpha; // coefficients of thermal expansion
    static bool Wlocal_stress, Wlocal_stress_disp, Wlocal_stress_thic;
  public:
    ShellMaterialType1(doublereal *_coef, doublereal *_aframe, doublereal _rhoh, doublereal _h = 0.,
                       doublereal _Ta = 0.)
      : coef(_coef), aframe(_aframe), rhoh(_rhoh), h(_h), Ta(_Ta), Alpha(_coef+36) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                 doublereal *staten = 0, doublereal *statenp = 0);
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                      doublereal *staten = 0, doublereal *statenp = 0);
    doublereal* GetCoefOfConstitutiveLaw() { return coef.data(); }
    doublereal GetShellThickness();
    doublereal GetAreaDensity() { return rhoh; }
    doublereal GetAmbientTemperature() { return Ta; }
    
    void GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh,
                                                   doublereal *eframe, int gp, doublereal temp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh,
                                                        doublereal *, int);
    void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                        doublereal *eframe, int gp);
};

//     ---------------------------------------------- 
//     COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES 
//     ---------------------------------------------- 
template<typename doublereal>
class ShellMaterialTypes2And3 : public ShellMaterial<doublereal>
{
    int nlayer;         // number of layers of the composite element
    Eigen::Map<Eigen::Matrix<doublereal,12,Eigen::Dynamic> > mtlayer; // material properties of each layer
    bool couple;        // type of constitutive law: true  -> bending-membrane coupling
                        //                           false -> no bending-membrane coupling
    doublereal *aframe; // arbitrary 3x3 frame of the constitutive law
    doublereal h;       // shell thickness
    doublereal rhoh;    // area density (mass per unit surface area)
    doublereal Ta;      // ambient temperature
    doublereal nsm;     // non-structural mass per unit surface area

  public:
    ShellMaterialTypes2And3(int _nlayer, doublereal *_mtlayer, bool _couple, doublereal *_aframe, doublereal _Ta = 0.,
                            doublereal _nsm = 0.);

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                 doublereal *staten = 0, doublereal *statenp = 0);
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                      doublereal *staten = 0, doublereal *statenp = 0);
    doublereal GetShellThickness() { return h; }
    doublereal GetAreaDensity() { return rhoh; }
    doublereal GetAmbientTemperature() { return Ta; }

    void GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh,
                                                   doublereal *eframe, int gp, doublereal temp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh,
                                                        doublereal *, int);
    void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                        doublereal *eframe, int gp);
};

//     ------------------------------------------------ 
//     ISOTROPIC LINEAR ELASTIC - J2 PLASTIC MATERIAL   
//     ------------------------------------------------ 
template<typename doublereal, typename localmaterial>
class ShellMaterialType4 : public ShellMaterial<doublereal>
{
    doublereal h;   // shell thickness
    doublereal nu;  // Poisson's ratio
    doublereal rho; // density
    localmaterial **mat;
    int nlayer;     // number of material points through the thickness of the shell
    int maxgus;     // number of material points over the area of the shell, per layer
  public:
    ShellMaterialType4(doublereal _h, doublereal _nu, doublereal _rho, localmaterial *_mat, int _nlayer, int _maxgus)
      : h(_h), nu(_nu), rho(_rho), nlayer(_nlayer), maxgus(_maxgus) {
      mat = new localmaterial * [_nlayer*_maxgus];
      for (int i = 0; i < _nlayer*_maxgus; ++i) mat[i] = _mat->Clone();
    }
    ~ShellMaterialType4() {
      for(int i = 0; i < nlayer*maxgus; ++i) delete mat[i];
      delete [] mat;
    }

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                 doublereal *staten = 0, doublereal *statenp = 0);
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                      doublereal *staten = 0, doublereal *statenp = 0);
    doublereal GetShellThickness() { return h; }
    doublereal GetAreaDensity() { return rho*h; }
    doublereal GetAmbientTemperature() { return 0.; }
    int GetNumStates() { return nlayer*maxgus*7; }
    int GetNumLocalStates() { return 7; }
    void SetState(doublereal *state);
    void GetState(doublereal *state);
    void UpdateState(doublereal *Upsilon, doublereal *staten, doublereal *statenp, int gp, doublereal temp, doublereal dt = 0);
    std::vector<doublereal> GetLocalPlasticStrain(int nd, doublereal z, doublereal *statenp = 0);
    std::vector<doublereal> GetLocalBackStress(int nd, doublereal z, doublereal *statenp = 0);
    doublereal GetLocalEquivalentPlasticStrain(int nd, doublereal z, doublereal *statenp = 0);
    doublereal GetLocalDamage(int nd, doublereal z, doublereal *statenp = 0);
    doublereal GetDissipatedEnergy(int gp, doublereal *statenp = 0);
    bool CheckFailure(doublereal *statenp = 0);

    void GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh,
                                                   doublereal *eframe, int gp, doublereal temp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh,
                                                        doublereal *, int);
    void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                        doublereal *eframe, int gp);
};

//     --------------------------------------------------- 
//     ORTHOTROPIC MATERIAL WITH KNOWN CONSTITUTIVE MATRIX 
//     --------------------------------------------------- 
template<typename doublereal>
class ShellMaterialType5 : public ShellMaterial<doublereal>
{
    Eigen::Map<Eigen::Matrix<doublereal,6,6,Eigen::RowMajor> > coef; // coefficients of the constitutive law
    doublereal *aframe; // arbitrary 3x3 frame of the constitutive law
    doublereal rhoh;    // area density (mass per unit surface area)
    doublereal h;       // shell thickness
    doublereal Ta;      // ambient temperature
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Alpha; // coefficients of thermal expansion
    static bool Wlocal_stress, Wlocal_stress_disp, Wlocal_stress_thic;
  public:
    ShellMaterialType5(doublereal *_coef, doublereal *_aframe, doublereal _rhoh, doublereal _h,
                       doublereal _Ta = 0.)
      : coef(_coef), aframe(_aframe), rhoh(_rhoh), h(_h), Ta(_Ta), Alpha(_coef+36) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                 doublereal *staten = 0, doublereal *statenp = 0);
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                      doublereal *staten = 0, doublereal *statenp = 0);
    doublereal* GetCoefOfConstitutiveLaw() { return coef.data(); }
    doublereal GetShellThickness() { return h; }
    doublereal GetAreaDensity() { return rhoh; }
    doublereal GetAmbientTemperature() { return Ta; }

    void GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh,
                                                   doublereal *eframe, int gp, doublereal temp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh,
                                                        doublereal *, int);
    void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                        doublereal *eframe, int gp);
};

//     ------------------------------------------------ 
//     PLANE STRESS MATERIAL   
//     ------------------------------------------------ 
template<typename doublereal, typename localmaterial>
class ShellMaterialType6 : public ShellMaterial<doublereal>
{
    doublereal h;   // shell thickness
    doublereal rho; // density
    localmaterial *mat;
    int nlayer;     // number of material points through the thickness of the shell
    int maxgus;     // number of material points over the area of the shell, per layer
  public:
    ShellMaterialType6(localmaterial *_mat, int _nlayer, int _maxgus)
      : h(_mat->getThickness()), rho(_mat->getDensity()), mat(_mat), nlayer(_nlayer), maxgus(_maxgus) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp, doublereal dt = 0,
                                 doublereal *staten = 0, doublereal *statenp = 0);
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp, doublereal dt = 0, 
                                      doublereal *staten = 0, doublereal *statenp = 0);
    doublereal GetShellThickness() { return h; }
    doublereal GetAreaDensity() { return rho*h; }
    doublereal GetAmbientTemperature() { return mat->getReferenceTemperature(); }
    int GetNumStates() { return nlayer*maxgus*mat->getNumStates(); }
    int GetNumLocalStates() { return mat->getNumStates(); }
    void UpdateState(doublereal *Upsilon, doublereal *staten, doublereal *statenp, int gp, doublereal temp, doublereal dt = 0);
    std::vector<doublereal> GetLocalPlasticStrain(int nd, doublereal z, doublereal *statenp = 0);
    std::vector<doublereal> GetLocalBackStress(int nd, doublereal z, doublereal *statenp = 0);
    doublereal GetLocalEquivalentPlasticStrain(int nd, doublereal z, doublereal *statenp = 0);
    doublereal GetLocalDamage(int nd, doublereal z, doublereal *statenp = 0);
    doublereal GetDissipatedEnergy(int gp, doublereal *statenp = 0);
    bool CheckFailure(doublereal *statenp = 0);

    void GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh,
                                                   doublereal *eframe, int gp, doublereal temp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh,
                                                        doublereal *, int);
    void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                        doublereal *eframe, int gp);
};

#endif

#endif
