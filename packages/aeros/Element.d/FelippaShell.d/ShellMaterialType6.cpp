#ifdef USE_EIGEN3
#include <cmath>
#include <stdexcept>
#include <vector>
#include <Eigen/Core>
#include <Math.d/TTensor.h>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <iostream>

extern int quietFlag;

template<typename doublereal, typename localmaterial>
void
ShellMaterialType6<doublereal,localmaterial>
::GetConstitutiveResponse(doublereal *_Upsilon, doublereal *_Sigma, doublereal *_D,
                          doublereal*, int point, doublereal temp, doublereal dt,
                          doublereal *staten, doublereal *statenp)
{
    // Local variables 
    int ilayer;
    doublereal z;
    Eigen::Matrix<doublereal,3,3> Ch;
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon), Sigma(_Sigma);
    Eigen::Map<Eigen::Matrix<doublereal,6,6> > D(_D);

    SymTensor<double,2> en, enp, stress;
    Eigen::Map<Eigen::Matrix<doublereal,3,1> > sigmah(&stress[0]), epsilon(&enp[0]);
    SymTensor<SymTensor<double,2>,2> *tm = (_D) ? new SymTensor<SymTensor<double,2>,2>() : NULL;

    Eigen::Block< Eigen::Map<Eigen::Matrix<doublereal,6,6> >,3,3>
      Dm = D.template topLeftCorner<3,3>(),     Dmb = D.template topRightCorner<3,3>(),
      Dbm = D.template bottomLeftCorner<3,3>(), Db = D.template bottomRightCorner<3,3>();
    Eigen::VectorBlock<Eigen::Map<Eigen::Matrix<doublereal,6,1> >,3>
      e = Upsilon.template head<3>(), chi = Upsilon.template tail<3>();
    Eigen::VectorBlock<Eigen::Map<Eigen::Matrix<doublereal,6,1> >,3>
      N = Sigma.template head<3>(), M = Sigma.template tail<3>();

// .....CLEAR THE CONSTITUTIVE MATRIX AND GENERALIZED STRESS VECTOR 

    if(_D) D.setZero();
    Sigma.setZero();

//     -------------------------------------------------- 
//       (NUMERICAL INTEGRATION THROUGH THE THICKNESS)    
//     -------------------------------------------------- 

    // 5 point Gauss-Legendre rule
    doublereal nodes[5] = { -0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664 };
    doublereal weights[5] = { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 };

    for (ilayer = 0; ilayer < nlayer; ++ilayer) {

// .....[z] COORDINATE AT THE THRU-THICKNESS GAUSS POINT

        z = nodes[ilayer]*h/2;

// .....COMPUTE THE LOCAL STRAINS [epsilon] = {epsilonxx,epsilonyy,gammaxy} ON THE SPECIFIED SURFACE

        epsilon = e + z*chi;
        epsilon[2] /= 2; // epsilonxy = gammaxy/2

// .....CALCULATE THE LOCAL STRESSES AND CONSISTENT TANGENT ELASTOPLASTIC MODULUS 

        double *stateni = (staten) ? staten+(nlayer*point+ilayer)*mat->getNumStates() : 0;
        double *statenpi = (statenp) ? statenp+(nlayer*point+ilayer)*mat->getNumStates() : 0;

        if(_D) {

            mat->integrate(&stress, tm, en, enp, stateni, statenpi, temp, NULL, dt);

            Ch << (*tm)[0][0], (*tm)[0][1], (*tm)[0][2],
                  (*tm)[1][0], (*tm)[1][1], (*tm)[1][2],
                  (*tm)[2][0], (*tm)[2][1], (*tm)[2][2];

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING 

            Db += weights[ilayer] * 1./2 * z * z * Ch;

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE 

            Dm += weights[ilayer] * 1./2 * Ch;

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE

            Dbm += weights[ilayer] * 1./2 * z * Ch.transpose();

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING

            Dmb += weights[ilayer] * 1./2 * z * Ch;

        }
        else {

          mat->integrate(&stress, en, enp, stateni, statenpi, temp, NULL, dt);

        }

// .....ASSEMBLE THE GENERALIZED "STRESSES"

        N += weights[ilayer] * 1./2 * sigmah;
        M += weights[ilayer] * 1./2 * z * sigmah;

    }

    if(_D) delete tm;
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType6<doublereal,localmaterial>
::GetLocalConstitutiveResponse(doublereal *_Upsilon, doublereal *sigma, doublereal z,
                               doublereal*, int nd, doublereal temp, doublereal dt,
                               doublereal *staten, doublereal *statenp)
{
    // Local variables 
    int ilayer;
    doublereal epszz;
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon);
    Eigen::VectorBlock< Eigen::Map<Eigen::Matrix<doublereal,6,1> > >
      e = Upsilon.head(3), chi = Upsilon.tail(3);

    SymTensor<double,2> en, enp, stress;
    Eigen::Map<Eigen::Matrix<doublereal,3,1> > epsilon(&enp[0]);

    if(z < 0) ilayer = 0;       // lower surface
    else if(z == 0) ilayer = 1; // median surface
    else ilayer = 2;            // upper surface

// .....COMPUTE THE LOCAL STRAINS [epsilon] = {epsilonxx,epsilonyy,gammaxy} ON THE SPECIFIED SURFACE

    epsilon = e + z * chi;
    epsilon[2] /= 2; // epsilonxy = gammaxy/2

// .....CALCULATE THE LOCAL STRESSES [sigma] = {sigmaxx,sigmayy,sigmaxy} ON THE SPECIFIED SURFACE

    double *stateni = (staten) ? staten+(nlayer*nd+ilayer)*mat->getNumStates() : 0;
    double *statenpi = (statenp) ? statenp+(nlayer*nd+ilayer)*mat->getNumStates() : 0;

    mat->integrate(&stress, en, enp, stateni, statenpi, temp, NULL, dt);

    sigma[0] = stress[0]/h; // xx
    sigma[1] = stress[1]/h; // yy
    sigma[2] = stress[2]/h; // xy
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType6<doublereal,localmaterial>
::UpdateState(doublereal *_Upsilon, doublereal *staten, doublereal *statenp, int point, doublereal temp, doublereal dt)
{
    // Local variables 
    int ilayer;
    doublereal z;
    SymTensor<double,2> en, enp, stress;
    Eigen::Map<Eigen::Matrix<doublereal,3,1> > epsilon(&enp[0]);
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon);

    Eigen::VectorBlock<Eigen::Map<Eigen::Matrix<doublereal,6,1> >,3>
      e = Upsilon.template head<3>(), chi = Upsilon.template tail<3>();

    // hardcoded 5 point Gauss-Legendre rule;
    doublereal nodes[5] = { -0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664 };

    for (ilayer = 0; ilayer < nlayer; ++ilayer) {

// .....[z] COORDINATE AT THE THRU-THICKNESS GAUSS POINT (OR SURFACE) AND STRAINS

        if(nlayer == 5) {
          z = nodes[ilayer]*h/2;
        }
        else {
           if(ilayer == 0) z = -h/2;     // lower surface
           else if(ilayer == 1) z = 0;   // median surface
           else z = h/2;                 // upper surface
        }

// .....COMPUTE THE LOCAL STRAINS [epsilon] = {epsilonxx,epsilonyy,gammaxy} ON THE SPECIFIED SURFACE

        epsilon = e + z * chi;
        epsilon[2] /= 2; // epsilonxy = gammaxy/2

// .....CALCULATE THE LOCAL STRESSES [sigma] = {sigmaxx,sigmayy,sigmaxy} AND UPDATE THE LOCAL STATE (INTERNAL VARIABLES)

        double *stateni = (staten) ? staten+ilayer*mat->getNumStates() : 0;
        double *statenpi = (statenp) ? statenp+ilayer*mat->getNumStates() : 0;

        mat->integrate(&stress, en, enp, stateni, statenpi, temp, NULL, dt);
    }
}

template<typename doublereal, typename localmaterial>
std::vector<doublereal>
ShellMaterialType6<doublereal,localmaterial>
::GetLocalPlasticStrain(int nd, doublereal z, doublereal *statenp)
{
  // return a copy of the 3-vector of plastic strain
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  double *statenpi = (statenp) ? statenp+(nlayer*nd+ilayer)*mat->getNumStates() : 0;
  Tensor_d0s2_Ss12 plasticstrain;

  if(mat->getPlasticStrain(statenpi, &plasticstrain)) {
    std::vector<doublereal> ret(3);
    ret[0] = plasticstrain[0]; // xx
    ret[1] = plasticstrain[3]; // yy
    ret[2] = plasticstrain[1]; // xy
    return ret;
  }
  else return std::vector<doublereal>();
}

template<typename doublereal, typename localmaterial>
std::vector<doublereal>
ShellMaterialType6<doublereal,localmaterial>
::GetLocalBackStress(int nd, doublereal z, doublereal *statenp)
{
  // return a copy of the 3-vector of back-stress
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  double *statenpi = (statenp) ? statenp+(nlayer*nd+ilayer)*mat->getNumStates() : 0;
  Tensor_d0s2_Ss12 backstress;
  
  if(mat->getBackStress(statenpi, &backstress)) {
    std::vector<doublereal> ret(3);
    ret[0] = backstress[0]; // xx
    ret[1] = backstress[3]; // yy
    ret[2] = backstress[1]; // xy
    return ret;
  }
  else return std::vector<doublereal>();
}

template<typename doublereal, typename localmaterial>
doublereal
ShellMaterialType6<doublereal,localmaterial>
::GetLocalEquivalentPlasticStrain(int nd, doublereal z, doublereal *statenp)
{
  // return the equivalent plastic strain
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  double *statenpi = (statenp) ? statenp+(nlayer*nd+ilayer)*mat->getNumStates() : 0;

  return mat->getEquivPlasticStrain(statenpi);
}

template<typename doublereal, typename localmaterial>
doublereal
ShellMaterialType6<doublereal,localmaterial>
::GetLocalDamage(int nd, doublereal z, doublereal *statenp)
{
  // return the scalar damage
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  double *statenpi = (statenp) ? statenp+(nlayer*nd+ilayer)*mat->getNumStates() : 0;

  return mat->getDamage(statenpi);
}

template<typename doublereal, typename localmaterial>
doublereal
ShellMaterialType6<doublereal,localmaterial>
::GetDissipatedEnergy(int point, doublereal *statenp)
{
    doublereal D = 0;
    int ilayer;

//     -------------------------------------------------- 
//       (NUMERICAL INTEGRATION THROUGH THE THICKNESS)    
//     -------------------------------------------------- 

    // 5 point Gauss-Legendre rule
    doublereal nodes[5] = { -0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664 };
    doublereal weights[5] = { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 };

    for (ilayer = 0; ilayer < nlayer; ++ilayer) {

      double *statenpi = (statenp) ? statenp+(nlayer*point+ilayer)*mat->getNumStates() : 0;

      D += weights[ilayer]*h/2*mat->getDissipatedEnergy(statenpi);

    }

    return D;
}

template<typename doublereal, typename localmaterial>
bool
ShellMaterialType6<doublereal,localmaterial>
::CheckFailure(doublereal *statenp)
{
  for(int i = 0; i < nlayer*maxgus; ++i) {
    if(mat->getDamage(statenp) < 1) return false;
    statenp += mat->getNumStates();
  }
  return true;
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType6<doublereal,localmaterial>
::GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D, doublereal *eframe, int gp)
{
  fprintf(stderr," *** ERROR: Stiffness w.r.t. displacement sensitivity output is \n"
                 "            not implemented for shell element types 15 and 1515 \n"
                 "            with a numerically enforced and integrated plane-   \n"
                 "            stress constitutive law.\n");
  exit(-1);
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType6<doublereal,localmaterial>
::GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z, doublereal *eframe, int gp)
{
  fprintf(stderr," *** ERROR: Local stress w.r.t. displacement sensitivity output \n"
                 "            is not implemented for shell element types 15 and   \n"
                 "            1515 with a numerically enforced and integrated     \n"
                 "            plane-stress constitutive law.                      \n");
  exit(-1);
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType6<doublereal,localmaterial>
::GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh, doublereal *, int, doublereal temp)
{
  fprintf(stderr," *** ERROR: Stiffness w.r.t. thickness sensitivity output is not\n"
                 "            implemented for shell element types 15 and 1515 with\n"
                 "            a numerically enforced and integrated plane-stress  \n"
                 "            constitutive law.                                   \n");
  exit(-1);
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType6<doublereal,localmaterial>
::GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh, doublereal *, int)
{
  fprintf(stderr," *** ERROR: Local stress w.r.t. displacement sensitivity output \n"
                 "            is not implemented for shell element types 15 and   \n"
                 "            1515 with a numerically enforced and integrated     \n"
                 "            plane-stress constitutive law.                      \n");
  exit(-1);
}

#include <Element.d/NonLinearity.d/NLMaterial.h>
template
void
ShellMaterialType6<double,NLMaterial>
::GetConstitutiveResponse(double *, double *, double *, double *, int, double, double, double *, double *);

template
void
ShellMaterialType6<double,NLMaterial>
::GetLocalConstitutiveResponse(double *, double *, double, double *, int, double, double, double *, double *);

template
void
ShellMaterialType6<double,NLMaterial>
::UpdateState(double *, double *, double *, int, double, double);

template
std::vector<double>
ShellMaterialType6<double,NLMaterial>
::GetLocalPlasticStrain(int, double, double *);

template
std::vector<double>
ShellMaterialType6<double,NLMaterial>
::GetLocalBackStress(int, double, double *);

template
double
ShellMaterialType6<double,NLMaterial>
::GetLocalEquivalentPlasticStrain(int, double, double *);

template
double
ShellMaterialType6<double,NLMaterial>
::GetLocalDamage(int, double, double *);

template
double
ShellMaterialType6<double,NLMaterial>
::GetDissipatedEnergy(int, double *);

template
bool
ShellMaterialType6<double,NLMaterial>
::CheckFailure(double *);

template
void
ShellMaterialType6<double,NLMaterial>
::GetConstitutiveResponseSensitivityWRTthic(double *, double *, double *, double *, int, double);

template
void
ShellMaterialType6<double,NLMaterial>
::GetConstitutiveResponseSensitivityWRTdisp(double *, double *, double *, double *, int);

template
void
ShellMaterialType6<double,NLMaterial>
::GetLocalConstitutiveResponseSensitivityWRTthic(double *, double *, double, double *, int);

template
void
ShellMaterialType6<double,NLMaterial>
::GetLocalConstitutiveResponseSensitivityWRTdisp(double *, double *, double, double *, int);
#endif
