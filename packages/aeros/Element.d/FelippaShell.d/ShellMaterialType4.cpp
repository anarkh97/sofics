#ifdef USE_EIGEN3
#include <cmath>
#include <stdexcept>
#include <vector>
#include <Eigen/Core>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <iostream>

extern int quietFlag;

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::GetConstitutiveResponse(doublereal *_Upsilon, doublereal *_Sigma, doublereal *_D,
                          doublereal *, int point, doublereal temp, doublereal dt,
                          doublereal *, doublereal *)
{
  // Local variables 
  int ilayer;
  doublereal z;
  Eigen::Matrix<doublereal,3,3> C;
  Eigen::Matrix<doublereal,3,1> sigma, epsilon;
  Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon), Sigma(_Sigma);
  Eigen::Map<Eigen::Matrix<doublereal,6,6> > D(_D);
  std::vector<doublereal> F(9);
  std::vector<doublereal> _sigma(9);
  std::vector<doublereal> *__C = (_D) ? new std::vector<doublereal>(9) : NULL;
  bool UpdateFlag = (_D) ? false : true; // by convention

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
    F[2] = F[5] = F[6] = F[7] = 0;

//     -------------------------------------------------- 
//       (NUMERICAL INTEGRATION THROUGH THE THICKNESS)    
//     -------------------------------------------------- 

    // 5 point Gauss-Legendre rule
    doublereal nodes[5] = { -0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664 };
    doublereal weights[5] = { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 };

    for (ilayer = 0; ilayer < nlayer; ++ilayer) {

// .....[z] COORDINATE AT THE THRU-THICKNESS GAUSS POINT AND STRAINS

        z = nodes[ilayer]*h/2;

        epsilon = e + z*chi;

// .....CALCULATE THE LOCAL STRESSES AND CONSISTENT TANGENT ELASTOPLASTIC MODULUS 

        F[0] = 1.0 + epsilon[0]; // xx
        F[1] = 0.5 * epsilon[2]; // xy
        F[3] = 0.5 * epsilon[2]; // yx
        F[4] = 1.0 + epsilon[1]; // yy
        F[8] = 1.0 - nu/(1-nu)*(epsilon[0]+epsilon[1]); // zz

        if(! mat[nlayer*point+ilayer]->ComputeElastoPlasticConstitutiveResponse(F, &_sigma, __C, UpdateFlag, dt) )
          throw std::runtime_error("ShellMaterialType4::GetConstitutiveResponse failed\n");

        sigma << _sigma[0], _sigma[4], _sigma[1];

        if(_D) {

            C << (*__C)[0], (*__C)[1], (*__C)[2],
                 (*__C)[3], (*__C)[4], (*__C)[5],
                 (*__C)[6], (*__C)[7], (*__C)[8];

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING 

            Db += weights[ilayer] * h/2 * z * z * C;

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE 

            Dm += weights[ilayer]  * h/2 * C;

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
            Dbm += weights[ilayer]  * h/2 * z * C.transpose();

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
            Dmb += weights[ilayer] * h/2 * z * C;

        }

// .....ASSEMBLE THE GENERALIZED "STRESSES"

        N += weights[ilayer] * h/2 * sigma;
        M += weights[ilayer] * h/2 * z * sigma;

    }

    if(_D) delete __C;
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::GetLocalConstitutiveResponse(doublereal *_Upsilon, doublereal *sigma, doublereal z,
                               doublereal*, int nd, doublereal temp, doublereal dt,
                               doublereal*, doublereal *)
{
    // Local variables 
    int ilayer;
    doublereal epszz;
    Eigen::Matrix<doublereal,3,1> epsilon;
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon);
    Eigen::VectorBlock< Eigen::Map<Eigen::Matrix<doublereal,6,1> > >
      e = Upsilon.head(3), chi = Upsilon.tail(3);

    if(z < 0) ilayer = 0;       // lower surface
    else if(z == 0) ilayer = 1; // median surface
    else ilayer = 2;            // upper surface

// .....COMPUTE THE LOCAL STRAINS [epsilon] = {epsilonxx,epsilonyy,gammaxy} ON THE SPECIFIED SURFACE

    epsilon = e + z * chi;

// .....CALCULATE THE LOCAL STRESSES [sigma] = {sigmaxx,sigmayy,sigmaxy} ON THE SPECIFIED SURFACE

    std::vector<doublereal> F(9);
    std::vector<doublereal> _sigma(9);
    F[0] = 1.0 + epsilon[0]; // xx
    F[1] = 0.5 * epsilon[2]; // xy
    F[2] = 0.0;              // xz
    F[3] = 0.5 * epsilon[2]; // yx
    F[4] = 1.0 + epsilon[1]; // yy
    F[5] = 0.0;              // yz
    F[6] = 0.0;              // zx
    F[7] = 0.0;              // zy
    F[8] = 1.0 - nu/(1-nu)*(epsilon[0]+epsilon[1]); // zz

    mat[nlayer*nd+ilayer]->ComputeElastoPlasticConstitutiveResponse(F, &_sigma, NULL, true, dt);

    sigma[0] = _sigma[0]; // xx
    sigma[1] = _sigma[4]; // yy
    sigma[2] = _sigma[1]; // xy
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::SetState(doublereal *state)
{
  std::vector<doublereal> PlasticStrain;
  std::vector<doublereal> BackStress;
  int l = 0;
  for(int i = 0; i < maxgus; ++i) {
    for(int j = 0; j < nlayer; ++j) {
      // set the internal variables
      // these should be the converged values at the beginning of the time step  
      PlasticStrain.clear();
      for (int k = 0; k < 3; ++k) PlasticStrain.push_back(state[l++]);
      mat[nlayer*i+j]->SetMaterialPlasticStrain(PlasticStrain);
      BackStress.clear();
      for (int k = 0; k < 3; ++k) BackStress.push_back(state[l++]);
      mat[nlayer*i+j]->SetMaterialBackStress(BackStress);
      mat[nlayer*i+j]->SetMaterialEquivalentPlasticStrain(state[l++]);
    }
  }
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::GetState(doublereal *state)
{
  int l = 0;
  for(int i = 0; i < maxgus; ++i) {
    for(int j = 0; j < nlayer; ++j) {
      // get the internal variables
      const std::vector<doublereal> &PlasticStrain = mat[nlayer*i+j]->GetMaterialPlasticStrain();
      state[l++] = PlasticStrain[0];
      state[l++] = PlasticStrain[1];
      state[l++] = PlasticStrain[2];
      const std::vector<doublereal> &BackStress = mat[nlayer*i+j]->GetMaterialBackStress();
      state[l++] = BackStress[0];
      state[l++] = BackStress[1];
      state[l++] = BackStress[2];
      state[l++] = mat[nlayer*i+j]->GetMaterialEquivalentPlasticStrain();
    }
  }
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::UpdateState(doublereal *_Upsilon, doublereal *, doublereal *state, int point, doublereal temp, doublereal dt)
{
  // Local variables 
  int ilayer;
  doublereal z;
  Eigen::Matrix<doublereal,3,1> epsilon;
  Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon);
  std::vector<doublereal> F(9);
  std::vector<doublereal> sigma(9);
  std::vector<doublereal> PlasticStrain;
  std::vector<doublereal> BackStress;

  Eigen::VectorBlock<Eigen::Map<Eigen::Matrix<doublereal,6,1> >,3>
    e = Upsilon.template head<3>(), chi = Upsilon.template tail<3>();

    F[2] = F[5] = F[6] = F[7] = 0;

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

        epsilon = e + z*chi;

// .....CALCULATE THE LOCAL STRESSES AND UPDATE THE LOCAL STATE (INTERNAL VARIABLES)

        F[0] = 1.0 + epsilon[0]; // xx
        F[1] = 0.5 * epsilon[2]; // xy
        F[3] = 0.5 * epsilon[2]; // yx
        F[4] = 1.0 + epsilon[1]; // yy
        F[8] = 1.0 - nu/(1-nu)*(epsilon[0]+epsilon[1]); // zz

        mat[nlayer*point+ilayer]->ComputeElastoPlasticConstitutiveResponse(F, &sigma, NULL, true, dt);

// .....MAKE A COPY OF THE LOCAL STATE (INTERNAL VARAIBLES)

        PlasticStrain = mat[nlayer*point+ilayer]->GetMaterialPlasticStrain();
        for (int j = 0; j < 3; ++j) state[7*ilayer+j] = PlasticStrain[j];
        BackStress = mat[nlayer*point+ilayer]->GetMaterialBackStress();
        for (int j = 0; j < 3; ++j) state[7*ilayer+3+j] = BackStress[j];
        state[7*ilayer+6] = mat[nlayer*point+ilayer]->GetMaterialEquivalentPlasticStrain();
    }
}

template<typename doublereal, typename localmaterial>
std::vector<doublereal>
ShellMaterialType4<doublereal,localmaterial>
::GetLocalPlasticStrain(int nd, doublereal z, doublereal *)
{
  // return a copy of the 3-vector of plastic strain
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  return mat[nlayer*nd+ilayer]->GetMaterialPlasticStrain();
}

template<typename doublereal, typename localmaterial>
std::vector<doublereal>
ShellMaterialType4<doublereal,localmaterial>
::GetLocalBackStress(int nd, doublereal z, doublereal *)
{
  // return a copy of the 3-vector of back-stress
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  return mat[nlayer*nd+ilayer]->GetMaterialBackStress();
}

template<typename doublereal, typename localmaterial>
doublereal
ShellMaterialType4<doublereal,localmaterial>
::GetLocalEquivalentPlasticStrain(int nd, doublereal z, doublereal *)
{
  // return the equivalent plastic strain
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  return mat[nlayer*nd+ilayer]->GetMaterialEquivalentPlasticStrain();
}

template<typename doublereal, typename localmaterial>
doublereal
ShellMaterialType4<doublereal,localmaterial>
::GetLocalDamage(int nd, doublereal z, doublereal *)
{
  // return the scalar damage
  // at node nd, and layer determined by z, as follows:
  int ilayer;
  if(z < 0) ilayer = 0;       // lower surface
  else if(z == 0) ilayer = 1; // median surface
  else ilayer = 2;            // upper surface

  return (mat[nlayer*nd+ilayer]->GetMaterialEquivalentPlasticStrain() >= mat[nlayer*nd+ilayer]->GetEquivalentPlasticStrainAtFailure()) ? 1 : 0;
}

template<typename doublereal, typename localmaterial>
doublereal
ShellMaterialType4<doublereal,localmaterial>
::GetDissipatedEnergy(int point, doublereal *)
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

      D += weights[ilayer]*h/2*mat[nlayer*point+ilayer]->GetDissipatedEnergy();

    }

    return D;
}

template<typename doublereal, typename localmaterial>
bool
ShellMaterialType4<doublereal,localmaterial>
::CheckFailure(doublereal *)
{
  for(int i = 0; i < nlayer*maxgus; ++i) {
    if(mat[i]->GetMaterialEquivalentPlasticStrain() < mat[i]->GetEquivalentPlasticStrainAtFailure()) return false;
  }
  return true;
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D, doublereal *eframe, int gp)
{
  fprintf(stderr," *** ERROR: Stiffness w.r.t. displacement sensitivity output is \n"
                 "            not implemented for shell element types 15 and 1515 \n"
                 "            with an elasto-plastic constitutive law.            \n");
  exit(-1);
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z, doublereal *eframe, int gp)
{
  fprintf(stderr," *** ERROR: Local stress w.r.t. displacement sensitivity output \n"
                 "            is not implemented for shell element types 15 and   \n"
                 "            1515 with an elasto-plastic constitutive law.       \n");
  exit(-1);
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::GetConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dSigmadh, doublereal *dDdh, doublereal *, int, doublereal temp)
{
  fprintf(stderr," *** ERROR: Stiffness w.r.t. thickness sensitivity output is not\n"
                 "            implemented for shell element types 15 and 1515 with\n"
                 "            an elasto-plastic constitutive law.                 \n");
  exit(-1);
}

template<typename doublereal, typename localmaterial>
void
ShellMaterialType4<doublereal,localmaterial>
::GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *Upsilon, doublereal *dsigmadh, doublereal dzdh, doublereal *, int)
{
  fprintf(stderr," *** ERROR: Local stress w.r.t. displacement sensitivity output \n"
                 "            is not implemented for shell element types 15 and   \n"
                 "            1515 with an elasto-plastic constitutive law.       \n");
  exit(-1);
}

#include <Material.d/IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h>
template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetConstitutiveResponse(double *, double *, double *, double *, int, double, double, double *, double *);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetLocalConstitutiveResponse(double *, double *, double, double *, int, double, double, double *, double *);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::SetState(double *);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetState(double *);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::UpdateState(double *, double *, double *, int, double, double);

template
std::vector<double>
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetLocalPlasticStrain(int, double, double *);

template
std::vector<double>
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetLocalBackStress(int, double, double *);

template
double
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetLocalEquivalentPlasticStrain(int, double, double *);

template
double
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetLocalDamage(int, double, double *);

template
double
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetDissipatedEnergy(int, double *);

template
bool
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::CheckFailure(double *);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetConstitutiveResponseSensitivityWRTthic(double *, double *, double *, double *, int, double);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetConstitutiveResponseSensitivityWRTdisp(double *, double *, double *, double *, int);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetLocalConstitutiveResponseSensitivityWRTthic(double *, double *, double, double *, int);

template
void
ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>
::GetLocalConstitutiveResponseSensitivityWRTdisp(double *, double *, double, double *, int);
#endif
