#ifndef _SHELLMATERIALTYPE1_CPP_
#define _SHELLMATERIALTYPE1_CPP_

#ifdef USE_EIGEN3
#include <cmath>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

template<typename doublereal>
bool ShellMaterialType1<doublereal>::Wlocal_stress = true;
template<typename doublereal>
bool ShellMaterialType1<doublereal>::Wlocal_stress_disp = true;
template<typename doublereal>
bool ShellMaterialType1<doublereal>::Wlocal_stress_thic = true;

template<typename doublereal>
void
ShellMaterialType1<doublereal>::GetConstitutiveResponse(doublereal *_Upsilon, doublereal *_Sigma, doublereal *_D,
                                                        doublereal *eframe, int, doublereal temp, doublereal,
                                                        doublereal *, doublereal *)
{
  // Initialized data 
  doublereal zero = 0.;

  // Local variables 
  Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon), Sigma(_Sigma);
  doublereal *data = (_D == NULL) ? new doublereal[36] : _D;
  Eigen::Map<Eigen::Matrix<doublereal,6,6> > D(data);
  Eigen::Block< Eigen::Map<Eigen::Matrix<doublereal,6,6> > >
    Dm = D.topLeftCorner(3,3),     Dmb = D.topRightCorner(3,3),
    Dbm = D.bottomLeftCorner(3,3), Db = D.bottomRightCorner(3,3);
  Eigen::Matrix<doublereal,3,3> invT;

// ==================================================================== 
//                                                                      
//     Perform =   Assembles the 6 by 6 constitutive matrix and
//     ---------   computes the generalized stress resultants          
//                                                                      
//                                                                      
//     Input/Output =                                                   
//     --------------                                                   
//     Upsilon <input>  generalized strains
//     Sigma   <output> generalized stresses
//     D       <output> 3 by 3 constitutive matrix                      
//     eframe  <input>  element level 3x3 frame                         
//     temp    <input>  temperature
//                                                                      
//                                                                      
//     Computations =                                                   
//     --------------                                                   
//                             [ [D_m]   [D_mb] ]                       
//     [Constitutive_Matrix] = [                ]                       
//            6 by 6           [ [D_bm]  [D_b]  ]                       
//                                                                      
//     where "b" and "m" stand for bending and membrane, respectively:  
//                                                                      
//     The constitutive matrix [D_b] relates the element's curvatures  
//     to the bending moments and the constitutive matrix [D_m]        
//     relates the element's normal efforts to the strains. Similarly,  
//     the constitutive matrices [D_bm] and [D_mb] couple the bending   
//     and membrane effects:                                            
//                                                                      
//     [ M_x  ]            [ k_x  ]                                     
//     [ M_y  ] = [D_b]  * [ k_y  ]                                     
//     [ M_xy ]            [ k_xy ]                                     
//                                                                      
//     [ N_x  ]            [ e_x  ]                                     
//     [ N_y  ] = [D_m]  * [ e_y  ]                                     
//     [ N_xy ]            [ e_xy ]                                     
//                                                                      
//     [ M_x  ]            [ e_x  ]                                     
//     [ M_y  ] = [D_bm] * [ e_y  ]                                     
//     [ M_xy ]            [ e_xy ]                                     
//                                                                      
//     [ N_x  ]            [ k_x  ]                                     
//     [ N_y  ] = [D_mb] * [ k_y  ]                                     
//     [ N_xy ]            [ k_xy ]                                     
//                                                                      
//     The symmetry constraints are:                                    
//                                                                      
//     [D_b]^T  = [D_b]                                                
//     [D_m]^T  = [D_m]                                                
//     [D_bm]^T = [D_mb]                                                
//                                                                      
//     The assembly of these matrices is performed according to the     
//     type of constitutive law.     
//     In the following, the assembly of each one of the four matrices  
//     is briefly summarized.                                           
//                                                                      
//     2.  Constitutive Law of type-1:                                  
//     - - - - - - - - - - - - - - - -                                  
//     The coefficients [d_ij] for i=1 to 6 and j=1 to 6 are given in   
//     the output. They are stored in a vector of length 36 according   
//     to the following convention:                                     
//                                                                      
//     [  d_11  d_12  d_13  d_14  d_15  d_16  ]                         
//     [  d_12  d_22  d_23  d_24  d_25  d_26  ]                         
//     [  d_13  d_23  d_33  d_34  d_35  d_36  ]                         
//     [  d_14  d_24  d_33  d_44  d_45  d_46  ]                         
//     [  d_15  d_25  d_33  d_44  d_55  d_56  ]                         
//     [  d_16  d_26  d_33  d_44  d_55  d_66  ]                         
//                                                                      
//     is stored in the vector [coef] of size 36 at the following       
//     location:                                                        
//                                                                      
//     [   01    02    03    04    05    06   ]                         
//     [   07    08    09    10    11    12   ]                         
//     [   13    14    15    16    17    18   ]                         
//     [   19    20    21    22    23    24   ]                         
//     [   25    26    27    28    29    30   ]                         
//     [   31    32    33    34    35    36   ]                         
//                                                                      
//     Therefore, the entry on the ith row and jth column is stored at  
//     position number 6*(i-1)+j in the vector [coef].                  
//                                                                      
//     2.1  Pure Bending:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//              [ d_11 d_12 d_13 ]                    [ 22  23  24 ]    
//     [D_b]  = [ d_12 d_22 d_23 ] found at positions [ 28  29  30 ]    
//              [ d_13 d_23 d_33 ]                    [ 34  35  36 ]    
//                                                                      
//     2.2 Pure Membrane:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//              [ d_44 d_42 d_43 ]                    [ 01  02  03 ]    
//     [D_m]  = [ d_45 d_55 d_56 ] found at positions [ 07  08  09 ]    
//              [ d_43 d_56 d_66 ]                    [ 13  14  15 ]    
//                                                                      
//     2.3 Coupling Bending-Membrane ("BM"):                            
//     - - - - - - - - - - - - - - - - - - -                            
//                                                                      
//              [ d_14 d_15 d_16 ]                    [ 19  20  21 ]    
//     [D_bm] = [ d_24 d_25 d_26 ] found at positions [ 25  26  27 ]    
//              [ d_34 d_35 d_36 ]                    [ 31  32  33 ]    
//                                                                      
//     2.4 Coupling Membrane-Bending ("MB"):                            
//     - - - - - - - - - - - - - - - - - - -                            
//                                                                      
//              [ d_41 d_42 d_43 ]                    [ 04  05  06 ]    
//     [D_mb] = [ d_51 d_52 d_53 ] found at positions [ 10  11  12 ]    
//              [ d_61 d_62 d_63 ]                    [ 16  17  18 ]    
//                                                                      
//                                                                      
//     Caution =   It is assumed that the element has a constant        
//     ---------   thickness so that no numerical interpolation is      
//                 required. It is also assumed that the symmetry of    
//                 the [D] matrix has been checked when its 36 entries  
//                 are inputed.                                         
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9th, 1995                                             
// Version = 2.0                                                        
// Comment =                                                            
// ==================================================================== 

// .....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES 
// .....THE STRAINS [ex], [ey] AND [exy] TO THE STRESSES [sx], [sy] 
// .....AND [sxy] IN THE TRIANGULAR COORDINATE SYSTEM {x;y}: 
// .....[D] = [invT] * [D'] * [invT]^t 

    invT = this->andesinvt(eframe, aframe, zero);

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING 

    Db = invT * coef.bottomRightCorner(3,3) * invT.transpose();

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE

    Dm = invT * coef.topLeftCorner(3,3) * invT.transpose();

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE

    Dbm = -invT * coef.bottomLeftCorner(3,3) * invT.transpose();

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING

    Dmb = -invT * coef.topRightCorner(3,3) * invT.transpose();

// .....COMPUTE THE GENERALIZED "STRESSES"
//
    Sigma = D*(Upsilon - (temp-Ta)*Alpha);


    if(_D == NULL) delete [] data;
}

template<typename doublereal>
void
ShellMaterialType1<doublereal>::GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                                          doublereal *eframe, int gp)
{
    for(int i=0; i<18; ++i) { GetConstitutiveResponse(dUpsilondu, dSigmadu, D, eframe, gp, Ta); dUpsilondu += 6; dSigmadu += 6; }
}

template<typename doublereal>
void
ShellMaterialType1<doublereal>::GetConstitutiveResponseSensitivityWRTthic(doublereal *, doublereal *_dSigmadh, doublereal *_dDdh,
                                                                          doublereal *, int, doublereal)
{
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > dSigmadh(_dSigmadh);
    dSigmadh.setZero();

    if(_dDdh) {
      Eigen::Map<Eigen::Matrix<doublereal,6,6> > D(_dDdh);
      D.setZero();
    }
}

template<typename doublereal>
doublereal
ShellMaterialType1<doublereal>::GetShellThickness()
{
// .....INITIALIZE THE THICKNESS FOR A TYPE-1 CONSTITUTIVE LAW 
// .....TAKE THE DEFAULT THICKNESS ASSUMED CONSTANT 
// .....IF ZERO, ESTIMATE THE CONSTANT THICKNESS USING THE 
// .....FIRST COEFFICIENTS OF EXTENTIONAL AND BENDING STIFFNESSES 

  using std::sqrt;

  if (h == 0.) {
    if (coef(0,0) == 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellMaterialType1::getShellThickness ***\n"
          "*** The first coefficient of extentional stiffness is    ***\n"
          "*** equal to zero. Cannot estimate the thickness.        ***\n");
    }
    doublereal appxh2 = 3 * coef(3,3) / coef(0,0);
    if (appxh2 <= 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellMaterialType1::getShellThickness ***\n"
          "*** The ratio Between the first coefficient of bending   ***\n"
          "*** stiffness and the first coefficient of extensional   ***\n"
          "*** is negative or zero. Cannot take the square root and ***\n"
          "*** estimate the shell thickness.                        ***\n");
    }
    h = sqrt(appxh2);
  }

  return h;
}

template<typename doublereal>
void
ShellMaterialType1<doublereal>::GetLocalConstitutiveResponse(doublereal *, doublereal *sigma, doublereal,
                                                             doublereal *, int, doublereal, doublereal,
                                                             doublereal *, doublereal *)
{
  sigma[0] = sigma[1] = sigma[2] = 0.0;
  if(Wlocal_stress) {
    fprintf(stderr," *** WARNING: Local stress output is not available for shell      \n"
                   "              element types 15 and 1515 with a COEF-type composite\n"
                   "              constitutive law.                                   \n");
    Wlocal_stress = false;
  }
}

template<typename doublereal>
void
ShellMaterialType1<doublereal>::GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *, doublereal *dsigmadu, doublereal,
                                                                               doublereal *, int)
{
  for(int i=0; i<3*18; ++i) dsigmadu[i] = 0;
  if(Wlocal_stress_disp) {
    fprintf(stderr," *** WARNING: Local stress w.r.t. displacement sensitivity output \n"
                   "              is not available for shell element types 15 and 1515\n"
                   "              with a COEF-type composite constitutive law.        \n");
    Wlocal_stress_disp = false;
  }
}

template<typename doublereal>
void
ShellMaterialType1<doublereal>::GetLocalConstitutiveResponseSensitivityWRTthic(doublereal *, doublereal *dsigmadh,
                                                                               doublereal, doublereal *, int)
{
  dsigmadh[0] = dsigmadh[1] = dsigmadh[2] = 0;
  if(Wlocal_stress_thic) {
    fprintf(stderr," *** WARNING: Local stress w.r.t. displacement sensitivity output \n"
                   "              is not available for shell element types 15 and 1515\n"
                   "              15/1515 with a COEF-type composite constitutive law.\n");
    Wlocal_stress_thic = false;
  }
}

template
double* 
ShellMaterialType1<double>
::GetCoefOfConstitutiveLaw();

template
void
ShellMaterialType1<double>
::GetConstitutiveResponse(double *, double *, double *, double *, int, double, double, double *, double *);

template
double 
ShellMaterialType1<double>
::GetShellThickness();

template
void
ShellMaterialType1<double>
::GetLocalConstitutiveResponse(double *, double *, double, double *, int, double, double, double *, double *);

template
void
ShellMaterialType1<double>
::GetConstitutiveResponseSensitivityWRTthic(double *, double *, double *, double *, int, double);

template
void
ShellMaterialType1<double>
::GetConstitutiveResponseSensitivityWRTdisp(double *, double *, double *, double *, int);

template
void
ShellMaterialType1<double>
::GetLocalConstitutiveResponseSensitivityWRTthic(double *, double *, double, double *, int);

template
void
ShellMaterialType1<double>
::GetLocalConstitutiveResponseSensitivityWRTdisp(double *, double *, double, double *, int);
#endif

#endif
