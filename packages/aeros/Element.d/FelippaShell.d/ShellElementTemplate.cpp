#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTTEMPLATE_CPP_
#define _SHELLELEMENTTEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <Utils.d/SolverInfo.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

extern SolverInfo &solInfo;

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andescrd(int elm, doublereal *_x, doublereal *_y, doublereal *_z,
           doublereal *_eframe, doublereal *_xlp, 
           doublereal *_ylp, doublereal *_zlp, doublereal &area)
{
  // Builtin functions 
  using std::sqrt;
  using std::abs;

  // Local variables 
  int i, j;
  doublereal x21, y21, z21, x13, y13, z13, x32, y32, z32, signedarea,
             projection, side21length, side32length;
  Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y), z(_z), xlp(_xlp), ylp(_ylp), zlp(_zlp);
  Eigen::Matrix<doublereal,3,3> xyz; xyz << x[0], x[1], x[2],
                                            y[0], y[1], y[2],
                                            z[0], z[1], z[2];
  Eigen::Matrix<doublereal,3,1> side21, side13, side32, cg, xp, yp, zp;
  Eigen::Map<Eigen::Matrix<doublereal,3,3> > eframe(_eframe);

// ==================================================================== 
//                                                                      
//     Perform =    This subroutine computes basic quantities needed    
//     ---------    for the assembly of the basic and higher order      
//                  stiffness matrices.                       
//                                                                      
//                                                                      
//     Inputs/Outputs =                                                 
//     ----------------                                                 
//     ELM    <input>   finite element number                           
//     X      <input>   nodal coordinates in the X-direction            
//     Y      <input>   nodal coordinates in the Y-direction            
//     Z      <input>   nodal coordinates in the Z-direction            
//     EFRAME <output>  element frame                                   
//     XLP    <output>  triangular coordinates in the X-direction       
//     YLP    <output>  triangular coordinates in the Y-direction       
//     ZLP    <output>  triangular coordinates in the Z-direction       
//     AREA   <output>  element area
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9th, 1994                                             
// Version = 1.0                                                        
// Comment =                                                            
// ==================================================================== 

// .....COMPUTE THE NODAL COORDINATE DIFFERENCES 

    side21 = xyz.col(1)-xyz.col(0);
    side13 = xyz.col(0)-xyz.col(2);
    side32 = xyz.col(2)-xyz.col(1);

// .....COMPUTE THE LENGTH OF SIDE 2-1 

    side21length = side21.norm();

// .....CHECK IF LENGTH 2-1 IS DIFFERENT FROM ZERO 

    if (side21length == 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andescrd ***\n"
          "*** Side between nodes 1 and 2 Has zero length    ***\n"
          "*** Check coordinates and FE topology             ***\n");
    }

// .....COMPUTE THE LENGTH OF SIDE 3-2 

    side32length = side32.norm();

// .....CHECK IF LENGTH 3-2 IS DIFFERENT FROM ZERO 

    if (side32length == 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andescrd ***\n"
          "*** Side between nodes 2 and 3 Has zero length    ***\n"
          "*** Check coordinates and FE topology             ***\n");
    }

// .....COMPUTE THE DISTANCE OF THE OPPOSING NODE 3 TO SIDE 2-1 

    projection = abs(side21.dot(side32))/side21length;

// .....GET THE AREA OF THE TRIANGLE 

    signedarea = side32length * side32length - projection * projection;

    if (signedarea <= 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andescrd ***\n"
          "*** The area is negative or zero                  ***\n"
          "*** Check coordinates and FE topology             ***\n");
    }

    area = side21length * .5 * sqrt(signedarea);

// .....COMPUTE THE DIRECTION COSINES OF THE LOCAL SYSTEM 
// .....DIRECTION [X] IS DIRECTED PARALLEL TO THE SIDE 2-1 
// .....DIRECTION [Z] IS THE EXTERNAL NORMAL (COUNTERCLOCKWISE) 
// .....DIRECTION [Y] IS COMPUTED AS [Z] x [X] (TENSORIAL PRODUCT) 

    xp = side21.normalized();
    zp = (side21.cross(side32)).normalized();
    yp = (zp.cross(xp)).normalized();

// .....COMPUTE THE COORDINATES FOR THE CENTER OF GRAVITY 

    cg = xyz.rowwise().sum()/3;

// .....CONSTRUCT THE ELEMENT FRAME 

    eframe << xp, yp, zp;

// .....COMPUTE THE LOCAL COORDINATES 

    Eigen::Matrix<doublereal,3,3> lp = (xyz-cg.replicate(1,3)).transpose()*eframe;
    xlp = lp.col(0);
    ylp = lp.col(1);
    zlp = lp.col(2); 
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesfrm(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *aframe,
           doublereal *_cframe)
{
  // Local variables 
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area;

  Eigen::Matrix<doublereal,3,3> eframe;
  Eigen::Map<Eigen::Matrix<doublereal,3,3,Eigen::RowMajor> > cframe(_cframe);

// ================================================================== 
//                                                                    
//     Perform =    This subroutine will form the composite frame     
//     ---------    for the 3D 3-node ANDES-EFF shell element.        
//                                                                    
//                                                                    
//     Inputs/Outputs =                                               
//     ----------------                                               
//     ELM     <input>  finite element number                         
//     X       <input>  nodal coordinates in the X-direction          
//     Y       <input>  nodal coordinates in the Y-direction          
//     Z       <input>  nodal coordinates in the Z-direction          
//     AFRAME  <input>  reference composite frame                     
//     CFRAME  <output> composite frame (projected onto element)      
//                                                                    
// ================================================================== 
// Authors = Francois M. Hemez and Philip J. S. Avery                                       
// Date    = September 13, 2011                                             
// Version = 3.0                                                      
// ================================================================== 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, x, y, z, eframe.data(), xlp, ylp, zlp, area);

    // andesinvt returns the transformation from element to fibre coordinate system
    // note: eframe.transpose()*v transforms v from global to element coordinate system
    // return the matrix that transforms from global to fibre coordinate system, i.e.:
    cframe = ShellMaterial<doublereal>::andesinvt(eframe.data(), aframe, 0.)*eframe.transpose(); // XXX is this correct?
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesgf(int elm, doublereal *_x, doublereal *_y, doublereal *_z, doublereal *_gravityForce,
          doublereal *gamma, int gravflg, doublereal rhoh)
{
  doublereal grvfor[3];
  doublereal totmas;
  andesms(elm, _x, _y, _z, gamma, grvfor, totmas, rhoh);
  
  Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y), z(_z);
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > gravityForce(_gravityForce);
  grvfor[0] /= 3.0;
  grvfor[1] /= 3.0;
  grvfor[2] /= 3.0;

  doublereal mx[3],my[3],mz[3];
  int i;
  for(i=0; i<3; ++i) {
    mx[i]=0.0;
    my[i]=0.0;
    mz[i]=0.0;
  }

  // Lumped
  if(gravflg == 0) {

  } 
  // Consistent or lumped with fixed end moments.  Compute treating shell as 3 beams.
  else {

    Eigen::Matrix<doublereal,3,1> T1,T2,T3;
    // Vector 1 from Node 1->2
    T1[0] = x[1] - x[0];
    T1[1] = y[1] - y[0];
    T1[2] = z[1] - z[0];
    T1.normalize();
    // Vector 2 from Node 1->3
    T2[0] = x[2] - x[0];
    T2[1] = y[2] - y[0];
    T2[2] = z[2] - z[0];
    T2.normalize();
    // Local Z-axis as cross between V1 and V2
    T3 = T1.cross(T2);
    T3.normalize();

    int beam, beamnode[3][2];
    beamnode[0][0] = 0;
    beamnode[0][1] = 1;
    beamnode[1][0] = 0;
    beamnode[1][1] = 2;
    beamnode[2][0] = 1;
    beamnode[2][1] = 2;

    for(beam=0; beam<3; ++beam) {
      doublereal length, dx, dy, dz, localg[3];
      int n1, n2;
      n1 = beamnode[beam][0];
      n2 = beamnode[beam][1];
      dx = x[n2] - x[n1];
      dy = y[n2] - y[n1];
      dz = z[n2] - z[n1];
      length = sqrt(dx*dx + dy*dy + dz*dz);
      // Local X-axis from Node 1->2
      T1[0] = x[n2] - x[n1];
      T1[1] = y[n2] - y[n1];
      T1[2] = z[n2] - z[n1];
      T1.normalize();
      // Local Y-axis as cross between Z and X
      T2 = T3.cross(T1); 
      T2.normalize();

      for(i=0; i<3; ++i)
        localg[i] = 0.0;
      for(i=0; i<3; ++i) {
        localg[0] += T1[i]*grvfor[i];
        localg[1] += T2[i]*grvfor[i];
        localg[2] += T3[i]*grvfor[i];
      }
      doublereal lmy,lmz;
      if (gravflg == 2) { // consistent
        lmy = -localg[2]*length/12.0;
        lmz = localg[1]*length/12.0;
      }
      else { // lumped with fixed-end moments
        lmy = -localg[2]*length/16.0;
        lmz = localg[1]*length/16.0;
      }
      mx[n1] += ((T2[0]*lmy) + (T3[0]*lmz));
      my[n1] += ((T2[1]*lmy) + (T3[1]*lmz));
      mz[n1] += ((T2[2]*lmy) + (T3[2]*lmz));
      mx[n2] -= ((T2[0]*lmy) + (T3[0]*lmz));
      my[n2] -= ((T2[1]*lmy) + (T3[1]*lmz));
      mz[n2] -= ((T2[2]*lmy) + (T3[2]*lmz));
    }
  }

  // set gravity force
  gravityForce[0]  = grvfor[0];
  gravityForce[1]  = grvfor[1];
  gravityForce[2]  = grvfor[2];
  gravityForce[3]  = mx[0];
  gravityForce[4]  = my[0];
  gravityForce[5]  = mz[0];
  gravityForce[6]  = grvfor[0];
  gravityForce[7]  = grvfor[1];
  gravityForce[8]  = grvfor[2];
  gravityForce[9]  = mx[1];
  gravityForce[10] = my[1];
  gravityForce[11] = mz[1];
  gravityForce[12] = grvfor[0];
  gravityForce[13] = grvfor[1];
  gravityForce[14] = grvfor[2];
  gravityForce[15] = mx[2];
  gravityForce[16] = my[2];
  gravityForce[17] = mz[2];
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesmm(int elm, doublereal *x, doublereal *y, doublereal *z, 
          doublereal *_emass, doublereal rhoh, int mflg)
{
  // Reference:
  // "Parametrized variational principles in dynamics applied
  //  to the optimization of dynamic models of plates",
  // F. J. Brito Castro, C. Militello, C. A. Felippa
  // Computational Mechanics 20 (1997) 285-294
  // Notes: this is the LMM (lumped mass matrix) which is formed by row lumping of the BCMM
  //        either with or without mass augmentation

  // Builtin functions 
  using std::sqrt;
  using std::abs;

  // Local variables 
  int i, j, i1, i2, i3;
  doublereal twicearea2, x21, x13, y13, z13, x32, y32, z32, y21, z21,
    ix, iy, iz, rlb, bpr, rlr, area, dist[3], mass0, mass1, 
    mass2, mass3, alpha;

  Eigen::Map<Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> > emass(_emass); 

// ==================================================================== 
//                                                                      
//     Performs =   This subroutine will form the elemental mass        
//     ----------   matrix of the 3D 3-node ANDES-EFF shell element.      
//                  Lumping is assumed here.                            
//                                                                      
//                                                                      
//     Inputs/Outputs =                                                 
//     ----------------                                                 
//     X        <input>   nodal coordinates in the X-direction          
//     Y        <input>   nodal coordinates in the Y-direction          
//     Z        <input>   nodal coordinates in the Z-direction          
//                                                                      
//     Computations =                                                   
//     ---------------                                                  
//                                                                      
//     The lumped mass matrix [M] is equal to:                          
//                                                                      
//           [ mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0 ]   
//     [M] = [ 0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3]   
//                                                                      
//     with the following ordering of local degrees of freedom:         
//                                                                      
//                                [     U_x1 ]                          
//                                [     U_y1 ]                          
//                                [     U_z1 ]                          
//                                [ theta_x1 ]                          
//                                [ theta_y1 ]                          
//                                [ theta_z1 ]                          
//                                [     U_x2 ]                          
//                                [     U_y2 ]                          
//                                [     U_z2 ]                          
//     Ordering_of_Local_DOFs  =  [ theta_x2 ]                          
//                                [ theta_y2 ]                          
//                                [ theta_z2 ]                          
//                                [     U_x3 ]                          
//                                [     U_y3 ]                          
//                                [     U_z3 ]                          
//                                [ theta_x3 ]                          
//                                [ theta_y3 ]                          
//                                [ theta_z3 ]                          
//                                                                      
//     No rotation of local-to-global basis is implemented since the    
//     mass matrix [M] is formed of 3 by 3 blocks proportional to the   
//     identity. For mflg = 0, the lumping factors are equal to:                      
//                                                                      
//     [mt] = [rhoh] * [A]        /    3.0                              
//     [m1] = [rhoh] * [A] * [Ix] / 1260.0                              
//     [m2] = [rhoh] * [A] * [Iy] / 1260.0                              
//     [m3] = [rhoh] * [A] * [Iz] / 1260.0                              
//
//     For mflg = 1, the lumping factors are equal to:
//
//     [mt] = [rhoh] * [A]       /  3.0                              
//     [m1] = [rhoh] * [A] * [A] / 24.0                              
//     [m2] = [rhoh] * [A] * [A] / 24.0                              
//     [m3] = [rhoh] * [A] * [A] / 24.0 
//                                                                      
//     where:                                                           
//                                                                      
//     [rhoh]                area density                               
//     [A]                   area                                       
//     [Ix], [Iy] and [Iz]   equivalent pseudo-moments of inertia       
//                                                                      
//     Caution =                                                        
//     ---------                                                        
//     The finite element is assumed to have constant thickness so      
//     that no numerical interpolation is required.                     
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9, 1995                                               
// Version = 2.0                                                        
// ==================================================================== 

// .....CLEAR THE OUTPUT MASS MATRIX

    emass.setZero();

//     -------------------------- 
//     CHECKS AND INITIALIZATIONS 
//     -------------------------- 

// .....COMPUTE THE DISTANCE BETWEEN X-, Y- AND Z- NODAL COORDINATES

    x21 = x[1] - x[0];
    y21 = y[1] - y[0];
    z21 = z[1] - z[0];

    x32 = x[2] - x[1];
    y32 = y[2] - y[1];
    z32 = z[2] - z[1];

    x13 = x[0] - x[2];
    y13 = y[0] - y[2];
    z13 = z[0] - z[2];

// .....COMPUTE THE DISTANCE BETWEEN NODES 1-2, 2-3 AND 3-1

    dist[0] = sqrt(x21 * x21 + y21 * y21 + z21 * z21);
    dist[1] = sqrt(x32 * x32 + y32 * y32 + z32 * z32);
    dist[2] = sqrt(x13 * x13 + y13 * y13 + z13 * z13);

// .....COMPUTE THE LENGTH OF SIDE 1-2

    rlr = sqrt(x21 * x21 + y21 * y21 + z21 * z21);

// .....CHECK FOR ZERO-SIDE LENGTH

    if (rlr == 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andesmm ***\n"
          "*** Side between nodes 1 and 2 has zero length   ***\n"
          "*** Check coordinates and FE topology            ***\n");
    }

// .....COMPUTE THE DISTANCE OF THE OPPOSING NODE (3) TO THAT SIDE (1-2)

    rlb = sqrt(x32 * x32 + y32 * y32 + z32 * z32);
    bpr = abs(x21 * x32 + y21 * y32 + z21 * z32) / rlr;

// .....COMPUTE THE SQUARE OF TWICE THE TRIANGLE'S AREA

    twicearea2 = rlb * rlb - bpr * bpr;

// .....CHECK IF THE TRIANGLE'S AREA IS POSITIVE

    if (twicearea2 <= 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andesmm ***\n"
          "*** The area is negative or zero                 ***\n"
          "*** Check coordinates and FE topology            ***\n");
    }

// .....COMPUTE THE AREA OF THE TRIANGLE

    area = rlr * .5 * sqrt(twicearea2);

// .....COMPUTE THE THREE PSEUDO MOMENTS OF INERTIA

    ix = dist[0] * dist[0] + dist[2] * dist[2];
    iy = dist[0] * dist[0] + dist[1] * dist[1];
    iz = dist[1] * dist[1] + dist[2] * dist[2];

// .....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM

    mass0 = rhoh * area / 3.;
    if(mflg == 0) {
        mass1 = rhoh * area * ix / 1260.;
        mass2 = rhoh * area * iy / 1260.;
        mass3 = rhoh * area * iz / 1260.;
    }
    else {
        alpha = area / 8.;
        mass1 = mass0 * alpha;
        mass2 = mass0 * alpha;
        mass3 = mass0 * alpha;
    }
//     ------------------------------------- 
//     ASSEMBLY OF THE ELEMENTAL MASS MATRIX 
//     ------------------------------------- 

// .....FORM THE LUMPED ELEMENT MASS MATRIX 

    for (i = 0; i < 3; ++i) {
        i2 = i + 6;
        i3 = i + 12;
        emass(i, i) = mass0;
        emass(i2, i2) = mass0;
        emass(i3, i3) = mass0;
    }

    for (i = 0; i < 3; ++i) {
        i1 = i + 3;
        i2 = i + 9;
        i3 = i + 15;
        emass(i1, i1) = mass1;
        emass(i2, i2) = mass2;
        emass(i3, i3) = mass3;
    }

}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesms(int elm, doublereal *x, doublereal *y, doublereal *z, 
          doublereal *gamma, doublereal *grvfor, doublereal &totmas,
          doublereal rhoh)
{
  // Builtin functions 
  using std::sqrt;
  using std::abs;

  // Local variables 
  doublereal twicearea2, x21, x13, y13, z13, x32, y32, z32, y21, z21,
    rlb, bpr, rlr, area, dist[3], mass0;

// ==================================================================== 
//                                                                      
//     Performs =   This subroutine will compute the body force due to        
//     ----------   gravity and/or the total mass of the element.      
//                  Lumping is assumed here.                            
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9, 1995                                               
// Version = 2.0                                                        
// ==================================================================== 

//     -------------------------- 
//     CHECKS AND INITIALIZATIONS 
//     -------------------------- 

// .....COMPUTE THE DISTANCE BETWEEN X-, Y- AND Z- NODAL COORDINATES

    x21 = x[1] - x[0];
    y21 = y[1] - y[0];
    z21 = z[1] - z[0];

    x32 = x[2] - x[1];
    y32 = y[2] - y[1];
    z32 = z[2] - z[1];

    x13 = x[0] - x[2];
    y13 = y[0] - y[2];
    z13 = z[0] - z[2];

// .....COMPUTE THE DISTANCE BETWEEN NODES 1-2, 2-3 AND 3-1

    dist[0] = sqrt(x21 * x21 + y21 * y21 + z21 * z21);
    dist[1] = sqrt(x32 * x32 + y32 * y32 + z32 * z32);
    dist[2] = sqrt(x13 * x13 + y13 * y13 + z13 * z13);

// .....COMPUTE THE LENGTH OF SIDE 1-2

    rlr = sqrt(x21 * x21 + y21 * y21 + z21 * z21);

// .....CHECK FOR ZERO-SIDE LENGTH

    if (rlr == 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andesms ***\n"
          "*** Side between nodes 1 and 2 has zero length   ***\n"
          "*** Check coordinates and FE topology            ***\n");
    }

// .....COMPUTE THE DISTANCE OF THE OPPOSING NODE (3) TO THAT SIDE (1-2)

    rlb = sqrt(x32 * x32 + y32 * y32 + z32 * z32);
    bpr = abs(x21 * x32 + y21 * y32 + z21 * z32) / rlr;

// .....COMPUTE THE SQUARE OF TWICE THE TRIANGLE'S AREA

    twicearea2 = rlb * rlb - bpr * bpr;

// .....CHECK IF THE TRIANGLE'S AREA IS POSITIVE

    if (twicearea2 <= 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andesms ***\n"
          "*** The area is negative or zero                 ***\n"
          "*** Check coordinates and FE topology            ***\n");
    }

// .....COMPUTE THE AREA OF THE TRIANGLE

    area = rlr * .5 * sqrt(twicearea2);

// .....FORM THE MASS COEFFICIENT

    mass0 = rhoh * area / 3.;

// ..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED

    if (gamma) {
        grvfor[0] = mass0 * 3. * gamma[0];
        grvfor[1] = mass0 * 3. * gamma[1];
        grvfor[2] = mass0 * 3. * gamma[2];
    }

// .... COMPUTE THE TOTAL MASS

    totmas = mass0 * 3.;

}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesstf(int elm, doublereal *_estiff, doublereal *_fint, doublereal nu,
           doublereal *x, doublereal *y, doublereal *z, doublereal *_v,
           int ctyp, ShellMaterial<doublereal> *gpmat, int flag,
           int tflg, doublereal *ndtemps, doublereal dt, doublereal *staten,
           doublereal *statenp)
{
  // Initialized data 
  bool debug = false;
  doublereal clr = solInfo.andes_clr;
  doublereal cqr = solInfo.andes_cqr;
  doublereal betab = solInfo.andes_betab;
  doublereal alpha = solInfo.andes_alpha; // using 0 here seems to work better for thin shell
  doublereal betam = solInfo.andes_betam; // using 0 here seems to work better for very thin shell

  // Local variables 
  int i, j;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area;
  doublereal temp;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;

  Eigen::Map< Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> > K(_estiff);
  Eigen::Map< Eigen::Matrix<doublereal,18,1> > F(_fint);
  Eigen::Matrix<doublereal,6,1> Upsilon, Sigma;
  Eigen::Matrix<doublereal,6,6> D;
  Eigen::Matrix<doublereal,3,3> eframe;
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);

  // Some convenient definitions 
  Eigen::Block< Eigen::Map<Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> >,9,9 >
    Km = K.template topLeftCorner<9,9>(),     Kmb = K.template topRightCorner<9,9>(),
    Kbm = K.template bottomLeftCorner<9,9>(), Kb = K.template bottomRightCorner<9,9>();
  Eigen::VectorBlock< Eigen::Map<Eigen::Matrix<doublereal,18,1> >,9 >
    Fm = F.template head<9>(), Fb = F.template tail<9>();
  Eigen::Block< Eigen::Matrix<doublereal,6,6>,3,3 > 
    Dm = D.template topLeftCorner<3,3>(),     Dmb = D.template topRightCorner<3,3>(),
    Dbm = D.template bottomLeftCorner<3,3>(), Db = D.template bottomRightCorner<3,3>();
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1>,3 >
    e = Upsilon.template head<3>(), chi = Upsilon.template tail<3>();
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1>,3 >
    N = Sigma.template head<3>(), M = Sigma.template tail<3>();

// ================================================================== 
//                                                                    
//     Perform =    This subroutine will form the element stiffness   
//     ---------    matrix and/or element internal force vector for   
//                  the 3D 3-node ANDES-EFF shell element.            
//                                                                    
//                                                                    
//     Inputs/Outputs =                                               
//     ----------------                                               
//     ELM     <input>  finite element number                         
//     ESTIFF  <output> element stiffness matrix                      
//     NU      <input>  Poisson's ratio                               
//     X       <input>  nodal coordinates in the X-direction          
//     Y       <input>  nodal coordinates in the Y-direction          
//     Z       <input>  nodal coordinates in the Z-direction          
//     V       <input>  nodal generalized displacements               
//     CTYP    <input>  type of constitutive law              
//     FLAG    <input>  int specifying whether to return              
//                      transformed element stiffness matrix or       
//                      global element stiffness matrix               
//     NDTEMPS <input>  nodal temperatures                            
//     TFLG    <input>  int specifying whether to use constant        
//                      or linear interpolation for temperature       
//                                                                    
//                                                                    
//     Computations =                                                 
//     --------------                                                 
//     This subroutine evaluates the stiffness matrix for the 18      
//     degrees of freedom 3-node composite triangle obtained as a     
//     combination of the Assumed Quadratic Rotations bending         
//     triangle plus the EFF membrane with driling degrees of freedom     
//     developed by Militello, Felippa et al. For documentation, see     
//     "The first ANDES elements: 9-dof plate bending triangles",
//      Militello & Felippa, Comput. Methods Appl. Mech. Engrg,
//      Vol. 93 (1991) pp. 217-246
//     "Membrane triangles with corner drilling freedoms I. The EFF element"
//      Alvin, de la Fuente, Haugen & Felippa,
//      Finite Elements in Analysis and Design Vol. 12 (1992) pp. 163-187
//     
//                                                                    
//     The original version of the ANS shell element has been         
//     generalized here to the case of a composite shell element      
//     with complete bending-membrane coupling and/or layered
//     plastic shell with basic-higher order coupling. Five types of         
//     constitutive laws have been used:                    
//        type-0: isotropic linear                              
//        type-1: constitutive coefficients are given                 
//        type-2: properties of each layer are given and no coupling  
//                   is assumed between bending and membrane          
//        type-3: properties of each layer are given and coupling     
//                   between bending and membrane is assumed          
//        type-4: isotropic nonlinear
//                                                                    
//     The stiffness matrix [K] is assembled as the combination of    
//     the basic stiffness and the higher order stiffness matrices, see
//     "Plastic buckling and collapse of thin shell structures,
//      using layered plastic modeling and co-rotational ANDES
//      finite element", Dal Cortivo, Felippa, Bavestrello & Silva,
//      Comput. Methods Appl. Mech. Engrg, Vol. 198 (2009) pp. 785-798.
//                                                                    
//                                                                    
//     Caution =                                                      
//     ---------                                                      
//     The finite element is assumed to have a constant thickness     
//     so that no numerical interpolation is required.
//                                                                    
// ================================================================== 
// Authors = Francois M. Hemez and Philip J. S. Avery                                       
// Date    = September 13, 2011                                             
// Version = 3.0                                                      
// ================================================================== 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, x, y, z, eframe.data(), xlp, ylp, zlp, area);

// .....CONSTRUCT THE PERMUTATION MATRIX 

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> P(indices);

// .....ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL 
//      FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT) 
//      AND APPLY PERMUTATION to {M,B} ordering

    if(flag == 1) {

        for(i = 0; i < 18; i += 3)
            vd.segment(i,3) = eframe.transpose()*v.segment(i,3);

        vd = P.transpose()*vd;
    }
    else {

        // for Corotational Formulation, v is already local
        vd = P.transpose()*v; 

    }

    // Note: betam = 0.32 is max(1/2*(1-4*nu^2),0.01) assuming nu to be 0.3
    // Reference: "Membrane triangles with corner drilling freedoms III. Implementation and performance evaluation"
    //             Carlos A. Felippa and Scott Alexander
    // It's not clear what to do when nu is not available (e.g. for type 1 composite) or varies through the thickness
    // (e.g. for types 2 and 3 composite)

// .....FORM THE BASIC CURVATURE-TO-DISPLACEMENT MATRIX
 
    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

// .....FORM THE BASIC EXTENSION-TO-DISPLACEMENT MATRIX

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

    if(tflg == 0) {
        // if ndtemps (nodal temperatures) is not null, then the mean temperature is used at each gauss point.
        temp = (ndtemps) ? doublereal((ndtemps[0]+ndtemps[1]+ndtemps[2])/3) : gpmat->GetAmbientTemperature();
    }
    else {
        // if ndtemps (nodal temperatures) is not null, then the temperature at each gauss point is interpolated.
        if(!ndtemps) temp = gpmat->GetAmbientTemperature();
    }

    if(_estiff) K.setZero();
    if(_fint) F.setZero();
    doublereal zeta[3][3] = { { 0.,.5,.5 }, { .5,0.,.5 }, { .5,.5,0. } }; // triangular coordinates of gauss integration points
    doublereal weight[3] = { 1/3., 1/3., 1/3. };
    for(i = 0; i < 3; ++i) {

// .....FORM HIGHER ORDER INTEGRATED CURVATURE-TO-DISPLACEMENT MATRIX
//      AND THE ELEMENT CURVATURE VECTOR 

        Bb = 1/area*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
        chi = Bb*vd.tail(9);

// .....FORM THE HIGHER ORDER INTEGRATED EXTENSION-TO-DISPLACEMENT MATRIX
//      AND THE ELEMENT EXTENSION VECTOR

        Bm = 1/area*Lm.transpose() + Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
        e = Bm*vd.head(9);

// .....GET THE TANGENT CONSTITUTIVE MATRIX [D = {Dm,Dmb;Dbm,Db}] 
//      AND THE GENERALIZED STRESSES [Sigma = {N,M}]

        if(i == 0 || ctyp == 4 || (ndtemps && tflg != 0) || _fint) {
          if(ndtemps && tflg != 0) 
            temp = zeta[i][0]*ndtemps[0] + zeta[i][1]*ndtemps[1] + zeta[i][2]*ndtemps[2];
          doublereal *_D = (_estiff) ? D.data() : NULL;
          gpmat->GetConstitutiveResponse(Upsilon.data(), Sigma.data(), _D, eframe.data(), i, temp, dt,
                                         staten, statenp);
        }

        if(_estiff) {

// .....FORM STIFFNESS FOR PURE BENDING

            Kb.noalias() += (area*weight[i])*(Bb.transpose()*Db*Bb);

// .....FORM STIFFNESS FOR PURE MEMBRANE 

            Km.noalias() += (area*weight[i])*(Bm.transpose()*Dm*Bm);

            if(ctyp != 0 && ctyp != 2) {

// .....FORM STIFFNESS FOR BENDING-MEMBRANE COUPLING

              Kbm.noalias() += (area*weight[i])*(Bb.transpose()*Dbm*Bm);

// .....FORM STIFFNESS FOR MEMBRANE-BENDING COUPLING

              Kmb.noalias() += (area*weight[i])*(Bm.transpose()*Dmb*Bb);

            }

// .....CHECK THE CONSTITUTIVE MATRIX (USED FOR DEBUGGING ONLY)
#ifdef DEBUG_SHELL_ELEMENT_TEMPLATE
            if(debug && (i == 0 || ctyp == 4)) {
                std::cerr << "Here are the eigenvalues of the constitutive matrix (element " << elm 
                          << ", gauss point " << i << ") :\n"
                          << D.eigenvalues().transpose() << std::endl;
            }
#endif
        }

        if(_fint) {

// .....FORM THE INTERNAL FORCE FOR BENDING

            Fb.noalias() += (area*weight[i])*Bb.transpose()*M;

// .....FORM THE INTERNAL FORCE FOR MEMBRANE

            Fm.noalias() += (area*weight[i])*Bm.transpose()*N;

        }

    }

// .....APPLY PERMUTATION 

    if(_estiff) K = P*K*P.transpose();
    if(_fint)   F = P*F;

// .....ROTATE ELEMENT STIFFNESS AND/OR FORCE TO GLOBAL COORDINATES 
//      (only if flag equals 1)

    if (flag == 1 && _estiff) {
        for(i = 0; i < 18; i += 3)
            for(j = 0; j < 18; j += 3)
                 K.template block<3,3>(i,j) = eframe*K.template block<3,3>(i,j)*eframe.transpose();
    }

    if (flag == 1 && _fint) {
        for(i = 0; i < 18; i += 3)
            F.template segment<3>(i) = eframe*F.template segment<3>(i);
    }

// .....CHECK THE POSITIVITY OF THE OUTPUT STIFFNESS MATRIX 
// .....(USED FOR DEBUGGING ONLY) 
#ifdef DEBUG_SHELL_ELEMENT_TEMPLATE
    if(debug && _estiff) {
        std::cerr << "Here are the eigenvalues of the stiffness matrix (element " << elm << ") :\n"
                  << K.eigenvalues().transpose() << std::endl;
    }
#endif
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesstfWRTthic(int elm, doublereal *_destiffdh, doublereal *_dfintdh, doublereal nu,
                  doublereal *x, doublereal *y, doublereal *z, doublereal *_v,
                  int ctyp, ShellMaterial<doublereal> *gpmat, int flag,
                  int tflg, doublereal *ndtemps)
{
  // Initialized data 
  doublereal clr = solInfo.andes_clr;
  doublereal cqr = solInfo.andes_cqr;
  doublereal betab = solInfo.andes_betab;
  doublereal alpha = solInfo.andes_alpha;
  doublereal betam = solInfo.andes_betam;

  // Local variables 
  int i, j;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area;
  doublereal temp;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;

  Eigen::Map< Eigen::Matrix<doublereal,18,1> > dFdh(_dfintdh);
  Eigen::Map< Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> > dKdh(_destiffdh);
  Eigen::Matrix<doublereal,6,1> Upsilon, dSigmadh;
  Eigen::Matrix<doublereal,6,6> D;
  Eigen::Matrix<doublereal,3,3> eframe;
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);

  // Some convenient definitions 
  Eigen::Block< Eigen::Map<Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> >,9,9 >
    Km = dKdh.template topLeftCorner<9,9>(),     Kmb = dKdh.template topRightCorner<9,9>(),
    Kbm = dKdh.template bottomLeftCorner<9,9>(), Kb = dKdh.template bottomRightCorner<9,9>();
  Eigen::VectorBlock< Eigen::Map<Eigen::Matrix<doublereal,18,1> >,9 >
    Fm = dFdh.template head<9>(), Fb = dFdh.template tail<9>();
  Eigen::Block< Eigen::Matrix<doublereal,6,6>,3,3 >
    Dm = D.template topLeftCorner<3,3>(),     Dmb = D.template topRightCorner<3,3>(),
    Dbm = D.template bottomLeftCorner<3,3>(), Db = D.template bottomRightCorner<3,3>();
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1>,3 >
    e = Upsilon.template head<3>(), chi = Upsilon.template tail<3>();
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1>,3 >
    N = dSigmadh.template head<3>(), M = dSigmadh.template tail<3>();

// ================================================================== 
//                                                                    
//     Perform =    This subroutine will form the sensitivity of the  
//     ---------    element stiffness matrix and/or internal force
//                  vector wrt thickness for the 3D 3-node ANDES-EFF
//                  shell element.
//                                                                    
//                                                                    
//     Inputs/Outputs =                                               
//     ----------------                                               
//     ELM       <input>  finite element number                         
//     DESTIFFDH <output> sensitivity of element stiffness matrix wrt   
//                        thickness                                     
//     DFINTDH   <output> sensitivity of element internal force 
//                        vector wrt thickness
//     NU        <input>  Poisson's ratio                               
//     X         <input>  nodal coordinates in the X-direction          
//     Y         <input>  nodal coordinates in the Y-direction          
//     Z         <input>  nodal coordinates in the Z-direction          
//     V         <input>  nodal generalized displacements
//     CTYP      <input>  type of constitutive law                   
//     FLAG      <input>  int specifying whether to return              
//                        transformed element stiffness sensitivity     
//                        matrix or global element stiffness            
//                        sensitivity matrix                            
//     NDTEMPS   <input>  nodal temperatures                            
//     TFLG      <input>  int specifying whether to use constant        
//                        or linear interpolation for temperature       
//                                                                    
//                                                                    
// ================================================================== 
// Author  = Youngsoo Choi                                            
// Date    = January 17, 2014                                             
// Version = 1.0                                                      
// ================================================================== 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, x, y, z, eframe.data(), xlp, ylp, zlp, area);

// .....CONSTRUCT THE PERMUTATION MATRIX 

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> P(indices);

// .....ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL 
//      FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT) 
//      AND APPLY PERMUTATION to {M,B} ordering

    if(flag == 1) {

        for(i = 0; i < 18; i += 3)
            vd.segment(i,3) = eframe.transpose()*v.segment(i,3);

        vd = P.transpose()*vd;
    }
    else {

        // for Corotational Formulation, v is already local
        vd = P.transpose()*v;

    }

// .....FORM THE BASIC CURVATURE-TO-DISPLACEMENT MATRIX
 
    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

// .....FORM THE BASIC EXTENSION-TO-DISPLACEMENT MATRIX

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

    if(tflg == 0) {
        // if ndtemps (nodal temperatures) is not null, then the mean temperature is used at each gauss point.
        temp = (ndtemps) ? doublereal((ndtemps[0]+ndtemps[1]+ndtemps[2])/3) : gpmat->GetAmbientTemperature();
    }
    else {
        // if ndtemps (nodal temperatures) is not null, then the temperature at each gauss point is interpolated.
        if(!ndtemps) temp = gpmat->GetAmbientTemperature();
    }

    if(_destiffdh) dKdh.setZero();
    if(_dfintdh) dFdh.setZero();
    doublereal zeta[3][3] = { { 0.,.5,.5 }, { .5,0.,.5 }, { .5,.5,0. } }; // triangular coordinates of gauss integration points
    doublereal weight[3] = { 1/3., 1/3., 1/3. };
    for(i = 0; i < 3; ++i) {

// .....FORM HIGHER ORDER INTEGRATED CURVATURE-TO-DISPLACEMENT MATRIX
//      AND THE ELEMENT CURVATURE VECTOR 

        Bb = 1/area*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
        chi = Bb*vd.tail(9);

// .....FORM THE HIGHER ORDER INTEGRATED EXTENSION-TO-DISPLACEMENT MATRIX
//      AND THE ELEMENT EXTENSION VECTOR

        Bm = 1/area*Lm.transpose() + Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
        e = Bm*vd.head(9);

// .....GET THE TANGENT CONSTITUTIVE SENSITIVITY MATRIX [D = {Dm,Dmb;Dbm,Db}] 

        if(i == 0 || ctyp == 4 || (ndtemps && tflg != 0) || _dfintdh) {
          if(ndtemps && tflg != 0)
            temp = zeta[i][0]*ndtemps[0] + zeta[i][1]*ndtemps[1] + zeta[i][2]*ndtemps[2];
          doublereal *_D = (_destiffdh) ? D.data() : NULL;
          gpmat->GetConstitutiveResponseSensitivityWRTthic(Upsilon.data(), dSigmadh.data(), _D, eframe.data(), i, temp);
        }

        if(_destiffdh) {

// .....FORM STIFFNESS SENSITIVITY FOR PURE BENDING

            Kb.noalias() += (area*weight[i])*(Bb.transpose()*Db*Bb);

// .....FORM STIFFNESS SENSITIVITY FOR PURE MEMBRANE 

            Km.noalias() += (area*weight[i])*(Bm.transpose()*Dm*Bm);

            if(ctyp != 0 && ctyp != 2) {

// .....FORM STIFFNESS SENSITIVITY FOR BENDING-MEMBRANE COUPLING

                Kbm.noalias() += (area*weight[i])*(Bb.transpose()*Dbm*Bm);

// .....FORM STIFFNESS SENSITIVITY FOR MEMBRANE-BENDING COUPLING

                Kmb.noalias() += (area*weight[i])*(Bm.transpose()*Dmb*Bb);

            }
        }

        if(_dfintdh) {

// .....FORM THE INTERNAL FORCE SENSITIVITY FOR BENDING

            Fb.noalias() += (area*weight[i])*Bb.transpose()*M;

// .....FORM THE INTERNAL FORCE SENSITIVITY FOR MEMBRANE

            Fm.noalias() += (area*weight[i])*Bm.transpose()*N;

        }

    }

// .....APPLY PERMUTATION 

    if(_destiffdh) dKdh = P*dKdh*P.transpose();
    if(_dfintdh)   dFdh = P*dFdh;

// .....ROTATE ELEMENT STIFFNESS AND/OR FORCE SENSITIVITIES TO GLOBAL COORDINATES 
//      (only if flag equals 1)

    if (flag == 1 && _destiffdh) {
        for(i = 0; i < 18; i += 3)
            for(j = 0; j < 18; j += 3)
                 dKdh.template block<3,3>(i,j) = eframe*dKdh.template block<3,3>(i,j)*eframe.transpose();
    }

    if (flag == 1 && _dfintdh) {
        for(i = 0; i < 18; i += 3)
            dFdh.template segment<3>(i) = eframe*dFdh.template segment<3>(i);
    }
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesvms(int elm, int maxstr, doublereal nu, 
           doublereal *X, doublereal *Y, doublereal *Z,
           doublereal *_v, doublereal *_stress,
           int ctyp, ShellMaterial<doublereal> *nmat, 
           int strainflg, int surface, int sflg,
           doublereal *ndtemps, int flag, doublereal *staten,
           doublereal *statenp)
{
  // Initialized data 
  doublereal clr = solInfo.andes_clr;
  doublereal cqr = solInfo.andes_cqr;
  doublereal betab = solInfo.andes_betab;
  doublereal alpha = solInfo.andes_alpha;
  doublereal betam = solInfo.andes_betam;

  // Local variables 
  int i, j;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area, h;
  doublereal str[6];
  doublereal z, epszz;
  doublereal temp;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;
  Eigen::Matrix<doublereal,3,3> eframe, gframe = Eigen::Matrix<doublereal,3,3>::Identity();
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
  Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > stress(_stress,maxstr,3);
  Eigen::Matrix<doublereal,3,1> sigma, epsilon;
  Eigen::Matrix<doublereal,6,1> Upsilon, Sigma;

  // Some convenient definitions 
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    e = Upsilon.head(3), chi = Upsilon.tail(3);
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    N = Sigma.head(3), M = Sigma.tail(3);


// ==================================================================== 
//                                                                      
//     -----------------                                                
//     V A R I A B L E S                                                
//     -----------------                                                
//                                                                      
//     elm      <input>   Finite Element Number                         
//     maxstr   <input>   Maximum Number of Stresses                    
//     nu       <input>   Poisson's Ratio (for an Isotropic Element)    
//     X        <input>   X- Nodal Coordinates                          
//     Y        <input>   Y- Nodal Coordinates                          
//     Z        <input>   Z- Nodal Coordinates                          
//     V        <input>   Nodal Generalized Displacements               
//     stress   <output>  Stresses (Von Mises Stress) of the Element    
//     ctyp     <input>   Type of Constitutive Law (0, 1, 2, 3, 4, or 5)   
//                                                                      
// ==================================================================== 
// Author   = Francois M. Hemez                                         
// Date     = June 10th, 1995                                           
// Version  = 2.0                                                       
// Modified = K. H. Pierson                                             
// Date     = April 11, 1997                                            
// Reason   = Added stress calculations for sigmaxx, sigmayy, sigmaxy   
//            and von mises stress at top, median and bottom surfaces   
//            Also added strain calculations for epsilonxx, epsilonyy,  
//            epsilonzz, epsilonxy and an equivalent strain at top,     
//            median and bottom surfaces.                               
// ==================================================================== 

    h = nmat->GetShellThickness();

//     ---------------------------------- 
//     STEP 1                             
//     COMPUTE THE TRIANGULAR COORDINATES 
//     ---------------------------------- 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, X, Y, Z, eframe.data(), xlp, ylp, zlp, area);

//     --------------------------------------------------- 
//     STEP 2                                              
//     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT 
//     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)    
//     --------------------------------------------------- 

    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

//     ------------------------------------------------- 
//     STEP 3                                            
//     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT  
//     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX) 
//     ------------------------------------------------- 

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

//     ------------------------------------------- 
//     STEP 4                                      
//     ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL 
//     FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT)   
//     ------------------------------------------- 

// .....CONSTRUCT THE PERMUTATION MATRIX
    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> P(indices);

    if(flag == 1) {

        for(i = 0; i < 18; i += 3)
            vd.segment(i,3) = eframe.transpose()*v.segment(i,3);

        vd = P.transpose()*vd;
    }
    else {
        // v is already local
        vd = P.transpose()*v;

    }

//     ----------------------------------------------------- 
//     STEP 5                                                
//     COMPUTE THE ELEMENTAL EXTENSION AND CURVATURE VECTORS 
//     ----------------------------------------------------- 

// .....COMPUTE THE Z- COORDINATE OF THE SELECTED SURFACE

    if(surface == 1) z = h/2;     // upper surface
    else if(surface == 2) z = 0;  // median surface
    else z = -h/2;                // lower surface

    // compute stresses and strains at the nodes
    doublereal zeta[3][3] = { { 1.,0.,0. }, { 0.,1.,0. }, { 0.,0.,1. } }; // triangular coordinates of nodes
    for(i = 0; i < 3; ++i) {

        if(sflg == 0) {
// .....ELEMENTAL CURVATURE COMPUTATION

            chi = (1/area)*Lb.transpose()*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION

            e = (1/area)*Lm.transpose()*vd.head(9);
        }
        else {
// .....ELEMENTAL CURVATURE COMPUTATION (including now the higher order contribution)

            Bb = (1/area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
            chi = Bb*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION (including now the higher order contribution)

            Bm = (1/area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
            e = Bm*vd.head(9);
        }

// .....NODAL TEMPERATURE
        temp = (ndtemps) ? ndtemps[i] : nmat->GetAmbientTemperature();

//     -------------------------------------------------
//     STEP 6
//     COMPUTE THE STRESS OR STRAIN OR HISTORY VARIABLES
//     -------------------------------------------------

        switch(strainflg) {

          case 1 : {

// .....COMPUTE THE LOCAL STRAINS ON THE SPECIFIED SURFACE

            epsilon = e + z * chi;

            // TODO nu is not necessarily defined, and should be obtained from the material
            // furthermore this equation is only correct for isotropic plane stress material
            epszz = -nu / (1. - nu) * (epsilon[0] + epsilon[1]);

// .....CALCULATE VON MISES EQUIVALENT STRAIN

            stress(6, i) = equivstr(epsilon[0], epsilon[1], epszz, 0.5*epsilon[2]);

// .....ROTATE LOCAL STRAINS TO GLOBAL AND CONVERT SHEAR STRAINS TO ENGINEERING SHEAR STRAINS

            str[0] = epsilon[0];
            str[1] = epsilon[1];
            str[2] = epszz;
            str[3] = 0.5*epsilon[2];
            str[4] = 0.;
            str[5] = 0.;
            transform(eframe.data(), gframe.data(), str);
            for (j = 0; j < 3; ++j)
                stress(j, i) = str[j];
            for (j = 3; j < 6; ++j) 
                stress(j, i) = 2*str[j];

          } break;

          default :  
          case 0 : {

            if((sflg == 0 && (ctyp == 0 || ctyp == 2 || ctyp == 3)) || ctyp == 1 || ctyp == 5) {

// .....COMPUTE THE GENERALIZED STRESSES [Sigma = {N,M}] WHICH ARE
// .....FORCE AND MOMENT PER UNIT LENGTH

                if(i == 0 || sflg != 0)
                    nmat->GetConstitutiveResponse(Upsilon.data(), Sigma.data(), NULL, eframe.data(), i, temp);

                if (surface == 1) {

// .....ESTIMATE THE LOCAL STRESSES ON THE UPPER SURFACE

                    sigma = N/h + 6*M/(h*h); 

                }

                else if (surface == 2) {

// .....ESTIMATE THE LOCAL STRESSES ON THE MEDIAN SURFACE

                    sigma = N/h;

                }

                else if (surface == 3) {

// .....ESTIMATE THE LOCAL STRESSES ON THE LOWER SURFACE

                    sigma = N/h - 6*M/(h*h);

                }
            }
            else {

// .....COMPUTE THE LOCAL STRESSES ON THE SPECIFIED SURFACE
                nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i, temp, 0, staten, statenp);

            }

// .....CALCULATE VON MISES EQUIVALENT STRESS

            stress(6, i) = equivstr(sigma[0], sigma[1], 0, sigma[2]);

// .....ROTATE LOCAL STRESSES TO GLOBAL

            str[0] = sigma[0];
            str[1] = sigma[1];
            str[2] = 0.;
            str[3] = sigma[2];
            str[4] = 0.;
            str[5] = 0.;
            transform(eframe.data(), gframe.data(), str);
            for (j = 0; j < 6; ++j)
                stress(j, i) = str[j];

          } break;

          case 2 : {

// .....COMPUTE THE EQUIVALENT PLASTIC STRAIN FOR ELASTO-PLASTIC MATERIALS

            if(ctyp == 4) {
              nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i, temp, 0, staten, statenp);
              stress(0, i) = nmat->GetLocalEquivalentPlasticStrain(i, z, statenp);
            }
            else {
              stress(0, i) = 0;
            }

          } break;

          case 3 : {

// .....COMPUTE THE BACKSTRESS FOR ELASTO-PLASTIC MATERIALS

            if(ctyp == 4) {
              nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i, temp, 0, staten, statenp);
              std::vector<doublereal> sigma = nmat->GetLocalBackStress(i, z, statenp);

// .....ROTATE LOCAL STRESSES TO GLOBAL

              str[0] = sigma[0];
              str[1] = sigma[1];
              str[2] = 0.;
              str[3] = sigma[2];
              str[4] = 0.;
              str[5] = 0.;
              transform(eframe.data(), gframe.data(), str);
              for (j = 0; j < 6; ++j)
                stress(j, i) = str[j];
  
            }
            else {
              for (j = 0; j < 6; ++j)
                stress(j, i) = 0;
            }

          } break;

          case 4 : {

// .....COMPUTE THE PLASTIC STRAIN TENSOR FOR ELASTO-PLASTIC MATERIALS

            if(ctyp == 4) {
              nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i, temp, 0, staten, statenp);
              std::vector<doublereal> epsilon = nmat->GetLocalPlasticStrain(i, z, statenp);

// .....ROTATE LOCAL STRAINS TO GLOBAL AND CONVERT SHEAR STRAINS TO ENGINEERING SHEAR STRAINS

              str[0] = epsilon[0];
              str[1] = epsilon[1];
              str[2] = -(epsilon[0]+epsilon[1]);
              str[3] = 0.5*epsilon[2];
              str[4] = 0.;
              str[5] = 0.;
              transform(eframe.data(), gframe.data(), str);
              for (j = 0; j < 3; ++j)
                stress(j, i) = str[j];
              for (j = 3; j < 6; ++j)
                stress(j, i) = 2*str[j];

            } 
            else {
              for (j = 0; j < 6; ++j)
                stress(j, i) = 0;
            }
          
          } break;

          case 5 : {

// .....COMPUTE THE SCALAR DAMAGE FOR ELASTO-PLASTIC MATERIALS
            if(ctyp == 4) {
              nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i, temp, 0, staten, statenp);
              stress(0, i) = nmat->GetLocalDamage(i, z, statenp);
            }
            else {
              stress(0, i) = 0;
            }

          } break;

        }
    }
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesvmsWRTdisp(int elm, doublereal nu, 
                  doublereal *X, doublereal *Y, doublereal *Z,
                  doublereal *_v, doublereal *_vmsWRTdisp,
                  int ctyp, ShellMaterial<doublereal> *nmat,
                  int surface, int sflg, doublereal *ndtemps)
{
  // Initialized data 
  doublereal clr = solInfo.andes_clr;
  doublereal cqr = solInfo.andes_cqr;
  doublereal betab = solInfo.andes_betab;
  doublereal alpha = solInfo.andes_alpha;
  doublereal betam = solInfo.andes_betam;

  // Local variables 
  int i, j;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area, h;
  doublereal str[6];
  doublereal z, epszz;
  doublereal temp;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,6,18> B;
  Eigen::Matrix<doublereal,3,3> eframe;
  Eigen::Map<Eigen::Matrix<doublereal,3,18> > vmsWRTdisp(_vmsWRTdisp);
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
  Eigen::Matrix<doublereal,3,18> dsigmadu;
  Eigen::Matrix<doublereal,3,1> sigma, epsilon;
  Eigen::Matrix<doublereal,6,18> dUpsilondu;
  Eigen::Matrix<doublereal,6,1> Upsilon, Sigma;
  Eigen::Matrix<doublereal,6,18> dSigmadu;
  dSigmadu.setZero();

  // Some convenient definitions 
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    e = Upsilon.head(3), chi = Upsilon.tail(3);
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    N = Sigma.head(3), M = Sigma.tail(3);
  Eigen::Block< Eigen::Matrix<doublereal,6,18>,3,9 >
    Bm = B.template topLeftCorner<3,9>(),     Bmb = B.template topRightCorner<3,9>(),
    Bbm = B.template bottomLeftCorner<3,9>(), Bb = B.template bottomRightCorner<3,9>();
  Bmb.setZero(); Bbm.setZero();

// ==================================================================== 
//                                                                      
//     -----------------                                                
//     V A R I A B L E S                                                
//     -----------------                                                
//                                                                      
//     elm        <input>   Finite Element Number                         
//     nu         <input>   Poisson's Ratio (for an Isotropic Element)    
//     globalX    <input>   X- Nodal Coordinates                          
//     globalY    <input>   Y- Nodal Coordinates                          
//     globalZ    <input>   Z- Nodal Coordinates                          
//     globalU    <input>   Global Displacements at the Nodal Joints      
//     vmsWRTdisp <output>  Derivative of Von Mises Stress w.r.t        
//                          displacement                                
//     ctyp       <input>   Type of Constitutive Law (0, 1, 2, 3, 4, or 5)      
//                                                                      
// ==================================================================== 
// Author  = Youngsoo Choi                                              
// Date    = January 17, 2014                                           
// Version = 1.0                                                        
// ==================================================================== 

    h = nmat->GetShellThickness();

//     ---------------------------------- 
//     STEP 1                             
//     COMPUTE THE TRIANGULAR COORDINATES 
//     ---------------------------------- 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, X, Y, Z, eframe.data(), xlp, ylp, zlp, area);

//     --------------------------------------------------- 
//     STEP 2                                              
//     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT 
//     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)    
//     --------------------------------------------------- 

    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

//     ------------------------------------------------- 
//     STEP 3                                            
//     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT  
//     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX) 
//     ------------------------------------------------- 

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

//     ------------------------------------------- 
//     STEP 4                                      
//     ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL 
//     FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT)   
//     ------------------------------------------- 

    for(i = 0; i < 18; i += 3)
        vd.segment(i,3) = eframe.transpose()*v.segment(i,3);

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> P(indices);

    vd = P.transpose()*vd;

// .....COMPUTE THE Z- COORDINATE OF THE SELECTED SURFACE

    if(surface == 1) z = h/2;     // upper surface
    else if(surface == 2) z = 0;  // median surface
    else z = -h/2;                // lower surface

    doublereal zeta[3][3] = { { 1.,0.,0. }, { 0.,1.,0. }, { 0.,0.,1. } }; // triangular coordinates of nodes
    for(i = 0; i < 3; ++i) {

//     ----------------------------------------------------- 
//     STEP 5                                                
//     COMPUTE THE ELEMENTAL EXTENSION AND CURVATURE VECTORS 
//     AND THEIR SENSITIVITIES                               
//     ----------------------------------------------------- 

        if(sflg == 0) {
            Bb = (1./area)*Lb.transpose();
            Bm = (1./area)*Lm.transpose();
        }
        else {
            Bb = (1/area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
            Bm = (1/area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
        }

        chi = Bb*vd.tail(9);
        e = Bm*vd.head(9);

        dUpsilondu = B*P.transpose();
        for(int j=0; j<6; j+=3) for(int k=0; k<18; k+=3) dUpsilondu.template block<3,3>(j,k) *= eframe.transpose();

// .....NODAL TEMPERATURE
        temp = (ndtemps) ? ndtemps[i] : nmat->GetAmbientTemperature();

//     -------------------------------------------------
//     STEP 6
//     COMPUTE THE VONMISES STRESS SENSITIVITY
//     -------------------------------------------------

        if((sflg == 0 && (ctyp == 0 || ctyp == 2 || ctyp == 3)) || ctyp == 1 || ctyp == 5) {

// .....COMPUTE THE GENERALIZED STRESSES [Sigma = {N,M}] WHICH ARE
// .....FORCE AND MOMENT PER UNIT LENGTH, AND THEIR SENSITIVITIES

            if(i == 0 || sflg != 0) {
                nmat->GetConstitutiveResponse(Upsilon.data(), Sigma.data(), NULL, eframe.data(), i, temp);
                nmat->GetConstitutiveResponseSensitivityWRTdisp(dUpsilondu.data(), dSigmadu.data(), NULL, eframe.data(), i);
            }

            if (surface == 1) {

// .....ESTIMATE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE UPPER SURFACE

                sigma = N/h + 6*M/(h*h); 
                dsigmadu = 1./h*dSigmadu.topRows(3) + 6./(h*h)*dSigmadu.bottomRows(3);
            }

            else if (surface == 2) {

// .....ESTIMATE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE MEDIAN SURFACE

                sigma = N/h;
                dsigmadu = 1./h*dSigmadu.topRows(3);
            }

            else if (surface == 3) {

// .....ESTIMATE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE LOWER SURFACE

                sigma = N/h - 6*M/(h*h);
                dsigmadu = 1./h*dSigmadu.topRows(3) - 6./(h*h)*dSigmadu.bottomRows(3);
            }
        }
        else {

// .....COMPUTE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE SPECIFIED SURFACE

            nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i, temp);
            nmat->GetLocalConstitutiveResponseSensitivityWRTdisp(dUpsilondu.data(), dsigmadu.data(), z, eframe.data(), i);
        }

// .....CALCULATE DERIVATIVE OF VON MISES EQUIVALENT STRESS WRT NODAL DISPLACEMENTS

        doublereal vms = equivstr(sigma[0], sigma[1], 0, sigma[2]);
        vmsWRTdisp.row(i) = equivstrSensitivityWRTdisp(vms, sigma[0], sigma[1], 0, sigma[2], dsigmadu);
    }

}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesvmsWRTthic(int elm, doublereal nu, 
                  doublereal *X, doublereal *Y, doublereal *Z,
                  doublereal *_v, doublereal *_vmsWRTthic,
                  int ctyp, ShellMaterial<doublereal> *nmat,
                  int surface, int sflg, doublereal *ndtemps)
{
  // Initialized data 
  doublereal clr = solInfo.andes_clr;
  doublereal cqr = solInfo.andes_cqr;
  doublereal betab = solInfo.andes_betab;
  doublereal alpha = solInfo.andes_alpha;
  doublereal betam = solInfo.andes_betam;

  // Local variables 
  int i, j;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area, h;
  doublereal str[6];
  doublereal z, epszz, dzdh;
  doublereal temp;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;
  Eigen::Matrix<doublereal,3,3> eframe, gframe = Eigen::Matrix<doublereal,3,3>::Identity();
  Eigen::Map<Eigen::Matrix<doublereal,3,1> > vmsWRTthic(_vmsWRTthic);
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
  Eigen::Matrix<doublereal,18,18> Eframe;
  Eigen::Matrix<doublereal,3,1> sigma, epsilon, dsigmadh;
  Eigen::Matrix<doublereal,6,1> Upsilon, Sigma;
  Eigen::Matrix<doublereal,6,1> dSigmadh;

  // Some convenient definitions 
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    e = Upsilon.head(3), chi = Upsilon.tail(3);
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    N = Sigma.head(3), M = Sigma.tail(3);


// ==================================================================== 
//                                                                      
//     -----------------                                                
//     V A R I A B L E S                                                
//     -----------------                                                
//                                                                      
//     elm        <input>   Finite Element Number                       
//     nu         <input>   Poisson's Ratio (for an Isotropic Element)  
//     X          <input>   X- Nodal Coordinates                        
//     Y          <input>   Y- Nodal Coordinates                        
//     Z          <input>   Z- Nodal Coordinates                        
//     V          <input>   Nodal Generalized Displacements             
//     vmsWRTthic <output>  Derivative of Von Mises Stress w.r.t        
//                          thickness                                   
//     ctyp       <input>   Type of Constitutive Law (0, 1, 2, 3, 4, or 5) 
//                                                                      
// ==================================================================== 
// Author  = Youngsoo Choi                                              
// Date    = January 17, 2014                                           
// Version = 1.0                                                        
// ====================================================================  

    h = nmat->GetShellThickness();

//     ---------------------------------- 
//     STEP 1                             
//     COMPUTE THE TRIANGULAR COORDINATES 
//     ---------------------------------- 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, X, Y, Z, eframe.data(), xlp, ylp, zlp, area);

//     --------------------------------------------------- 
//     STEP 2                                              
//     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT 
//     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)    
//     --------------------------------------------------- 

    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

//     ------------------------------------------------- 
//     STEP 3                                            
//     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT  
//     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX) 
//     ------------------------------------------------- 

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

//     ------------------------------------------- 
//     STEP 4                                      
//     ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL 
//     FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT)   
//     ------------------------------------------- 

    for(i = 0; i < 18; i += 3)
        vd.segment(i,3) = eframe.transpose()*v.segment(i,3);

//     ----------------------------------------------------- 
//     STEP 5                                                
//     COMPUTE THE ELEMENTAL EXTENSION AND CURVATURE VECTORS 
//     ----------------------------------------------------- 

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> P(indices);

    vd = P.transpose()*vd;

// .....COMPUTE THE Z- COORDINATE OF THE SELECTED SURFACE

    if(surface == 1) { z = h/2; dzdh = 1./2.; } // upper surface
    else if(surface == 2) { z = 0;  dzdh = 0; } // median surface
    else { z = -h/2; dzdh = -1./2.; }           // lower surface

    // compute stresses and strains at the nodes
    doublereal zeta[3][3] = { { 1.,0.,0. }, { 0.,1.,0. }, { 0.,0.,1. } }; // triangular coordinates of nodes
    for(i = 0; i < 3; ++i) {

        if(sflg == 0) {
// .....ELEMENTAL CURVATURE COMPUTATION

            chi = (1./area)*Lb.transpose()*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION

            e = (1./area)*Lm.transpose()*vd.head(9);
        }
        else {
// .....ELEMENTAL CURVATURE COMPUTATION (including now the higher order contribution)

            Bb = (1./area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
            chi = Bb*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION (including now the higher order contribution)

            Bm = (1./area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
            e = Bm*vd.head(9);
        }

// .....NODAL TEMPERATURE
        temp = (ndtemps) ? ndtemps[i] : nmat->GetAmbientTemperature();

//     -------------------------------------------------
//     STEP 6
//     COMPUTE THE VON MISES STRESS SENSITIVITY
//     -------------------------------------------------

        if((sflg == 0 && (ctyp == 0 || ctyp == 2 || ctyp == 3)) || ctyp == 1 || ctyp == 5) {

// .....COMPUTE THE GENERALIZED STRESSES [Sigma = {N,M}] WHICH ARE
// .....FORCE AND MOMENT PER UNIT LENGTH, AND THEIR SENSITIVITIES

            if(i == 0 || sflg != 0) {
                nmat->GetConstitutiveResponse(Upsilon.data(), Sigma.data(), NULL, eframe.data(), i, temp);
                nmat->GetConstitutiveResponseSensitivityWRTthic(Upsilon.data(), dSigmadh.data(), NULL, eframe.data(), i, temp);
            }

            if (surface == 1) {

// .....ESTIMATE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE UPPER SURFACE

                sigma = N/h + 6*M/(h*h); 
                dsigmadh = 1./h*dSigmadh.head(3) + 6./(h*h)*dSigmadh.tail(3)
                         - N/(h*h) - 12*M/(h*h*h);
            }

            else if (surface == 2) {

// .....ESTIMATE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE MEDIAN SURFACE

                sigma = N/h;
                dsigmadh = 1./h*dSigmadh.head(3) - N/(h*h);
            }

            else if (surface == 3) {

// .....ESTIMATE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE LOWER SURFACE

                sigma = N/h - 6*M/(h*h);
                dsigmadh = 1./h*dSigmadh.head(3) - 6./(h*h)*dSigmadh.tail(3)
                         - N/(h*h) + 12*M/(h*h*h);
            }
        }
        else {

// .....COMPUTE THE LOCAL STRESSES AND THEIR SENSITIVITIES ON THE SPECIFIED SURFACE

            nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i, temp);
            dsigmadh.setZero();
            nmat->GetLocalConstitutiveResponseSensitivityWRTthic(Upsilon.data(), dsigmadh.data(), dzdh, eframe.data(), i);

        }

// .....CALCULATE DERIVATIVE OF VON MISES EQUIVALENT STRESS WRT THICKNESS

        doublereal vms = equivstr(sigma[0], sigma[1], 0, sigma[2]);
        vmsWRTthic[i] = equivstrSensitivityWRTthic(vms, sigma[0], sigma[1], 0, sigma[2], dsigmadh);
    }

}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
doublereal
ShellElementTemplate<doublereal,Membrane,Bending>
::equivstr(doublereal sxx, doublereal syy, doublereal szz, doublereal sxy)
{
    // Builtin functions 
    using std::sqrt;

    // Local variables 
    doublereal s0, dsxx, dsyy, dszz, eq;

// ... COMPUTE MEAN HYDROSTATIC STRESS OR STRAIN 

    s0 = (sxx + syy + szz) / 3.;

// ... COMPUTE DEVIATORIC STRESSES OR STRAINS 

    dsxx = sxx - s0;
    dsyy = syy - s0;
    dszz = szz - s0;

// ... COMPUTE EQUIVALENT STRESS OR STRAIN 

    eq = (dsxx * dsxx + dsyy * dsyy + dszz * dszz) / 2. + sxy * sxy;
    eq = sqrt(eq * 3.);

    return eq;
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
Eigen::Matrix<doublereal,1,18>
ShellElementTemplate<doublereal,Membrane,Bending>
::equivstrSensitivityWRTdisp(doublereal vms, doublereal sxx, doublereal syy, doublereal szz,
                             doublereal sxy, Eigen::Matrix<doublereal,3,18> &dsigmadu)
{
    if(vms == 0) {
      Eigen::Matrix<doublereal,1,18> a;
      a.setZero();
      return a;
    }

    // Builtin functions 
    using std::sqrt;

    // Local variables 
    doublereal s0, dsxx, dsyy, dszz, eq;

// ... COMPUTE MEAN HYDROSTATIC STRESS OR STRAIN 

    s0 = (sxx + syy + szz) / 3.;

// ... COMPUTE DEVIATORIC STRESSES OR STRAINS 

    dsxx = sxx - s0;
    dsyy = syy - s0;
    dszz = szz - s0;

// ... COMPUTE EQUIVALENT STRESS OR STRAIN 

    Eigen::Matrix<doublereal,3,18> dsdu;
    Eigen::Matrix<doublereal,3,3> D;
    D << 2./3., -1./3., 0.,
        -1./3.,  2./3., 0.,
        -1./3., -1./3., 0.;
    dsdu = D*dsigmadu;

    return 3*dsxx/(2*vms)*dsdu.row(0) + 
           3*dsyy/(2*vms)*dsdu.row(1) + 
           3*dszz/(2*vms)*dsdu.row(2) + 
           3*sxy/vms*dsigmadu.row(2);
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
doublereal
ShellElementTemplate<doublereal,Membrane,Bending>
::equivstrSensitivityWRTthic(doublereal vms, doublereal sxx, doublereal syy, doublereal szz,
                             doublereal sxy, Eigen::Matrix<doublereal,3,1> &dsigmadh)
{
    // Builtin functions 
    using std::sqrt;

    // Local variables 
    doublereal s0, dsxx, dsyy, dszz, eq;

// ... COMPUTE MEAN HYDROSTATIC STRESS OR STRAIN 

    s0 = (sxx + syy + szz) / 3.;

// ... COMPUTE DEVIATORIC STRESSES OR STRAINS 

    dsxx = sxx - s0;
    dsyy = syy - s0;
    dszz = szz - s0;

// ... COMPUTE EQUIVALENT STRESS OR STRAIN 

    Eigen::Matrix<doublereal,3,1> dsdh;
    Eigen::Matrix<doublereal,3,3> D;
    D << 2./3., -1./3., 0.,
        -1./3.,  2./3., 0.,
        -1./3., -1./3., 0.;
    dsdh = D*dsigmadh;

    return 3*dsxx/(2*vms)*dsdh[0] + 
           3*dsyy/(2*vms)*dsdh[1] + 
           3*dszz/(2*vms)*dsdh[2] + 
           3*sxy/vms*dsigmadh[2];
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::transform(doublereal *_lframe, doublereal *_gframe, doublereal *_str)
{
    // Local variables
    doublereal l1, l2, l3, m1, m2, m3, n1, n2, n3;
    Eigen::Map<Eigen::Matrix<doublereal,3,3> > lframe(_lframe), gframe(_gframe);
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > str(_str);
    Eigen::Matrix<doublereal,6,6> T;


// **********************************************************************

// Purpose: to form a transformation matrix from local coordinates to 
//          global coordinates 

// input variables: 
//      xl = x local unit vector = first column of lframe
//      yl = y local unit vector = second column of lframe
//      zl = z local unit vector = third column of lframe
//      xg = x global unit vector = first column of gframe
//      yg = y global unit vector = second column of gframe
//      zg = z global unit vector = third column of gframe
//      str = stress/strain 6x1 vector 
//            sigmaxx, sigmayy, sigmazz, sigma12, sigma23, sigma13 

// local variables: 
//      l1 = direction cosine between xl and xg 
//      l2 = direction cosine between xl and yg 
//      l3 = direction cosine between xl and zg 
//      m1 = direction cosine between yl and xg 
//      m2 = direction cosine between yl and yg 
//      m3 = direction cosine between yl and zg 
//      n1 = direction cosine between zl and xg 
//      n2 = direction cosine between zl and yg 
//      n3 = direction cosine between zl and zg 
//      T  = transformation matrix from local to global 

// **********************************************************************

// Compute direction cosines 
/*
    l1 = xg[0] * xl[0] + xg[1] * xl[1] + xg[2] * xl[2];
    l2 = yg[0] * xl[0] + yg[1] * xl[1] + yg[2] * xl[2];
    l3 = zg[0] * xl[0] + zg[1] * xl[1] + zg[2] * xl[2];
    m1 = xg[0] * yl[0] + xg[1] * yl[1] + xg[2] * yl[2];
    m2 = yg[0] * yl[0] + yg[1] * yl[1] + yg[2] * yl[2];
    m3 = zg[0] * yl[0] + zg[1] * yl[1] + zg[2] * yl[2];
    n1 = xg[0] * zl[0] + xg[1] * zl[1] + xg[2] * zl[2];
    n2 = yg[0] * zl[0] + yg[1] * zl[1] + yg[2] * zl[2];
    n3 = zg[0] * zl[0] + zg[1] * zl[1] + zg[2] * zl[2];
*/
    Eigen::Matrix<doublereal,3,3> dc = gframe.transpose()*lframe;

    l1 = dc(0,0);
    l2 = dc(1,0);
    l3 = dc(2,0);
    m1 = dc(0,1);
    m2 = dc(1,1);
    m3 = dc(2,1);
    n1 = dc(0,2);
    n2 = dc(1,2);
    n3 = dc(2,2);

// Construct the 6x6 transformation matrix 

    T(0, 0) = l1 * l1;
    T(0, 1) = m1 * m1;
    T(0, 2) = n1 * n1;
    T(0, 3) = l1 * 2 * m1;
    T(0, 4) = m1 * 2 * n1;
    T(0, 5) = n1 * 2 * l1;

    T(1, 0) = l2 * l2;
    T(1, 1) = m2 * m2;
    T(1, 2) = n2 * n2;
    T(1, 3) = l2 * 2 * m2;
    T(1, 4) = m2 * 2 * n2;
    T(1, 5) = n2 * 2 * l2;

    T(2, 0) = l3 * l3;
    T(2, 1) = m3 * m3;
    T(2, 2) = n3 * n3;
    T(2, 3) = l3 * 2 * m3;
    T(2, 4) = m3 * 2 * n3;
    T(2, 5) = n3 * 2 * l3;

    T(3, 0) = l1 * l2;
    T(3, 1) = m1 * m2;
    T(3, 2) = n1 * n2;
    T(3, 3) = l1 * m2 + l2 * m1;
    T(3, 4) = m1 * n2 + m2 * n1;
    T(3, 5) = n1 * l2 + n2 * l1;

    T(4, 0) = l2 * l3;
    T(4, 1) = m2 * m3;
    T(4, 2) = n2 * n3;
    T(4, 3) = l2 * m3 + l3 * m2;
    T(4, 4) = m2 * n3 + m3 * n2;
    T(4, 5) = n2 * l3 + n3 * l2;

    T(5, 0) = l3 * l1;
    T(5, 1) = m3 * m1;
    T(5, 2) = n3 * n1;
    T(5, 3) = l3 * m1 + l1 * m3;
    T(5, 4) = m3 * n1 + m1 * n3;
    T(5, 5) = n3 * l1 + n1 * l3;

// Perform the multiplication {str'} = T{str} 

    str = T*str;
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesups(int elm, doublereal *staten, doublereal *statenp, doublereal *X, doublereal *Y, doublereal *Z,
           doublereal *_v, ShellMaterial<doublereal> *gpmat, ShellMaterial<doublereal> *nmat,
           int sflg, int tflg, doublereal *ndtemps, doublereal dt)
{
  // Initialized data 
  doublereal clr = solInfo.andes_clr;
  doublereal cqr = solInfo.andes_cqr;
  doublereal betab = solInfo.andes_betab;
  doublereal alpha = solInfo.andes_alpha;
  doublereal betam = solInfo.andes_betam;

  // Local variables 
  int i;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area;
  doublereal temp;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;
  Eigen::Matrix<doublereal,3,3> eframe;
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
  Eigen::Matrix<doublereal,6,1> Upsilon;

  // Some convenient definitions 
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    e = Upsilon.head(3), chi = Upsilon.tail(3);

// ==================================================================== 
//                                                                      
//     -----------------                                                
//     V A R I A B L E S                                                
//     -----------------                                                
//                                                                      
//     elm      <input>   Finite Element Number                         
//     staten   <input>   Material States at t_n             
//     statenp  <output>  Material States at t_{n+1}
//     X        <input>   X- Nodal Coordinates                          
//     Y        <input>   Y- Nodal Coordinates                          
//     Z        <input>   Z- Nodal Coordinates                          
//     _v       <input>   Global Displacements at the Nodal Joints      
//                                                                      
// ==================================================================== 
// Author   = Philip J. S. Avery                                       
// Date     = September 13, 2011        
// Version  = 1.0                                                        
// Comment  =   
// ==================================================================== 

//     ---------------------------------- 
//     STEP 1                             
//     COMPUTE THE TRIANGULAR COORDINATES 
//     ---------------------------------- 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, X, Y, Z, eframe.data(), xlp, ylp, zlp, area);

//     --------------------------------------------------- 
//     STEP 2                                              
//     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT 
//     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)    
//     --------------------------------------------------- 

    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

//     ------------------------------------------------- 
//     STEP 3                                            
//     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT  
//     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX) 
//     ------------------------------------------------- 

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

//     ----------------------------------------------------- 
//     STEP 5                                                
//     COMPUTE THE ELEMENTAL EXTENSION AND CURVATURE VECTORS 
//     ----------------------------------------------------- 

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> P(indices);

    vd = P.transpose()*v; // note: v is already local in this case

    if(tflg == 0) {
        // if ndtemps (nodal temperatures) is not null, then the mean temperature is used at each gauss point.
        temp = (ndtemps) ? doublereal((ndtemps[0]+ndtemps[1]+ndtemps[2])/3) : gpmat->GetAmbientTemperature();
    }
    else {
        // if ndtemps (nodal temperatures) is not null, then the temperature at each gauss point is interpolated.
        if(!ndtemps) temp = gpmat->GetAmbientTemperature();
    }

    // compute updated material state at the gauss points
    doublereal zeta[3][3] = { { 0.,.5,.5 }, { .5,0.,.5 }, { .5,.5,0. } }; // triangular coordinates of gauss integration points
    for(i = 0; i < 3; ++i) {

// .....ELEMENTAL CURVATURE COMPUTATION

        Bb = (1/area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
        chi = Bb*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION

        Bm = (1/area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
        e = Bm*vd.head(9);

// .....COMPUTE THE UPDATED MATERIAL STATE AT GAUSS POINT i

        if(ndtemps && tflg != 0)
            temp = zeta[i][0]*ndtemps[0] + zeta[i][1]*ndtemps[1] + zeta[i][2]*ndtemps[2];
        gpmat->UpdateState(Upsilon.data(), staten, statenp, i, temp, dt);
        staten += 5*gpmat->GetNumLocalStates(); // 5 is the number of layers per gauss point
        statenp += 5*gpmat->GetNumLocalStates();

    }

    // compute updated material state at the nodes
    doublereal eta[3][3] = { { 1.,0.,0. }, { 0.,1.,0. }, { 0.,0.,1. } }; // triangular coordinates of nodes
    for(i = 0; i < 3; ++i) {

        if(sflg == 0) {
// .....ELEMENTAL CURVATURE COMPUTATION

            chi = (1/area)*Lb.transpose()*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION

            e = (1/area)*Lm.transpose()*vd.head(9);
        }
        else {
// .....ELEMENTAL CURVATURE COMPUTATION (including now the higher order contribution)

            Bb = (1/area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, eta[i]);
            chi = Bb*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION (including now the higher order contribution)

            Bm = (1/area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, eta[i]);
            e = Bm*vd.head(9);
        }

// .....COMPUTE THE UPDATED MATERIAL STATE AT NODE i

        temp = (ndtemps) ? ndtemps[i] : gpmat->GetAmbientTemperature();
        nmat->UpdateState(Upsilon.data(), staten, statenp, i, temp, dt);
        staten += 3*nmat->GetNumLocalStates(); // 3 is the number of layers per node
        statenp += 3*nmat->GetNumLocalStates();

    }

}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesare(int elm, doublereal *_x, doublereal *_y, doublereal *_z,
           doublereal &area)
{
  // Builtin functions 
  using std::sqrt;
  using std::abs;

  // Local variables 
  int i, j;
  doublereal x21, y21, z21, x13, y13, z13, x32, y32, z32, signedarea,
             projection, side21length, side32length;
  Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y), z(_z);
  Eigen::Matrix<doublereal,3,3> xyz; xyz << x[0], x[1], x[2],
                                            y[0], y[1], y[2],
                                            z[0], z[1], z[2];
  Eigen::Matrix<doublereal,3,1> side21, side13, side32;

// ==================================================================== 
//                                                                      
//     Perform =    This subroutine computes the element area   
//     --------- 
//                                                                      
//                                                                      
//     Inputs/Outputs =                                                 
//     ----------------                                                 
//     ELM    <input>   finite element number                           
//     X      <input>   nodal coordinates in the X-direction            
//     Y      <input>   nodal coordinates in the Y-direction            
//     Z      <input>   nodal coordinates in the Z-direction            
//     AREA   <output>  element area                                    
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9th, 1994                                             
// Version = 1.0                                                        
// Comment =                                                            
// ==================================================================== 

// .....COMPUTE THE NODAL COORDINATE DIFFERENCES 

    side21 = xyz.col(1)-xyz.col(0);
    side13 = xyz.col(0)-xyz.col(2);
    side32 = xyz.col(2)-xyz.col(1);

// .....COMPUTE THE LENGTH OF SIDE 2-1 

    side21length = side21.norm();

// .....CHECK IF LENGTH 2-1 IS DIFFERENT FROM ZERO 

    if (side21length == 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andesare ***\n"
          "*** Side between nodes 1 and 2 Has zero length    ***\n"
          "*** Check coordinates and FE topology             ***\n");
    }

// .....COMPUTE THE LENGTH OF SIDE 3-2 

    side32length = side32.norm();

// .....CHECK IF LENGTH 3-2 IS DIFFERENT FROM ZERO 

    if (side32length == 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andesare ***\n"
          "*** Side between nodes 2 and 3 Has zero length    ***\n"
          "*** Check coordinates and FE topology             ***\n");
    }

// .....COMPUTE THE DISTANCE OF THE OPPOSING NODE 3 TO SIDE 2-1 

    projection = abs(side21.dot(side32))/side21length;

// .....GET THE AREA OF THE TRIANGLE 

    signedarea = side32length * side32length - projection * projection;

    if (signedarea <= 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in ShellElementTemplate::andesare ***\n"
          "*** The area is negative or zero                  ***\n"
          "*** Check coordinates and FE topology             ***\n");
    }

    area = side21length * .5 * sqrt(signedarea);
}

//(AN) Following routine now takes in statenp and passes it to the material dissipated energy function.

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesden(int elm, doublereal *X, doublereal *Y, doublereal *Z,
           ShellMaterial<doublereal> *gpmat, doublereal *statenp, doublereal &D)
{
  // Local variables
  int i;
  doublereal area;

// ==================================================================== 
//                                                                      
//     Perform =    This subroutine computes the cumulative energy   
//     ---------    dissipated by the element
//                                                                      
//                                                                      
//     Inputs/Outputs =                                                 
//     ----------------                                                 
//     ELM    <input>   finite element number                           
//     X      <input>   nodal coordinates in the X-direction            
//     Y      <input>   nodal coordinates in the Y-direction            
//     Z      <input>   nodal coordinates in the Z-direction            
//     D      <output>  cumulative dissipated energy                                   
//                                                                      
// ==================================================================== 
// Author  = Philip J. S. Avery                                          
// Date    = December 7th, 2013                                            
// Version = 1.0                                                        
// Comment =                                                            
// ==================================================================== 

// .....GET THE AREA OF THE TRIANGLE

    andesare(elm, X, Y, Z, area);

// .....INTEGRATE OVER THE AREA OF THE TRIANGLE

    D = 0.0;
    doublereal zeta[3][3] = { { 0.,.5,.5 }, { .5,0.,.5 }, { .5,.5,0. } }; // triangular coordinates of gauss integration points
    doublereal weight[3] = { 1/3., 1/3., 1/3. };
    for(i = 0; i < 3; ++i) {
      D += area*weight[i]*gpmat->GetDissipatedEnergy(i, statenp);
    }
}

#include <Element.d/Function.d/SpaceDerivatives.h>
#include <unsupported/Eigen/NumericalDiff>
#include <Element.d/FelippaShell.d/ShellElementGravityForceWRTNodalCoordinateSensitivity.h>

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesgfWRTcoord(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *_J,
                  doublereal *gamma, int gravflg, doublereal rhoh)
{
  const int senMethod = 1;
  doublereal eps = 1e-6;

  Eigen::Array<double,4,1> dconst;
  dconst << gamma[0], gamma[1], gamma[2], rhoh;

  Eigen::Array<int,1,1> iconst;
  iconst[0] = gravflg;

  Eigen::Matrix<double,9,1> q;
  q << x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2];

  Eigen::Map<Eigen::Matrix<double,18,9> > J(_J);

  switch(senMethod) {
    default:
#ifndef DEBUG_SHELL_ELEMENT_TEMPLATE
    case 1: { // automatic differentiation
      Simo::Jacobian<double,ShellElementGravityForceWRTNodalCoordinateSensitivity> dfdx(dconst,iconst);
      J = dfdx(q, 0);
    } break;
#endif
    case 2 : { // finite difference approximation 
      Simo::SpatialView<double,ShellElementGravityForceWRTNodalCoordinateSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementGravityForceWRTNodalCoordinateSensitivity>,Eigen::Central> nd(sv, eps);
      Eigen::Matrix<double,18,9> jac;
      nd.df(q, jac);
      J = jac;
    } break;
  }
}

#include <Element.d/FelippaShell.d/ShellElementMassWRTNodalCoordinateSensitivity.h>

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesmsWRTcoord(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *_J,
                  doublereal rhoh)
{
  const int senMethod = 1;
  doublereal eps = 1e-6;

  Eigen::Array<double,1,1> dconst;
  dconst << rhoh;

  Eigen::Array<int,0,1> iconst;

  Eigen::Matrix<double,9,1> q;
  q << x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2];

  Eigen::Map<Eigen::Matrix<double,1,9> > J(_J);

  switch(senMethod) {
    default:
#ifndef DEBUG_SHELL_ELEMENT_TEMPLATE
    case 1: { // automatic differentiation
      Simo::Jacobian<double,ShellElementMassWRTNodalCoordinateSensitivity> dfdx(dconst,iconst);
      J = dfdx(q, 0);
    } break;
#endif
    case 2 : { // finite difference approximation 
      Simo::SpatialView<double,ShellElementMassWRTNodalCoordinateSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementMassWRTNodalCoordinateSensitivity>,Eigen::Central> nd(sv, eps);
      Eigen::Matrix<double,1,9> jac;
      nd.df(q, jac);
      J = jac;
    } break;
  }
}

#include <Element.d/FelippaShell.d/ShellElementStiffnessWRTNodalCoordinateSensitivity.h>
#include <Element.d/FelippaShell.d/ShellMaterial.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType1.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType5.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesstfWRTcoord(int elm, doublereal *_destiffdx[9], doublereal E, doublereal nu, doublereal rho,
                   doublereal eh, doublereal Ta, doublereal W, doublereal *cFrame, doublereal *x,
                   doublereal *y, doublereal *z, int ctyp, doublereal *coefs, int flag, int tflg,
                   doublereal *ndtemps)
{
  const int senMethod = 1;
  doublereal eps = 1e-6;

  Eigen::Array<double,60,1> dconst;
  dconst[0] = E;
  dconst[1] = nu;
  dconst[2] = rho;
  dconst[3] = eh;
  if(ctyp == 1 || ctyp == 5) {
    dconst.segment<9>(4) = Eigen::Map<Eigen::Matrix<double,9,1> >(cFrame);
    dconst.segment<42>(13) = Eigen::Map<Eigen::Matrix<double,42,1> >(coefs);
  }
  dconst[55] = Ta;
  dconst[56] = W;
  if(ndtemps)
    dconst.segment<3>(57) = Eigen::Map<Eigen::Matrix<double,3,1> >(ndtemps);
  else
    dconst.segment<3>(57).setConstant(Ta);

  Eigen::Array<int,2,1> iconst;
  iconst[0] = ctyp;
  iconst[1] = tflg;

  Eigen::Matrix<double,9,1> q;
  q << x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2];

  Eigen::Array<Eigen::Matrix<double,18,18>,1,9> J;

  switch(senMethod) {
    default:
#ifndef DEBUG_SHELL_ELEMENT_TEMPLATE
    case 1: { // automatic differentiation
      Simo::FirstPartialSpaceDerivatives<double,ShellElementStiffnessWRTNodalCoordinateSensitivity> dfdx(dconst,iconst);
      J = dfdx(q, 0);
      for(int i=0; i<9; ++i) for(int j=0; j<324; ++j) _destiffdx[i][j] = J[i].data()[j];
    } break;
#endif
    case 2 : { // finite difference approximation 
      Simo::SpatialView<double,ShellElementStiffnessWRTNodalCoordinateSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementStiffnessWRTNodalCoordinateSensitivity>,Eigen::Central> nd(sv, eps);
      Eigen::Matrix<double,324,9> jac;
      nd.df(q, jac);
      for(int i=0; i<9; ++i) for(int j=0; j<324; ++j) _destiffdx[i][j] = jac.col(i)[j];
    } break;
  }
}

#include <Element.d/FelippaShell.d/ShellElementStressWRTNodalCoordinateSensitivity.h>

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesvmsWRTcoord(int elm, doublereal E, doublereal nu, doublereal rho, doublereal eh,
                   doublereal Ta, doublereal W, doublereal *cFrame, doublereal *x,
                   doublereal *y, doublereal *z, doublereal *v, doublereal *_J, int ctyp,
                   doublereal *coefs, int surface, int sflg, doublereal *ndtemps)
{
  const int senMethod = 1;
  doublereal eps = 1e-6;

  Eigen::Array<double,78,1> dconst;
  dconst.segment<18>(0) = Eigen::Map<Eigen::Matrix<double,18,1> >(v);
  dconst[18] = E;
  dconst[19] = nu;
  dconst[20] = rho;
  dconst[21] = eh;
  if(ctyp == 1 || ctyp == 5) {
    dconst.segment<9>(22) = Eigen::Map<Eigen::Matrix<double,9,1> >(cFrame);
    dconst.segment<42>(31) = Eigen::Map<Eigen::Matrix<double,42,1> >(coefs); 
  }
  dconst[73] = Ta;
  dconst[74] = W;
  if(ndtemps)
    dconst.segment<3>(75) = Eigen::Map<Eigen::Matrix<double,3,1> >(ndtemps);
  else 
    dconst.segment<3>(75).setConstant(Ta);

  Eigen::Array<int,3,1> iconst;
  iconst[0] = surface;
  iconst[1] = ctyp;
  iconst[2] = sflg;

  Eigen::Matrix<double,9,1> q;
  q << x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2];

  Eigen::Map<Eigen::Matrix<double,3,9> > J(_J);

  switch(senMethod) {
    default:
#ifndef DEBUG_SHELL_ELEMENT_TEMPLATE
    case 1: { // automatic differentiation
      Simo::Jacobian<double,ShellElementStressWRTNodalCoordinateSensitivity> dfdx(dconst,iconst);
      J = dfdx(q, 0);
      for(int i=0; i<3; ++i)
        for(int j=0; j<9; ++j)
          if(std::isnan(J(i,j))) J(i,j) = 0;
    } break;
#endif
    case 2 : { // finite difference approximation 
      Simo::SpatialView<double,ShellElementStressWRTNodalCoordinateSensitivity> sv(dconst,iconst,0.);
      Eigen::NumericalDiff<Simo::SpatialView<double,ShellElementStressWRTNodalCoordinateSensitivity>,Eigen::Central> nd(sv, eps);
      Eigen::Matrix<double,3,9> jac;
      nd.df(q, jac);
      J = jac;
    } break;
  }
}

#endif
#endif
