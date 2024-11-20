#ifdef USE_EIGEN3
#ifndef _SHEARPANELTEMPLATE_CPP_
#define _SHEARPANELTEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/Shear.d/ShearPanelTemplate.hpp>
#include <Element.d/Quad4.d/QuadElementTemplate.cpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal>
void
ShearPanelTemplate<doublereal>
::spstress(doublereal *_xg, doublereal *_yg, doublereal *_zg, doublereal *_v,
           doublereal G, doublereal E, doublereal F1, doublereal F2,
           doublereal *_stress, doublereal *_strain,
           doublereal &vmssig, doublereal &vmseps)
{ 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Gregrory W. Brown
C January, 2000
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     SPSTRESS calculates stresses and strains for the
C     four-node shear panel element.
C
C            call spstress (xg,yg,zg,v,G,E,F1,F2,
C    $                      stress,strain,vmssig,vmseps)
C
C     where the input arguments are
C
C       xg        (4 x 1) array of global x coordinates
C       yg        (4 x 1) array of global y coordinates
C       zg        (4 x 1) array of global y coordinates
C       v         (12x 1) array of node dispalcements
C                 arranged, vx1,vy1,vz1,...,vz4
C       G         Shear Modulus of element
C       E         Youngs Modulus of element
C       F1        Extensional factor 1
C       F2        Extensional factor 2
C
C     The outputs are:
C
C       stress    (6x4) corner node stresses arranged
C                  (sigxx, sigyy, sigzz, tauxy, tauyz, tauxz)
C       strain    (6x4) corner node strains 
C       vmssig    Vommises stress
C       vmseps    Vommises strain
C
*/
      using std::abs;
      using std::sqrt;

      Eigen::Map<Eigen::Matrix<doublereal,6,4> > stress(_stress), strain(_strain);
      Eigen::Map<Eigen::Matrix<doublereal,4,1> > xg(_xg), yg(_yg), zg(_zg);
      Eigen::Map<Eigen::Matrix<doublereal,12,1> > v(_v);

//                   L O C A L   V A R I A B L E S

// Local Coordinates
      Eigen::Matrix<doublereal,4,1> x,y,z,zdl;
      Eigen::Matrix<doublereal,8,1> vl;
// Working Vectors
      Eigen::Matrix<doublereal,3,4> node, local;
      Eigen::Matrix<doublereal,3,1> v1,v2,v3,v4;
      doublereal len;
// Consitutive Matrix
      Eigen::Matrix<doublereal,3,3> c;
// For sands2 routine
      Eigen::Matrix<doublereal,7,4> quadstress, quadstrain;
      int four, three, one, seven;
      bool falseflag;
// Effective Extentional Areas
      Eigen::Matrix<doublereal,4,1> area;
// Rod Coordinates and Stiffness
      Eigen::Matrix<doublereal,6,4> rodstrain;
      doublereal epsrod;
      Eigen::Matrix<doublereal,2,1> xr, yr, vrx, vry;
      Eigen::Matrix<int, 2,4> rnod; 
      doublereal l1,l2,q1,q2,dq;
// For vms routine
      Eigen::Matrix<doublereal,7,4> vmsstress, vmsstrain;
// Output Control
      bool output;
// Other Variables
      int i, j, k;
      doublereal avg;

// Set Output Status
      output = false;

/*
C     COMPUTE LOCAL COORDINATE SYSTEM 
C  (works best if element is plane)
C
C      Y                             
C     /|\
C      |  3            2             
C      |  /------------|             
C      | /             |             
C      |/              |             
C      /---------------|-----> X     
C      0               1             
C
C... store nodal coordinates in "node"
C... shift origin to node 0
*/
      for(i = 0; i<4; ++i) {
        node(0,i) = xg(i)-xg(0);
        node(1,i) = yg(i)-yg(0);
        node(2,i) = zg(i)-zg(0);
      }

      Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, " ");
//... Vector from 0->1 (Local X axis = v1)
      v1(0) = node(0,1) - node(0,0);
      v1(1) = node(1,1) - node(1,0);
      v1(2) = node(2,1) - node(2,0);
      v1.normalize();

//... Vector from 3->1
      v2(0) = node(0,1) - node(0,3);
      v2(1) = node(1,1) - node(1,3);
      v2(2) = node(2,1) - node(2,3);
      v2.normalize();

//... Vector from 0->2
      v3(0) = node(0,2) - node(0,0);
      v3(1) = node(1,2) - node(1,0);
      v3(2) = node(2,2) - node(2,0);
      v3.normalize();

//... Define perpendicular as cross between v2 and v3
      v4 = v2.cross(v3);
      v4.normalize();
      len = v4.norm();
      if (len == 0.0) { 
        std::cerr << " ... ERROR: Element has zero area.\n";
      }

//... Cross v4 and v1 (Local Y axis = v2)
      v2 = v4.cross(v1);
      v2.normalize();

//... Cross v1 and v2 (Local Z axis = v3)
      v3 = v1.cross(v2);
      v3.normalize();

//... Compute Local Coordinates: store in "local"

      for(i = 0; i < 4; ++i) {
        x[i] = 0.0;
        y[i] = 0.0;
        z[i] = 0.0;
        zdl[i] = 0.0;
        for(j = 0; j < 3; ++j) {
          local(j,i) = 0.0;
        }
      }
      for(i = 0; i < 8; ++i) {
        vl[i] = 0.0;
      }

//... v1,v2,v3 form transformation matrix
//... local = T*node

      for(i = 0; i < 4; ++i) {
        for(j = 0; j < 3; ++j) {
          local(0,i) = local(0,i) + v1[j]*node(j,i);
          local(1,i) = local(1,i) + v2[j]*node(j,i);
          local(2,i) = local(2,i) + v3[j]*node(j,i);
        }
      }

//... Transform displacements
//... vl = T*v

      for(i = 0; i<4; ++i) {
        for(j = 0; j<3; ++j) {
          vl[2*i]   = vl[2*i]   + v1[j]*v[(3*i)+j];
          vl[2*i+1] = vl[2*i+1] + v2[j]*v[(3*i)+j];
          zdl[i]    = zdl[i]    + v3[j]*v[(3*i)+j];
        }
      }

//... Store Local Coordinates in x,y,z

      for(i = 0; i < 4; ++i) {
        x(i) = local(0,i);
        y(i) = local(1,i);
        z(i) = local(2,i);
      }

//... Check amount of out-of-plane behavior
//    Use x coordinate of node 2 as representative length

      for(i = 0; i<4; ++i) {
        len = abs(local(2,i)/local(0,1));
        if (len > 0.1) {
          std::cerr << " ... Error: Element has warp > 10%\n";
        }
      }

//... Set Parameters for sands2.f

//... Create consitutive matrix relevant to Shear element
//    The only non-zero entry is for the shear c(3,3) = G

      for(i = 0; i<3; ++i) {
        for(j = 0; j<3; ++j) {
          c(i,j) = 0.0;
        }
      }
      c(2,2) = G;

      char escm[7] = "EXTRA";

      seven = 7;
      four  = 4;
      three = 3;
      one   = 1;
      falseflag = false;
      quadstress.setZero();
      quadstrain.setZero();


//... Get Stress for Shear Portion
      this->sands2(escm,x.data(),y.data(),c.data(),vl.data(),quadstress.data(),quadstrain.data(),
                   four,seven,one,one,falseflag,falseflag,0,0,0);


//... Compute Extensional Stresses
//    This is defined as rods along the edges

      avg = 1.0;
      for(i = 0; i<4; ++i) {
        area(i) = 0.0;
      }
      if (F1 > 0) {
        avg = avg+1.0;
        area[0] = 1.;
        area[1] = 1.;
      }
      if (F2 > 0) {
        avg = avg+1.0;
        area[2] = 1.;
        area[3] = 1.;
      }

// Node definition of each rod
// Side 0-1
      rnod(0,0) = 0;
      rnod(1,0) = 1;
// Side 2-3
      rnod(0,1) = 2;
      rnod(1,1) = 3;
// Side 0-3
      rnod(0,2) = 0;
      rnod(1,2) = 3;
// Side 1-2
      rnod(0,3) = 1;
      rnod(1,3) = 2;

      rodstrain.setZero();
/*
// Begin Loop Over 4 Rods
      for(k = 0; k < 4; ++k) {
        if (area[k] != 0.) {

// Get nodes of rod
// and put origin at 1st node
          for(i = 0; i < 2; ++i) {
            xr[i] = x[rnod(i,k)] - x[rnod(0,k)];
            yr[i] = y[rnod(i,k)] - y[rnod(0,k)];
          }

// Length of Rod
          len =  xr[1]*xr[1] + yr[1]*yr[1];
          if (len < 0) {
            std::cerr << " ... Error: Rod has zero or negative length.\n";
            exit(-1);
          }
          len = abs(sqrt(len));

// Get deformations of rod
          for(i = 0; i < 2; ++i) {
            vrx[i] = vl[2*rnod(i,k)];
            vry[i] = vl[2*rnod(i,k)+1];
          }

// Stress Computation copied from sands1.f by P.R. Stern

//.... Calculate the defromation difference

          l1 = xr[1]/len;
          l2 = yr[1]/len;

          q1 = l1*vrx[0] + l2*vry[0];
          q2 = l1*vrx[1] + l2*vry[1];

          dq = q2-q1;

//.... Calculate the Cauchy Strains
          epsrod = dq/len;

//.... Add to strain matrix for each node
//     Transformation Matrix is  T = |  cos  sin |
//                                   | -sin  cos |
//     strain = T^t * | epsrod  0 | *T
//                    |   0     0 |   
// cos is stored in l1
// sin is stored in l1

          rodstrain(0,rnod(0,k)) = rodstrain(0,rnod(0,k)) + epsrod*l1*l1;
          rodstrain(1,rnod(0,k)) = rodstrain(1,rnod(0,k)) + epsrod*l2*l2;
          rodstrain(3,rnod(0,k)) = rodstrain(3,rnod(0,k)) + epsrod*l1*l2;
          rodstrain(0,rnod(1,k)) = rodstrain(0,rnod(1,k)) + epsrod*l1*l1;
          rodstrain(1,rnod(1,k)) = rodstrain(1,rnod(1,k)) + epsrod*l2*l2;
          rodstrain(3,rnod(1,k)) = rodstrain(3,rnod(1,k)) + epsrod*l1*l2;

        }

// End Loop on Rods
      }

// Compute Stresses From Strains
// Add Stresses/Strains to Total from Shear

      for(i = 0; i < 6; ++i) {
        for(j = 0; j < 4; ++j) {
          quadstress(i,j) =  quadstress(i,j) + E*rodstrain(i,j);
          quadstrain(i,j) =  quadstrain(i,j) + rodstrain(i,j);
        }
      }

// Average Stresses
      for(i = 0; i < 6; ++i) {
        for(j = 0; j < 4; ++j) {
          quadstress(i,j) =  quadstress(i,j)/avg;
          quadstrain(i,j) =  quadstrain(i,j)/avg;
        }
      }
*/

// Transform to Global Frame

//     stress = T^t stress_local T

// each part of quadstress holds [sigxx,sigyy,tauxy]
// each part of quadstress holds [sigxx,sigyy,sigzz,tauxy,tauyz,tauxz]

      for(j = 0; j < 4; ++j) {
        stress(0,j) = quadstress(0,j)*v1[0]*v1[0] + quadstress(3,j)*v2[0]*v1[0] + quadstress(3,j)*v1[0]*v2[0] + quadstress(1,j)*v2[0]*v2[0]; // Global XX Stress
        stress(1,j) = quadstress(0,j)*v1[1]*v1[1] + quadstress(3,j)*v2[1]*v1[1] + quadstress(3,j)*v1[1]*v2[1] + quadstress(1,j)*v2[1]*v2[1]; // Global YY Stress
        stress(2,j) = quadstress(0,j)*v1[2]*v1[2] + quadstress(3,j)*v2[2]*v1[2] + quadstress(3,j)*v1[2]*v2[2] + quadstress(1,j)*v2[2]*v2[2]; // Global ZZ Stress
        stress(3,j) = quadstress(0,j)*v1[0]*v1[1] + quadstress(3,j)*v2[0]*v1[1] + quadstress(3,j)*v1[0]*v2[1] + quadstress(1,j)*v2[0]*v2[1]; // Global XY Stress
        stress(4,j) = quadstress(0,j)*v1[1]*v1[2] + quadstress(3,j)*v2[1]*v1[2] + quadstress(3,j)*v1[1]*v2[2] + quadstress(1,j)*v2[1]*v2[2]; // Global YZ Stress
        stress(5,j) = quadstress(0,j)*v1(0)*v1(2) + quadstress(3,j)*v2[0]*v1[2] + quadstress(3,j)*v1[0]*v2[2] + quadstress(1,j)*v2[0]*v2[2]; // Global XZ Stress

        strain(0,j) = quadstrain(0,j)*v1[0]*v1[0] + quadstrain(3,j)*v2[0]*v1[0] + quadstrain(3,j)*v1[0]*v2[0] + quadstrain(1,j)*v2[0]*v2[0]; // Global XX Strain
        strain(1,j) = quadstrain(0,j)*v1[1]*v1[1] + quadstrain(3,j)*v2[1]*v1[1] + quadstrain(3,j)*v1[1]*v2[1] + quadstrain(1,j)*v2[1]*v2[1]; // Global YY Strain
        strain(2,j) = quadstrain(0,j)*v1[2]*v1[2] + quadstrain(3,j)*v2[2]*v1[2] + quadstrain(3,j)*v1[2]*v2[2] + quadstrain(1,j)*v2[2]*v2[2]; // Global ZZ Strain
        strain(3,j) = quadstrain(0,j)*v1[0]*v1[1] + quadstrain(3,j)*v2[0]*v1[1] + quadstrain(3,j)*v1[0]*v2[1] + quadstrain(1,j)*v2[0]*v2[1]; // Global XY Strain
        strain(4,j) = quadstrain(0,j)*v1[1]*v1[2] + quadstrain(3,j)*v2[1]*v1[2] + quadstrain(3,j)*v1[1]*v2[2] + quadstrain(1,j)*v2[1]*v2[2]; // Global YZ Strain
        strain(5,j) = quadstrain(0,j)*v1(0)*v1(2) + quadstrain(3,j)*v2[0]*v1[2] + quadstrain(3,j)*v1[0]*v2[2] + quadstrain(1,j)*v2[0]*v2[2]; // Global XZ Strain
      }

//... Compute the Von Mises Stress

      for(j = 0; j < 4; ++j) {
        vmsstress(6,j) = 0;
        vmsstrain(6,j) = 0;
        for(i = 0; i < 6; ++i) {
          vmsstress(i,j) = stress(i,j);
          vmsstrain(i,j) = strain(i,j);
        }
      }
      this->vmelmv(vmsstress.data(),four,seven,one,one,four);
      this->strainvm(vmsstrain.data(),four,seven,one,four);
      vmssig = vmsstress(6,1);
      vmseps = vmsstrain(6,1);
}

template<typename doublereal>
void
ShearPanelTemplate<doublereal>
::vmssWRTdisp(doublereal *_xg, doublereal *_yg, doublereal *_zg, doublereal *_v,
              doublereal G, doublereal E, doublereal *_vmsWRTdisp, doublereal &vmssig)
{ 
/*
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     vmssWRTdisp calculates sensitivity of von Mises stress
C     wrt displacements on four-node shear panel element.
C
C       xg        (4 x 1) array of global x coordinates
C       yg        (4 x 1) array of global y coordinates
C       zg        (4 x 1) array of global y coordinates
C       v         (12x 1) array of node dispalcements
C                 arranged, vx1,vy1,vz1,...,vz4
C       G         Shear Modulus of element
C       E         Youngs Modulus of element
C
C     The outputs are:
C
C       stress    (6x4) corner node stresses arranged
C                  (sigxx, sigyy, sigzz, tauxy, tauyz, tauxz)
C       vmssig    Vommises stress
*/
      using std::abs;
      using std::sqrt;

      Eigen::Map<Eigen::Matrix<doublereal,4,12> > vmsWRTdisp(_vmsWRTdisp);
      Eigen::Map<Eigen::Matrix<doublereal,4,1> > xg(_xg), yg(_yg), zg(_zg);
      Eigen::Map<Eigen::Matrix<doublereal,12,1> > v(_v);

//                   L O C A L   V A R I A B L E S

// Local Coordinates
      Eigen::Matrix<doublereal,4,1> x,y,z;
      Eigen::Matrix<doublereal,8,1> vl;
// Working Vectors
      Eigen::Matrix<doublereal,3,4> node, local;
      Eigen::Matrix<doublereal,3,1> v1, v2,v3,v4;
      doublereal len;
// Consitutive Matrix
      Eigen::Matrix<doublereal,3,3> c;
// For sands2 routine
      Eigen::Matrix<doublereal,7,4> quadstress, quadstrain;
      int four, three, one, seven;
      bool falseflag;
// Effective Extentional Areas
      Eigen::Matrix<doublereal,4,1> area;
// Rod Coordinates and Stiffness
      doublereal epsrod;
      Eigen::Matrix<doublereal,2,1> xr, yr, vrx, vry;
      doublereal l1,l2,q1,q2,dq;
      Eigen::Matrix<doublereal,6,4> stress;
// For vms routine
      Eigen::Matrix<doublereal,7,4> vmsstress, vmsstrain;
// Output Control
      bool output;
// Other Variables
      int i, j, k;
      doublereal avg;

// Set Output Status
      output = false;

/*
C     COMPUTE LOCAL COORDINATE SYSTEM 
C  (works best if element is plane)
C
C      Y                             
C     /|\                            
C      |  3            2             
C      |  /------------|             
C      | /             |             
C      |/              |             
C      /---------------|-----> X     
C      0               1             
C
C... store nodal coordinates in "node"
C... shift origin to node 0
*/
      for(i = 0; i<4; ++i) {
        node(0,i) = xg(i)-xg(0);
        node(1,i) = yg(i)-yg(0);
        node(2,i) = zg(i)-zg(0);
      }

      Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, " ");
//... Vector from 0->1 (Local X axis = v1)
      v1(0) = node(0,1) - node(0,0);
      v1(1) = node(1,1) - node(1,0);
      v1(2) = node(2,1) - node(2,0);
      v1.normalize();

//... Vector from 3->1
      v2(0) = node(0,1) - node(0,3);
      v2(1) = node(1,1) - node(1,3);
      v2(2) = node(2,1) - node(2,3);
      v2.normalize();

//... Vector from 0->2
      v3(0) = node(0,2) - node(0,0);
      v3(1) = node(1,2) - node(1,0);
      v3(2) = node(2,2) - node(2,0);
      v3.normalize();

//... Define perpendicular as cross between v2 and v3
      v4 = v2.cross(v3);
      v4.normalize();
      len = v4.norm();
      if (len == 0.0) { 
        std::cerr << " ... ERROR: Element has zero area.\n";
      }


//... Cross v4 and v1 (Local Y axis = v2)
      v2 = v4.cross(v1);
      v2.normalize();

//... Cross v1 and v2 (Local Z axis = v3)
      v3 = v1.cross(v2);
      v3.normalize();

      std::cerr << "printing v1, v2, v3\n" << std::endl;
      std::cerr << v1 << std::endl << v2 << std::endl << v3 << std::endl;

//... Compute Local Coordinates: store in "local"

      for(i = 0; i < 4; ++i) {
        x[i] = 0.0;
        y[i] = 0.0;
        z[i] = 0.0;
        for(j = 0; j < 3; ++j) {
          local(j,i) = 0.0;
        }
      }
      for(i = 0; i < 8; ++i) {
        vl[i] = 0.0;
      }

//... v1,v2,v3 form transformation matrix
//... local = T*node

      for(i = 0; i < 4; ++i) {
        for(j = 0; j < 3; ++j) {
          local(0,i) = local(0,i) + v1[j]*node(j,i);
          local(1,i) = local(1,i) + v2[j]*node(j,i);
          local(2,i) = local(2,i) + v3[j]*node(j,i);
        }
      }

//... Transform displacements
//... vl = T*v

      for(i = 0; i<4; ++i) {
        for(j = 0; j<3; ++j) {
          vl[2*i]   = vl[2*i]   + v1[j]*v[(3*i)+j];
          vl[2*i+1] = vl[2*i+1] + v2[j]*v[(3*i)+j];
        }
      }
      std::cerr << "printint vl\n" << vl << std::endl;

//... Construct T
      Eigen::Matrix<doublereal, 8, 12> T;
      T.setZero();
      Eigen::Matrix<doublereal, 1, 3> zero13;
      zero13.setZero();
      T << v1.transpose(), zero13, zero13, zero13,
           v2.transpose(), zero13, zero13, zero13,
           zero13, v1.transpose(), zero13, zero13,
           zero13, v2.transpose(), zero13, zero13,
           zero13, zero13, v1.transpose(), zero13,
           zero13, zero13, v2.transpose(), zero13,
           zero13, zero13, zero13, v1.transpose(),
           zero13, zero13, zero13, v2.transpose();

      std::cerr << "printing T\n" << T << std::endl;

//... Store Local Coordinates in x,y,z

      for(i = 0; i < 4; ++i) {
        x(i) = local(0,i);
        y(i) = local(1,i);
        z(i) = local(2,i);
      }

//... Check amount of out-of-plane behavior
//    Use x coordinate of node 2 as representative length

      for(i = 0; i<4; ++i) {
        len = abs(local(2,i)/local(0,1));
        if (len > 0.1) {
          std::cerr << " ... Error: Element has warp > 10%\n";
        }
      }

//... Set Parameters for sands2.f

//... Create consitutive matrix relevant to Shear element
//    The only non-zero entry is for the shear c(3,3) = G

      for(i = 0; i<3; ++i) {
        for(j = 0; j<3; ++j) {
          c(i,j) = 0.0;
        }
      }
      c(2,2) = G;

      char escm[7] = "EXTRA";

      seven = 7;
      four  = 4;
      three = 3;
      one   = 1;
      falseflag = false;
      quadstress.setZero();
      quadstrain.setZero();

//... Get Stress for Shear Portion
      this->sands2(escm,x.data(),y.data(),c.data(),vl.data(),quadstress.data(),quadstrain.data(),
                   four,seven,one,one,falseflag,falseflag,0,0,0);
      Eigen::Matrix<doublereal,12,8> DquadstressDdisp;
      Eigen::Matrix<doublereal,4,8> DvmsDdisp2;
      DquadstressDdisp.setZero();
      DvmsDdisp2.setZero();
      this->vms2WRTdisp(escm,x.data(),y.data(),c.data(),vl.data(),DvmsDdisp2.data(),DquadstressDdisp.data(),
                        four,seven,one,one,true,true,0,0,0);
      std::cerr << "printing DvmsDdisp2\n" << DvmsDdisp2 << std::endl;
      std::cerr << "printing DquadstressDdisp\n" << DquadstressDdisp << std::endl;

// Transform to Global Frame

//     stress = T^t stress_local T

// each part of quadstress holds [sigxx,sigyy,tauxy]
// each part of quadstress holds [sigxx,sigyy,sigzz,tauxy,tauyz,tauxz]

      for(j = 0; j < 4; ++j) {
        stress(0,j) = quadstress(0,j)*v1[0]*v1[0] + quadstress(3,j)*v2[0]*v1[0] + quadstress(3,j)*v1[0]*v2[0] + quadstress(1,j)*v2[0]*v2[0]; // Global XX Stress
        stress(1,j) = quadstress(0,j)*v1[1]*v1[1] + quadstress(3,j)*v2[1]*v1[1] + quadstress(3,j)*v1[1]*v2[1] + quadstress(1,j)*v2[1]*v2[1]; // Global YY Stress
        stress(2,j) = quadstress(0,j)*v1[2]*v1[2] + quadstress(3,j)*v2[2]*v1[2] + quadstress(3,j)*v1[2]*v2[2] + quadstress(1,j)*v2[2]*v2[2]; // Global ZZ Stress
        stress(3,j) = quadstress(0,j)*v1[0]*v1[1] + quadstress(3,j)*v2[0]*v1[1] + quadstress(3,j)*v1[0]*v2[1] + quadstress(1,j)*v2[0]*v2[1]; // Global XY Stress
        stress(4,j) = quadstress(0,j)*v1[1]*v1[2] + quadstress(3,j)*v2[1]*v1[2] + quadstress(3,j)*v1[1]*v2[2] + quadstress(1,j)*v2[1]*v2[2]; // Global YZ Stress
        stress(5,j) = quadstress(0,j)*v1(0)*v1(2) + quadstress(3,j)*v2[0]*v1[2] + quadstress(3,j)*v1[0]*v2[2] + quadstress(1,j)*v2[0]*v2[2]; // Global XZ Stress
      }
 
      std::cerr << "printing stress\n" << stress << std::endl;

      Eigen::Matrix<doublereal, 24, 12> DstressDquadstress;
      DstressDquadstress.setZero(); 
      for(j = 0; j < 4; ++j) {
        DstressDquadstress(6*j,3*j) = v1[0]*v1[0];
        DstressDquadstress(6*j,3*j+1) = v2[0]*v2[0];
        DstressDquadstress(6*j,3*j+2) = 2*v2[0]*v1[0]; //

        DstressDquadstress(6*j+1,3*j) = v1[1]*v1[1];
        DstressDquadstress(6*j+1,3*j+1) = v2[1]*v2[1];
        DstressDquadstress(6*j+1,3*j+2) = 2*v2[1]*v1[1]; //

        DstressDquadstress(6*j+2,3*j) = v1[2]*v1[2];
        DstressDquadstress(6*j+2,3*j+1) = v2[2]*v2[2];
        DstressDquadstress(6*j+2,3*j+2) = 2*v2[2]*v1[2]; //

        DstressDquadstress(6*j+3,3*j) = v1[0]*v1[1];
        DstressDquadstress(6*j+3,3*j+1) = v2[0]*v2[1];
        DstressDquadstress(6*j+3,3*j+2) = v2[0]*v1[1] + v1[0]*v2[1]; //

        DstressDquadstress(6*j+4,3*j) = v1[1]*v1[2];
        DstressDquadstress(6*j+4,3*j+1) = v2[1]*v2[2];
        DstressDquadstress(6*j+4,3*j+2) = v2[1]*v1[2] + v1[1]*v2[2]; //

        DstressDquadstress(6*j+5,3*j) = v1[0]*v1[2];
        DstressDquadstress(6*j+5,3*j+1) = v2[0]*v2[2];
        DstressDquadstress(6*j+5,3*j+2) = v2[0]*v1[2] + v1[0]*v2[2]; //
      }
/*
      DstressDquadstress(0,0) = v1[0]*v1[0];
      DstressDquadstress(0,1) = v2[0]*v2[0];
      DstressDquadstress(0,2) = 2*v2[0]*v1[0]; //

      DstressDquadstress(1,0) = v1[1]*v1[1];
      DstressDquadstress(1,1) = v2[1]*v2[1];
      DstressDquadstress(1,2) = 2*v2[1]*v1[1]; //

      DstressDquadstress(2,0) = v1[2]*v1[2];
      DstressDquadstress(2,1) = v2[2]*v2[2];
      DstressDquadstress(2,2) = 2*v2[2]*v1[2]; //

      DstressDquadstress(3,0) = v1[0]*v1[1];
      DstressDquadstress(3,1) = v2[0]*v2[1];
      DstressDquadstress(3,2) = v2[0]*v1[1] + v1[0]*v2[1]; //

      DstressDquadstress(4,0) = v1[1]*v1[2];
      DstressDquadstress(4,1) = v2[1]*v2[2];
      DstressDquadstress(4,2) = v2[1]*v1[2] + v1[1]*v2[2]; //

      DstressDquadstress(5,0) = v1[0]*v1[2];
      DstressDquadstress(5,1) = v2[0]*v2[2];
      DstressDquadstress(5,2) = v2[0]*v1[2] + v1[0]*v2[2]; //
*/
      std::cerr << "printing DstressDquadstress\n" << DstressDquadstress << std::endl;
//... Compute the Von Mises Stress

      for(j = 0; j < 4; ++j) {
        vmsstress(6,j) = 0;
        for(i = 0; i < 6; ++i) {
          vmsstress(i,j) = stress(i,j);
        }
      }
      this->vmelmv(vmsstress.data(),four,seven,one,one,four);
      vmssig = vmsstress(6,0);
      std::cerr << vmsstress(6,0) << " "  << vmsstress(6,1) <<  " " << vmsstress(6,2) << " " << vmsstress(6,3) << std::endl;

      Eigen::Matrix<doublereal,4,24> dvmsdStress;
      dvmsdStress.setZero();
/*      dvmsdStress(3,0) = dvmsdStress(2,0) = dvmsdStress(1,0) = dvmsdStress(0,0) = (2.*stress(0,0)-stress(1,0)-stress(2,0))/(2.*vmsstress(6,0));
      dvmsdStress(3,1) = dvmsdStress(2,1) = dvmsdStress(1,1) = dvmsdStress(0,1) = (2.*stress(1,0)-stress(0,0)-stress(2,0))/(2.*vmsstress(6,0));
      dvmsdStress(3,2) = dvmsdStress(2,2) = dvmsdStress(1,2) = dvmsdStress(0,2) = (2.*stress(2,0)-stress(1,0)-stress(0,0))/(2.*vmsstress(6,0));
      dvmsdStress(3,3) = dvmsdStress(2,3) = dvmsdStress(1,3) = dvmsdStress(0,3) = (3.*stress(3,0))/vmsstress(6,0);
      dvmsdStress(3,4) = dvmsdStress(2,4) = dvmsdStress(1,4) = dvmsdStress(0,4) = (3.*stress(4,0))/vmsstress(6,0);
      dvmsdStress(3,5) = dvmsdStress(2,5) = dvmsdStress(1,5) = dvmsdStress(0,5) = (3.*stress(5,0))/vmsstress(6,0); */
      for(int n=0; n<4; ++n) {
        dvmsdStress(n,6*n) = (2.*stress(0,n)-stress(1,n)-stress(2,n))/(2.*vmsstress(6,n));
        dvmsdStress(n,6*n+1) = (2.*stress(1,n)-stress(0,n)-stress(2,n))/(2.*vmsstress(6,n));
        dvmsdStress(n,6*n+2) = (2.*stress(2,n)-stress(1,n)-stress(0,n))/(2.*vmsstress(6,n));
        dvmsdStress(n,6*n+3) = (3.*stress(3,n))/vmsstress(6,n);
        dvmsdStress(n,6*n+4) = (3.*stress(4,n))/vmsstress(6,n);
        dvmsdStress(n,6*n+5) = (3.*stress(5,n))/vmsstress(6,n);
      }

      std::cerr << "printing dvmsdStress\n" << dvmsdStress << std::endl;

      vmsWRTdisp = (dvmsdStress * DstressDquadstress) * (DquadstressDdisp * T);
      std::cerr << "printing vmsWRTdisp\n" << vmsWRTdisp << std::endl;
 
}

#endif
#endif
