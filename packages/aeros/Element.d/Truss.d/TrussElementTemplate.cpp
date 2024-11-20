#ifdef USE_EIGEN3
#ifndef _TRUSSELEMENTSEMITEMPLATE_CPP_
#define _TRUSSELEMENTSEMITEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/Truss.d/TrussElementTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal>
void
TrussElementTemplate<doublereal>
::stiffnessMatrix(doublereal *_stiffness, doublereal *_x,doublereal *_y, doublereal *_z,
                  doublereal E, doublereal A, doublereal preload)
{
       using std::sqrt;

       Eigen::Map<Eigen::Matrix<doublereal,2,1> > x(_x), y(_y), z(_z);
       Eigen::Map<Eigen::Matrix<doublereal,6,6> > stiffness(_stiffness);

       doublereal dx = x[1] - x[0];
       doublereal dy = y[1] - y[0];
       doublereal dz = z[1] - z[0];

       doublereal length = sqrt( dx*dx + dy*dy + dz*dz );
       doublereal c1[3];

       c1[0] = dx/length;
       c1[1] = dy/length;
       c1[2] = dz/length;

       doublereal elementK = E*A/length;
       int i,j;
       for(i=0; i < 3; ++i) {
         for(j=0; j < 3; ++j) {
            stiffness(i,j)     = elementK*c1[i]*c1[j];
            stiffness(i+3,j+3) = elementK*c1[i]*c1[j];
            stiffness(i+3,j)   = -stiffness(i,j);
            stiffness(i,j+3)   = -stiffness(i,j);
         }
       }

       if(preload != 0.0) {
         for(i=0; i < 3; ++i) {
            stiffness(i,i)     += preload/length;
            stiffness(i+3,i+3) += preload/length;
            stiffness(i+3,i)   = -stiffness(i,i);
            stiffness(i,i+3)   = -stiffness(i,i);
          }
       }
}


template<typename doublereal>
void
TrussElementTemplate<doublereal>
::vmst(doublereal *_stress, doublereal *_x,doublereal *_y, doublereal *_z, doublereal *_elDisp,
       doublereal E, doublereal A, doublereal W, doublereal Ta, doublereal preload, 
       int strInd, int surface, doublereal *_ndTemps, doublereal ylayer, doublereal zlayer, int avgnum)
{
   using std::sqrt;

#ifndef SALINAS
   if(strInd != 6) {
     std::cerr << " ... Error: strInd is not 6. exiting!\n";
     exit(-1);
   }

   Eigen::Map<Eigen::Matrix<doublereal,2,1> > x(_x), y(_y), z(_z);
   Eigen::Map<Eigen::Matrix<doublereal,2,1> > stress(_stress);
   Eigen::Map<Eigen::Matrix<doublereal,2,1> > ndTemps(_ndTemps);
   Eigen::Map<Eigen::Matrix<doublereal,6,1> > elDisp(_elDisp);

   doublereal dx = x[1] - x[0];
   doublereal dy = y[1] - y[0];
   doublereal dz = z[1] - z[0];

   doublereal length = sqrt(dx*dx + dy*dy + dz*dz);

   // scale dx, dy, and dz by the length
   dx /= length;
   dy /= length;
   dz /= length;

   // Compute the change in length of the element
   doublereal dq = dx*(elDisp[3]-elDisp[0])
                 + dy*(elDisp[4]-elDisp[1])
                 + dz*(elDisp[5]-elDisp[2]);

   // Compute axial strain
   doublereal exx = dq/length;

    switch (avgnum) {

      case 0:
      {
           // Compute axial force
           doublereal f = A*E*exx;
	
           // Add Preload
           f  += preload;

           // Compute thermal force
           doublereal coefficient = E*A*W;
           doublereal Tref = Ta;
 
           doublereal fth1(0); 
           doublereal fth2(0); 
           if(_ndTemps) {
              fth1 = coefficient*(ndTemps[0]-Tref);
              fth2 = coefficient*(ndTemps[1]-Tref);
           }

           // compute stresses 
	         doublereal elForce[2]={0.0,0.0};
           elForce[0] = -f + fth1;
           elForce[1] =  f - fth2;
           stress[0] = sqrt(elForce[0]*elForce[0]);
           stress[1] = sqrt(elForce[1]*elForce[1]);
           break;
      }

      case 1:
      { 
        // Compute von Mises stress resultant
        doublereal f = A*E*exx;
        f += preload;

        // Compute thermal force
        doublereal coefficient = E*A*W;
        doublereal Tref = Ta;
        doublereal fth1(0);
        doublereal fth2(0);
        if(_ndTemps) {
          fth1 = coefficient*(ndTemps[0]-Tref);
          fth2 = coefficient*(ndTemps[1]-Tref);
        }
        stress[0] = sqrt((-f + fth1)*(-f + fth1));
        stress[1] = sqrt(( f - fth1)*( f - fth1));      

        break;
      }

      case 2:
      {
        stress[0] = 0.0;
        stress[1] = 0.0;
        break;
      }

      default:
        std::cerr << "avgnum = " << avgnum << " is not a valid number\n";
    }
#endif
}

#endif
#endif
