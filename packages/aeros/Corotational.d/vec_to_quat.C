#include <cmath>
#include <Corotational.d/utilities.h>
#include <iostream>
#include <limits>

const double epsilon = std::numeric_limits<double>::epsilon(); //1.0e-15;

typedef double Quat[4];
int EM_To_Q(double v[3], Quat q, int reparam);

void vec_to_quat( double rvec[3], double q[4] )
/*****************************************************************
 *  Compute the quaterion representation from a rotation vector
 *
 *  Input:
 *  rvec:    rotation vector
 *
 *  Output:
 *  q :      quaterion representation of the rotation ( Euler parameters )
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
      double coef;

      // compute norm of rotation vector rvec
      double th = std::sqrt( rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2] );

      q[0] = std::cos( 0.5*th );

      if ( th < epsilon )
        coef = 0.5;
      else
        coef = std::sin(0.5*th)/th;

      q[1] = coef*rvec[0];
      q[2] = coef*rvec[1];
      q[3] = coef*rvec[2];
         
      // Compute norm of quaterion q
      double norm = std::sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
      if(norm == 0) {
        std::cerr << "WARNING: Skipping quaternion normalize( null )\n";
      }
      else {
        // Normalize quaterion q
        q[0] /= norm;
        q[1] /= norm;
        q[2] /= norm;
        q[3] /= norm;
      }

      //if(q[0] < 0.0) std::cerr << "WARNING: negative quaternion in vec_to_quat\n";
      //if(q[0] < 0.0) for(int l=0; l<4; ++l) q[l] = -q[l]; // PJSA DEBUG

/*
  Quat quat;
  EM_To_Q(rvec, quat, 1);
  //std::cerr << "2. q = " << quat[0] << " " << quat[1] << " " << quat[2] << " " << quat[3] << std::endl;
  q[0] = quat[3];
  q[1] = quat[0];
  q[2] = quat[1];
  q[3] = quat[2];
*/  
}

void dvec_to_quat( double rvec[3],double drvec[3], double q[4], double dq[4] )
/*****************************************************************
 *  Compute the quaterion representation from a rotation vector
 *
 *  Input:
 *  rvec:    rotation vector
 *
 *  Output:
 *  q :      quaterion representation of the rotation ( Euler parameters )
 *
 *  Coded by Bjorn Haugen
 *****************************************************************/
{
      double  coef;
      double dcoef;
      double dth;
      
      // compute norm of rotation vector rvec
      double th = sqrt( rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2] );
      if (th > 0){
         dth= (rvec[0]*drvec[0] + rvec[1]*drvec[1] + rvec[2]*drvec[2])/th;
	}
      else{
         dth=0;
         //dth = sqrt( drvec[0]*drvec[0] + 
	 //            drvec[1]*drvec[1] + drvec[2]*drvec[2] );
      }           
      
      q[0] = cos( 0.5*th );
      dq[0]=-0.5*dth*sin(0.5*th);

      if ( th < epsilon ){
        coef = 0.5;
       dcoef = 0.0;}
      else {
        coef = sin(0.5*th)/th;
        dcoef= 0.5*dth*cos(0.5*th)/th-dth*sin(0.5*th)/(th*th);}

      q[1] = coef*rvec[0];
     dq[1] = dcoef*rvec[0]+coef*drvec[0]; 
      q[2] = coef*rvec[1];
     dq[2] = dcoef*rvec[1]+coef*drvec[1];
      q[3] = coef*rvec[2];
     dq[3] = dcoef*rvec[2]+coef*drvec[2];  
     
      // Compute norm of quaterion q
     double norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
     double dnorm =( dq[0]*q[0] + dq[1]*q[1] + dq[2]*q[2] + dq[3]*q[3])/norm;

      // Normalize quaterion q and dq
      
      dq[0]= (dq[0]*norm-q[0]*dnorm)/(norm*norm);
      dq[1]= (dq[1]*norm-q[1]*dnorm)/(norm*norm);
      dq[2]= (dq[2]*norm-q[2]*dnorm)/(norm*norm);
      dq[3]= (dq[3]*norm-q[3]*dnorm)/(norm*norm);
      
      q[0] /= norm;     
      q[1] /= norm;      
      q[2] /= norm;     
      q[3] /= norm;
}
