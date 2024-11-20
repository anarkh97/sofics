#include <Corotational.d/utilities.h>
#include <Math.d/FullSquareMatrix.h>

void tran_stiff(  FullSquareMatrix &stiff, double t0n[3][3])
/*****************************************************************
 *
 *  Purpose:
 *     Transform the stiffness matrix to global system
 *     or local predefined system if this is defined
 *
 *  Method:           
 *            stiff_global = T'*stiff*T
 *
 *  Input:
 *     stiff : local coordinate stiffness matrix
 *     t0n   : transformation matrix T to global system
 *
 *  Output:
 *     stiff : now with respect to global coordinate system
 *
 *  Coded by: Kendall H. Pierson
 *******************************************************************/
{
	double prod1[12];
	int i,j,k;
//
// premultiplication by rotation matrix
//
	
	for(j=0; j<12; ++j) {
	  for(k=0; k<3; ++k) {
            prod1[k]   = 0.0;
            prod1[k+3] = 0.0;
            prod1[k+6] = 0.0;
            prod1[k+9] = 0.0;
	    for(i=0; i<3; ++i) {
              prod1[k]   += t0n[k][i]*stiff[i][j];
              prod1[k+3] += t0n[k][i]*stiff[i+3][j];
              prod1[k+6] += t0n[k][i]*stiff[i+6][j];
              prod1[k+9] += t0n[k][i]*stiff[i+9][j];
	    }  
	  }
	  for(k=0; k<12; ++k)
	    stiff[k][j] = prod1[k];
	}
//
// postmultiplication by rotation matrix transposed
//
	for(j=0; j<12; ++j) {
	  for(k=0; k<3; ++k) {
            prod1[k]   = 0.0;
            prod1[k+3] = 0.0;
            prod1[k+6] = 0.0;
            prod1[k+9] = 0.0;
	    for(i=0; i<3; ++i) {
              prod1[k]   += stiff[j][i]*t0n[k][i];
              prod1[k+3] += stiff[j][i+3]*t0n[k][i];
              prod1[k+6] += stiff[j][i+6]*t0n[k][i];
              prod1[k+9] += stiff[j][i+9]*t0n[k][i];
	    }  
	  }
	  for(k=0; k<12; ++k)
	    stiff[j][k] = prod1[k];
	}

}
