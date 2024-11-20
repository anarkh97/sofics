#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Utils.d/linkfc.h> 
#include <Utils.d/MyComplex.h>

/**************************************************
 *  Function PSTRESS                              *
 *                                                *
 *  Computes 3 principal stresses from 3x3        *
 *  full stress matrix.                           *
 *  Stress matrix is assumed symmetric and        *
 *  function takes vector of 6 independent        *
 *  components as input, in the following order   *
 *  [sigxx,sigyy,sigzz,tauxy,tauyz,tauxz]         *
 *                                                *
 *  Output is vector of length 3 holding          *
 *  [sig1,sig2,sig3]                              *
 *  where sig1 >= sig2 >= sig3                    *
 *                                                *
 *  //Direct computation of principals            *
 *    Uses FORTRAN routine dsyev.f
 *                                                *
 *  Gregory W. Brown                              *
 *  June, 2000                                    *
 **************************************************/

template<class Scalar>
void pstress(Scalar stress[6], Scalar pstress[3], Scalar *pdirections)
{
  //  ... BUILD STRESS MATRIX
  Scalar smatrix[3][3]; // stress matrix
  smatrix[0][0] = stress[0];
  smatrix[1][1] = stress[1];
  smatrix[2][2] = stress[2];
  smatrix[0][1] = stress[3];
  smatrix[1][2] = stress[4];
  smatrix[0][2] = stress[5];
  smatrix[1][0] = smatrix[0][1];
  smatrix[2][1] = smatrix[1][2];
  smatrix[2][0] = smatrix[0][2];

  //  ... GET EIGEN VALUES USING DSYEV
  Scalar singval[3];
  Scalar work[24];
  int info = 0;

  char jobz = (pdirections) ? 'V' : 'N';  // PJSA 3-24-05: flag to compute eigenvectors if required
#if defined(_OPENMP)
  // note: some implementations of dsyev are apparently not threadsafe
  #pragma omp critical
#endif
  Tdsyev(jobz, 'U', 3, (Scalar*)smatrix, 3, singval, work, 24, info);

  //  ... CHECK RETURN STATUS OF DSVDC
  if(info <  0)
    fprintf(stderr," *** ERROR: Illegal value in argument #%d.\n",info);
  if(info >  0)
    fprintf(stderr," *** ERROR: Eigenvalues did not converge with %d iter.\n",info);

  // eigenvalues are sorted in increasing order, so need to reverse direction
  pstress[0] = singval[2];
  pstress[1] = singval[1];
  pstress[2] = singval[0];
  if(pdirections) {
    for(int i=0;i<3;++i) 
      for(int j=0;j<3;++j) pdirections[3*i+j] = smatrix[2-i][j];
  }
/*
  for (int i=0;i<3;++i) {
    std::cerr << "principal stress/strain #" << i+1 << " = " << pstress[i] << std::endl;
    if(pdirections) std::cerr << "principal direction #" << i+1 << " = " << pdirections[3*i+0] 
                              << " " << pdirections[3*i+1] << " " << pdirections[3*i+2] << std::endl;
  }
  std::cerr << std::endl;
*/
  return;
}
