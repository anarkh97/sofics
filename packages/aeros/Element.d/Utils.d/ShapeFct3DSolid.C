//#ifdef LOOP_UNROLL_4

// HB (05-21-05): for given reference shape fct derivatives dShape & nodes coordinates (X.Y,Z), 
//                compute the jacobian & shape fcts derivative w.r.t to the global coordinate system.
//                assuming an iso-parametric formulation.
double
computeDShapeFct3DSolid(double (*dShape)[3], int nnodes, double* X, double* Y, double* Z, double (*DShape)[3]=0)
{
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0; xd2 = 0.0; xd3 = 0.0;
  yd1 = 0.0; yd2 = 0.0; yd3 = 0.0;
  zd1 = 0.0; zd2 = 0.0; zd3 = 0.0;

#ifdef LOOP_UNROLL_4
  fprintf(stderr," *** ERROR: In computeDShapeFct3DSolid.C: Loop enrolling implementation NOT supported. Abort.\n");
  exit(-1);
  int nloops = nnodes/4;
  int i;
  for(int i=0; i<nnodes; i+=4) {
    xd1 += dShape[i][0]*X[i] + dShape[i+1][0]*X[i+1] + dShape[i+2][0]*X[i+2] + dShape[i+3][0]*X[i+3];
    xd2 += dShape[i][1]*X[i] + dShape[i+1][1]*X[i+1] + dShape[i+2][1]*X[i+2] + dShape[i+3][1]*X[i+3];
    xd3 += dShape[i][2]*X[i] + dShape[i+1][2]*X[i+1] + dShape[i+2][2]*X[i+2] + dShape[i+3][2]*X[i+3];
    yd1 += dShape[i][0]*Y[i] + dShape[i+1][0]*Y[i+1] + dShape[i+2][0]*Y[i+2] + dShape[i+3][0]*Y[i+3];
    yd2 += dShape[i][1]*Y[i] + dShape[i+1][1]*Y[i+1] + dShape[i+2][1]*Y[i+2] + dShape[i+3][1]*Y[i+3];
    yd3 += dShape[i][2]*Y[i] + dShape[i+1][2]*Y[i+1] + dShape[i+2][2]*Y[i+2] + dShape[i+3][2]*Y[i+3];
    zd1 += dShape[i][0]*Z[i] + dShape[i+1][0]*Z[i+1] + dShape[i+2][0]*Z[i+2] + dShape[i+3][0]*Z[i+3];
    zd2 += dShape[i][1]*Z[i] + dShape[i+1][1]*Z[i+1] + dShape[i+2][1]*Z[i+2] + dShape[i+3][1]*Z[i+3];
    zd3 += dShape[i][2]*Z[i] + dShape[i+1][2]*Z[i+1] + dShape[i+2][2]*Z[i+2] + dShape[i+3][2]*Z[i+3];
  }
  for(int i-=4; i<nnodes; i++) {
    xd1 += dShape[i][0]*X[i]; xd2 += dShape[i][1]*X[i]; xd3 += dShape[i][2]*X[i];
    yd1 += dShape[i][0]*Y[i]; yd2 += dShape[i][1]*Y[i]; yd3 += dShape[i][2]*Y[i];
    zd1 += dShape[i][0]*Z[i]; zd2 += dShape[i][1]*Z[i]; zd3 += dShape[i][2]*Z[i];
  }
#else
   for(int i=0; i<nnodes; i++) {
    xd1 += dShape[i][0]*X[i]; xd2 += dShape[i][1]*X[i]; xd3 += dShape[i][2]*X[i];
    yd1 += dShape[i][0]*Y[i]; yd2 += dShape[i][1]*Y[i]; yd3 += dShape[i][2]*Y[i];
    zd1 += dShape[i][0]*Z[i]; zd2 += dShape[i][1]*Z[i]; zd3 += dShape[i][2]*Z[i];
  }
#endif  

  double a11 = yd2*zd3 - yd3*zd2; double a12 = yd3*zd1 - yd1*zd3; double a13 = yd1*zd2 - yd2*zd1;
  double a21 = xd3*zd2 - xd2*zd3; double a22 = xd1*zd3 - zd1*xd3; double a23 = xd2*zd1 - xd1*zd2;
  double a31 = xd2*yd3 - xd3*yd2; double a32 = yd1*xd3 - xd1*yd3; double a33 = xd1*yd2 - yd1*xd2;
                                                                                                                             
  /* -------> DETERMINANT OF THE JACOBIAN MATRIX <--- */
                                                                                                                             
  double J = xd1*a11 + yd1*a21 + zd1*a31;
  if(J == 0.0) {
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computeDShapeFct3DSolid.C ROUTINE.\n");
  }
  if(J < 0.0) {
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computeDShapeFct3DSolid.C ROUTINE.\n");
  }
  double cdet = 1.0/J;
                                                                                                                             
  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- */
  if(DShape){
#ifdef LOOP_UNROLL_4
    fprintf(stderr," *** ERROR: In computeDShapeFct3DSolid.C: Loop enrolling implementation NOT supported. Abort.\n");
    exit(-1);
    int nloops = nnodes/4;
    for(int i=0; i<nloops; i+=4) {
      DShape[i  ][0] = cdet * ( a11*dShape[i  ][0] + a12*dShape[i  ][1] + a13*dShape[i  ][2] );
      DShape[i  ][1] = cdet * ( a21*dShape[i  ][0] + a22*dShape[i  ][1] + a23*dShape[i  ][2] );
      DShape[i  ][2] = cdet * ( a31*dShape[i  ][0] + a32*dShape[i  ][1] + a33*dShape[i  ][2] );
                                                                                                            
      DShape[i+1][0] = cdet * ( a11*dShape[i+1][0] + a12*dShape[i+1][1] + a13*dShape[i+1][2] );
      DShape[i+1][1] = cdet * ( a21*dShape[i+1][0] + a22*dShape[i+1][1] + a23*dShape[i+1][2] );
      DShape[i+1][2] = cdet * ( a31*dShape[i+1][0] + a32*dShape[i+1][1] + a33*dShape[i+1][2] );
                                                                                                            
      DShape[i+2][0] = cdet * ( a11*dShape[i+2][0] + a12*dShape[i+2][1] + a13*dShape[i+2][2] );
      DShape[i+2][1] = cdet * ( a21*dShape[i+2][0] + a22*dShape[i+2][1] + a23*dShape[i+2][2] );
      DShape[i+2][2] = cdet * ( a31*dShape[i+2][0] + a32*dShape[i+2][1] + a33*dShape[i+2][2] );
                                                                                                            
      DShape[i+3][0] = cdet * ( a11*dShape[i+3][0] + a12*dShape[i+3][1] + a13*dShape[i+3][2] );
      DShape[i+3][1] = cdet * ( a21*dShape[i+3][0] + a22*dShape[i+3][1] + a23*dShape[i+3][2] );
      DShape[i+3][2] = cdet * ( a31*dShape[i+3][0] + a32*dShape[i+3][1] + a33*dShape[i+3][2] );
    } 
#else
    for(int i=0; i<nnodes; i++) {
      DShape[i  ][0] = cdet * ( a11*dShape[i  ][0] + a12*dShape[i  ][1] + a13*dShape[i  ][2] );
      DShape[i  ][1] = cdet * ( a21*dShape[i  ][0] + a22*dShape[i  ][1] + a23*dShape[i  ][2] );
      DShape[i  ][2] = cdet * ( a31*dShape[i  ][0] + a32*dShape[i  ][1] + a33*dShape[i  ][2] );
     }
#endif  
  }
  return(J);
}
