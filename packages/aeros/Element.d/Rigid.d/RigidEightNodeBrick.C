#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidEightNodeBrick.h>
#include <Element.d/Rigid.d/RigidTwoNodeTruss.h>

extern "C" {
void _FORTRAN(hxgaus)(int &, int &, int &, int &, int &, int &,
                      double &,  double &, double &, double &);
void _FORTRAN(h8shpe)(double &, double &, double &, double *, double *,
                      double *, double *, double *, double *, double *,
                      double &);
void _FORTRAN(lgauss)(const int &, int &, double *, double *);
}

double Hexa8ShapeFct(double Shape[8], double DShape[8][3], double m[3], double X[8], double Y[8], double Z[8]);
void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);

RigidEightNodeBrick::RigidEightNodeBrick(int *_nn)
 : SuperElement(true)
{
  nnodes = 8;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 18;
  subElems = new Element * [nSubElems];
  int indices[18][2] = {{0,1},
                        {1,2},
                        {2,3},
                        {4,5},
                        {5,6},
                        {6,7},
                        {0,2},
                        {4,6},
                        {0,3},
                        {4,7},
                        {0,4},
                        {1,5},
                        {2,6},
                        {3,7},
                        {0,5},
                        {1,6},
                        {2,7},
                        {0,7}};
  for(int i = 0; i < nSubElems; ++i)
    subElems[i] = new RigidTwoNodeTruss(indices[i]);
}

double
RigidEightNodeBrick::getMass(const CoordSet& cs) const
{
  if(prop == NULL || prop->rho == 0) return 0.0;

  const int nnodes = 8;
  const int numgauss = 2;

  double X[8], Y[8], Z[8];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double m[3], Shape[8], DShape[8][3];
  double wx, wy, wz, dOmega;
  double volume = 0.0;
  for(int i = 1; i <= numgauss; i++) {
    _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
    for(int j = 1; j <= numgauss; j++) {
      _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
      for(int k = 1; k <= numgauss; k++) {
        _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);
        dOmega = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
        volume += fabs(dOmega)*wx*wy*wz;
      }
    }
  }

  return prop->rho*volume;
}

void
RigidEightNodeBrick::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                                     Vector& gravityForce, int gravflg, GeomState *geomState)
{
  if(prop == NULL || prop->rho == 0) {
    gravityForce.zero();
    return;
  }

  const int nnodes = 8;

  // Lumped
  if (gravflg != 2) {

    double totmas = getMass(cs);

    double m[24*24];
    FullSquareMatrix result = massMatrix(cs, m, 1);
    std::vector<double> factors;
    lumpMatrix(result, factors);

    // divvy up the total body force using same ratio as the corresponding diagonal of the lumped mass matrix to the total mass
    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = totmas*gravityAcceleration[0]*(3.0*factors[3*i+0]);
      gravityForce[3*i+1] = totmas*gravityAcceleration[1]*(3.0*factors[3*i+1]);
      gravityForce[3*i+2] = totmas*gravityAcceleration[2]*(3.0*factors[3*i+2]);
    }
  }
  // Consistent
  else {

    int numgauss = 2;

    double X[8], Y[8], Z[8];
    cs.getCoordinates(nn, nnodes, X, Y, Z);

    double lforce[8];
    for(int i=0; i<nnodes; ++i) lforce[i] = 0.0;

    for (int pt1 = 1; pt1 <= numgauss; pt1++)  {
      for (int pt2 = 1; pt2 <= numgauss; pt2++)  {
        for (int pt3 = 1; pt3 <= numgauss; pt3++)  {
          // get gauss point
          double xi, eta, mu, wt;
          _FORTRAN(hxgaus)(numgauss, pt1, numgauss, pt2, numgauss, pt3, xi, eta, mu, wt);

          //compute shape functions
          double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];
          double detJ; // det of jacobian

          _FORTRAN(h8shpe)(xi, eta, mu, X, Y, Z,
                           shapeFunc, shapeGradX, shapeGradY, shapeGradZ, detJ);
          for (int i = 0; i < 8; ++i)
            lforce[i] += wt*shapeFunc[i]*detJ;
        }
      }
    }
    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*prop->rho*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*prop->rho*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*prop->rho*gravityAcceleration[2];
    }
  }
}

FullSquareMatrix
RigidEightNodeBrick::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  const int nnodes= 8;
  const int ndofs = 24;
  const int numgauss= 2;

  FullSquareMatrix M(ndofs,mel);

  if(prop == NULL || prop->rho == 0) {
    M.zero();
    return M;
  }

  if(cmflg) { // consistent mass matrix
    M.zero();
    int ls[24] = {0,3,6,9,12,15,18,21,
                  1,4,7,10,13,16,19,22,
                  2,5,8,11,14,17,20,23};

    double X[8], Y[8], Z[8];
    cs.getCoordinates(nn, nnodes, X, Y, Z);

    // integration: loop over Gauss pts
    double m[3], Shape[8], DShape[8][3];
    double wx, wy, wz, w;
    double dOmega; // det of jacobian
#ifdef CHECK_JACOBIAN
    int jSign = 0;
#endif
    for(int i = 1; i <= numgauss; i++) {
      _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
      for(int j = 1; j <= numgauss; j++) {
        _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
        for(int k = 1; k <= numgauss; k++) {
          _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);
          dOmega = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
          checkJacobian(&dOmega, &jSign, "RigidEightNodeBrick::massMatrix");
#endif
          w = fabs(dOmega)*wx*wy*wz*prop->rho;
          addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
        }
      }
    }
  }
  else { // Lumped mass matrix
    fprintf(stderr," *** In RigidEightNodeBrick::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return(M);
}

#endif
