#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidTwoNodeTruss.h>

RigidTwoNodeTruss::RigidTwoNodeTruss(int* _nn)
  : ConstantDistanceConstraint(_nn)
{
}

RigidTwoNodeTrussWithMass::RigidTwoNodeTrussWithMass(int* _nn)
  : ConstantDistanceConstraint(_nn)
{
}

double
RigidTwoNodeTrussWithMass::getMass(const CoordSet& cs) const
{
  if(prop == NULL || prop->rho == 0 || prop->A == 0)
    return 0;

  auto &nd1 = cs.getNode( nn[0] );
  auto &nd2 = cs.getNode( nn[1] );

  double x[2], y[2], z[2];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

  double dx = x[1] - x[0];
  double dy = y[1] - y[0];
  double dz = z[1] - z[0];

  double length = sqrt( dx*dx + dy*dy + dz*dz );
  double   mass = (length*prop->A*prop->rho);

  return mass;
}

void
RigidTwoNodeTrussWithMass::getGravityForce(CoordSet& cs, double *gravityAcceleration, Vector &gravityForce,
                                           int gravflg, GeomState *geomState)
{
  if(prop == NULL || prop->rho == 0 || prop->A == 0) {
    gravityForce.zero();
    return;
  }

  const double massPerNode = 0.5 * getMass(cs);

  const double fx = massPerNode * gravityAcceleration[0];
  const double fy = massPerNode * gravityAcceleration[1];
  const double fz = massPerNode * gravityAcceleration[2];

  gravityForce[0] = fx;
  gravityForce[1] = fy;
  gravityForce[2] = fz;
  gravityForce[3] = fx;
  gravityForce[4] = fy;
  gravityForce[5] = fz;
}

FullSquareMatrix
RigidTwoNodeTrussWithMass::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  FullSquareMatrix elementMassMatrix(6, mel);
  elementMassMatrix.zero();

  if(prop == NULL || prop->rho == 0 || prop->A == 0) {
     return elementMassMatrix;
  }

  const double mass = getMass(cs);

  if (cmflg) {
    // Consistent mass matrix
    const double outDiagMass = mass / 6.0;
    const double diagMass = 2.0 * outDiagMass;

    for (int i = 0; i < 6; ++i) {
      elementMassMatrix[i][i] = diagMass;
    }
    for (int i = 0; i < 3; ++i) {
      const int j = i + 3;
      elementMassMatrix[i][j] = outDiagMass;
      elementMassMatrix[j][i] = outDiagMass;
    }
  }
  else {
    // Lumped mass matrix
    const double massPerNode = 0.5 * mass;
    for (int i = 0; i < 6; ++i) {
      elementMassMatrix[i][i] = massPerNode;
    }
  }

  return elementMassMatrix;
}

#endif
