#if defined(USE_EIGEN3)
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <Element.d/DiscreteMass.d/DiscreteMass6Dof.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Element.d/Function.d/Rotation.d/IncrementalRotationVector.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Corotational.d/utilities.h>
#include <Utils.d/dofset.h>

#include <Eigen/Geometry>

DiscreteMass6Dof::DiscreteMass6Dof(int* nodenums)
 : C0(NULL)
{
  nn[0] = nodenums[0];
  f0.setZero();
}

DiscreteMass6Dof::~DiscreteMass6Dof()
{
  if(C0) delete C0;
}

void
DiscreteMass6Dof::setFrame(EFrame *elemframe)
{
  C0 = new Eigen::Matrix3d;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      C0->coeffRef(i,j) = (*elemframe)[i][j];
}

void
DiscreteMass6Dof::renum(const int *table)
{
  if(nn[0] > -1) nn[0] = table[nn[0]];
}

void
DiscreteMass6Dof::renum(EleRenumMap& table)
{
  if(nn[0] > -1) nn[0] = table[nn[0]];
}

int*
DiscreteMass6Dof::nodes(int* p) const
{
  if(p == 0) p = new int[1];
  p[0] = nn[0];
  return p;
}

int *
DiscreteMass6Dof::dofs(DofSetArray &dsa, int *p) const
{
  if(p == 0) p = new int[6];
  dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p);
  return p;
}

void
DiscreteMass6Dof::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, 1, DofSet::XYZdisp | DofSet::XYZrot);
}

double
DiscreteMass6Dof::getMass(const CoordSet& cs) const
{
  return prop->rho;
}

FullSquareMatrix
DiscreteMass6Dof::massMatrix(const CoordSet &cs, double *mdata, int cmflg) const
{
  // ref: "Rigid Body Dynamics of Mechanisms 1 Theoretical Basis" by Hubert Hahn, Section 4.2.4.1
  double m = prop->rho;
  Eigen::Matrix3d J;
  Eigen::Vector3d c;
  if(!C0) {
    // in this case the components of the inertia tensor and offset vector are in the global frame
    J << prop->Ixx, prop->Ixy, prop->Ixz,
         prop->Ixy, prop->Iyy, prop->Iyz,
         prop->Ixz, prop->Iyz, prop->Izz;
    c << prop->cx, prop->cy, prop->cz;
  }
  else {
    // in this case the components of the inertia tensor and offset vector are in the element frame
    Eigen::Matrix3d I;
    I << prop->Ixx, prop->Ixy, prop->Ixz,
         prop->Ixy, prop->Iyy, prop->Iyz,
         prop->Ixz, prop->Iyz, prop->Izz;
    J = C0->transpose()*I*(*C0);
    Eigen::Vector3d b;
    b << prop->cx, prop->cy, prop->cz;
    c = C0->transpose()*b;
  }

  Eigen::Matrix3d C;
  C << 0,        -prop->cz,  prop->cy,
         prop->cz,  0,        -prop->cx,
        -prop->cy,  prop->cx,  0;

  Eigen::Map<Eigen::Matrix<double,6,6,Eigen::RowMajor> > M(mdata);
  M << m*Eigen::Matrix3d::Identity(), -m*C,
       -m*C.transpose(),              J - m*C*C;

  return FullSquareMatrix(6, mdata);
}

FullSquareMatrix
DiscreteMass6Dof::stiffness(const CoordSet &cs, double *k, int flg) const
{
  FullSquareMatrix ret(6,k);
  ret.zero();
  return ret;
}

void
DiscreteMass6Dof::getGravityForce(CoordSet&, double *g, Vector &_f, int, GeomState *gs)
{
  Eigen::Vector3d c;
  f0 << g[0]*prop->rho, g[1]*prop->rho, g[2]*prop->rho;
  if(!C0) {
    // in this case the components of c are given in the global frame
    c << prop->cx, prop->cy, prop->cz;
  }
  else {
    // in this case the components of f0 and c are given in the element frame
    Eigen::Vector3d b;
    b << prop->cx, prop->cy, prop->cz;
    c = C0->transpose()*b;
  }

  Eigen::Map<Eigen::Matrix<double,6,1> > f(_f.data());
  f.head<3>() = f0;
  f.tail<3>() = c.cross(f0);
}

void
DiscreteMass6Dof::computePressureForce(CoordSet&, Vector& elPressureForce, GeomState *, int, double)
{
  elPressureForce.zero();
}

void
DiscreteMass6Dof::getStiffAndForce(GeomState& curState, CoordSet& c0, FullSquareMatrix& K,
                                   double* f, double dt, double t)
{
  getStiffAndForce(NULL, curState, c0, K, f, dt, t);
}

void
DiscreteMass6Dof::getStiffAndForce(GeomState *, GeomState &curState, CoordSet &,
                                   FullSquareMatrix &_K, double *_f, double, double)
{
  Eigen::Vector3d c;
  if(!C0) {
    // in this case the components of c are given in the global frame
    c << prop->cx, prop->cy, prop->cz;
  }
  else {
    // in this case the components of f0 and c are given in the element frame
    Eigen::Vector3d b;
    b << prop->cx, prop->cy, prop->cz;
    c = C0->transpose()*b;
  }

  // compute the discrepancy between the gravity force in the undeformed and the current configurations
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&curState[nn[0]].R[0][0]);
  Eigen::Map<Eigen::Matrix<double,6,1> > f(_f);
  f.head<3>() = Eigen::Vector3d::Zero();
  f.tail<3>() = -(R*c).cross(f0) + c.cross(f0);

  // and the corresponding load stiffness matrix
  Eigen::Map<Eigen::Matrix<double,6,6,Eigen::RowMajor> > K(_K.data());

  K.setZero();
  Eigen::Matrix3d ab = (R*c)*f0.transpose();

  K(3,3) = ab(1,1) + ab(2,2);
  K(4,4) = ab(0,0) + ab(2,2);
  K(5,5) = ab(0,0) + ab(1,1);

  K(3,4) = K(4,3) = -0.5*(ab(0,1) + ab(1,0));
  K(3,5) = K(5,3) = -0.5*(ab(0,2) + ab(2,0));
  K(4,5) = K(5,4) = -0.5*(ab(1,2) + ab(2,1));
}

void
DiscreteMass6Dof::getInternalForce(GeomState *, GeomState &curState, CoordSet &,
                                   FullSquareMatrix&, double *_f, double, double)
{
  Eigen::Vector3d c;
  if(!C0) {
    // in this case the components of c are given in the global frame
    c << prop->cx, prop->cy, prop->cz;
  }
  else {
    // in this case the components of f0 and c are given in the element frame
    Eigen::Vector3d b;
    b << prop->cx, prop->cy, prop->cz;
    c = C0->transpose()*b;
  }

  // compute the discrepancy between the gravity force in the undeformed and the current configurations
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&curState[nn[0]].R[0][0]);
  Eigen::Map<Eigen::Matrix<double,6,1> > f(_f);
  f.head<3>() = Eigen::Vector3d::Zero();
  f.tail<3>() = -(R*c).cross(f0) + c.cross(f0);
}

double
DiscreteMass6Dof::getElementEnergy(GeomState &curState, CoordSet &c0)
{
  Eigen::Vector3d c;
  if(!C0) {
    // in this case the components of c are given in the global frame
    c << prop->cx, prop->cy, prop->cz;
  }
  else {
    // in this case the components of f0 and c are given in the element frame
    Eigen::Vector3d b;
    b << prop->cx, prop->cy, prop->cz;
    c = C0->transpose()*b;
  }

  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&curState[nn[0]].R[0][0]);
  Eigen::Map<Eigen::Vector3d> theta(curState[nn[0]].theta);
  Eigen::Vector3d u;
  u << curState[nn[0]].x - c0[nn[0]]->x,
       curState[nn[0]].y - c0[nn[0]]->y,
       curState[nn[0]].z - c0[nn[0]]->z;

  // compute the discrepancy between the scalar potential of the gravity force and
  // the work done by the constant part of the gravity force computed in getGravityForce
  Eigen::Vector3d u0 = u+R*c-c;
  return -u0.dot(f0) + u.dot(f0) + theta.dot(c.cross(f0));
}

void
DiscreteMass6Dof::getInertialStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0,
                                           FullSquareMatrix &_K, double *_f, double dt, double t,
                                           double beta, double gamma, double alphaf, double alpham)
{
  double m = prop->rho;
  Eigen::Matrix3d J;
  Eigen::Vector3d c;
  if(!C0) {
    // in this case the components of the inertia tensor and offset vector are in the global frame
    J << prop->Ixx, prop->Ixy, prop->Ixz,
         prop->Ixy, prop->Iyy, prop->Iyz,
         prop->Ixz, prop->Iyz, prop->Izz;
    c << prop->cx, prop->cy, prop->cz;
  }
  else {
    // in this case the components of the inertia tensor and offset vector are in the element frame
    Eigen::Matrix3d I;
    I << prop->Ixx, prop->Ixy, prop->Ixz,
         prop->Ixy, prop->Iyy, prop->Iyz,
         prop->Ixz, prop->Iyz, prop->Izz;
    J = C0->transpose()*I*(*C0);
    Eigen::Vector3d b;
    b << prop->cx, prop->cy, prop->cz;
    c = C0->transpose()*b;
  }

  Eigen::Matrix3d C, J0;
  C << 0,        -prop->cz,  prop->cy,
         prop->cz,  0,        -prop->cx,
        -prop->cy,  prop->cx,  0;
  J0 = J - m*C*C;

  Eigen::Map<Eigen::Matrix<double,6,1> > f(_f);
  Eigen::Map<Eigen::Matrix<double,6,6,Eigen::RowMajor> > K(_K.data());
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R(&c1[nn[0]].R[0][0]);

  if(beta == 0 || t == 0) {
    Eigen::Map<Eigen::Vector3d> Omega(&c1[nn[0]].v[3]);

    f.head<3>() = m*R*Omega.cross(Omega.cross(c));
    f.tail<3>() = R*Omega.cross(J0*Omega);

    K.setZero();

    return;
  }

  if(c1.getNumRotationDof(nn[0]) == 2) {
    std::cerr << " *** WARNING: DiscreteMass6Dof::getInertialStiffAndForce with 1 constrained rotation dof is not implemented\n";
    exit(-1);
  }

  Eigen::Vector3d u, u_n, F, G, P, Q;
  Eigen::Matrix3d E;
  Eigen::Matrix<double,6,1> a, v, inc_displacement;

  Eigen::Map<Eigen::Matrix<double,6,1> > a_n((*refState)[nn[0]].a), v_n((*refState)[nn[0]].v);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > R_n(&(*refState)[nn[0]].R[0][0]);
  u_n << (*refState)[nn[0]].x, (*refState)[nn[0]].y, (*refState)[nn[0]].z;
  u   << c1[nn[0]].x, c1[nn[0]].y, c1[nn[0]].z;

  Eigen::VectorBlock<Eigen::Matrix<double,6,1>,3> udot  = v.head<3>(), Omega = v.tail<3>(),
                                                  uddot = a.head<3>(), Alpha = a.tail<3>();

  double s1 = gamma/(dt*beta);
  double s2 = (1-alpham)/(dt*dt*beta*(1-alphaf));

  inc_displacement.head<3>() = u - u_n;
  mat_to_vec<double>(R_n.transpose()*R, inc_displacement.tail<3>());

  v = s1*inc_displacement + (1-(1-alphaf)*gamma/beta)*v_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*a_n;
  a = s2*inc_displacement - (1-alpham)/(dt*beta)*v_n + ((alpham-1)/(2*beta)+1)*a_n;

  E = R*C*R.transpose();
  F = m*E*uddot + R*(J0*Alpha + Omega.cross(J0*Omega));
  G = R*(C*Alpha + Omega.cross(C*Omega));
  P = J0*Omega;
  Q = C*Omega;

  Eigen::Array<double,18,1> sconst;
  sconst << (*refState)[nn[0]].R[0][0], (*refState)[nn[0]].R[0][1], (*refState)[nn[0]].R[0][2],
            (*refState)[nn[0]].R[1][0], (*refState)[nn[0]].R[1][1], (*refState)[nn[0]].R[1][2],
            (*refState)[nn[0]].R[2][0], (*refState)[nn[0]].R[2][1], (*refState)[nn[0]].R[2][2],
            c1[nn[0]].R[0][0], c1[nn[0]].R[0][1], c1[nn[0]].R[0][2],
            c1[nn[0]].R[1][0], c1[nn[0]].R[1][1], c1[nn[0]].R[1][2],
            c1[nn[0]].R[2][0], c1[nn[0]].R[2][1], c1[nn[0]].R[2][2];
  Simo::Jacobian<double,Simo::IncrementalRotationVector> dPsidq(sconst, Eigen::Array<int,0,1>::Zero());
  Eigen::Matrix3d D = dPsidq(Eigen::Vector3d::Zero(), 0.);

  // compute inertial forces and moments, minus M*a
  f.head<3>() = m*(-R*(c.cross(Alpha) - Omega.cross(Omega.cross(c))) + c.cross(Alpha));
  f.tail<3>() = R*(-m*(R.transpose()*uddot).cross(c) + J0*Alpha + Omega.cross(J0*Omega)) + m*C.transpose()*uddot - J0*Alpha;

  // compute the jacobian of the inertial forces and moments, minus s2*M
  K.topLeftCorner<3,3>()     = Eigen::Matrix3d::Zero();
  K.topRightCorner<3,3>()    = m*(skew(G) - R*(s2*C + s1*(skew(Omega)*C - skew(Q)))*D + s2*C);
  K.bottomLeftCorner<3,3>()  = s2*m*(E + C.transpose());
  K.bottomRightCorner<3,3>() = -0.5*skew(F) + m*E*skew(uddot) + R*(s2*J0 + s1*(skew(Omega)*J0 - skew(P)))*D - s2*J0;
}

void
DiscreteMass6Dof::getNLVonMises(Vector& stress, Vector& weight,
                                GeomState&, CoordSet&, int)
{
  stress.zero();
  weight.zero();
}

#endif
