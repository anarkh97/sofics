#include <Math.d/matrix.h> 
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <cstdio>

LinearStrain linearStrain;
GreenLagrangeStrain greenLagrangeStrain;
LogarithmicStrain logarithmicStrain;
PrincipalStretches principalStretches;
LogarithmicPrincipalStretches logarithmicPrincipalStretches;
ElasticLogarithmicPrincipalStretches elasticLogarithmicPrincipalStretches;
DeformationGradient deformationGradient;

Tensor *
LinearStrain::getTMInstance()
{
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}

Tensor *
LinearStrain::getStressInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LinearStrain::getStrainInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LinearStrain::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}

Tensor *
LinearStrain::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_null *DB = new Tensor_d2s2_Sd12s34_null(numdofs);
  return DB;
}

Tensor *
LinearStrain::getCacheInstance()
{
  return NULL;
}

void 
LinearStrain::getEBandDB(Tensor &_e, Tensor & bB, Tensor &DB,
                         const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *, double *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(bB);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(gradU + tgradU);
  enonsym.convertToSym(e);

  B.setZero();
  B.addSymPart(dgradUdqk);
}

void
LinearStrain::getEandB(Tensor &_e, Tensor & __B, 
                       const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *, double *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(gradU + tgradU);
  enonsym.convertToSym(e);

  B.setZero();
  B.addSymPart(dgradUdqk);
}

void 
LinearStrain::getE(Tensor &_e, Tensor &_gradU, Tensor *, double *)
{
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);
  Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(gradU + tgradU);
  enonsym.convertToSym(e);
}

void
LinearStrain::transformStress(Tensor &_stress, Tensor *, Tensor_d0s2_Ss12 &S)
{
  // do nothing: transformation is only applied for finite-strain materials
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(_stress);
  S = stress;
}

Tensor *
GreenLagrangeStrain::getTMInstance()
{
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}

Tensor *
GreenLagrangeStrain::getStressInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
GreenLagrangeStrain::getStrainInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
GreenLagrangeStrain::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}

Tensor *
GreenLagrangeStrain::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_sparse *DB = new Tensor_d2s2_Sd12s34_sparse(numdofs);
  return DB;
}

Tensor *
GreenLagrangeStrain::getCacheInstance()
{
  // (AN) Cache required for crushable foam
  Tensor_d1s2_full *cache = new Tensor_d1s2_full(2);
  return cache;
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void 
GreenLagrangeStrain::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_sparse & DB = static_cast<Tensor_d2s2_Sd12s34_sparse &>(_DB);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  // e = 1/2(gradU^t + gradU + gradU^t|gradU)
  // de/dq = 1/2(dgradUdq^t + dgradUdq + ...)
  Tensor_d0s2 temp1;
  temp1 = tgradU|gradU;

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(tgradU + (gradU + temp1));
  enonsym.convertToSym(e);

  int size = B.getSize();
  Tensor_d1s2_full temp2(size);
  temp2 = tgradU | dgradUdqk;

  B.assignSymPart(dgradUdqk, temp2);

  dgradUdqk.getSymSquare(DB);

  // save F to cache
#ifdef USE_EIGEN3
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();
  // second item updated later.
  static_cast<Tensor_d1s2_full &>(*cache)[0] = static_cast<Tensor_d1s2_full &>(*cache)[1] = F;
#else
  std::cerr << "ERROR: GreenLagrangeStrain::getEBandDB can not calculate deformation gradient without EIGEN.\n";
  exit(-1);
#endif
}

void
GreenLagrangeStrain::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  // e = 1/2(gradU^t + gradU + gradU^t|gradU)
  // de/dq = 1/2(dgradUdq^t + dgradUdq + ...)
  Tensor_d0s2 temp1;
  temp1 = tgradU|gradU;

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(tgradU + (gradU + temp1));
  enonsym.convertToSym(e);

  int size = B.getSize();
  Tensor_d1s2_full temp2(size);
  temp2 = tgradU | dgradUdqk;

  B.assignSymPart(dgradUdqk, temp2);
  // save F to cache
#ifdef USE_EIGEN3
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();
  // second item updated later.
  static_cast<Tensor_d1s2_full &>(*cache)[0] = static_cast<Tensor_d1s2_full &>(*cache)[1] = F;
#else
  std::cerr << "ERROR: GreenLagrangeStrain::getEandB can not calculate deformation gradient without EIGEN.\n";
  exit(-1);
#endif
}

void 
GreenLagrangeStrain::getE(Tensor &_e, Tensor &_gradU, Tensor *cache, double *)
{
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);
  Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);

  // e = 1/2(gradU^t + gradU + gradU^t|gradU)
  Tensor_d0s2 tgradU;
  gradU.getTranspose(tgradU);

  Tensor_d0s2 temp;
  temp = tgradU|gradU;

  Tensor_d0s2 enonsym;
  enonsym = (1/2.)*(tgradU + (gradU + temp));
  enonsym.convertToSym(e);

  // save F to cache
#ifdef USE_EIGEN3
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();
  // second item updated later.
  static_cast<Tensor_d1s2_full &>(*cache)[0] = static_cast<Tensor_d1s2_full &>(*cache)[1] = F;
#else
  std::cerr << "ERROR: GreenLagrangeStrain::getE can not calculate deformation gradient without EIGEN.\n";
  exit(-1);
#endif
}

void
GreenLagrangeStrain::transformStress(Tensor &_stress, Tensor *cache, Tensor_d0s2_Ss12 &S)
{
  // do nothing: stress is already PK2 in this case

  //Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(_stress);
  //S = stress;

  // AN: Temporarily output set to cauchy stress
  Tensor_d0s2_Ss12 &stress = static_cast<Tensor_d0s2_Ss12 &>(_stress);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F = static_cast<Tensor_d1s2_full &>(*cache)[0].matrix();
  double J = F.determinant();
  if(J==0) {
    fprintf(stderr, "***ERROR: Deformation gradient has zero determinant.\n");
    exit(-1);
  }

  Eigen::Matrix3d Spk;
  stress.assignTo(Spk);
  Eigen::Matrix3d sigma = (1/J)*F*Spk*F.transpose();
  S = sigma;
}

void
GreenLagrangeStrain::getL(const Tensor &_gradUnp, const Tensor &_gradUn, Tensor *cache, double dt)
{
  const Tensor_d0s2 &gradUn = static_cast<const Tensor_d0s2 &>(_gradUn);
  const Tensor_d0s2 &gradUnp = static_cast<const Tensor_d0s2 &>(_gradUnp);

#ifdef USE_EIGEN3
  Eigen::Matrix3d F = gradUnp.matrix() + Eigen::Matrix3d::Identity();
  Eigen::Matrix3d dUn, dUnp;
  gradUn.assignTo(dUn);
  gradUnp.assignTo(dUnp);
  
  Eigen::Matrix3d Finv, Fdot, L;
  if(dt!=0) Fdot = (dUnp - dUn)/dt;
  else {
    fprintf(stderr, "ERROR: Incorrect time step size provided (dt=%e).\n", dt);
    exit(-1);
  }
  Finv = F.inverse();
  L = Fdot*Finv;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = L;
#else
  std::cerr << "ERROR: GreenLagrangeStrain::getL can not calculate velocity gradient tensor without EIGEN.\n";
  exit(-1);
#endif
}

Tensor *
LogarithmicStrain::getTMInstance()
{
  Tensor_d0s4_Ss12s34 *s = new Tensor_d0s4_Ss12s34;
  return s;
}

Tensor *
LogarithmicStrain::getStressInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LogarithmicStrain::getStrainInstance()
{
  Tensor_d0s2_Ss12 *s = new Tensor_d0s2_Ss12;
  return s;
}

Tensor *
LogarithmicStrain::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23 *B = new Tensor_d1s2_Ss23(numdofs);
  return B;
}

Tensor *
LogarithmicStrain::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_dense *DB = new Tensor_d2s2_Sd12s34_dense(numdofs);
  return DB;
}

Tensor *
LogarithmicStrain::getCacheInstance()
{
  Tensor_d1s2_full *cache = new Tensor_d1s2_full(2);
  return cache;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicStrain::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense & DB = static_cast<Tensor_d2s2_Sd12s34_dense &>(_DB);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix<double,3,3>::Identity();

  if(F.isIdentity()) {
    //note: the eigenvalue decomposition is apparently not differentiable in this case.
    //using LinearStrain for consistency with linear elasticity at small strains
    e = 0.5*(gradU.matrix() + gradU.matrix().transpose());
    for(int k=0; k<numdofs; ++k) {
      B[k] = 0.5*(dGradUdq[k] + dGradUdq[k].transpose());
    }
    DB.setZero();
    static_cast<Tensor_d1s2_full &>(*cache)[0] = static_cast<Tensor_d1s2_full &>(*cache)[1] = Eigen::Matrix3d::Identity();
  }
  else {
    // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3> > dec((F.transpose()*F).eval());
    const Eigen::Vector3d &l = dec.eigenvalues();
    const Eigen::Matrix3d &N = dec.eigenvectors();

    // logarithmic (Hencky) strain, Lagrangean description E = log(U) = 0.5*N*log(λ²)*N^T
    Eigen::Vector3d lnl = l.array().log();
    Eigen::Matrix3d Nlnl = N*lnl.asDiagonal();
    e = 0.5*Nlnl*N.transpose();

    // Moore-Penrose pseudo inverses of C - λ²_i*I
    double tol = std::numeric_limits<double>::epsilon();
    Eigen::Array<Eigen::Matrix3d,3,1> mpinverse;
    for(int i=0; i<3; ++i) {
      Eigen::Vector3d singularValues = l - Eigen::Vector3d::Constant(l[i]);
      Eigen::Vector3d invertedSingularVals;
      for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
      mpinverse[i] = N * invertedSingularVals.asDiagonal() * N.transpose();
    }

    // some more precomputation...
    Eigen::Array<double,3,1> linv = l.array().inverse();
    Eigen::Array<double,3,1> linv2 = linv.square();

    // allocate memory for intermediate derivatives
    Eigen::Array<Eigen::Array<double,3,1>,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
    Eigen::Array<Eigen::Vector3d,Eigen::Dynamic,1> dlnldq(numdofs);
    Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dCdq(numdofs), d2Cdqkdq(numdofs), dNdq(numdofs), d2Ndqkdq(numdofs), d2Edqkdq(numdofs), tmp2(numdofs);
    Eigen::Matrix3d tmp1,tmp3,toto,dEdqk;

    for(int k=0; k<numdofs; ++k) {

      // first derivative of C with respect to q_k
      dCdq[k] = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

      // second derivatives of C with respect to q_k,j
      for(int j=0; j<=k; ++j) {
        if((j-k)%3 == 0) {
          Eigen::Matrix3d tmp = dGradUdq[j].transpose()*dGradUdq[k];
          d2Cdqkdq[j] = tmp + tmp.transpose();
        }
      }

      Eigen::Matrix3d dCdqkN = dCdq[k]*N;

      // first derivative of λ² with respect to q_k
      dldq[k] = (N.transpose()*dCdqkN).diagonal();

      // first derivatve of N with respect to q_k
      for(int i=0; i<3; ++i) dNdq[k].col(i) = -mpinverse[i]*dCdqkN.col(i);

      for(int j=0; j<=k; ++j) {
        Eigen::Matrix<double,3,3> dCdqjdNdqk = dCdq[j]*dNdq[k];
        if((j-k)%3 == 0) {
          Eigen::Matrix3d d2CdqkdqjN = d2Cdqkdq[j]*N;

          // second derivatives of λ² with respect to q_k,j for non-zero d2Cdqkdq[j]
          d2ldqkdq[j] = (N.transpose()*(d2CdqkdqjN + 2*dCdqjdNdqk)).diagonal();

          // second derivatives of N with respect to q_k,j 
          toto = dCdqjdNdqk + dCdq[k]*dNdq[j] + d2CdqkdqjN
                 - dNdq[k]*dldq[j].matrix().asDiagonal() - dNdq[j]*dldq[k].matrix().asDiagonal();
        }
        else {
          // second derivatives of λ² with respect to q_k,j for zero d2Cdqkdq[j]
          d2ldqkdq[j] = (N.transpose()*(2*dCdqjdNdqk)).diagonal();

          toto = dCdqjdNdqk + dCdq[k]*dNdq[j] 
                 - dNdq[k]*dldq[j].matrix().asDiagonal() - dNdq[j]*dldq[k].matrix().asDiagonal();
        }
        // second derivatives of N with respect to q_k,j
        d2Ndqkdq[j] = -N*(dNdq[j].transpose()*dNdq[k]).diagonal().asDiagonal();
        for(int i=0; i<3; ++i)
          d2Ndqkdq[j].col(i) -= mpinverse[i]*toto.col(i);
      }
      // first derivative of E with respect to q_k
      dlnldq[k] = linv*dldq[k];
      tmp1 = dNdq[k]*lnl.asDiagonal();
      tmp2[k] = dlnldq[k].asDiagonal()*N.transpose();
      dEdqk.triangularView<Eigen::Upper>() = 0.5*(dNdq[k]*Nlnl.transpose()
                    + N*dlnldq[k].asDiagonal()*N.transpose()
                      + Nlnl*dNdq[k].transpose());
      for(int j=0; j<=k; ++j) {
        // second derivative of E with respect to q_k,j
        Eigen::Vector3d d2lnldqkdqj = linv*d2ldqkdq[j] - linv2*dldq[j]*dldq[k];
        tmp3 = d2Ndqkdq[j]*Nlnl.transpose() + dNdq[k]*tmp2[j] + dNdq[j]*(tmp1.transpose() + tmp2[k]);
        d2Edqkdq[j].triangularView<Eigen::Upper>() = 0.5*(tmp3 + tmp3.transpose()
                             + N*d2lnldqkdqj.asDiagonal()*N.transpose());
      }

      B[k] = dEdqk;
      for(int j=0; j<=k; ++j) DB[j*(2*numdofs-j-1)/2+k] = d2Edqkdq[j];
    }

    // store copies of N and λ² in cache
    static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
    static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();
  }

#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicStrain::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23 & B = static_cast<Tensor_d1s2_Ss23 &>(__B);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  int numdofs = dgradUdqk.getSize();

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix<double,3,3>::Identity();

  if(F.isIdentity()) {
    //note: the eigenvalue decomposition is apparently not differentiable in this case.
    //using LinearStrain for consistency with linear elasticity at small strains
    e = 0.5*(gradU.matrix() + gradU.matrix().transpose());
    for(int k=0; k<numdofs; ++k) {
      B[k] = 0.5*(dGradUdq[k] + dGradUdq[k].transpose());
    }
    static_cast<Tensor_d1s2_full &>(*cache)[0] = static_cast<Tensor_d1s2_full &>(*cache)[1] = Eigen::Matrix3d::Identity();
  }
  else {
    // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3> > dec((F.transpose()*F).eval());
    const Eigen::Vector3d &l = dec.eigenvalues();
    const Eigen::Matrix3d &N = dec.eigenvectors();

    // logarithmic (Hencky) strain, Lagrangean description E = log(U) = 0.5*N*log(λ²)*N^T
    Eigen::Vector3d lnl = l.array().log();
    Eigen::Matrix3d ylnl = N*lnl.asDiagonal();
    e = 0.5*ylnl*N.transpose();

    // Moore-Penrose pseudo inverses of C - λ²_i*I
    double tol = std::numeric_limits<double>::epsilon();
    Eigen::Array<Eigen::Matrix3d,3,1> mpinverse;
    for(int i=0; i<3; ++i) {
      Eigen::Vector3d singularValues = l - Eigen::Vector3d::Constant(l[i]);
      Eigen::Vector3d invertedSingularVals;
      for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
      mpinverse[i] = N * invertedSingularVals.asDiagonal() * N.transpose();
    }

    // some more precomputation...
    Eigen::Array<double,3,1> linv = l.array().inverse();

    // allocate memory for intermediate derivatives
    Eigen::Array<double,3,1> dldqk, dlnldqk;
    Eigen::Matrix3d dCdqk, dCdqky, dydqk, dEdqk;

    for(int k=0; k<numdofs; ++k) {

      // first derivative of C with respect to q_k
      dCdqk = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

      dCdqky = dCdqk*N;

      // first derivative of λ² with respect to q_k
      dldqk = (N.transpose()*dCdqky).diagonal();

      // first derivatve of y with respect to q_k
      for(int i=0; i<3; ++i) dydqk.col(i) = -mpinverse[i]*dCdqky.col(i);

      // first derivative of E with respect to q_k
      dlnldqk = linv*dldqk;
      dEdqk.triangularView<Eigen::Upper>() = 0.5*(dydqk*ylnl.transpose()
                    + N*dlnldqk.matrix().asDiagonal()*N.transpose()
                      + ylnl*dydqk.transpose());

      B[k] = dEdqk;
    }

    // store copies of N and λ² in cache
    static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
    static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();
  }

#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicStrain::getE(Tensor &_e, Tensor &_gradU, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  Tensor_d0s2_Ss12 & e = static_cast<Tensor_d0s2_Ss12 &>(_e);

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();

  // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &N = dec.eigenvectors();

  // logarithmic (Hencky) strain, Lagrangean description E = log(U) = 0.5*N*log(λ²)*N^T
  e = 0.5*(N * l.array().log().matrix().asDiagonal() * N.transpose());

  // store copies of N and λ² in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();

#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicStrain::transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Eigen::Matrix3d T; static_cast<Tensor_d0s2_Ss12 &>(stress).assignTo(T);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > N = static_cast<Tensor_d1s2_full &>(*cache)[0].matrix();
  Eigen::Array3d l = static_cast<Tensor_d1s2_full &>(*cache)[1].matrix().diagonal();

  // 2nd Piola-Kirchhoff stress S = F^{-1}*τ*F^{-t} = F^{-1}*R*T*R^{t}*F^{-t} = U^{-1}*T*U^{-t}
  // note: T is the rotated Kirchhoff stress which is conjugate to the Lagrangean Hencky strain E = log(U), i.e. T = ∂ω(E)/∂E
  //       and U = N*diag(1/λ)*N^T  is the right stretch tensor corresponding to the polar decomposition F = R*U
  Eigen::Matrix3d Uinv = N*(1/l.sqrt()).matrix().asDiagonal()*N.transpose();
  S = Uinv*T*Uinv.transpose();
  
#else
  std::cerr << " *** ERROR: LogarithmicStrain requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

Tensor *
PrincipalStretches::getTMInstance()
{
  Tensor_d0s4_Ss12s34_diag *s = new Tensor_d0s4_Ss12s34_diag;
  return s;
}

Tensor *
PrincipalStretches::getStressInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
PrincipalStretches::getStrainInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
PrincipalStretches::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23_diag *B = new Tensor_d1s2_Ss23_diag(numdofs);
  return B;
}

Tensor *
PrincipalStretches::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_dense_diag *DB = new Tensor_d2s2_Sd12s34_dense_diag(numdofs);
  return DB;
}

Tensor *
PrincipalStretches::getCacheInstance()
{
  Tensor_d1s2_full *cache = new Tensor_d1s2_full(2);
  return cache;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
PrincipalStretches::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense_diag & DB = static_cast<Tensor_d2s2_Sd12s34_dense_diag &>(_DB);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();
  if(F.isIdentity()) F += std::numeric_limits<double>::epsilon()*Eigen::Matrix3d::Random(); // XXX

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &N = dec.eigenvectors();

  // principal stretches λ
  Eigen::Array<double,3,1> sqrtl = l.array().sqrt();
  e = sqrtl;

  // Moore-Penrose pseudo inverses of C - λ²_i*I
  double tol = std::numeric_limits<double>::epsilon();
  Eigen::Array<Eigen::Matrix<double,3,3>,3,1> mpinverse;
  for(int i=0; i<3; ++i) {
    Eigen::Vector3d singularValues = l - Eigen::Vector3d::Constant(l[i]);
    Eigen::Vector3d invertedSingularVals;
    for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
    mpinverse[i] = N * invertedSingularVals.asDiagonal() * N.transpose();
  }

  // first and second derivatives of λ with respect to λ²
  Eigen::Array<double,3,1> dsqrtldl = 0.5/sqrtl;
  Eigen::Array<double,3,1> d2sqrtldl2 = -0.5*dsqrtldl/l.array();

  // allocate memory for intermediate derivatives
  Eigen::Array<Eigen::Array<double,3,1>,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dCdq(numdofs);
  Eigen::Matrix3d dCdqkN, d2Cdqkdqj, dNdqk;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C with respect to q_k
    dCdq[k] = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of λ² with respect to q_k
    dCdqkN = dCdq[k]*N;
    dldq[k] = (N.transpose()*dCdqkN).diagonal();

    // first derivative of N with respect to q_k
    for(int i=0; i<3; ++i) dNdqk.col(i) = -mpinverse[i]*dCdqkN.col(i);

    for(int j=0; j<=k; ++j) {
      if((j-k)%3 == 0) {
        // second derivatives of C with respect to q_k,j
        d2Cdqkdqj = dGradUdq[j].transpose()*dGradUdq[k] + (dGradUdq[j].transpose()*dGradUdq[k]).transpose();

        // second derivatives of λ² with respect to q_k,j for non-zero d2Cdqkdqj
        d2ldqkdq[j] = (N.transpose()*(d2Cdqkdqj*N + 2*dCdq[j]*dNdqk)).diagonal();
      }
      else {
        // second derivatives of λ² with respect to q_k,j for zero d2Cdqkdqj
        d2ldqkdq[j] = (N.transpose()*(2*dCdq[j]*dNdqk)).diagonal();
      }
    }

    // first derivative of λ with respect to q_k
    B[k] = (dsqrtldl*dldq[k]).eval();
    for(int j=0; j<=k; ++j) {
      // second derivative of λ with respect to q_k,j
      DB[j*(2*numdofs-j-1)/2+k] = (dsqrtldl*d2ldqkdq[j] + d2sqrtldl2*dldq[j]*dldq[k]).eval();
    }
  }

  // store copies of N and λ² in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();

#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
PrincipalStretches::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &N = dec.eigenvectors();

  // principal stretches λ
  Eigen::Array<double,3,1> sqrtl = l.array().sqrt();
  e = sqrtl;

  // first derivative of λ with respect to λ²
  Eigen::Array<double,3,1> dsqrtldl = 0.5/sqrtl;

  // allocate memory for intermediate derivatives
  Eigen::Array<double,3,1> dldqk;
  Eigen::Matrix3d dCdqk;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C with respect to q_k
    dCdqk = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of λ² of C with respect to q_k
    dldqk = (N.transpose()*dCdqk*N).diagonal();

    // first derivative of λ with respect to q_k
    B[k] = (dsqrtldl*dldqk).eval();
  }

  // store copies of N and λ² in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();

#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
PrincipalStretches::getE(Tensor &_e, Tensor &_gradU, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();

  // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &N = dec.eigenvectors();

  // principal stretches λ
  e = l.array().sqrt().eval();

  // store copies of N and λ² in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();

#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
PrincipalStretches::transformStress(Tensor &_stress, Tensor *cache, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Array<double,3,1> > stress(&static_cast<Tensor_d0s2_Ss12_diag &>(_stress)[0]);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > N = static_cast<Tensor_d1s2_full &>(*cache)[0].matrix();
  Eigen::Array3d l = static_cast<Tensor_d1s2_full &>(*cache)[1].matrix().diagonal();

  // 2nd Piola-Kirchhoff stress S = N*diag(β/λ²)*N^T
  // note: the principal Kirchhoff stress β = λ*∂ω/∂λ where ω is the strain energy density and ∂ω/∂λ is the stress measure conjugate to λ
  S = N*(stress/l.sqrt()).matrix().asDiagonal()*N.transpose();
#else
  std::cerr << " *** ERROR: PrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

Tensor *
LogarithmicPrincipalStretches::getTMInstance()
{
  Tensor_d0s4_Ss12s34_diag *s = new Tensor_d0s4_Ss12s34_diag;
  return s;
}

Tensor *
LogarithmicPrincipalStretches::getStressInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
LogarithmicPrincipalStretches::getStrainInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
LogarithmicPrincipalStretches::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23_diag *B = new Tensor_d1s2_Ss23_diag(numdofs);
  return B;
}

Tensor *
LogarithmicPrincipalStretches::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_dense_diag *DB = new Tensor_d2s2_Sd12s34_dense_diag(numdofs);
  return DB;
}

Tensor *
LogarithmicPrincipalStretches::getCacheInstance()
{
  Tensor_d1s2_full *cache = new Tensor_d1s2_full(2);
  return cache;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicPrincipalStretches::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense_diag & DB = static_cast<Tensor_d2s2_Sd12s34_dense_diag &>(_DB);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();
  if(F.isIdentity()) F += std::numeric_limits<double>::epsilon()*Eigen::Matrix3d::Random(); // XXX

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &N = dec.eigenvectors();

  // principal stretches λ
  Eigen::Array<double,3,1> sqrtl = l.array().sqrt();

  // logarithmic principal stretches ε = ln(λ)
  e = sqrtl.log().eval();

  // Moore-Penrose pseudo inverses of C - λ²_i*I
  double tol = std::numeric_limits<double>::epsilon();
  Eigen::Array<Eigen::Matrix<double,3,3>,3,1> mpinverse;
  for(int i=0; i<3; ++i) {
    Eigen::Vector3d singularValues = l - Eigen::Vector3d::Constant(l[i]);
    Eigen::Vector3d invertedSingularVals;
    for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
    mpinverse[i] = N * invertedSingularVals.asDiagonal() * N.transpose();
  }

  // first and second derivatives of ε with respect to λ²
  Eigen::Array<double,3,1> dedl = 0.5/l.array();
  Eigen::Array<double,3,1> d2edl2 = -dedl/l.array();

  // allocate memory for intermediate derivatives
  Eigen::Array<Eigen::Array<double,3,1>,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dCdq(numdofs);
  Eigen::Matrix3d dCdqkN, dNdqk, d2Cdqkdqj;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C with respect to q_k
    dCdq[k] = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of λ² with respect to q_k
    dCdqkN = dCdq[k]*N;
    dldq[k] = (N.transpose()*dCdqkN).diagonal();

    // first derivative of N with respect to q_k
    for(int i=0; i<3; ++i) dNdqk.col(i) = -mpinverse[i]*dCdqkN.col(i);

    for(int j=0; j<=k; ++j) {
      if((j-k)%3 == 0) {
        // second derivatives of C with respect to q_k,j
        d2Cdqkdqj = dGradUdq[j].transpose()*dGradUdq[k] + (dGradUdq[j].transpose()*dGradUdq[k]).transpose();

        // second derivatives of λ² with respect to q_k,j for non-zero d2Cdqkdqj
        d2ldqkdq[j] = (N.transpose()*(d2Cdqkdqj*N + 2*dCdq[j]*dNdqk)).diagonal();
      }
      else {
        // second derivatives of λ² of C with respect to q_k,j for zero d2Cdqkdqj
        d2ldqkdq[j] = (N.transpose()*(2*dCdq[j]*dNdqk)).diagonal();
      }
    }

    // first derivative of ε with respect to q_k
    B[k] = (dedl*dldq[k]).eval();
    for(int j=0; j<=k; ++j) {
      // second derivative of ε with respect to q_k,j
      DB[j*(2*numdofs-j-1)/2+k] = (dedl*d2ldqkdq[j] + d2edl2*dldq[j]*dldq[k]).eval();
    }
  }

  // store copies of N and λ² in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();

#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
LogarithmicPrincipalStretches::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  int numdofs = dgradUdqk.getSize();

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &N = dec.eigenvectors();
  
  // principal stretches λ
  Eigen::Array<double,3,1> sqrtl = l.array().sqrt();

  // logarithmic principal stretches ε = ln(λ)
  e = sqrtl.log().eval();

  // first derivative of ε with respect to λ²
  Eigen::Array<double,3,1> dedl = 0.5/l.array();

  // allocate memory for intermediate derivatives
  Eigen::Array<double,3,1> dldqk;
  Eigen::Matrix3d dCdqk;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of C with respect to q_k
    dCdqk = F.transpose()*dGradUdq[k] + (F.transpose()*dGradUdq[k]).transpose();

    // first derivative of λ² with respect to q_k
    dldqk = (N.transpose()*dCdqk*N).diagonal();

    // first derivative of ε with respect to q_k
    B[k] = (dedl*dldqk).eval();
  }

  // store copies of N and λ² in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();

#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicPrincipalStretches::getE(Tensor &_e, Tensor &_gradU, Tensor *cache, double *)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();

  // eigenvalues λ² and eigenvectors N of right Cauchy-Green strain tensor C = F^T*F
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F.transpose()*F).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &N = dec.eigenvectors();

  // principal stretches λ
  Eigen::Array<double,3,1> sqrtl = l.array().sqrt();

  // logarithmic principal stretches ε = ln(λ)
  e = sqrtl.log().eval();

  // store copies of N and λ² in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = N;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = l.asDiagonal();

#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
LogarithmicPrincipalStretches::transformStress(Tensor &_stress, Tensor *cache, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Array<double,3,1> > beta(&static_cast<Tensor_d0s2_Ss12_diag &>(_stress)[0]);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > N = static_cast<Tensor_d1s2_full &>(*cache)[0].matrix();
  Eigen::Array3d l = static_cast<Tensor_d1s2_full &>(*cache)[1].matrix().diagonal();

  // 2nd Piola-Kirchhoff stress S = N*diag(β/λ²)*N^T
  S = N*(beta/l).matrix().asDiagonal()*N.transpose();

#else
  std::cerr << " *** ERROR: LogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

Tensor *
ElasticLogarithmicPrincipalStretches::getTMInstance()
{
  Tensor_d0s4_Ss12s34_diag *s = new Tensor_d0s4_Ss12s34_diag;
  return s;
}

Tensor *
ElasticLogarithmicPrincipalStretches::getStressInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
ElasticLogarithmicPrincipalStretches::getStrainInstance()
{
  Tensor_d0s2_Ss12_diag *s = new Tensor_d0s2_Ss12_diag;
  return s;
}

Tensor *
ElasticLogarithmicPrincipalStretches::getBInstance(int numdofs)
{
  Tensor_d1s2_Ss23_diag *B = new Tensor_d1s2_Ss23_diag(numdofs);
  return B;
}

Tensor *
ElasticLogarithmicPrincipalStretches::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_dense_diag *DB = new Tensor_d2s2_Sd12s34_dense_diag(numdofs);
  return DB;
}

Tensor *
ElasticLogarithmicPrincipalStretches::getCacheInstance()
{
  Tensor_d1s2_full *cache = new Tensor_d1s2_full(2);
  return cache;
}

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
ElasticLogarithmicPrincipalStretches::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk,
                                                 Tensor *cache, double *state)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d2s2_Sd12s34_dense_diag & DB = static_cast<Tensor_d2s2_Sd12s34_dense_diag &>(_DB);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  // plastic right Cauchy-Green strain tensor Cp = Ep+I
  Eigen::Matrix3d Cp;
  Cp << state[0]+1, state[1],   state[2],
        state[1],   state[3]+1, state[4],
        state[2],   state[4],   state[5]+1;

  int numdofs = dgradUdqk.getSize();

  // deformaton gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();
  if(F.isIdentity()) F += std::numeric_limits<double>::epsilon()*Eigen::Matrix3d::Random(); // XXX

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // eigenvalues λ² and eigenvectors n of elastic left Cauchy-Green strain tensor be = F*Cp^{-1}*F^T
  Eigen::Matrix3d Cpinv = Cp.inverse();
  Eigen::Matrix3d FCpinv = F*Cpinv;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((FCpinv*F.transpose()).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &n = dec.eigenvectors();

  // elastic logarithmic principal stretches εe = ln(λ)
  e = l.array().sqrt().log().eval();

  // Moore-Penrose pseudo inverses of be - λ²_i*I
  Eigen::Array<Eigen::Matrix<double,3,3>,3,1> mpinverse;
  double tol = std::numeric_limits<double>::epsilon();
  for(int i=0; i<3; ++i) {
    Eigen::Vector3d singularValues = l - Eigen::Vector3d::Constant(l[i]);
    Eigen::Vector3d invertedSingularVals;
    for(int j=0; j<3; ++j) invertedSingularVals[j] = (fabs(singularValues[j]) < tol) ? 0 : 1/singularValues[j];
    mpinverse[i] = n * invertedSingularVals.asDiagonal() * n.transpose();
  }

  // first and second derivatives of εe with respect to λ²
  Eigen::Array<double,3,1> dedl = 0.5/l.array();
  Eigen::Array<double,3,1> d2edl2 = -dedl/l.array();

  // allocate memory for intermediate derivatives
  Eigen::Array<Eigen::Array<double,3,1>,Eigen::Dynamic,1> dldq(numdofs), d2ldqkdq(numdofs);
  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dbedq(numdofs);
  Eigen::Matrix3d dbedqkn, dndqk;
  Eigen::Matrix<double,3,3> d2bedqkdqj;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of be with respect to q_k
    dbedq[k] = FCpinv*dGradUdq[k].transpose() + (FCpinv*dGradUdq[k].transpose()).transpose();

    // first derivative of λ² with respect to q_k
    dbedqkn = dbedq[k]*n;
    dldq[k] = (n.transpose()*dbedqkn).diagonal();

    // first derivative of n with respect to q_k
    for(int i=0; i<3; ++i) dndqk.col(i) = -mpinverse[i]*dbedqkn.col(i);

    for(int j=0; j<=k; ++j) {
      // second derivatives of be with respect to q_k,j
      d2bedqkdqj = dGradUdq[j]*Cpinv*dGradUdq[k].transpose() + (dGradUdq[j]*Cpinv*dGradUdq[k].transpose()).transpose();

      // second derivatives of λ² with respect to q_k,j
      d2ldqkdq[j] = (n.transpose()*(d2bedqkdqj*n + 2*dbedq[j]*dndqk)).diagonal();
    }

    // first derivative of εe with respect to q_k
    B[k] = (dedl*dldq[k]).eval();
    for(int j=0; j<=k; ++j) {
      // second derivative of εe with respect to q_k,j
      DB[j*(2*numdofs-j-1)/2+k] = (dedl*d2ldqkdq[j] + d2edl2*dldq[j]*dldq[k]).eval();
    }
  }

  // store copies of F and n in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = F;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = n;

#else
  std::cerr << " *** ERROR: ElasticLogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

#ifdef USE_EIGEN3
#if EIGEN_GNUC_AT_LEAST(4,7)
__attribute__((flatten))
#endif
#endif
void
ElasticLogarithmicPrincipalStretches::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk,
                                               Tensor *cache, double *state)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_Ss23_diag & B = static_cast<Tensor_d1s2_Ss23_diag &>(__B);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  // plastic right Cauchy-Green strain tensor Cp = Ep+I
  Eigen::Matrix3d Cp;
  Cp << state[0]+1, state[1],   state[2],
        state[1],   state[3]+1, state[4],
        state[2],   state[4],   state[5]+1;

  int numdofs = dgradUdqk.getSize();

  // deformation gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix()+Eigen::Matrix3d::Identity();

  Eigen::Array<Eigen::Matrix3d,Eigen::Dynamic,1> dGradUdq(numdofs);
  dgradUdqk.assignTo(dGradUdq);

  // eigenvalues (λ²) and eigenvectors (n) of elastic left Cauchy-Green strain tensor be = F*Cp^{-1}*F^T
  Eigen::Matrix3d FCpinv = F*Cp.inverse();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((FCpinv*F.transpose()).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &n = dec.eigenvectors();

  // elastic logarithmic principal stretches εe = ln(λ)
  e = l.array().sqrt().log().eval();

  // first derivative of εe with respect to λ²
  Eigen::Array<double,3,1> dedl = 0.5/l.array();

  // allocate memory for intermediate derivatives
  Eigen::Array<double,3,1> dldqk;
  Eigen::Matrix<double,3,3> dbedqk;

  for(int k=0; k<numdofs; ++k) {

    // first derivative of be with respect to q_k
    dbedqk = FCpinv*dGradUdq[k].transpose() + (FCpinv*dGradUdq[k].transpose()).transpose();

    // first derivative of λ² with respect to q_k
    dldqk = (n.adjoint()*dbedqk*n).diagonal();

    // first derivative of εe with respect to q_k
    B[k] = (dedl*dldqk).eval();
  }

  // store copies of F and n in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = F;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = n;

#else
  std::cerr << " *** ERROR: ElasticLogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
ElasticLogarithmicPrincipalStretches::getE(Tensor &_e, Tensor &_gradU, Tensor *cache, double *state)
{
#ifdef USE_EIGEN3
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  Tensor_d0s2_Ss12_diag & e = static_cast<Tensor_d0s2_Ss12_diag &>(_e);

  // plastic right Cauchy-Green strain tensor Cp = Ep+I
  Eigen::Matrix3d Cp;
  Cp << state[0]+1, state[1],   state[2],
        state[1],   state[3]+1, state[4],
        state[2],   state[4],   state[5]+1;

  // deformation gradient F = ∇u+I
  Eigen::Matrix3d F = gradU.matrix() + Eigen::Matrix3d::Identity();

  // eigenvalues λ² and eigenvectors n of elastic left Cauchy-Green strain tensor be = F*Cp^{-1}*F^T
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> dec((F*Cp.inverse()*F.transpose()).eval());
  const Eigen::Vector3d &l = dec.eigenvalues();
  const Eigen::Matrix3d &n = dec.eigenvectors();

  // elastic logarithmic principal stretches εe = ln(λ)
  e = l.array().sqrt().log().eval();

  // store copies of F and n in cache
  static_cast<Tensor_d1s2_full &>(*cache)[0] = F;
  static_cast<Tensor_d1s2_full &>(*cache)[1] = n;

#else
  std::cerr << " *** ERROR: ElasticLogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

void
ElasticLogarithmicPrincipalStretches::transformStress(Tensor &_stress, Tensor *cache, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Vector3d> beta(&static_cast<Tensor_d0s2_Ss12_diag &>(_stress)[0]);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F = static_cast<Tensor_d1s2_full &>(*cache)[0].matrix();
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > n = static_cast<Tensor_d1s2_full &>(*cache)[1].matrix();

  // 2nd Piola-Kirchhoff stress S = F^{-1}*n*diag(β)*n^T*F^{-T}
  S = F.inverse()*n*beta.asDiagonal()*n.transpose()*F.inverse().transpose();

#else
  std::cerr << " *** ERROR: ElasticLogarithmicPrincipalStretches requires AERO-S configured with Eigen library. Exiting..." << std::endl;
  exit(-1);
#endif
}

Tensor *
DeformationGradient::getTMInstance()
{
  Tensor_d0s4 *s = new Tensor_d0s4;
  return s;
}

Tensor *
DeformationGradient::getStressInstance()
{
  Tensor_d0s2 *s = new Tensor_d0s2;
  return s;
}

Tensor *
DeformationGradient::getStrainInstance()
{
  Tensor_d0s2 *s = new Tensor_d0s2;
  return s;
}

Tensor *
DeformationGradient::getBInstance(int numdofs)
{
  Tensor_d1s2_full *B = new Tensor_d1s2_full(numdofs);
  return B;
}

Tensor *
DeformationGradient::getDBInstance(int numdofs)
{
  Tensor_d2s2_Sd12s34_null *DB = new Tensor_d2s2_Sd12s34_null(numdofs);
  return DB;
}

Tensor *
DeformationGradient::getCacheInstance()
{
  Tensor_d0s2 *cache = new Tensor_d0s2;
  return cache;
}

void 
DeformationGradient::getEBandDB(Tensor &_e, Tensor &__B, Tensor &_DB, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_full & B = static_cast<Tensor_d1s2_full &>(__B);
  
  Tensor_d0s2 & e = static_cast<Tensor_d0s2 &>(_e);

  // deformation gradient F = ∇u+I
  Tensor_d0s2 identity;
  identity[0] = identity[4] = identity[8] = 1;
  static_cast<Tensor_d0s2 &>(*cache) = e = gradU + identity;

  // first derivative of F with respect to q_k
  B = dgradUdqk;
}

void 
DeformationGradient::getEandB(Tensor &_e, Tensor &__B, const Tensor &_gradU, const Tensor &_dgradUdqk, Tensor *cache, double *)
{
  const Tensor_d0s2 & gradU = static_cast<const Tensor_d0s2 &>(_gradU);
  const Tensor_d1s2_sparse & dgradUdqk = static_cast<const Tensor_d1s2_sparse &>(_dgradUdqk);
  Tensor_d1s2_full & B = static_cast<Tensor_d1s2_full &>(__B);
  
  Tensor_d0s2 & e = static_cast<Tensor_d0s2 &>(_e);

  // deformation gradient F = ∇u+I
  Tensor_d0s2 identity;
  identity[0] = identity[4] = identity[8] = 1;
  static_cast<Tensor_d0s2 &>(*cache) = e = gradU + identity;

  // first derivative of F with respect to q_k
  B = dgradUdqk;
}

void
DeformationGradient::getE(Tensor &_e, Tensor &_gradU, Tensor *cache, double *)
{
  Tensor_d0s2 & e = static_cast<Tensor_d0s2 &>(_e);
  Tensor_d0s2 & gradU = static_cast<Tensor_d0s2 &>(_gradU);

  // deformation gradient F = ∇u+I
  Tensor_d0s2 identity;
  identity[0] = identity[4] = identity[8] = 1;
  static_cast<Tensor_d0s2 &>(*cache) = e = gradU + identity;
}

void
DeformationGradient::transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S)
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > P = static_cast<Tensor_d0s2 &>(stress).matrix();
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F = static_cast<Tensor_d0s2 &>(*cache).matrix();

  // 2nd Piola-Kirchhoff stress S = F^{-1}*P
  S = F.inverse()*P;

#else
  std::cerr << "ERROR: DeformationGradient::transformStress is not implemented\n";
  S.setZero();
#endif
}

