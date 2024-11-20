#ifndef _STRAINEVALUATOR_H_
#define _STRAINEVALUATOR_H_
#include <Utils.d/NodeSpaceArray.h>

//Computes the geometrical part of the stiffness matrix 

class StrainEvaluator
{
  public:
    virtual Tensor *getTMInstance() = 0;
    virtual Tensor *getStressInstance() = 0;
    virtual Tensor *getStrainInstance() = 0;
    virtual Tensor *getBInstance(int numdofs) = 0;
    virtual Tensor *getDBInstance(int numdofs) = 0;
    virtual Tensor *getCacheInstance() = 0;
    virtual void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                            Tensor *cache, double *state=0) = 0;
    virtual void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                          Tensor *cache, double *state=0) = 0;
    virtual void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0) = 0;
    virtual void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S) = 0; // returns PK2 stress for finite-strain materials
    virtual void getL(const Tensor &gradUnp, const Tensor &gradUn, Tensor *cache, double dt=0) { // returns velocity gradient tensor
      //fprintf(stderr, " ***WARNING: Trying to calculate velocity gradient tensor for unsupported strain measure!\n");
    }
    virtual bool isStrainCompatibleWithUpdLag() { return false; };
};

class LinearStrain : public StrainEvaluator
{
  // to be used when the appropriate strain measure is the
  // infintesimal strain tensor e = 0.5*(gradU^t + gradU)
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getCacheInstance();
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                    Tensor *cache, double *state=0);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                  Tensor *cache, double *state=0);
    void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0);
    void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S);
};

class GreenLagrangeStrain : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the
  // Green-Lagrange strain tensor E = 0.5*(gradU^t + gradU + gradU^T gradU) which is a symmetric rank 2 tensor
  // Constitutive models based on this strain measure should return 
  // the second Piola-Kirchoff stress tensor S (rank 2, symmetric) and
  // the second (material) elasticity tensor M (rank 4, major and minor symmetries)
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getCacheInstance();
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                    Tensor *cache, double *state=0);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                  Tensor *cache, double *state=0);
    void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0);
    void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S);
    void getL(const Tensor &gradUnp, const Tensor &gradUn, Tensor *cache, double dt=0) override;
};

class LogarithmicStrain : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the
  // Hencky strain tensor H = log(sqrt(F^T*F)) which is a symmetric rank 2 tensor
  // Constitutive models based on this strain measure should return 
  // the rotated Kirchoff stress tensor T (rank 2, symmetric)
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getCacheInstance();
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                    Tensor *cache, double *state=0);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                  Tensor *cache, double *state=0);
    void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0);
    void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S);
};

class PrincipalStretches : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the principal stretches
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getCacheInstance();
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                    Tensor *cache, double *state=0);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                  Tensor *cache, double *state=0);
    void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0);
    void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S);
};

class LogarithmicPrincipalStretches : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the logarithmic principal stretches
  // Constitutive models based on this strain measure should return 
  // the principal Kirchoff stresses beta
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getCacheInstance();
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                    Tensor *cache, double *state=0);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                  Tensor *cache, double *state=0);
    void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0);
    void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S);
};

class ElasticLogarithmicPrincipalStretches : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the elastic logarithmic principal stretches
  // Constitutive models based on this strain measure should return 
  // the principal Kirchoff stresses beta
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getCacheInstance(int numdofs);
    Tensor *getCacheInstance();
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                    Tensor *cache, double *state=0);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                  Tensor *cache, double *state=0);
    void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0);
    void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S);
};

class DeformationGradient : public StrainEvaluator
{
  // To be used when the appropriate strain measure is the
  // deformation gradient tensor F = I + gradU, which is a nonsymmetric rank 2 tensor
  // Constitutive models based on this strain measure should return 
  // the first Piola-Kirchoff stress tensor P (rank 2, nonsymmetric)
  // and the first elasticity tensor A (rank 4, major symmetries only)
  public:
    Tensor *getTMInstance();
    Tensor *getStressInstance();
    Tensor *getStrainInstance();
    Tensor *getBInstance(int numdofs);
    Tensor *getDBInstance(int numdofs);
    Tensor *getCacheInstance();
    void getEBandDB(Tensor &e, Tensor &B, Tensor &DB, const Tensor &gradU, const Tensor &dgradUdqk,
                    Tensor *cache, double *state=0);
    void getEandB(Tensor &e, Tensor &B, const Tensor &gradU, const Tensor &dgradUdqk,
                  Tensor *cache, double *state=0);
    void getE(Tensor &e, Tensor &gradU, Tensor *cache, double *state=0);
    void transformStress(Tensor &stress, Tensor *cache, Tensor_d0s2_Ss12 &S);
};


template<class TensorTypes>
class GenStrainEvaluator
{
  public:
    virtual void getEBandDB(typename TensorTypes::StrainTensor &e, 
                            typename TensorTypes::BTensor &B, 
                            typename TensorTypes::DBTensor &DB,
                            typename TensorTypes::GradUTensor &gradU, 
                            typename TensorTypes::GradUDerivTensor &dgradUdqk) = 0;
    virtual void getEandB(typename TensorTypes::StrainTensor &e,
                          typename TensorTypes::BTensor &B,
                          typename TensorTypes::GradUTensor &gradU,
                          typename TensorTypes::GradUDerivTensor &dgradUdqk) = 0;
    virtual void getE(typename TensorTypes::StrainTensor &e, typename TensorTypes::GradUTensor &gradU) = 0;
    virtual bool isNonLinear() { return false; }
};

#endif
