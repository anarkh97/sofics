#ifndef _GAUSSINTGELEM_H_
#define _GAUSSINTGELEM_H_
#include <alloca.h>
#include <Math.d/FullSquareMatrix.h>
#include <Element.d/NLElement.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/ShapeFunction.h>
#include <Utils.d/NodeSpaceArray.h>
#include <Math.d/TTensor.h>
#include <Utils.d/dbg_alloca.h>

class StrainEvaluator;
template <class TT> class GenStrainEvaluator;

class GaussIntgElement : public MatNLElement
{
protected:

    double *tframe;
    virtual int getNumGaussPoints() const = 0;
    virtual void getGaussPointAndWeight(int, double *, double &) const = 0;
    virtual void getLocalNodalCoords(int, double *) = 0;
    virtual ShapeFunction *getShapeFunction() const = 0;
    virtual StrainEvaluator *getStrainEvaluator() const = 0;
    virtual NLMaterial *getMaterial() const = 0;

public:
    GaussIntgElement() : tframe(NULL) {}
    virtual ~GaussIntgElement() { if(tframe) delete [] tframe; }
    void getStiffAndForce(Node *nodes, double *disp,
			  double *state, FullSquareMatrix &kTan,
			  double *force) override;
    FullSquareMatrix stiffness(const CoordSet& cs, double *k, int flg=1) const override;
    FullSquareMatrix massMatrix(const CoordSet& cs, double *m, int flg=1) const override;
    void updateStates(Node *node, double *state, double *un, double *unp, double *temps, double dt) override;
    void integrate(Node *nodes, double *dispn, double *staten,
		   double *dispnp, double *statenp,
		   FullSquareMatrix &kTan, double *force, double dt,
		   double *temps) override;
    void integrate(Node *nodes, double *dispn, double *staten,
		   double *dispnp, double *statenp, double *force, double dt,
		   double *temps) override;
    int numStates() override {
 	int nGP = getNumGaussPoints();
	NLMaterial *mat = getMaterial();
	int nsGP = (mat) ? mat->getNumStates() : 0;
	return nGP*nsGP;
    }
    void initStates(double *) override;
    // the following functions return result for postprocessing at every node
    void getStrainTens(Node *nodes, double *dispnp, double (*result)[9], int avgnum) override;
    void getVonMisesStrain(Node *nodes, double *dispnp, double *result, int avgnum) override;
    void getStressTens(Node *nodes, double *dispn, double *staten,
		       double *dispnp, double *statenp, double (*result)[9], int avgnum,
		       double *temps) override;
    void getVonMisesStress(Node *nodes, double *dispn, double *staten,
			   double *dispnp, double *statenp, double *result, int avgnum,
			   double *temps) override;
    void getEquivPlasticStrain(double *statenp, double *result, int avgnum) override;
    void getBackStressTens(double *statenp, double (*result)[9], int avgnum) override;
    void getPlasticStrainTens(double *statenp, double (*result)[9], int avgnum) override;
    void getDamage(double *statenp, double *result, int avgnum) override;
    double getStrainEnergy(Node *nodes, double *dispnp, double *state, double *temps) override;
    double getDissipatedEnergy(Node *nodes, double *state) override;
    bool checkFailure(double *state) override;
};

template <class TensorTypes>
class GenGaussIntgElement : public MatNLElement
{
protected:
	// TODO Remove the mutable declarations. Computing a stiffness should not modify the element.
    mutable double *tframe;
    mutable double *tframe_inv;
    double *eframe;
	std::vector<double> preload;
	virtual int getNumGaussPoints() const = 0;
	virtual void getGaussPointAndWeight(int, double *, double &) const = 0;
    virtual void getLocalNodalCoords(int, double *) = 0;
	virtual GenShapeFunction<TensorTypes> *getShapeFunction() const = 0;
	virtual GenStrainEvaluator<TensorTypes> *getGenStrainEvaluator() const = 0;
	virtual const NLMaterial *getMaterial() const = 0;
    virtual NLMaterial *getLinearMaterial() const = 0;

public:
    GenGaussIntgElement() : tframe(NULL), tframe_inv(NULL), eframe(NULL) {}
    virtual ~GenGaussIntgElement() { if(tframe) delete [] tframe;
                                     if(tframe_inv) delete [] tframe_inv;
                                     if(eframe) delete [] eframe; }
	void getStiffAndForce(Node *nodes, double *disp,
						  double *state, FullSquareMatrix &kTan,
						  double *force) override;
	FullSquareMatrix  stiffness(const CoordSet& cs, double *k, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double *m, int flg=1) const override;
	void updateStates(Node *node, double *state, double *un, double *unp, double *temps, double dt) override;
	void integrate(Node *nodes, double *dispn, double *staten,
				   double *dispnp, double *statenp,
				   FullSquareMatrix &kTan, double *force, double dt, double *temps) override;
	void integrate(Node *nodes, double *dispn, double *staten,
				   double *dispnp, double *statenp,
				   double *force, double dt, double *temps) override;
	int numStates() override {
		int ngp = getNumGaussPoints();
		auto *mat = getMaterial();
		int nst = (mat) ? mat->getNumStates() : 0;
		return nst*ngp;
	}
	void initStates(double *) override;
	void getStrainTens(Node *nodes, double *dispnp, double (*result)[9], int avgnum) override;
	void getVonMisesStrain(Node *nodes, double *dispnp, double *result, int avgnum) override;
	void getStressTens(Node *nodes, double *dispn, double *staten,
					   double *dispnp, double *statenp, double (*result)[9], int avgnum,
					   double *temps) override;
	void getVonMisesStress(Node *nodes, double *dispn, double *staten,
						   double *dispnp, double *statenp, double *result,
						   int avgnum, double *temps) override;
	bool checkFailure(double *state) override;
};

template <class TensorTypes>
FullSquareMatrix
GenGaussIntgElement<TensorTypes>::stiffness(const CoordSet& cs, double *k, int) const
{
	int i;
	FullSquareMatrix kTan(numDofs(), k);

	Node *nodes = (Node *) dbg_alloca(numNodes()*sizeof(Node));
	int nnd = numNodes();
	int *ndn = (int *)dbg_alloca(nnd*sizeof(int));
	this->nodes(ndn);
	for(i = 0; i < nnd; ++i)
		nodes[i] = *cs[ndn[i]];
	int ndofs = numDofs();
	double *disp = (double *)dbg_alloca(ndofs*sizeof(double));
	for(i = 0; i < ndofs; ++i)
		disp[i] = 0;

	GenShapeFunction<TensorTypes> *shapeF = getShapeFunction();

	// Obtain the material model
  NLMaterial *material = getLinearMaterial();

  // Obtain the linear strain function
  GenStrainEvaluator<TensorTypes> *strainEvaluator = material->getGenStrainEvaluator();

	// Obtain the storage for gradU ( 3xndim )
	typename TensorTypes::GradUTensor gradU;

	// Obtain the storage for dgradUdqk ( ndof x3x ndim )
	typename TensorTypes::GradUDerivTensor dgradUdqk;

	// NDofsx3x3x-> 6xNDofs
	typename TensorTypes::BTensor B;

	// NdofsxNdofsx3x3x -> 6xNdofsxNdofs
	typename TensorTypes::DBTensor DB;

	typename TensorTypes::DTensor D;
	typename TensorTypes::StressTensor e;
	typename TensorTypes::StressTensor s;
	typename TensorTypes::BTensor temp1;
	Tensor_d2s0 temp2(ndofs,false);
	Tensor_d2s0 temp3(ndofs,false);

	//fprintf(stderr,"Je suis dans stiffness\n");

	int ngp = getNumGaussPoints();

	kTan.zero();
	for(i = 0; i < ngp; i++) {
		double point[3], weight, jac;
		getGaussPointAndWeight(i, point, weight);

		StackVector dispVec(disp,ndofs);
		shapeF->getGlobalGrads(gradU, dgradUdqk, &jac, nodes, point, dispVec);

		strainEvaluator->getEBandDB(e, B, DB, gradU, dgradUdqk);

    if(tframe) transformEtoC(e, tframe);

		material->integrate(&s, &D, e, e, 0, 0, 0, 0);

		for(int j=0; j<preload.size(); ++j) s[j] += preload[j];

    if(tframe) {
      transformCtoE(s, tframe_inv);
      transformCtoE(D, tframe, tframe_inv);
    }

		temp1 = dblContractTransp(D,B);
		temp2 =   B||temp1;
		// We should test if the strain is non linear to call the following two lines
		if(strainEvaluator->isNonLinear()) {
			SimpleTensor<SimpleTensor<double,TensorTypes::ndofs>,TensorTypes::ndofs> kg;
			kg=DB || s;
			for(int j=0; j < TensorTypes::ndofs; ++j)
				for(int k=0; k < TensorTypes::ndofs; ++k)
					temp3[j*TensorTypes::ndofs+k] = kg[j][k];
			temp2 = temp3 + temp2;
		}
		temp3 = (weight * fabs(jac))*temp2;
		kTan += temp3;
	}

	return kTan;
}

template <class TensorType>
FullSquareMatrix
GenGaussIntgElement<TensorType>::massMatrix(const CoordSet&, double *m, int) const
{
	return FullSquareMatrix(numDofs(), m);
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::getStiffAndForce(Node *nodes, double *disp,
												  double *state, FullSquareMatrix &kTan,
												  double *force)
{
	int ndofs = numDofs();

	GenShapeFunction<TensorType> *shapeF = getShapeFunction();

	// Obtain the strain function. It can be linear or non-linear
	GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

	// Obtain the material model
	const NLMaterial *material = getMaterial();

	// Obtain the storage for gradU ( 3x3 )
	typename TensorType::GradUTensor gradU;

	// Obtain the storage for dgradUdqk ( ndof x3x3 )
	typename TensorType::GradUDerivTensor dgradUdqk;

	// NDofsx3x3x-> 6xNDofs
	typename TensorType::BTensor B;

	// NdofsxNdofsx3x3x -> 6xNdofsxNdofs but diagonal in dofs
	typename TensorType::DBTensor DB;

	typename TensorType::DTensor D;
	typename TensorType::StrainTensor e;
	typename TensorType::StressTensor s;

	SimpleTensor<double, TensorType::ndofs> temp0, nodeforce;
	typename TensorType::BTensor temp1;
	Tensor_d2s0 temp2(ndofs,false);
	Tensor_d2s0 temp3(ndofs,false);

	//fprintf(stderr,"Je suis dans getStiffAndForce\n");

	int i,j;
	int ngp = getNumGaussPoints();

	kTan.zero();

	for(i = 0; i < ngp; i++) {
		double point[3], weight, jac;
		getGaussPointAndWeight(i, point, weight);

		StackVector dispVec(disp,ndofs);

		shapeF->getGlobalGrads(gradU, dgradUdqk,  &jac, nodes, point, dispVec);
		strainEvaluator->getEBandDB(e, B, DB, gradU, dgradUdqk);

		material->integrate(&s, &D,  e, e, 0, 0, 0, 0);

		temp0 = s || B;
		temp0 = (weight*jac)*temp0;
		nodeforce = nodeforce + temp0;
		temp1 =  dblContractTransp(D, B);
		temp2 =  B||temp1;
		if(strainEvaluator->isNonLinear()) {
			SimpleTensor<SimpleTensor<double,TensorType::ndofs>,TensorType::ndofs> kg;
			kg=DB || s;
			for(int j=0; j < TensorType::ndofs; ++j)
				for(int k=0; k < TensorType::ndofs; ++k)
					temp3[j*TensorType::ndofs+k] = kg[j][k];
			temp3 = temp3 + temp2;
		}
		temp3 = (weight * fabs(jac))*temp3;
		kTan += temp3;
	}

	for(j = 0; j < ndofs; ++j) {
		force[j] = - nodeforce[j];
	}

	////////////////////////////////////////////////////////
	//Check convergence with approximate tangent stiffness//
	////////////////////////////////////////////////////////

	typename TensorType::StrainTensor eleft, eright, sleft, sright;
	typename TensorType::DBTensor DBleft, DBright;
	typename TensorType::BTensor Bleft, Bright;
	typename TensorType::GradUTensor gradUleft;
	typename TensorType::GradUTensor gradUright;
	typename TensorType::GradUDerivTensor dgradUdqkleft;
	typename TensorType::GradUDerivTensor dgradUdqkright;

	double *kt = new double[ndofs*ndofs];
	FullSquareMatrix kTant(ndofs, kt);

	int a,b;
	double dispsqNorm =0.0;
	for (b=0; b < ndofs; ++b)
		dispsqNorm+=disp[b]*disp[b];
	double epsilon = 1e-8; //(1e-8)*sqrt(dispsqNorm);
	//if(dispsqNorm == 0.0) epsilon = 1e-8;
	kTant.zero();
	fprintf(stderr, "epsilon = %e\n" , epsilon);
	double* dispLeft = new double[ndofs];
	double* dispRight = new double[ndofs];
	for (b=0; b < ndofs; ++b){
		dispLeft[b] = disp[b];
		dispRight[b] = disp[b];
	}

	for (b=0; b < ndofs; ++b) {

		if(disp[b] > 1 || disp[b] < -1)
			epsilon = 1e-8*disp[b];
		else
			epsilon = 1e-8;
		dispLeft[b] = disp[b]+epsilon;
		dispRight[b] = disp[b]-epsilon;
		SimpleTensor<double, TensorType::ndofs> nodeforceLeft, nodeforceRight;
		StackVector dispVecTestLeft(dispLeft,ndofs);
		StackVector dispVecTestRight(dispRight,ndofs);

		for(i = 0; i < ngp; i++) {
			double point[3], weight, jac;
			getGaussPointAndWeight(i, point, weight);

			shapeF->getGlobalGrads(gradUleft, dgradUdqkleft,  &jac, nodes, point, dispVecTestLeft);
			strainEvaluator->getEBandDB(eleft, Bleft, DBleft, gradUleft, dgradUdqkleft);
			material->integrate(&sleft, &D, eleft, eleft, 0, 0, 0, 0);
			temp0 = sleft || Bleft;
			temp0 = (weight*jac)*temp0;
			nodeforceLeft = nodeforceLeft + temp0;

			shapeF->getGlobalGrads(gradUright, dgradUdqkright,  &jac, nodes, point, dispVecTestRight);
			strainEvaluator->getEBandDB(eright, Bright, DBright, gradUright, dgradUdqkright);
			material->integrate(&sright, &D,  eright, eright, 0, 0, 0, 0);
			temp0 = sright || Bright;
			temp0 = (weight*jac)*temp0;
			nodeforceRight = nodeforceRight + temp0;
		}

		for (a=0; a < ndofs; ++a) {
			kTant[a][b] = (nodeforceLeft[a] - nodeforceRight[a])/(2*epsilon);
		}

		dispLeft[b] = disp[b];
		dispRight[b] = disp[b];
	}

	double diffsqNorm = 0.0;
	double sqNorm = 0.0;
	double relativeNorm;

	for (a=0; a < ndofs; ++a)
		for (b=0; b < ndofs; ++b) {
			diffsqNorm+=(kTan[a][b]-kTant[a][b])*(kTan[a][b]-kTant[a][b]);
			sqNorm+=kTant[a][b]*kTant[a][b];
		}
	relativeNorm = sqrt(diffsqNorm/sqNorm);
	fprintf(stderr, "Relative Norm = %e\n", relativeNorm);
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::updateStates(Node *nodes, double *state, double *un, double *unp, double *temps, double dt) {}
/*
void
GaussIntgElement::updateStates(Node *nodes, double *state, double *un,double *unp)
{
  StrainEvaluator *strainEvaluator = getStrainEvaluator();
  ShapeFunction *shapeF = getShapeFunction();
  NLMaterial *material = getMaterial();

  Tensor &gradUn = *shapeF->getGradUInstance();
  Tensor &gradUnp = *shapeF->getGradUInstance();

  Tensor &en = *strainEvaluator->getStrainInstance();
  Tensor &enp = *strainEvaluator->getStrainInstance();

  int ndofs = numDofs();
  int ngp  = getNumGaussPoints();
  int nstatepgp = material->getNumStates();

  double point[3], weight;

  StackVector dispn(un,ndofs);
  StackVector dispnp(unp,ndofs);

  for (int i=0; i<ngp; ++i) {
	getGaussPointAndWeight(i, point, weight);
	shapeF->getGradU(&gradUn, nodes, point, dispn);
	shapeF->getGradU(&gradUnp, nodes, point, dispnp);
	strainEvaluator->getE(en, gradUn);
	strainEvaluator->getE(enp, gradUnp);
	material->updateStates(en, enp, state + nstatepgp*i);
  };
}
*/

inline void
transformEtoC(Stress2D &enp, double *tframe)
{
  // transform stress or strain from element frame (EFRAME) to material frame (CFRAME)
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>> T(tframe);
  Eigen::Vector3d v;
  v << enp[0], enp[1], enp[2];
  v = T*v;
  enp[0] = v[0];
  enp[1] = v[1];
  enp[2] = v[2];
#else
  std::cerr << " *** ERROR: transformEtoC in GaussIntgElem.h requires AERO-S to be configured with the Eigen library. Exiting...\n";
  exit(-1);
#endif
}

inline void
transformCtoE(Stress2D &s, double *tframe_inv)
{
  // transform stress or strain from material frame (CFRAME) to element frame (EFRAME)
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>> Tinv(tframe_inv);
  Eigen::Vector3d v;
  v << s[0], s[1], s[2];
  v = Tinv*v;
  s[0] = v[0];
  s[1] = v[1];
  s[2] = v[2];
#else
  std::cerr << " *** ERROR: transformCtoE in GaussIntgElem.h requires AERO-S to be configured with the Eigen library. Exiting...\n";
  exit(-1);
#endif
}

inline void
transformCtoE(SymTensor<Stress2D,2> &Dnp, double *tframe, double *tframe_inv)
{
  // transform constitutive tangent from material frame (CFRAME) to element frame (EFRAME)
  Eigen::Matrix3d D;
  D << Dnp[0][0], Dnp[0][1], 2*Dnp[0][2],
       Dnp[1][0], Dnp[1][1], 2*Dnp[1][2],
       Dnp[2][0], Dnp[2][1], 2*Dnp[2][2];
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>> T(tframe), Tinv(tframe_inv);
  D = Tinv*D*T;
  D.col(2) /= 2.;
  for(int i = 0; i < 3; ++i) for(int j = 0; j < 3; ++j) Dnp[i][j] = D(i,j);
}

inline void
transformEtoG(Stress2D &s, double *eframe, double *data)
{
  // transform stress or strain from element frame (EFRAME) to global frame
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix3d> E(eframe);
  Eigen::Matrix3d Se;
  Se << s[0], s[2], 0,
        s[2], s[1], 0,
           0,    0, 0;
  Eigen::Map< Eigen::Matrix<double,3,3,Eigen::RowMajor> > Sg(data, 3, 3);
  Sg = E*Se*E.transpose();
#else
  std::cerr << " *** ERROR: transformEtoG in GaussIntgElem.h requires AERO-S to be configured with the Eigen library. Exiting...\n";
  exit(-1);
#endif
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::integrate(Node *nodes, double *dispn,  double *staten,
										   double *dispnp, double *statenp,
										   FullSquareMatrix &kTan,
                                           double *force, double dt, double *temps)
{
	int ndofs = numDofs();
	GenShapeFunction<TensorType> *shapeF = getShapeFunction();

	// Obtain the strain function. It can be linear or non-linear
	GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

	// Obtain the material model
	const NLMaterial *material = getMaterial();

	// Obtain the storage for gradU ( 3x3 )
	typename TensorType::GradUTensor gradUn;
	typename TensorType::GradUTensor gradUnp;
	// Obtain the storage for dgradUdqk ( ndof x3x3 )
	typename TensorType::GradUDerivTensor dgradUdqkn;
	typename TensorType::GradUDerivTensor dgradUdqknp;

	// NDofsx3x3x-> 6xNDofs
	typename TensorType::BTensor Bn;
	typename TensorType::BTensor Bnp;

	// NdofsxNdofsx3x3x -> 6xNdofsxNdofs but sparse
	typename TensorType::DBTensor DBn;
	typename TensorType::DBTensor DBnp;

	typename TensorType::DTensor Dnp;
	typename TensorType::StressTensor en;
	typename TensorType::StressTensor enp;
	typename TensorType::StressTensor s;

	SimpleTensor<double,TensorType::ndofs> temp0;
	SimpleTensor<double,TensorType::ndofs> nodeforce;
	typename TensorType::BTensor temp1;
	Tensor_d2s0 temp2(ndofs,false);
	Tensor_d2s0 temp3(ndofs,false);

	int i,j;
	int ngp = getNumGaussPoints();
	int nstatepgp = material->getNumStates();

	kTan.zero();
	nodeforce = 0;

	for(i = 0; i < ngp; i++) {

		double point[3], weight, jacn, jacnp, tempnp;
		//StackVector dispVecn(dispn,ndofs);
		StackVector dispVecnp(dispnp,ndofs);

		getGaussPointAndWeight(i, point, weight);

		//shapeF->getGlobalGrads(gradUn, dgradUdqkn,  &jacn, nodes, point, dispVecn);
		shapeF->getGlobalGrads(gradUnp, dgradUdqknp,  &jacnp, nodes, point, dispVecnp);

		//strainEvaluator->getEBandDB(en, Bn, DBn, gradUn, dgradUdqkn);
		strainEvaluator->getEBandDB(enp, Bnp, DBnp, gradUnp, dgradUdqknp);

		// compute temperature at integration point
		tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    //if(tframe) transformEtoC(en, tframe);
    if(tframe) transformEtoC(enp, tframe);

		material->integrate(&s, &Dnp, en, enp,
                        staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, 0, dt);

		for(int j=0; j<preload.size(); ++j) s[j] += preload[j]; // note: for membrane element preload should have units of
		// force per unit length (ie. prestress multiplied by thickness)
    if(tframe) {
      transformCtoE(s, tframe_inv);
      transformCtoE(Dnp, tframe, tframe_inv);
    }

		temp0 = s || Bnp;
		temp0 = (weight*jacnp)*temp0;
		nodeforce = nodeforce + temp0;
		temp1 = dblContractTransp(Dnp, Bnp);
		temp2 =   Bnp||temp1;
		if(strainEvaluator->isNonLinear()) {
			SimpleTensor<SimpleTensor<double,TensorType::ndofs>,TensorType::ndofs> kg;
			kg=DBnp || s;
			for(int j=0; j < TensorType::ndofs; ++j)
				for(int k=0; k < TensorType::ndofs; ++k)
					temp3[j*TensorType::ndofs+k] = kg[j][k];
			temp2 = temp3 + temp2;
		}
		temp3 = (weight * fabs(jacnp))*temp2;
		kTan += temp3;
	}

	for(j = 0; j < ndofs; ++j) {
		force[j] = - nodeforce[j];}
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::integrate(Node *nodes, double *dispn,  double *staten,
										   double *dispnp, double *statenp,
                                           double *force, double dt, double *temps)
{
	int ndofs = numDofs();
	GenShapeFunction<TensorType> *shapeF = getShapeFunction();

	// Obtain the strain function. It can be linear or non-linear
	GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

	// Obtain the material model
	const NLMaterial *material = getMaterial();

	// Obtain the storage for gradU ( 3x3 )
	typename TensorType::GradUTensor gradUn;
	typename TensorType::GradUTensor gradUnp;
	// Obtain the storage for dgradUdqk ( ndof x3x3 )
	typename TensorType::GradUDerivTensor dgradUdqkn;
	typename TensorType::GradUDerivTensor dgradUdqknp;

	// NDofsx3x3x-> 6xNDofs
	typename TensorType::BTensor Bn;
	typename TensorType::BTensor Bnp;

	typename TensorType::StressTensor en;
	typename TensorType::StressTensor enp;
	typename TensorType::StressTensor s;

	SimpleTensor<double,TensorType::ndofs> temp0;
	SimpleTensor<double,TensorType::ndofs> nodeforce;

	int i,j;
	int ngp = getNumGaussPoints();
	int nstatepgp = material->getNumStates();

	nodeforce = 0;

	for(i = 0; i < ngp; i++) {

		double point[3], weight, jacn, jacnp, tempnp;
		//StackVector dispVecn(dispn,ndofs);
		StackVector dispVecnp(dispnp,ndofs);

		getGaussPointAndWeight(i, point, weight);

		//shapeF->getGlobalGrads(gradUn, dgradUdqkn,  &jacn, nodes, point, dispVecn);
		shapeF->getGlobalGrads(gradUnp, dgradUdqknp,  &jacnp, nodes, point, dispVecnp);

		//strainEvaluator->getEandB(en, Bn, gradUn, dgradUdqkn);
		strainEvaluator->getEandB(enp, Bnp, gradUnp, dgradUdqknp);

		// compute temperature at integration point
		tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    //if(tframe) transformEtoC(en, tframe);
    if(tframe) transformEtoC(enp, tframe);

		material->integrate(&s, en, enp,
                        staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, 0, dt);

		for(int j=0; j<preload.size(); ++j) s[j] += preload[j]; // note: for membrane element preload should have units of
		// force per unit length (ie. prestress multiplied by thickness)
    if(tframe) transformCtoE(s, tframe_inv);

		temp0 = s || Bnp;
		temp0 = (weight*jacnp)*temp0;
		nodeforce = nodeforce + temp0;
	}

	for(j = 0; j < ndofs; ++j) {
		force[j] = - nodeforce[j];}
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::initStates(double *st)
{
	auto *material = getMaterial();
	int ninterns = material->getNumStates();
	int ngp = getNumGaussPoints();

	for(int i = 0; i < ngp; ++i)
		const_cast<NLMaterial *>(material)->initStates(st+i*ninterns);
}

inline void
copyTens(Stress2D *stens, double *svec)
{
	svec[0] = (*stens)[0]; // sxx
	svec[1] = svec[3] = (*stens)[2]; // sxy
	svec[2] = svec[6] = 0; // sxz
	svec[4] = (*stens)[1]; // syy
	svec[5] = svec[7] = 0; // syz
	svec[8] = 0; // szz
}

#ifdef USE_EIGEN3
inline void
copyTens(Stress2D *stens, Eigen::Matrix3d &smat)
{
	smat << (*stens)[0], (*stens)[2], 0,
			(*stens)[2], (*stens)[1], 0,
			0,           0, 0;
}
#endif

template <class TensorType>
void
GenGaussIntgElement<TensorType>::getStrainTens(Node *nodes, double *dispnp, double (*result)[9],
											   int avgnum)
{
	int ndofs = numDofs();
	GenShapeFunction<TensorType> *shapeF = getShapeFunction();

	// Obtain the strain function. It can be linear or non-linear
	GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

	// Obtain the storage for gradU ( 3x3 )
	typename TensorType::GradUTensor gradUnp;

	typename TensorType::StressTensor enp;

	int i,j;
	int ngp = getNumGaussPoints();

	// evaluate at the gauss points and extrapolate to the nodes
	double (*gpstrain)[9] = new double[getNumGaussPoints()][9];

	for(i = 0; i < ngp; i++) {

		double point[3], weight;

		StackVector dispVecnp(dispnp,ndofs);

		getGaussPointAndWeight(i, point, weight);

		shapeF->getGradU(gradUnp, nodes, point, dispVecnp);

		strainEvaluator->getE(enp, gradUnp);

		if(avgnum == -1) {
      transformEtoG(enp, eframe, result[i]);
		}
		else {
      transformEtoG(enp, eframe, gpstrain[i]);
		}
	}

	if(avgnum != -1) {
		//TODO extrapolate from gauss points to nodes
		//temporary fix implemented here is to return the average of all the gauss points for each node
		double average[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		for(int i = 0; i < getNumGaussPoints(); i++)
			for(int j = 0; j < 9; ++j)
				average[j] += gpstrain[i][j];
		for(int i = 0; i < 9; ++i)
			average[i] /= getNumGaussPoints();
		for(int i = 0; i < numNodes(); ++i)
			for(int j = 0; j < 9; ++j)
				result[i][j] = average[j];
	}

	delete [] gpstrain;
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::getVonMisesStrain(Node *nodes, double *dispnp, double *result,
												   int avgnum)
{
	int ndofs = numDofs();
	GenShapeFunction<TensorType> *shapeF = getShapeFunction();

	// Obtain the strain function. It can be linear or non-linear
	GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

	// Obtain the storage for gradU ( 3x3 )
	typename TensorType::GradUTensor gradUnp;

	typename TensorType::StressTensor enp;

	int i,j;
	int ngp = getNumGaussPoints();

	// evaluate at the gauss points and extrapolate to the nodes
	double *gpstrain = new double[getNumGaussPoints()];

	for(i = 0; i < ngp; i++) {

		double point[3], weight, jacn, jacnp, tempnp;
		StackVector dispVecnp(dispnp,ndofs);

		getGaussPointAndWeight(i, point, weight);

		shapeF->getGradU(gradUnp, nodes, point, dispVecnp);

		strainEvaluator->getE(enp, gradUnp);

    // note: it is not necessary to transform the strain from element to global frame because the von Mises strain is frame invariant

#ifdef USE_EIGEN3
		Eigen::Matrix3d M;
		copyTens(&enp, M);
		// compute the deviatoric stress/strain tensor and it's second invariant
		Eigen::Matrix3d dev = M - (M.trace()/3)*Eigen::Matrix3d::Identity();
		double J2 = 0.5*(dev*dev).trace();
		// compute the effective stress/strain
		if(avgnum == -1) result[i] = sqrt(3*J2); else gpstrain[i] = sqrt(3*J2);
#else
		if(avgnum == -1) result[i] = 0; else gpstrain[i] = 0;
#endif
	}

	if(avgnum != -1) {
		//TODO extrapolate from gauss points to nodes
		//temporary fix implemented here is to return the average of all the gauss points for each node
		double average = 0;
		for(int i = 0; i < getNumGaussPoints(); i++)
			average += gpstrain[i];
		average /= getNumGaussPoints();
		for(int i = 0; i < numNodes(); ++i)
			result[i] = average;
	}

	delete [] gpstrain;
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::getStressTens(Node *nodes, double *dispn, double *staten,
											   double *dispnp, double *statenp, double (*result)[9],
											   int avgnum, double *temps)
{
	int ndofs = numDofs();
	GenShapeFunction<TensorType> *shapeF = getShapeFunction();

	// Obtain the strain function. It can be linear or non-linear
	GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

	// Obtain the material model
	const NLMaterial *material = getMaterial();

	// Obtain the storage for gradU ( 3x3 )
	typename TensorType::GradUTensor gradUn;
	typename TensorType::GradUTensor gradUnp;
	// Obtain the storage for dgradUdqk ( ndof x3x3 )
	typename TensorType::GradUDerivTensor dgradUdqkn;
	typename TensorType::GradUDerivTensor dgradUdqknp;

	// NDofsx3x3x-> 6xNDofs
	typename TensorType::BTensor Bn;
	typename TensorType::BTensor Bnp;

	typename TensorType::StressTensor en;
	typename TensorType::StressTensor enp;
	typename TensorType::StressTensor s;

	int i,j;
	int ngp = getNumGaussPoints();
	int nstatepgp = material->getNumStates();

	// evaluate at the gauss points and extrapolate to the nodes
	double (*gpstress)[9] = new double[getNumGaussPoints()][9];
	double t = material->getThickness();

	for(i = 0; i < ngp; i++) {

		double point[3], weight, jacn, jacnp, tempnp;
		//StackVector dispVecn(dispn,ndofs);
		StackVector dispVecnp(dispnp,ndofs);

		getGaussPointAndWeight(i, point, weight);

		//shapeF->getGlobalGrads(gradUn, dgradUdqkn,  &jacn, nodes, point, dispVecn);
		shapeF->getGlobalGrads(gradUnp, dgradUdqknp,  &jacnp, nodes, point, dispVecnp);

		//strainEvaluator->getEandB(en, Bn, gradUn, dgradUdqkn);
		strainEvaluator->getEandB(enp, Bnp, gradUnp, dgradUdqknp);

		// compute temperature at integration point
		tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    //if(tframe) transformEtoC(en, tframe);
    if(tframe) transformEtoC(enp, tframe);

		material->integrate(&s, en, enp,
                        staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, 0);

		for(int j=0; j<preload.size(); ++j) s[j] += preload[j]; // note: for membrane element preload should have units of
		// force per unit length (ie. prestress multiplied by thickness)
    if(tframe) transformCtoE(s, tframe_inv);

		if(avgnum == -1) {
      transformEtoG(s, eframe, result[i]);
			if(t > 0) for(int j=0; j<9; ++j) result[i][j] /= t;
		}
		else {
      transformEtoG(s, eframe, gpstress[i]);
			if(t > 0) for(int j=0; j<9; ++j) gpstress[i][j] /= t;
		}
	}

	if(avgnum != -1) {
		//TODO extrapolate from gauss points to nodes
		//temporary fix implemented here is to return the average of all the gauss points for each node
		double average[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		for(int i = 0; i < getNumGaussPoints(); i++)
			for(int j = 0; j < 9; ++j)
				average[j] += gpstress[i][j];
		for(int i = 0; i < 9; ++i)
			average[i] /= getNumGaussPoints();
		for(int i = 0; i < numNodes(); ++i)
			for(int j = 0; j < 9; ++j)
				result[i][j] = average[j];
	}

	delete [] gpstress;
}

template <class TensorType>
void
GenGaussIntgElement<TensorType>::getVonMisesStress(Node *nodes, double *dispn, double *staten,
												   double *dispnp, double *statenp, double *result,
												   int avgnum, double *temps)
{
	int ndofs = numDofs();
	GenShapeFunction<TensorType> *shapeF = getShapeFunction();

	// Obtain the strain function. It can be linear or non-linear
	GenStrainEvaluator<TensorType> *strainEvaluator = getGenStrainEvaluator();

	// Obtain the material model
	const NLMaterial *material = getMaterial();

	// Obtain the storage for gradU ( 3x3 )
	typename TensorType::GradUTensor gradUn;
	typename TensorType::GradUTensor gradUnp;
	// Obtain the storage for dgradUdqk ( ndof x3x3 )
	typename TensorType::GradUDerivTensor dgradUdqkn;
	typename TensorType::GradUDerivTensor dgradUdqknp;

	// NDofsx3x3x-> 6xNDofs
	typename TensorType::BTensor Bn;
	typename TensorType::BTensor Bnp;

	typename TensorType::StressTensor en;
	typename TensorType::StressTensor enp;
	typename TensorType::StressTensor s;

	int i,j;
	int ngp = getNumGaussPoints();
	int nstatepgp = material->getNumStates();

	// evaluate at the gauss points and extrapolate to the nodes
	double *gpstress = new double[getNumGaussPoints()];
	double t = material->getThickness();

	for(i = 0; i < ngp; i++) {

		double point[3], weight, jacn, jacnp, tempnp;
		//StackVector dispVecn(dispn,ndofs);
		StackVector dispVecnp(dispnp,ndofs);

		getGaussPointAndWeight(i, point, weight);

		//shapeF->getGlobalGrads(gradUn, dgradUdqkn,  &jacn, nodes, point, dispVecn);
		shapeF->getGlobalGrads(gradUnp, dgradUdqknp,  &jacnp, nodes, point, dispVecnp);

		//strainEvaluator->getEandB(en, Bn, gradUn, dgradUdqkn);
		strainEvaluator->getEandB(enp, Bnp, gradUnp, dgradUdqknp);

		// compute temperature at integration point
		tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    //if(tframe) transformEtoC(en, tframe);
    if(tframe) transformEtoC(enp, tframe);

		material->integrate(&s, en, enp,
							staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, 0); // XXX note should call getStress followed by transform

		for(int j=0; j<preload.size(); ++j) s[j] += preload[j]; // note: for membrane element preload should have units of
		// force per unit length (ie. prestress multiplied by thickness)

    // note: it is not necessary to transform the stress from material to global frame because the von Mises stress is frame invariant

#ifdef USE_EIGEN3
		Eigen::Matrix3d M;
		copyTens(&s, M);
		if(t > 0) M /= t;
		// compute the deviatoric stress/strain tensor and it's second invariant
		Eigen::Matrix3d dev = M - (M.trace()/3)*Eigen::Matrix3d::Identity();
		double J2 = 0.5*(dev*dev).trace();
		// compute the effective stress/strain
		if(avgnum == -1) result[i] = sqrt(3*J2); else gpstress[i] = sqrt(3*J2);
#else
		if(avgnum == -1) result[i] = 0; else gpstress[i] = 0;
#endif
	}

	if(avgnum != -1) {
		//TODO extrapolate from gauss points to nodes
		//temporary fix implemented here is to return the average of all the gauss points for each node
		double average = 0;
		for(int i = 0; i < getNumGaussPoints(); i++)
			average += gpstress[i];
		average /= getNumGaussPoints();
		for(int i = 0; i < numNodes(); ++i)
			result[i] = average;
	}

	delete [] gpstress;
}

template <class TensorType>
bool
GenGaussIntgElement<TensorType>::checkFailure(double *statenp)
{
	const NLMaterial *material = getMaterial();
	for(int i = 0; i < getNumGaussPoints(); i++) {
		if(material->getDamage(statenp + i*material->getNumStates()) < 1) return false;
	}
	return true;
}

#endif
