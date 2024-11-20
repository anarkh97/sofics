#ifndef _OPMAKE_C_
#define _OPMAKE_C_
#include <typeinfo>
#include <algorithm>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#ifdef USE_EIGEN3
#include <Math.d/EiSparseMatrix.h>
#endif
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Driver.d/GeoSource.h>
#include <Solvers.d/MappedAssembledSolver.h>
#include <Solvers.d/PCGSolver.h>
#include <Solvers.d/BCGSolver.h>
#include <Solvers.d/CRSolver.h>
#include <Driver.d/Dynam.h>
#include <Sfem.d/Sfem.h>
#include <Sfem.d/SfemBlockMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/DiagMatrix.h>
#include <Feti.d/DistrVector.h>
#include <Hetero.d/FlExchange.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/MpcElement.d/MpcElement.h>
#include <Utils.d/MFTT.h>
#include <Utils.d/MathUtils.h>
#include <Utils.d/DistHelper.h>
#include <Solvers.d/SolverFactory.h>
#include <Control.d/ControlInterface.h>
#include <Driver.d/SubDomain.h>
#include <Rom.d/VecNodeDof6Conversion.h>
#include <Rom.d/FileNameInfo.h>
#include <Rom.d/BasisFileStream.h>

extern Sfem* sfem;
extern int verboseFlag;
#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

template<class Scalar>
void
Domain::makeSparseOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
                      double Ccoef, GenSparseMatrix<Scalar> *mat,
                      FullSquareMatrix *kelArray, FullSquareMatrix *melArray, FullSquareMatrix *celArray)
{
	if(matrixTimers) matrixTimers->memoryForm -= memoryUsed();
	int makeMass = (Mcoef != 0 || ops.M != 0 || ops.C != 0);

	// Rayleigh damping coefficients: C = alpha*M + beta*K
	double alphaDamp = sinfo.alphaDamp, alpha;
	double  betaDamp = sinfo.betaDamp, beta;
	// Structural damping coefficient
	double etaDamp = sinfo.etaDamp, eta;

	int size = sizeof(double)*maxNumDOFs*maxNumDOFs;

	double *karray = new double[std::max(maxNumDOFs*maxNumDOFs,maxNumDOFsFluid*maxNumDOFsFluid)];
	double *ikarray = new double[maxNumDOFs*maxNumDOFs];
	double *marray = new double[maxNumDOFs*maxNumDOFs];

	FullSquareMatrix kel, mel, cel;

	double mratio = geoSource->getMRatio(); // 0 for lumped, 1 for consistent mass matrix

	complex<double> *kcarray = new complex<double>[4*maxNumDOFs*maxNumDOFs];
	complex<double> *mcarray = new complex<double>[4*maxNumDOFs*maxNumDOFs];

	FullSquareMatrixC kcel(1, kcarray);
	FullSquareMatrixC mcel(1, mcarray);
	FullSquareMatrix ikel(1, ikarray);

	bool isShifted = geoSource->isShifted();
	bool isDamped = (alphaDamp != 0.0) || (betaDamp != 0.0) || packedEset.hasDamping();

	double *izarray = (isShifted || isDamped) ? static_cast<double *>(dbg_alloca(size)) : 0;
	FullSquareMatrix izel(1, izarray); // izel stores the imaginary part of impedence matrix zel
	double omega, omega2;

	int sizeC = sizeof(DComplex)*maxNumDOFs*maxNumDOFs;
	DComplex *karrayC = new DComplex[maxNumDOFs*maxNumDOFs];
	DComplex *marrayC = new DComplex[maxNumDOFs*maxNumDOFs];
	FullSquareMatrixC kelC(1, karrayC);
	FullSquareMatrixC melC(1, marrayC);

	if(sinfo.isCoupled) computeCoupledScaleFactors();

	GenSubDomain<Scalar> *subCast = dynamic_cast<GenSubDomain<Scalar>*>(this);
	if(!subCast && (sinfo.ATDARBFlag >= 0.0 || sinfo.ATDROBalpha != 0.0)) checkSommerTypeBC(this);
	bool mdds_flag = (mat && subCast && sinfo.solvercntl->type == SolverSelection::Direct && sinfo.solvercntl->subtype != 13); // multidomain direct solver, dont use for ROM

	bool zeroRot = (sinfo.zeroRot && sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0);
	int *dofType = (zeroRot) ? dsa->makeDofTypeArray() : 0;

	if(sinfo.farfield) { addSBoundNodes(); makeKss(this); } // for Farfield output (TODO check with Radek)

	int iele;

	// LOOP over elements except for sommer
	for(iele = 0; iele < numele; ++iele) {
		StructProp *prop = packedEset[iele]->getProperty();
		if(packedEset[iele]->isSommerElement()) continue;
		if(prop && (prop->E0!=0.0 || prop->mu0!=0.0)) continue; // rubber handled below
		bool isComplexF = (prop && prop->fp.PMLtype != 0);
		if(packedEset[iele]->isConstraintElement()) { alpha = beta = eta = 0; }
		else {
			alpha = (packedEset[iele]->isDamped()) ? prop->alphaDamp : alphaDamp;
			beta = (packedEset[iele]->isDamped()) ? prop->betaDamp : betaDamp;
			eta = (packedEset[iele]->isSDamped()) ? prop->etaDamp : etaDamp;
		}
		bool isSDamped = (eta != 0.0);
		complex<double> kappa2 = packedEset[iele]->helmCoefC();
		omega2 = geoSource->shiftVal();
		omega = sqrt(omega2);

		// Form element real and complex stiffness matrices in kel and kcel
		if(matrixTimers) matrixTimers->formTime -= getTime();
		if(isComplexF) {
			kcel = packedEset[iele]->stiffness(nodes, kcarray);
		}
		else {
			if(kelArray) kel.copy(kelArray[iele]);
			else kel = packedEset[iele]->stiffness(nodes, karray);
			this->densProjectStiffness(kel, iele);
			this->transformMatrix(kel, iele);
		}
		if(sinfo.isCoupled) {
			if(isStructureElement(iele)) kel *= cscale_factor2;
			else if(packedEset[iele]->isFsiElement()) kel *= cscale_factor;
		}
		if(matrixTimers) matrixTimers->formTime += getTime();

		// Assemble element real and complex stiffness matrices in ops.K
		if(matrixTimers) matrixTimers->assemble -= getTime();
		if(isComplexF || (imag(kappa2) != 0)) {
			if(ops.K) ops.K->add(kcel,(*allDOFs)[iele]);
		}
		else {
			if(ops.K) ops.K->add(kel,(*allDOFs)[iele]);
			if(!isShifted && ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]);
			if(!isShifted && ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);

			if (isShifted && isSDamped) {
				int dim = kel.dim();
				ikel.setSize(dim);
				for(int i = 0; i < dim; ++i) for(int j = 0; j < dim; ++j)
						ikel[i][j] = -eta*kel[i][j];
				if(isStructureElement(iele)) if(ops.K) ops.K->addImaginary(ikel,(*allDOFs)[iele]);
			}
			if(packedEset[iele]->isConstraintElement()) {
				if(sinfo.isNonLin() && Mcoef == 1 && Kcoef == 0 && Ccoef == 0 && sinfo.newmarkBeta != 0) {
					//note: now I am using the tangent stiffness from kelArray so initial accelerations
					//      will be correctly computed even in the case of non-zero IDISP.
					//kel.~FullSquareMatrix();
					//kel = packedEset[iele]->stiffness(nodes, karray);
					if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
						mat->add(kel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
					}
					else if(mat) mat->add(kel,(*allDOFs)[iele]);
				}
				else {
					if(ops.Msolver) ops.Msolver->add(kel,(*allDOFs)[iele]);
				}
			}
		}
		if(matrixTimers) matrixTimers->assemble += getTime();

		int i, j, dim = kel.dim();
		// Form element real and complex mass matrices in mel and mcel
		if(matrixTimers) matrixTimers->formTime -= getTime();
		if(makeMass || isShifted) {
			if(isComplexF) {
				mcel = packedEset[iele]->massMatrix(nodes, mcarray);
			}
			else {
				if(melArray) mel.copy(melArray[iele]);
				else mel = packedEset[iele]->massMatrix(nodes, marray, mratio);
				this->densProjectStiffness(mel, iele);
				if(!sinfo.isNonLin()) this->transformMatrix(mel, iele);
			}
			if(sinfo.isCoupled) {
				if(isStructureElement(iele)) mel *= cscale_factor2;
				else if(packedEset[iele]->isFsiElement()) mel *= coupledScaling;
			}
		}

		if(isShifted) {
			if(isDamped && isStructureElement(iele)) {
				omega = sqrt(omega2);
				if(packedEset[iele]->hasDamping()) {
					izel = packedEset[iele]->dampingMatrix(nodes, izarray);
					this->densProjectStiffness(izel, iele);
					for(i = 0; i < dim; ++i)
						for(j = 0; j < dim; ++j)
							izel[i][j] = -omega*(izel[i][j]+beta*kel[i][j] + alpha*mel[i][j]);
				}
				else if(isDamped) {
					izel.setSize(mel.dim());
					for(i = 0; i < dim; ++i)
						for(j = 0; j < dim; ++j)
							izel[i][j] = -omega*(beta*kel[i][j] + alpha*mel[i][j]);
				}
			}
			if(isComplexF) {
				int dim = kcel.dim();
				for(i = 0; i < dim; ++i)
					for(j = 0; j < dim; ++j)
						kcel[i][j] -= kappa2*mcel[i][j];
			}
			else if(imag(kappa2) != 0.0) {
				kcel.setSize(kel.dim());
				for(i = 0; i < dim; ++i)
					for(j = 0; j < dim; ++j)
						kcel[i][j] = complex<double>(kel[i][j])-kappa2*mel[i][j];
			}
			else {
				double o2 = (isFluidElement(iele)) ? real(kappa2) : omega2;
				for(i = 0; i < dim; ++i)
					for(j = 0; j < dim; ++j) {
						if(isDamped && isStructureElement(iele))
							izel[i][j] = -omega*(beta*kel[i][j] + alpha*mel[i][j]);
						kel[i][j] -= o2*mel[i][j];
					}
			}
			if(sinfo.doFreqSweep && isFluidElement(iele)) {
				if(isComplexF) {
					mcel *= kappa2/omega2;
				}
				else if (imag(kappa2) != 0) {
					mcel.setSize(mel.dim());
					for(i = 0; i < dim; ++i)
						for(j = 0; j < dim; ++j)
							mcel[i][j] = kappa2/omega2*mel[i][j];
				}
				else {
					mel *= real(kappa2)/omega2;
				}
			}
		}
		else {
			if(makeMass && sinfo.isCoupled && isFluidElement(iele)) mel /= (1500.0 * 1500);
			if(sinfo.acoustic == 1 && prop) {
				double c = prop->ss; //speed of sound
				mel /= (c*c);
			}
		}
		if(matrixTimers) matrixTimers->formTime += getTime();

		// Assemble element real and complex mass matrices in ops.M, ops.Muc, ops.Mcc and ops.Msolver
		if(matrixTimers) matrixTimers->assemble -= getTime();
		if(isComplexF || (imag(kappa2) != 0)) {
			if(ops.M)   ops.M->add(mcel,(*allDOFs)[iele]);
			if(ops.Muc) ops.Muc->add(mcel,(*allDOFs)[iele]);
			if(ops.Mcc) ops.Mcc->add(mcel,(*allDOFs)[iele]);
			if(ops.Msolver) ops.Msolver->add(mcel,(*allDOFs)[iele]);
		}
		else {
			if(ops.M)   ops.M->add(mel,(*allDOFs)[iele]);
			if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[iele]);
			if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[iele]);
			if(ops.Msolver) ops.Msolver->add(mel,(*allDOFs)[iele]);
		}
		if(matrixTimers) matrixTimers->assemble += getTime();

		// Form the element impedance matrix in kel and the element damping matrix in cel
		if(matrixTimers) matrixTimers->formTime -= getTime();
		if(!isShifted) {
			if(makeMass) {
				if(celArray) cel.copy(celArray[iele]);
				else if(cel.dim() != dim) cel.setSize(dim);
				for(i = 0; i < dim; ++i) {
					for(j = 0; j < dim; ++j) {
						double m  = mel[i][j];
						double k  = kel[i][j];
						if(!celArray) {
							cel[i][j] = (zeroRot && (dofType[ (*allDOFs)[iele][i] ] == 1
							                         || dofType[ (*allDOFs)[iele][j] ] == 1)) ? 0 : alpha*m + beta*k;
						}
						kel[i][j] = Kcoef*k + Ccoef*cel[i][j] + Mcoef*m;
					}
				}
			}
			else {
				for(i = 0; i < dim; ++i)
					for(j = 0; j < dim; ++j)
						kel[i][j] *= Kcoef;
			}
		}
		if(matrixTimers) matrixTimers->formTime += getTime();

		// Assemble the element impedance matrix kel in mat, ops.Kuc and ops.spp
		if(matrixTimers) matrixTimers->assemble -= getTime();
		if(isComplexF || (imag(kappa2) != 0)) {
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(kcel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(kcel,(*allDOFs)[iele]);
			if(isShifted && ops.Kuc) ops.Kuc->add(kcel,(*allDOFs)[iele]);
			if(isShifted && ops.Kcc) ops.Kcc->add(kcel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(kcel,(*allDOFs)[iele]);
		}
		else {
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(kel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(kel,(*allDOFs)[iele]);
			if(isShifted && ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]); // note: Kuc is [K-omega2*M]_{uc} for IMPE (TODO check eigen)
			if(isShifted && ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(kel,(*allDOFs)[iele]);
			if(Kss) Kss->add(kel,(*allDOFs)[iele].data()); // for farfield output (TODO: check with Radek)
			if(isShifted && isDamped && isStructureElement(iele)) {
				if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
					mat->addImaginary(izel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
				}
				else if(mat) mat->addImaginary(izel,(*allDOFs)[iele]);
				if(ops.Kuc) ops.Kuc->addImaginary(izel,(*allDOFs)[iele]);
				if(ops.Kcc) ops.Kcc->addImaginary(izel,(*allDOFs)[iele]);
				if(ops.spp) ops.spp->addImaginary(izel,(*allDOFs)[iele]);
			}
			if(isShifted && isSDamped && isStructureElement(iele)) {
				if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
					mat->addImaginary(ikel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
				}
				else if(mat) mat->addImaginary(ikel,(*allDOFs)[iele]);
				if(ops.Kuc) ops.Kuc->addImaginary(ikel,(*allDOFs)[iele]);
				if(ops.Kcc) ops.Kcc->addImaginary(ikel,(*allDOFs)[iele]);
				if(ops.spp) ops.spp->addImaginary(ikel,(*allDOFs)[iele]);
			}

		}

		// Assemble the element damping matrix cel in ops.C and ops.C_deriv
		if(isShifted) {
			if(isDamped && isStructureElement(iele)) {
				izel /= omega;
				if(ops.C_deriv && ops.C_deriv[0]) ops.C_deriv[0]->addImaginary(izel,(*allDOFs)[iele]);
				if(ops.Cuc_deriv && ops.Cuc_deriv[0]) ops.Cuc_deriv[0]->addImaginary(izel,(*allDOFs)[iele]);
			}
		}
		else {
			if(ops.C) ops.C->add(cel,(*allDOFs)[iele]);
			if(ops.Cuc) ops.Cuc->add(cel,(*allDOFs)[iele]);
			if(ops.Ccc) ops.Ccc->add(cel,(*allDOFs)[iele]);
		}
		if(matrixTimers) matrixTimers->assemble += getTime();
	}

	// now deal with complex elements
	if(typeid(Scalar)==typeid(DComplex)) {
		// we assume they have no damping
		for(iele = 0; iele < numele; ++iele) {
			StructProp *prop = packedEset[iele]->getProperty();
			if(packedEset[iele]->isSommerElement()) continue;
			if(prop && (prop->E0!=0.0 || prop->mu0!=0.0)) continue; // rubber handled below
			if(!packedEset[iele]->isComplex()) continue;

			if(matrixTimers) matrixTimers->formTime -= getTime();
			kelC = packedEset[iele]->complexStiffness(nodes, karrayC);
			this->densProjectStiffnessC(kelC, iele);
			if(sinfo.isCoupled) {  // coupled scaling
				if(isStructureElement(iele)) kelC *= cscale_factor2;
				else if(packedEset[iele]->isFsiElement()) kelC *= cscale_factor;
			}

			int i, j, dim = kel.dim();
			if(makeMass || isShifted) {
				melC = packedEset[iele]->complexMassMatrix(nodes, marrayC, mratio);
				this->densProjectStiffnessC(melC, iele);
				if(sinfo.isCoupled) {
					if(isStructureElement(iele)) melC *= cscale_factor2;
					else if(packedEset[iele]->isFsiElement()) melC *= coupledScaling;
				}
			}
			if(isShifted) {
				omega2 = (isFluidElement(iele)) ? packedEset[iele]->helmCoef() : geoSource->shiftVal();
				if(isStructureElement(iele)) {
					omega = sqrt(omega2);
				}
				for(i = 0; i < dim; ++i)
					for(j = 0; j < dim; ++j)
						kelC[i][j] -= omega2*melC[i][j];
				if(sinfo.doFreqSweep && isFluidElement(iele))
					melC *= packedEset[iele]->helmCoef()/geoSource->shiftVal();
			}
			else {
				if(makeMass && sinfo.isCoupled && isFluidElement(iele)) melC /= (1500.0 * 1500);
				if(sinfo.acoustic == 1 && prop) {
					double c = prop->ss;//speed of sound
					melC /= (c*c);
				}
			}
			if(matrixTimers) matrixTimers->formTime += getTime();

			if(matrixTimers) matrixTimers->assemble -= getTime();
			if(ops.K)   ops.K->add(kelC,(*allDOFs)[iele]);
			if(ops.Kuc) ops.Kuc->add(kelC,(*allDOFs)[iele]);
			if(ops.Kcc) ops.Kcc->add(kelC,(*allDOFs)[iele]);
			if(ops.M)   ops.M->add(melC,(*allDOFs)[iele]);
			if(ops.Muc) ops.Muc->add(melC,(*allDOFs)[iele]);
			if(ops.Mcc) ops.Mcc->add(melC,(*allDOFs)[iele]);
			if(ops.Msolver) ops.Msolver->add(melC,(*allDOFs)[iele]);
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(kelC,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(kelC,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(kelC,(*allDOFs)[iele]);
			if(matrixTimers) matrixTimers->assemble += getTime();
		}
	}

	// Acoustic rubber elements
	int N = sinfo.doFreqSweep ? solInfo().getSweepParams()->nFreqSweepRHS - 1:0; // number of derivatives
	complex<double> *krarray = new complex<double>[(N+3)*maxNumDOFs*maxNumDOFs];
	for(iele = 0; iele < numele; ++iele) {
		StructProp *prop = packedEset[iele]->getProperty();

		if(!prop || (prop && (prop->E0==0.0 || prop->mu0==0.0)) || !isShifted) continue;

		omega2 = geoSource->shiftVal();
		SPropContainer& sProps = geoSource->getStructProps();
		SPropContainer::iterator it = sProps.begin();
		bool found = false;
		int i_arubber = 0;
		while(it != sProps.end()) {
			StructProp* p = &(it->second);
			if((p->E0!=0.0 || p->mu0!=0.0)) {
				if (prop==p) {
					found = true;
					break;
				}
				i_arubber++;
			}
			it++;
		}
		if (!found) fprintf(stderr,"Error: rubber material not found.\n");

		if(matrixTimers) matrixTimers->formTime -= getTime();
		packedEset[iele]->aRubberStiffnessDerivs(nodes, krarray,N,sqrt(omega2));
		mel = packedEset[iele]->massMatrix(nodes, marray,1);
		int ndof = packedEset[iele]->numDofs();
		if(sinfo.isCoupled) if(isStructureElement(iele)) {  // coupled scaling
				for(int i=0;i<ndof*ndof;i++)  marray[i] *= cscale_factor2;
				for(int i=0;i<ndof*ndof*(N+3);i++)  krarray[i] *= cscale_factor2;
			}

		if(matrixTimers) matrixTimers->formTime += getTime();
		if(matrixTimers) matrixTimers->assemble -= getTime();
		if(ops.M)   ops.M->add(mel,(*allDOFs)[iele]);
		if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[iele]);
		if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[iele]);
		if(ops.Msolver) ops.Msolver->add(mel,(*allDOFs)[iele]);
		if (sinfo.doFreqSweep) for(int i=0;i<=N;i++) {
				FullSquareMatrixC kdel(ndof,krarray+ndof*ndof*(i+2));
				if (ops.K_deriv) if (ops.K_deriv[i]) ops.K_deriv[i]->add(kdel,(*allDOFs)[iele]);
				if (ops.Kuc_deriv) if (ops.Kuc_deriv[i]) ops.Kuc_deriv[i]->add(kdel,(*allDOFs)[iele]);
			}

		FullSquareMatrixC krel_m(ndof,krarray+ndof*ndof*0);
		FullSquareMatrixC krel_l(ndof,krarray+ndof*ndof*1);


		if (found & sinfo.doFreqSweep) {
			if(ops.K_arubber_m[i_arubber])
				ops.K_arubber_m[i_arubber]->add(krel_m,(*allDOFs)[iele]);
			if(ops.K_arubber_l[i_arubber])
				ops.K_arubber_l[i_arubber]->add(krel_l,(*allDOFs)[iele]);
			if(ops.Kuc_arubber_m) if(ops.Kuc_arubber_m[i_arubber])
					ops.Kuc_arubber_m[i_arubber]->add(krel_m,(*allDOFs)[iele]);
			if(ops.Kuc_arubber_l) if(ops.Kuc_arubber_l[i_arubber])
					ops.Kuc_arubber_l[i_arubber]->add(krel_l,(*allDOFs)[iele]);
		}

		FullSquareMatrixC krel(ndof,krarray+ndof*ndof*2);
		if(ops.K)   ops.K->add(krel,(*allDOFs)[iele]);
		for(int i=0;i<ndof*ndof;i++)  krarray[i+ndof*ndof*2] -= omega2*marray[i];
		if(ops.Kuc) ops.Kuc->add(krel,(*allDOFs)[iele]);
		if(ops.Kcc) ops.Kcc->add(krel,(*allDOFs)[iele]);

		if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
			mat->add(krel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
		}
		else if(mat) mat->add(krel,(*allDOFs)[iele]);
		if(ops.spp) ops.spp->add(krel,(*allDOFs)[iele]);
		if(matrixTimers) matrixTimers->assemble += getTime();
	}
	delete[] krarray;

	if(sinfo.ATDARBFlag >= 1.0) {
		getCurvatures3Daccurate(this);
	}
	int passage = 0;
	bool Complx = false;//Test for stabilising absorbing condition of high order - JF

	// LOOP over the sommer elements
	if((sinfo.ATDARBFlag == 1.5) && Complx) { // kel goes for DComplex, mat and ops.K(s) too
		for(iele = 0; iele < numele; ++iele) {
			if(!packedEset[iele]->isSommerElement() || sinfo.ATDARBFlag == -2.0) continue;
			//add massMatrix
			mel.zero();
			mel = packedEset[iele]->massMatrix(nodes, marray,mratio);
			if (passage==0) {
				int si = mel.dim();
				cel.setSize(si);
				passage=1;
			}
			if(ops.M)   ops.M->add(mel,(*allDOFs)[iele]);
			if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[iele]);
			if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[iele]);
			if(ops.Msolver) ops.Msolver->add(mel,(*allDOFs)[iele]);
			mel *= Mcoef;
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(mel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(mel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(mel,(*allDOFs)[iele]);

			//add dampingMatrix
			cel.zero();
			cel = packedEset[iele]->dampingMatrix(nodes, marray);
			if(ops.C) ops.C->add(cel,(*allDOFs)[iele]);
			if(ops.Cuc) ops.Cuc->add(cel,(*allDOFs)[iele]);
			cel *= Ccoef;
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(cel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(cel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(cel,(*allDOFs)[iele]);

			//add stiffness
			kel.zero();
			kel = packedEset[iele]->stiffness(nodes, karray);
			if(ops.K) ops.K->add(kel,(*allDOFs)[iele]);
			if(ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]);
			if(ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);
			kel *= Kcoef;
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(kel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(kel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(kel,(*allDOFs)[iele]);
			//add Imaginary part of stiffness
			kel.zero();
			kel = packedEset[iele]->imStiffness(nodes, karray);
			if(ops.K) ops.K->addImaginary(kel,(*allDOFs)[iele]);
			if(ops.Kuc) ops.Kuc->addImaginary(kel,(*allDOFs)[iele]);
			if(ops.Kcc) ops.Kcc->addImaginary(kel,(*allDOFs)[iele]);
			kel *= Kcoef;
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->addImaginary(kel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->addImaginary(kel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->addImaginary(kel,(*allDOFs)[iele]);
		}
	}
	else {
		for(iele = 0; iele < numele; ++iele) {
			if(!packedEset[iele]->isSommerElement() || sinfo.ATDARBFlag == -2.0) continue;

			// add massMatrix
			mel.zero();
			mel = packedEset[iele]->massMatrix(nodes, marray, mratio);
			if (passage==0) {
				int si = mel.dim();
				cel.setSize(si);
				passage=1;
			}
			if(ops.M)   ops.M->add(mel,(*allDOFs)[iele]);
			if(ops.Muc) ops.Muc->add(mel,(*allDOFs)[iele]);
			if(ops.Mcc) ops.Mcc->add(mel,(*allDOFs)[iele]);
			if(ops.Msolver) ops.Msolver->add(mel,(*allDOFs)[iele]);
			mel *= Mcoef;
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(mel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(mel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(mel,(*allDOFs)[iele]);

			// add dampingMatrix
			cel.zero();
			cel = packedEset[iele]->dampingMatrix(nodes, marray);
			if(ops.C) ops.C->add(cel,(*allDOFs)[iele]);
			if(ops.Cuc) ops.Cuc->add(cel,(*allDOFs)[iele]);
			cel *= Ccoef;
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(cel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(cel,(*allDOFs)[iele]);
			if(ops.spp) ops.spp->add(cel,(*allDOFs)[iele]);

			// add stiffness
			kel.zero();
			kel = packedEset[iele]->stiffness(nodes, karray);
			if(ops.K) ops.K->add(kel,(*allDOFs)[iele]);
			if(ops.Kuc) ops.Kuc->add(kel,(*allDOFs)[iele]);
			if(ops.Kcc) ops.Kcc->add(kel,(*allDOFs)[iele]);
			kel *= Kcoef;
			if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
				mat->add(kel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			}
			else if(mat) mat->add(kel,(*domain->getAllDOFs())[subCast->getGlElems()[iele]]);
			if(ops.spp) ops.spp->add(kel,(*allDOFs)[iele]);
		}
	}

	if(makeMass || isShifted) {
		// Add discrete mass contribution to the Mass Matrix
		// Three matrices need to be changed. Mass matrix itself,
		// Damping Matrix and K tilda.
		Scalar m;
		DMassData *current = firstDiMass;
		double omega2 = geoSource->shiftVal();
		double omega = sqrt(geoSource->shiftVal());
		while(current != 0) {
			// PJSA: modified all addDiscreteMass functions to accept dsa dof rather than cdsa dof (much safer)
			// TODO what about Kuc, Cuc?
			int dof = dsa->locate(current->node, (1 << current->dof));
			if(dof == -1)
			{ current = current->next; continue; }
			if(current->jdof > -1) { // PJSA 10-9-06 for off-diagonal mass terms eg. products of inertia I21, I31, I32
				int jdof = dsa->locate(current->node, (1 << current->jdof));
				if(jdof == -1)
				{ current = current->next; continue; }
				if(isShifted) {
					// note: for coupled problems we are assuming lumped mass is structural not fluid
					double mass = sinfo.isCoupled ? current->diMass*cscale_factor2 : current->diMass;
					double m_real = -omega2*mass;
					double m_imag = isDamped ? -omega*alphaDamp*mass : 0.0;
					ScalarTypes::initScalar(m, m_real, m_imag);
					if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
						mat->add(domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->dof)),
						         domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->jdof)), m);
					}
					else if(mat) mat->add(dof, jdof, m);
					if(ops.spp) ops.spp->add(dof, jdof, m);
					if(isDamped) {
						ScalarTypes::initScalar(m, 0.0, alphaDamp*mass);
						if(ops.C_deriv && ops.C_deriv[0]) ops.C_deriv[0]->add(dof, jdof, m);
					}
					if(ops.M) ops.M->add(dof, jdof, mass);
					if(ops.Msolver) ops.Msolver->add(dof, jdof, mass);
				}
				else {
					if(ops.M) ops.M->add(dof,jdof,current->diMass);
					if(ops.Msolver) ops.Msolver->add(dof,jdof,current->diMass);
					if(ops.C) ops.C->add(dof,jdof,alphaDamp*current->diMass);
					double mass = Mcoef*current->diMass;
					if(isDamped) mass += Ccoef*alphaDamp*current->diMass;
					if(mass != 0) {
						if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
							mat->add(domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->dof)),
							         domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->jdof)), mass);
						}
						else if(mat) mat->add(dof,jdof,mass);
						if(ops.spp) ops.spp->add(dof,jdof,mass);
					}
				}
			}
			else {
				if(isShifted) {  // add discrete mass contributions to global matrices for frequency response analysis
					double mass = sinfo.isCoupled ? current->diMass*cscale_factor2 : current->diMass;
					double m_real = -omega2*mass;
					double m_imag = isDamped ? -omega*alphaDamp*mass : 0.0;
					ScalarTypes::initScalar(m, m_real, m_imag);
					if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
						mat->addDiscreteMass(domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->dof)), m);
					}
					else if(mat) mat->addDiscreteMass(dof, m);
					if(ops.spp) ops.spp->addDiscreteMass(dof, m);
					if(isDamped) {
						ScalarTypes::initScalar(m, 0.0, alphaDamp*mass);
						if(ops.C_deriv && ops.C_deriv[0]) ops.C_deriv[0]->addDiscreteMass(dof, m);
					}
					if(ops.M) ops.M->addDiscreteMass(dof, mass);
					if(ops.Msolver) ops.Msolver->addDiscreteMass(dof, mass);
				}
				else {
					if(ops.M) ops.M->addDiscreteMass(dof, current->diMass);
					if(ops.Msolver) ops.Msolver->addDiscreteMass(dof, current->diMass);
					if(ops.C) ops.C->addDiscreteMass(dof, alphaDamp*current->diMass);
					double mass = Mcoef*current->diMass;
					if(isDamped) mass += Ccoef*alphaDamp*current->diMass;
					if(mass != 0) {
						if(mdds_flag) {
#if defined(_OPENMP)
#pragma omp critical
#endif
							mat->addDiscreteMass(domain->getDSA()->locate(subCast->getGlNodes()[current->node], (1 << current->dof)), mass);
						}
						else if(mat) mat->addDiscreteMass(dof, mass);
						if(ops.spp) ops.spp->addDiscreteMass(dof, mass);
					}
				}
			}
			current = current->next;
		}
	}

	if((sinfo.HEV) && (probType() == SolverInfo::Modal)) {
		Mff = new GenBLKSparseMatrix<double>(nodeToNodeFluid, dsaFluid, c_dsaFluid, sinfo.solvercntl->trbm, *sinfo.solvercntl, 0);
		Mff->zeroAll();

		for(iele = 0; iele < geoSource->numElemFluid(); ++iele) {
			//element matrix already premultiplied by fluid density
			kel = (*(geoSource->getPackedEsetFluid()))[iele]->stiffness(nodes,karray);
			Mff->add(kel,(*allDOFsFluid)[iele].data());
		}
		fprintf(stderr," ... Factoring Fluid Mass Matrix ...\n");
		//double trbmtemp = sinfo.solvercntl->trbm;
		//sinfo.setTrbm(1e-16);
		Mff->factor();
		//sinfo.setTrbm(trbmtemp);
		fprintf(stderr," ... Fluid Mass Matrix Factored...\n");

		if(sinfo.addedMass == 1) {

			//get fluid density from first element since HOMOGENEOUS fluid is assumed
			double rhoFluid = ( *(geoSource->getPackedEsetFluid()) )[0]->getProperty()->rho;

			double** C_NZrows = getCMatrix();
			double** FCtinv_NZcol = new double* [domain->nuNonZero];

			int nnp = domain->numUnconFluid();

			for (int i=0; i < domain->nuNonZero; ++i)  {
				FCtinv_NZcol[i] = new double[nnp];
				for (int jp=0; jp < nnp; ++jp)  {
					C_NZrows[i][jp] = rhoFluid * C_NZrows[i][jp];
					FCtinv_NZcol[i][jp] = C_NZrows[i][jp];
				}
			}

			//obtain non zero parts of F^-1C^T
			Mff->reSolve(domain->nuNonZero,FCtinv_NZcol);

			//initialize the added mass operator (Ma_NZ only contains the non-zero parts)
			FullSquareMatrix Ma_NZ(domain->nuNonZero,0);
			Ma_NZ.zero();

			//obtain the added mass operator
			for (int i=0; i<domain->nuNonZero; ++i)  {
				for (int j=i; j<domain->nuNonZero; ++j)  {
					for (int k=0; k<domain->npNonZero[i]; ++k)  {
						int K = domain->pmap[i][k];
						Ma_NZ[i][j] += C_NZrows[i][K]*FCtinv_NZcol[j][K];
					}
				}
			}
			for (int i=0; i<domain->nuNonZero; ++i)  {
				for (int j=i; j<domain->nuNonZero; ++j)  {
					Ma_NZ[j][i] = Ma_NZ[i][j];
				}
			}

			fprintf(stderr," ... HEV Problem: Added mass obtained ...\n");

			//add added mass to structural mass for EVP!
			ops.M->add(Ma_NZ,domain->umap_add);
			fprintf(stderr," ... HEV Problem: Added mass contribution included in structural mass ...\n");

			if(isShifted) {
				omega2 = geoSource->shiftVal();
				Ma_NZ *= (-omega2*Kcoef);
				if (ops.K) ops.K->add(Ma_NZ,domain->umap_add);
				if (mat)   mat->add(Ma_NZ,domain->umap_add);
				if (ops.spp)   ops.spp->add(Ma_NZ,domain->umap_add);
				fprintf(stderr," ... HEV Problem: Added mass contribution included in shifted opertor ...\n");
			}
		}

	}

	delete [] karray;
	delete [] marray;
	delete [] ikarray;
	delete [] kcarray;
	delete [] mcarray;
	delete [] karrayC;
	delete [] marrayC;

	// TODO: assembleSommer and assembleATDROB need to be modified for mdds_flag == true
	if(sinfo.isAcousticHelm()) assembleSommer<Scalar>(mat, &ops);

	if(sinfo.ATDROBalpha != 0.0) assembleATDROB<Scalar>(mat, &ops,Kcoef);

	if(matrixTimers) matrixTimers->memoryForm += memoryUsed();
}

template<class Scalar>
GenDBSparseMatrix<Scalar> *
Domain::constructDBSparseMatrix(DofSetArray *dof_set_array, Connectivity *cn)
{
	if(dof_set_array == 0) dof_set_array = c_dsa;
	if(cn == 0)
		return new GenDBSparseMatrix<Scalar>(nodeToNode.get(), dsa, c_dsa);
	else
		return new GenDBSparseMatrix<Scalar>(cn, dsa, c_dsa);
}

#ifdef USE_EIGEN3
template<typename Scalar>
GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *
Domain::constructEiSparse(DofSetArray *c_dsa, Connectivity *nodeToNode, bool flag)
{
	if(c_dsa == nullptr)
		c_dsa = Domain::c_dsa;
	if(nodeToNode == nullptr)
		nodeToNode = Domain::nodeToNode.get();
	return new GenEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(nodeToNode, dsa, c_dsa, flag);
}
#endif

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass> *
Domain::constructEiSparseMatrix(DofSetArray *c_dsa, Connectivity *nodeToNode, bool flag)
{
#ifdef USE_EIGEN3
	if(c_dsa == nullptr)
		c_dsa = Domain::c_dsa;
	if(nodeToNode == nullptr)
		nodeToNode = Domain::nodeToNode.get();
	return new GenEiSparseMatrix<Scalar,SolverClass>(nodeToNode, dsa, c_dsa, flag);
#else
	std::cerr << "USE_EIGEN3 is not defined\n";
#endif
}

template<class Scalar>
GenNBSparseMatrix<Scalar> *
Domain::constructNBSparseMatrix()
{
	return new GenNBSparseMatrix<Scalar>(nodeToNode.get(), c_dsa);
}

template<class Scalar>
GenCuCSparse<Scalar> *
Domain::constructCuCSparse(DofSetArray *dof_set_array)
{
	if(dof_set_array == 0) dof_set_array = c_dsa;
	GenCuCSparse<Scalar> *sp = 0;
	if((dsa->size() - dof_set_array->size()) == 0)
		sp = 0;
	else
		sp = new GenCuCSparse<Scalar>(nodeToNode.get(), dsa, dof_set_array);
	return sp;
}

template<class Scalar>
GenCuCSparse<Scalar> *
Domain::constructCCSparse(DofSetArray *dof_set_array)
{
	if(dof_set_array == 0)
		dof_set_array = c_dsa;
	return new GenCuCSparse<Scalar>(nodeToNode.get(), dsa, dof_set_array->getConstrndNum(),
	                                dof_set_array->getConstrndNum());
}

template<class Scalar>
GenBLKSparseMatrix<Scalar> *
Domain::constructBLKSparseMatrix(DofSetArray *DSA, Rbm *rbm)
{
	if(DSA == 0) DSA = c_dsa;
	// use the correct zero pivot tolerance for GRBM or TRBM.
	double trbm = (rbm && geoSource->shiftVal() == 0) ? sinfo.solvercntl->trbm2 : sinfo.solvercntl->trbm;
	if(sinfo.newmarkBeta == 0.0 && (sinfo.acoustic || sinfo.inertiaLumping == 2)) {
		if(sinfo.ATDARBFlag >- 1.0 && sinfo.acoustic) {
			auto nToN = nodeToNode_sommer->modify();
			return new GenBLKSparseMatrix<Scalar>(&nToN, dsa, DSA, sinfo.solvercntl->trbm, *sinfo.solvercntl, rbm);
		}
		else {
			int numN=numNodes();
			// import a diagonal connectivity for the mass matrix
			Connectivity connForMas(numN);
			connForMas = connForMas.modify();
			GenBLKSparseMatrix<Scalar> * retur =
				new GenBLKSparseMatrix<Scalar>(&connForMas, dsa, DSA, sinfo.solvercntl->trbm, *sinfo.solvercntl, rbm);
			return retur;
		}
	}
	else {
		if(!sinfo.getDirectMPC())
			return new GenBLKSparseMatrix<Scalar>(nodeToNode.get(), dsa, DSA, sinfo.solvercntl->trbm, *sinfo.solvercntl, rbm);
		else {
			if(nodeToNodeDirect) delete nodeToNodeDirect;
			nodeToNodeDirect = prepDirectMPC();
			DOFMap *baseMap = new DOFMap[dsa->size()];
			DOFMap *eqMap = new DOFMap[DSA->size()];
			// TODO Examine when DSA can be different from c_dsa
			if(MpcDSA && sinfo.isNonLin()) delete MpcDSA;
			MpcDSA = makeMaps(dsa, c_dsa, baseMap, eqMap);
			typename WrapSparseMat<Scalar>::CtorData
				baseArg(nodeToNodeDirect, dsa, MpcDSA, sinfo.solvercntl->trbm, *sinfo.solvercntl, rbm);
			int nMappedEq = DSA->size();
			return new MappedAssembledSolver<WrapSparseMat<Scalar>, Scalar>(baseArg, dsa->size(), baseMap, nMappedEq, eqMap, c_dsa, rbm);
		}
	}
}

template<class Scalar>
void
Domain::buildPreSensitivities(AllSensitivities<Scalar> &allSens, Scalar *bcx)
{
	switch(sinfo.solvercntl->type) {
		default:
			fprintf(stderr," *** WARNING: Solver not Specified  ***\n");
		case SolverSelection::Direct:
			makePreSensitivities(allSens, bcx);
			break;
	}
}

template<class Scalar>
void
Domain::buildPostSensitivities(GenSolver<Scalar> *sysSolver,
                               GenSparseMatrix<Scalar> *K, GenSparseMatrix<Scalar> *spm,
                               AllSensitivities<Scalar> &allSens, GenVector<Scalar> *sol, Scalar *bcx, bool isDynam,
                               GeomState *refState, GeomState *geomState, Corotator **allCorot, bool isNonLin)
{
	switch(sinfo.solvercntl->type) {
		default:
			fprintf(stderr," *** WARNING: Solver not Specified  ***\n");
		case SolverSelection::Direct:
			makePostSensitivities(sysSolver, spm, allSens, sol, bcx, K, isDynam, refState, geomState, allCorot, isNonLin);
			break;
	}
}

template<class Scalar>
void
Domain::buildNLPostSensitivities(GenSolver<Scalar> *sysSolver, AllSensitivities<Scalar> &allSens,
                                 GeomState *refState, GeomState *geomState, Corotator **allCorot, bool isDynam)
{
	switch(sinfo.solvercntl->type) {
		default:
			fprintf(stderr," *** WARNING: Solver not Specified  ***\n");
		case SolverSelection::Direct:
			makeNLPostSensitivities(sysSolver, allSens, refState, geomState, allCorot, isDynam);
			break;
	}
}

template<class Scalar>
void
Domain::buildOps(AllOps<Scalar> &allOps, double Kcoef, double Mcoef, double Ccoef,
                 Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray,
                 FullSquareMatrix *celArray, bool factorize)
{
	if(matrixTimers) matrixTimers->memorySolve -= memoryUsed();

	if(allOps.sysSolver) delete allOps.sysSolver;
	GenSolver<Scalar> *systemSolver = 0;
	if(geoSource->shiftVal() != 0 || Mcoef != 0 || Ccoef != 0) rbm = 0; // PJSA: don't pass
	// geometric rbms to
	// solver in this case
	allOps.spm = 0;
	SfemBlockMatrix<Scalar> *sfbm = 0;
	int L = 1;
	int ndim = 0;
	if(sinfo.inpc || sinfo.noninpc) {
		L = sfem->getL();
		int n = domain->numUncon();
		int P = sfem->getP();
		ndim = sfem->getndim();
		int output_order = sfem->getoutput_order();
		if(sinfo.inpc) sfbm = new SfemBlockMatrix<Scalar>(L,n,P,ndim,output_order);
	}

	for(int i=0; i<ndim+1; ++i) {
		if(sinfo.inpc) domain->setNewProperties(i);
		if(sinfo.noninpc && i>0) break;
		switch(sinfo.solvercntl->type) {
			default:
				fprintf(stderr," *** WARNING: Solver not Specified %d  ***\n", static_cast<int>(sinfo.solvercntl->type));
			case SolverSelection::Direct:
				makeStaticOpsAndSolver<Scalar>(allOps, Kcoef, Mcoef, Ccoef,
				                               systemSolver, allOps.spm, rbm, kelArray, melArray, celArray); // also used for eigen
				break;
			case SolverSelection::Iterative:
				makeDynamicOpsAndSolver<Scalar>(allOps, Kcoef, Mcoef, Ccoef,
				                                systemSolver, allOps.spm, rbm, kelArray, melArray, celArray);
				break;
		}
		if(sinfo.inpc) {
			sfbm->setKi(allOps.spm,i);
			if(i==0 && sinfo.solvercntl->precond==3) {
				GenBLKSparseMatrix<Scalar> *prec_solver = constructBLKSparseMatrix<Scalar>(c_dsa, rbm);
				prec_solver->zeroAll();
				AllOps<Scalar> allOps_tmp;
				makeSparseOps<Scalar>(allOps_tmp,Kcoef,Mcoef,Ccoef,prec_solver,kelArray,melArray,celArray);
				prec_solver->factor();
				sfbm->setMeanSolver(prec_solver);
			}
			if(i==0) {
				allOps.rhs_inpc = new GenVector<Scalar>(domain->numUncon());
				domain->template buildRHSForce<Scalar>(*allOps.rhs_inpc, allOps.Kuc);
			}
			if(allOps.Kuc) allOps.Kuc->zeroAll();
		}
	}

	// Set allOps pointer to system solver
	if(!sinfo.inpc) allOps.sysSolver = systemSolver;
	else {
		switch(sinfo.solvercntl->iterType) {
			default:
			case 0 : {
				filePrint(stderr," ... CG Solver is Selected          ...\n");
				allOps.sysSolver = new GenPCGSolver<Scalar, GenVector<Scalar>, SfemBlockMatrix<Scalar> >(sfbm, sinfo.solvercntl->precond, sinfo.solvercntl->maxit,
				                                                                                         sinfo.solvercntl->tol, sinfo.solvercntl->maxvecsize,
				                                                                                         sinfo.solvercntl->verbose);
				break;
			}
			case 4: {
				filePrint(stderr," ... Bi-CG Solver is Selected       ...\n");
				allOps.sysSolver = new GenBCGSolver<Scalar, GenVector<Scalar>, SfemBlockMatrix<Scalar>, GenSolver<Scalar> >(sinfo.solvercntl->maxit, sinfo.solvercntl->tol, sfbm);
				break;
			}
			case 5: {
				filePrint(stderr," ... CR Solver is Selected          ...\n");
				allOps.sysSolver = new GenCRSolver<Scalar, GenVector<Scalar>, SfemBlockMatrix<Scalar>, GenSolver<Scalar> >(sinfo.solvercntl->maxit, sinfo.solvercntl->tol, sfbm);
				break;
			}
		}
	}

	/*if(sinfo.solvercntl->printMatLab && !dynamic_cast<GenSubDomain<Scalar>*>(this)) {
	  allOps.spm->printSparse(sinfo.solvercntl->printMatLabFile);
	}*/

	if(factorize)
	{
		// Time system matrix factorization
		if(matrixTimers) matrixTimers->factor -= getTime();

		if(systemSolver) {
			if(!verboseFlag && sinfo.isNonLin() && sinfo.isDynam()) systemSolver->setPrintNullity(false);
			systemSolver->factor();
		}

		if(matrixTimers) matrixTimers->factor += getTime();
	}

	if(matrixTimers) matrixTimers->memorySolve += memoryUsed();
}


template<class Scalar>
void
Domain::rebuildOps(AllOps<Scalar> &allOps, double Kcoef, double Mcoef, double Ccoef,
                   Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray,
                   FullSquareMatrix *celArray, bool factorize)
{
	allOps.spm->zeroAll();
	makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,allOps.spm,kelArray,melArray,celArray);

	if(factorize) allOps.sysSolver->factor();
}

template<class Scalar>
void
Domain::getSolverAndKuc(AllOps<Scalar> &allOps, FullSquareMatrix *kelArray, Rbm *rbm, bool factorize)
{
	// PJSA 10-5-04: new version, this one can pass all sys matricies back via allOps
	// not just K and Kuc. This is required for frequency sweep analysis

	// ... Call necessary Operator's constructors
	allOps.Kuc = constructCuCSparse<Scalar>();
	allOps.Kcc = constructCCSparse<Scalar>();
#ifdef USE_EIGEN3
	if(sinfo.sensitivity)
		allOps.K =
			constructEiSparseMatrix<Scalar, Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(
				c_dsa, nodeToNode.get(), false);
#endif

	// for freqency sweep: need M, Muc, C, Cuc
	bool isDamped = (sinfo.alphaDamp != 0.0) || (sinfo.betaDamp != 0.0) || packedEset.hasDamping();
	if(sinfo.doFreqSweep && (sinfo.getSweepParams()->nFreqSweepRHS > 1 || isDamped || sinfo.getSweepParams()->isAdaptSweep)) {
		//---- UH ----
		if(sinfo.getSweepParams()->freqSweepMethod == SweepParams::PadeLanczos ||
		   sinfo.getSweepParams()->freqSweepMethod == SweepParams::GalProjection ||
		   sinfo.getSweepParams()->freqSweepMethod == SweepParams::KrylovGalProjection ||
		   sinfo.getSweepParams()->freqSweepMethod == SweepParams::WCAWEGalProjection ||
		   sinfo.getSweepParams()->isAdaptSweep) {
			if (allOps.K)
				delete allOps.K;
			if (allOps.M) delete allOps.M;
			if (allOps.Muc) delete allOps.Muc;
			if (allOps.C_deriv) if (allOps.C_deriv[0]) delete allOps.C_deriv[0];

			if(sinfo.isCoupled) allOps.K = constructNBSparseMatrix<Scalar>(); // unsymmetric
			else allOps.K = constructDBSparseMatrix<Scalar>(); // symmetric
		}
		//---- UH ----
		if(sinfo.isCoupled) allOps.M = constructNBSparseMatrix<Scalar>();  // unsymmetric
		else allOps.M = constructDBSparseMatrix<Scalar>();  // symmetric
		allOps.Muc = constructCuCSparse<Scalar>();
		if(isDamped || (numSommer > 0)) {
			int numC_deriv;
			if((numSommer > 0) && ((sommerfeldType == 2) || (sommerfeldType == 4)))
				numC_deriv = sinfo.getSweepParams()->nFreqSweepRHS - 1;
			else
				numC_deriv = 1;
			allOps.C_deriv = new GenSparseMatrix<Scalar> * [sinfo.getSweepParams()->nFreqSweepRHS - 1];
			for(int n=0; n<numC_deriv; ++n)
				allOps.C_deriv[n] = constructDBSparseMatrix<Scalar>();
			for(int n=numC_deriv; n<sinfo.getSweepParams()->nFreqSweepRHS - 1; ++n)
				allOps.C_deriv[n] = 0;
			if(c_dsa->size() > 0 && (c_dsa->size() - dsa->size()) != 0) {
				allOps.Cuc_deriv = new GenSparseMatrix<Scalar> * [sinfo.getSweepParams()->nFreqSweepRHS - 1];
				for(int n=0; n<numC_deriv; ++n)
					allOps.Cuc_deriv[n] = constructCuCSparse<Scalar>();
				for(int n=numC_deriv; n<sinfo.getSweepParams()->nFreqSweepRHS - 1; ++n)
					allOps.Cuc_deriv[n] = 0;
			}

		}
	}
		// for rbmfilter: need M
	else if(sinfo.filterQ == 0 && (sinfo.filterFlags || sinfo.hzemFilterFlag || sinfo.slzemFilterFlag)) {
		allOps.M = constructDBSparseMatrix<Scalar>();
	}

	if(isDamped) {
		if (allOps.K_deriv) for(int i=0;i<=allOps.n_Kderiv;i++) if (allOps.K_deriv[i]) delete allOps.K_deriv[i];
		if (allOps.Kuc_deriv) for(int i=0;i<=allOps.n_Kderiv;i++) if (allOps.Kuc_deriv[i]) delete allOps.Kuc_deriv[i];
		int numK_deriv = sinfo.doFreqSweep ? sinfo.getSweepParams()->nFreqSweepRHS-1:0;
		allOps.n_Kderiv = numK_deriv;
		allOps.K_deriv = new GenSparseMatrix<Scalar> * [numK_deriv+1];
		for(int n=0; n<=numK_deriv; ++n)
			allOps.K_deriv[n] = constructDBSparseMatrix<Scalar>();
		if(c_dsa->size() > 0 && (c_dsa->size() - dsa->size()) != 0) {
			allOps.Kuc_deriv = new GenSparseMatrix<Scalar> * [numK_deriv+1];
			for(int n=0; n<=numK_deriv; ++n)
				allOps.Kuc_deriv[n] = constructCuCSparse<Scalar>();
		}

		if (allOps.K_arubber_l) for(int i=0;i<allOps.num_K_arubber;i++) {
				if (allOps.K_arubber_l[i]) delete allOps.K_arubber_l[i];
				if (allOps.K_arubber_m[i]) delete allOps.K_arubber_m[i];
			}
		if (allOps.Kuc_arubber_l) for(int i=0;i<allOps.num_K_arubber;i++) {
				if (allOps.Kuc_arubber_l[i]) delete allOps.Kuc_arubber_l[i];
				if (allOps.Kuc_arubber_m[i]) delete allOps.Kuc_arubber_m[i];
			}
		int num_K_arubber = geoSource->num_arubber;
		allOps.num_K_arubber= num_K_arubber;
		allOps.K_arubber_l = new GenSparseMatrix<Scalar> * [num_K_arubber];
		allOps.K_arubber_m = new GenSparseMatrix<Scalar> * [num_K_arubber];
		for(int n=0; n<num_K_arubber; ++n) {
			allOps.K_arubber_l[n] = constructDBSparseMatrix<Scalar>();
			allOps.K_arubber_m[n] = constructDBSparseMatrix<Scalar>();
		}
		if(c_dsa->size() > 0 && (c_dsa->size() - dsa->size()) != 0) {
			allOps.Kuc_arubber_l = new GenSparseMatrix<Scalar> * [num_K_arubber];
			allOps.Kuc_arubber_m = new GenSparseMatrix<Scalar> * [num_K_arubber];
			for(int n=0; n<num_K_arubber; ++n) {
				allOps.Kuc_arubber_l[n] = constructCuCSparse<Scalar>();
				allOps.Kuc_arubber_m[n] = constructCuCSparse<Scalar>();
			}
		}
	}

	// ... Build stiffness matrix K and Kuc, etc...
	buildOps<Scalar>(allOps, 1.0, 0.0, 0.0, rbm, kelArray, (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL, factorize);
}

template<class Scalar>
void
Domain::makeStaticOpsAndSolver(AllOps<Scalar> &allOps, double Kcoef, double Mcoef,
                               double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
                               Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray, FullSquareMatrix *celArray)
{
#ifdef USE_MPI
	FSCommunicator *com = new FSCommunicator(structCom);
#else
	FSCommunicator *com = 0;
#endif
	systemSolver = GenSolverFactory<Scalar>::getFactory()->createSolver(nodeToNode.get(), dsa, c_dsa, *sinfo.solvercntl, spm, rbm,
	                                                                    allOps.spp, allOps.prec, com);
	makeSparseOps<Scalar>(allOps,Kcoef,Mcoef,Ccoef,spm,kelArray,melArray,celArray);
}

template<class Scalar>
void
Domain::makeDynamicOpsAndSolver(AllOps<Scalar> &allOps, double Kcoef, double Mcoef,
                                double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
                                Rbm *rbm, FullSquareMatrix *kelArray, FullSquareMatrix *melArray, FullSquareMatrix *celArray)
{
#ifdef USE_MPI
	FSCommunicator *com = new FSCommunicator(structCom);
#else
	FSCommunicator *com = 0;
#endif
	systemSolver = GenSolverFactory<Scalar>::getFactory()->createSolver(nodeToNode.get(), dsa, c_dsa, *sinfo.solvercntl, spm, rbm,
	                                                                    allOps.spp, allOps.prec, com);
	makeSparseOps<Scalar>(allOps, Kcoef, Mcoef, Ccoef, spm, kelArray, melArray, celArray);
	if(allOps.prec) allOps.prec->factor();
}

template<class Scalar>
void
Domain::addGravityForce(GenVector<Scalar> &force)
{
	// ... ADD ELEMENT MASS CONTRIBUTION ...
	Vector elementGravityForce(maxNumDOFs);
	GenVector<Scalar> gravityForce(force.size(),0.0);
	int gravflg;
	for(int iele = 0; iele < numele; ++iele) {
		if(packedEset[iele]->getProperty() == 0) continue; // phantom element

		if(geoSource->consistentQFlag() && !(sinfo.isDynam() && packedEset[iele]->getMassType() == 0))
			gravflg = 2;                       // 2: consistent (for dynamics, consistent gravity should not be used if element only has a lumped mass matrix)
		else gravflg = geoSource->fixedEndM; // 1: lumped with fixed-end moments
		// 0: lumped without fixed-end moments

		elementGravityForce.zero();
		packedEset[iele]->getGravityForce(nodes, gravityAcceleration, elementGravityForce, gravflg);

		// transform vector from basic to DOF_FRM coordinates
		transformVector(elementGravityForce, iele);

		for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
			int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
			if(cn >= 0) {
				force[cn] += elementGravityForce[idof];
				gravityForce[cn] += elementGravityForce[idof];
			}
		}
	}

	// ... ADD DISCRETE MASS TO GRAVITY FORCE ...
	DMassData *current = firstDiMass;
	while(current != 0) {
		int thisDof = c_dsa->locate(current->node, (1 << current->dof));
		if(thisDof >= 0 ) {
			if(nodes[current->node]->cd == 0) {
				if(current->dof == 0)
					force[thisDof] += (current->diMass)*gravityAcceleration[0];
				if(current->dof == 1)
					force[thisDof] += (current->diMass)*gravityAcceleration[1];
				if(current->dof == 2)
					force[thisDof] += (current->diMass)*gravityAcceleration[2];
			}
			else {
				double g[3] = { gravityAcceleration[0], gravityAcceleration[1], gravityAcceleration[2] };
				transformVector(g, current->node, false);
				if(current->dof == 0)
					force[thisDof] += (current->diMass)*g[0];
				if(current->dof == 1)
					force[thisDof] += (current->diMass)*g[1];
				if(current->dof == 2)
					force[thisDof] += (current->diMass)*g[2];
			}
		}
		current = current->next;
	}
}

template<class Scalar>
void
Domain::addGravityForceSensitivity(GenVector<Scalar> &forceSen)
{
	// ... ADD ELEMENT MASS CONTRIBUTION ...
	Vector elementGravityForceSen(maxNumDOFs);
	int gravflg;
	for(int iele = 0; iele < numele; ++iele) {
		if(packedEset[iele]->getProperty() == 0) continue; // phantom element

		if(geoSource->consistentQFlag() && !(sinfo.isDynam() && packedEset[iele]->getMassType() == 0))
			gravflg = 2;                       // 2: consistent (for dynamics, consistent gravity should not be used if element only has a lumped mass matrix)
		else gravflg = geoSource->fixedEndM; // 1: lumped with fixed-end moments
		// 0: lumped without fixed-end moments

		elementGravityForceSen.zero();
		packedEset[iele]->getGravityForceThicknessSensitivity(nodes, gravityAcceleration, elementGravityForceSen, gravflg);

		// transform vector from basic to DOF_FRM coordinates
		transformVector(elementGravityForceSen, iele);

		for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
			int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
			if(cn >= 0)
				forceSen[cn] += elementGravityForceSen[idof];
		}
	}
}

template<class Scalar>
void
Domain::addPressureForce(GenVector<Scalar> &force, int which, double time)
{
	// which = 0: add only constant pressure
	//         1: add only varying dependent (MFTT/CONWEP)
	//         2: add both constant and varying pressure

	Vector elementPressureForce(maxNumDOFs);
	int cflg = 1; // NOW WE ALWAYS USE CONSISTENT PRESSURE
	PressureBCond *pbc;

	if(pressureFlag()) {
		for(int iele = 0; iele < numele; ++iele) {
			// If there is no pressure boundary condition defined for this element, skip it
			if((pbc = packedEset[iele]->getPressure()) == NULL) continue;

			// Get the force time table
			MFTTData *mftt = domain->getMFTT(pbc->loadsetid);

			// If the pressure is not to be included due to "which" setting, skip it
			bool PressureIsConstant = !((pbc->conwep && pbc->conwepswitch) || mftt);
			if((PressureIsConstant && which==1) || (!PressureIsConstant && which==0)) continue;

			// Compute the amplified pressure due to MFTT or constant load factor
			double loadFactor = (mftt && sinfo.isDynam()) ? mftt->getVal(std::max(time,0.0)) : domain->getLoadFactor(pbc->loadsetid);
			if(loadFactor == 0) continue;
			double p0 = pbc->val;
			pbc->val *= (loadFactor);

			// Compute element pressure force
			elementPressureForce.zero();
			packedEset[iele]->computePressureForce(nodes, elementPressureForce, (GeomState *) 0, cflg, time);
			pbc->val = p0;

			// Transform vector from basic to DOF_FRM coordinates
			transformVector(elementPressureForce, iele);

			// Assemble element pressure forces into domain force vector
			for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
				int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
				if(cn >= 0)
					force[cn] += elementPressureForce[idof];
			}
		}
	}

	for(int iele = 0; iele < numNeum; ++iele) {
		// If this is not pressure boundary condition, skip it
		if((pbc = neum[iele]->getPressure()) == NULL) continue;

		// Get the force time table
		MFTTData *mftt = domain->getMFTT(pbc->loadsetid);

		// If the pressure is not to be included due to "which" setting, skip it
		bool PressureIsConstant = !((pbc->conwep && pbc->conwepswitch) || mftt);
		if((PressureIsConstant && which==1) || (!PressureIsConstant && which==0)) continue;

		// Compute the amplified pressure due to MFTT or constant load factor
		double loadFactor = (mftt && sinfo.isDynam()) ? mftt->getVal(std::max(time,0.0)) : domain->getLoadFactor(pbc->loadsetid);
		if(loadFactor == 0) continue;
		double p0 = pbc->val;
		pbc->val *= loadFactor;

		// Compute structural element distributed Neumann force
		elementPressureForce.zero();
		neum[iele]->neumVector(nodes, elementPressureForce, 0, (GeomState*) 0, time);
		pbc->val = p0;

		// transform vector from basic to DOF_FRM coordinates
		transformNeumVector(elementPressureForce, iele);

		// Assemble element force vector into domain force vector
		int *dofs = neum[iele]->dofs(*c_dsa);
		for(int idof = 0; idof < neum[iele]->numDofs(); ++idof) {
			if(dofs[idof] >= 0)
				force[dofs[idof]] += elementPressureForce[idof];
		}
		delete [] dofs;
	}
}

template<class Scalar>
void
Domain::addAtddnbForce(GenVector<Scalar> &force, int which, double time)
{
	// Get the force time table. Note ATDDNB does not support the LOADSET_ID construct
	MFTTData *mftt = domain->getDefaultMFTT();

	// If the distributed Neumann force is not to be included due to "which" setting, return
	bool AtddnbIsConstant = !mftt;
	if((AtddnbIsConstant && which==1) || (!AtddnbIsConstant && which==0)) return;

	// compute the amplification factor, if applicable
	double loadFactor = (mftt) ? mftt->getVal(time) : 1.0;

	Vector elementAtddnbForce(maxNumDOFs);

	for(int iele = 0; iele < numNeum; ++iele) {

		// Compute acoustic element distributed Neumann force
		elementAtddnbForce.zero();
		neum[iele]->neumVector(nodes, elementAtddnbForce);
		elementAtddnbForce *= sinfo.ATDDNBVal;

		// Assemble element force vector into domain force vector
		int *dofs = neum[iele]->dofs(*c_dsa);
		for(int idof = 0; idof < neum[iele]->numDofs(); ++idof) {
			if(dofs[idof] >= 0)
				force[dofs[idof]] += loadFactor*elementAtddnbForce[idof];
		}
		delete [] dofs;
	}
}

template<class Scalar>
void
Domain::addAtdrobForce(GenVector<Scalar> &force, int which, double time)
{
	// Get the force time table. Note ATDROB does not support the LOADSET_ID construct
	MFTTData *mftt = domain->getDefaultMFTT();

	// If the distributed Robin boundary condition is not to be included due to "which" setting, skip it
	bool AtdrobIsConstant = !mftt;
	if((AtdrobIsConstant && which==1) || (!AtdrobIsConstant && which==0)) return;

	// compute the amplification factor, if applicable
	double loadFactor = (mftt) ? mftt->getVal(time) : 1.0;

	Vector elementAtdrobForce(maxNumDOFs);

	for(int iele = 0; iele < numScatter; ++iele) {

		// Compute acoustic element Robin distributed boundary condition
		elementAtdrobForce.zero();
		scatter[iele]->neumVector(nodes, elementAtdrobForce);
		elementAtdrobForce *= (sinfo.ATDROBVal/sinfo.ATDROBalpha);

		// Assemble element force vector into domain force vector
		int *dofs = scatter[iele]->dofs(*c_dsa);
		for(int idof = 0; idof < scatter[iele]->numDofs(); ++idof) {
			if(dofs[idof] >= 0)
				force[dofs[idof]] += loadFactor*elementAtdrobForce[idof];
		}
		delete [] dofs;
	}
}

template<class Scalar>
void
Domain::addThermalForce(GenVector<Scalar> &force)
{
	if(!temprcvd) initNodalTemperatures();
	Vector elementTemp(maxNumNodes);
	Vector elementThermalForce(maxNumDOFs);
	if(!elemToNode) elemToNode = new Connectivity(packedEset.asSet());

	for(int iele = 0; iele < numele;  ++iele) {
		// By convention phantom elements do not have thermal load
		if(packedEset[iele]->getProperty() == 0) continue;

		// Extract the element nodal temperatures from temprcvd and/or element property ambient temperature
		for(int inod = 0; inod < elemToNode->num(iele); ++inod) {
			double t = temprcvd[(*elemToNode)[iele][inod]];
			elementTemp[inod] = (t == defaultTemp) ? packedEset[iele]->getProperty()->Ta : t;
		}

		// Compute element thermal force in the global coordinates
		elementThermalForce.zero();
		packedEset[iele]->getThermalForce(nodes, elementTemp, elementThermalForce, 0);

		// transform vector from basic to DOF_FRM coordinates
		transformVector(elementThermalForce, iele);

		// Assemble element thermal forces into the force vector
		for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
			int dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
			if(dofNum >= 0)
				force[dofNum] += elementThermalForce[idof];
		}
	}
}

template<class Scalar>
void
Domain::addMpcRhs(GenVector<Scalar> &force, double t)
{
	Vector elementForce(maxNumDOFs);

	for(int iele = 0; iele < numele; ++iele) {
		// If there is element is not a MpcElement, skip it.
		if(!packedEset[iele]->isMpcElement()) continue; // this also works for superelements

		// Otherwise, compute element force due to mpc rhs
		packedEset[iele]->computePressureForce(nodes, elementForce, (GeomState *) 0, 0, t);

		// Assemble element pressure forces into domain force vector
		for(int idof = 0; idof < allDOFs->num(iele); ++idof) {
			int cn = c_dsa->getRCN((*allDOFs)[iele][idof]);
			if(cn >= 0)
				force[cn] += elementForce[idof];
		}
	}
}

template<class Scalar>
void
Domain::scaleDisp(Scalar *u)
{
	// PJSA: 9-22-06 this is for multi-point pade with coupled fluid-structure
	for(int inode = 0; inode < numnodes; ++inode) {
		int cdofs[6];
		c_dsa->number(inode, DofSet::XYZdisp | DofSet::XYZrot, cdofs);
		for(int jdof = 0; jdof<6; ++jdof)
			if(cdofs[jdof] >= 0) u[cdofs[jdof]] *= coupledScaling;
	}
}

template<class Scalar>
void
Domain::scaleInvDisp(Scalar *u)
{
	// RT: 9-10-09 this is for multi-point projection with coupled fluid-structure
	for(int inode = 0; inode < numnodes; ++inode) {
		int cdofs[6];
		c_dsa->number(inode, DofSet::XYZdisp | DofSet::XYZrot, cdofs);
		for(int jdof = 0; jdof<6; ++jdof)
			if(cdofs[jdof] >= 0) u[cdofs[jdof]] /= coupledScaling;
	}
}

template<class Scalar>
void
Domain::scaleDisp(Scalar *u, double alpha)
{
	for(int inode = 0; inode < numnodes; ++inode) {
		int cdofs[6];
		c_dsa->number(inode, DofSet::XYZdisp | DofSet::XYZrot, cdofs);
		for(int jdof = 0; jdof<6; ++jdof)
			if(cdofs[jdof] >= 0) u[cdofs[jdof]] *= alpha;
	}
}

template<class Scalar>
int Domain::mergeDistributedDisp(Scalar (*xyz)[11], Scalar *u, Scalar *bcx, Scalar (*xyz_loc)[11])
{
	// PJSA 9-22-06 u is already scaled
	int inode, nodeI;
	int realNode = -1;

	for(inode = 0; inode < numnodes; ++inode) {

		if(nodes[inode] == 0) continue;
		realNode++;
		nodeI = (outFlag) ? realNode : inode;

		int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);
		int xLoc1 =   dsa->locate(inode, DofSet::Xdisp);

		if (xLoc >= 0)
			xyz[nodeI][0] = u[xLoc];          // free
		else if (xLoc1 >= 0 && bcx)
			xyz[nodeI][0] = bcx[xLoc1];       // constrained
		else
			xyz[nodeI][0] = 0.0;

		int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);
		int yLoc1 =   dsa->locate(inode, DofSet::Ydisp);

		if (yLoc >= 0)
			xyz[nodeI][1] = u[yLoc];
		else if (yLoc1 >= 0 && bcx)
			xyz[nodeI][1] = bcx[yLoc1];
		else
			xyz[nodeI][1] = 0.0;

		int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);
		int zLoc1 =   dsa->locate(inode, DofSet::Zdisp);

		if (zLoc >= 0)
			xyz[nodeI][2] = u[zLoc];
		else if (zLoc1 >= 0 && bcx)
			xyz[nodeI][2] = bcx[zLoc1];
		else
			xyz[nodeI][2] = 0.0;

		int xRot  = c_dsa->locate(inode, DofSet::Xrot);
		int xRot1 =   dsa->locate(inode, DofSet::Xrot);

		if (xRot >= 0)
			xyz[nodeI][3] = u[xRot];
		else if(xRot1 >= 0 && bcx)
			xyz[nodeI][3] = bcx[xRot1];
		else
			xyz[nodeI][3] = 0.0;

		int yRot  = c_dsa->locate(inode, DofSet::Yrot);
		int yRot1 =   dsa->locate(inode, DofSet::Yrot);

		if (yRot >= 0)
			xyz[nodeI][4] = u[yRot];
		else if (yRot1 >= 0 && bcx)
			xyz[nodeI][4] = bcx[yRot1];
		else
			xyz[nodeI][4] = 0.0;

		int zRot  = c_dsa->locate(inode, DofSet::Zrot);
		int zRot1 =   dsa->locate(inode, DofSet::Zrot);

		if (zRot >= 0)
			xyz[nodeI][5] = u[zRot];
		else if (zRot1 >= 0 && bcx)
			xyz[nodeI][5] = bcx[zRot1];
		else
			xyz[nodeI][5] = 0.0;

		int xTemp  = c_dsa->locate(inode, DofSet::Temp);
		int xTemp1 =   dsa->locate(inode, DofSet::Temp);

		if (xTemp >= 0)
			xyz[nodeI][6] = u[xTemp];
		else if (xTemp1 >= 0 && bcx)
			xyz[nodeI][6] = bcx[xTemp1];
		else
			xyz[nodeI][6] = 0.0;

		int xHelm  = c_dsa->locate(inode, DofSet::Helm);
		int xHelm1 =   dsa->locate(inode, DofSet::Helm);

		if (xHelm >= 0)
			xyz[nodeI][7] = u[xHelm];
		else if (xHelm1 >= 0 && bcx)
			xyz[nodeI][7] = bcx[xHelm1];
		else
			xyz[nodeI][7] = 0.0;

		int xPot  = c_dsa->locate(inode, DofSet::Potential);
		int xPot1 =   dsa->locate(inode, DofSet::Potential);

		if (xPot >= 0)
			xyz[nodeI][10] = u[xPot];
		else if (xPot1 >= 0 && bcx)
			xyz[nodeI][10] = bcx[xPot1];
		else
			xyz[nodeI][10] = 0.0;

		// transform displacements and rotations (if present) from DOF_FRM to basic coordinates
		// and keep a copy of the original in xyz_loc
		if(!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
		   && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
			if(xyz_loc) for(int j=0; j<11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
			bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
			transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
		}

	}

	return ++realNode;
}

template<class Scalar>
void Domain::forceDistributedContinuity(Scalar *u, Scalar (*xyz)[11])//DofSet::max_known_nonL_dof
{
	int inode;
	int realNode = -1;
	for (inode = 0; inode < numnodes; ++inode){

		if(nodeToElem)
			if(nodeToElem->num(inode) < 0) continue;
		realNode++;
		int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);

		if (xLoc >= 0)
			u[xLoc] = xyz[inode][0];

		int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);

		if (yLoc >= 0)
			u[yLoc] = xyz[inode][1];

		int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);

		if (zLoc >= 0)
			u[zLoc] = xyz[inode][2];

		int xRot  = c_dsa->locate(inode, DofSet::Xrot);

		if (xRot >= 0)
			u[xRot] = xyz[inode][3];

		int yRot  = c_dsa->locate(inode, DofSet::Yrot);

		if (yRot >= 0)
			u[yRot] = xyz[inode][4];

		int zRot  = c_dsa->locate(inode, DofSet::Zrot);

		if (zRot >= 0)
			u[zRot] = xyz[inode][5];

		int xTemp  = c_dsa->locate(inode, DofSet::Temp);

		if (xTemp >= 0)
			u[xTemp] = xyz[inode][6];

		int xHelm  = c_dsa->locate(inode, DofSet::Helm);

		if (xHelm >= 0)
			u[xHelm] = xyz[inode][7];

		int xPot  = c_dsa->locate(inode, DofSet::Potential);

		if (xPot >= 0)
			u[xPot] = xyz[inode][10];
	}

//  return ++realNode;
}

template<class Scalar>
void
Domain::buildRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *kuc)
{
	double loadFactor; // load amplification factor used for combination load case

	if(! dynamic_cast<GenSubDomain<Scalar>*> (this))
		checkSommerTypeBC(this);
	// ... ZERO FORCE VECTOR INITIALLY
	force.zero();

	// ... COMPUTE EXTERNAL FORCE FROM REAL NEUMAN BC
	int i;
	for(i=0; i < numNeuman; ++i) {
		int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
		if(dof < 0) continue;
		switch(nbc[i].type) {
			case(BCond::Forces) : case(BCond::Flux) : case(BCond::Convection) : case(BCond::Hneu) : {
				double loadFactor = domain->getLoadFactor(nbc[i].loadsetid);
				ScalarTypes::addScalar(force[dof], loadFactor*nbc[i].val);
			} break;
			default :
				ScalarTypes::addScalar(force[dof], nbc[i].val);
		}
	}

	// ... COMPUTE EXTERNAL FORCE FROM COMPLEX NEUMAN BC
	for(i=0; i < numComplexNeuman; ++i) {
		int dof  = c_dsa->locate(cnbc[i].nnum, (1 << cnbc[i].dofnum));
		if(dof < 0) continue;
/*
    switch(cnbc[i].type) {
      case(BCond::Forces) : case(BCond::Flux) : case(BCond::Convection) : case(BCond::Hneu) : {
        double loadFactor = domain->getLoadFactor(cnbc[i].loadsetid);
        ScalarTypes::addScalar(force[dof], loadFactor*cnbc[i].reval, loadFactor*cnbc[i].imval);
      } break;
      default : */
		ScalarTypes::addScalar(force[dof], cnbc[i].reval, cnbc[i].imval);
		//}
	}

	if (implicitFlag) {
		int i, iele;
		double *direction = getWaveDirection();
		ComplexVector elementNeumanScatterForce(this->maxNumDOFs,0.0);
		int* edofs = (int*) dbg_alloca(this->maxNumDOFs*4*sizeof(int));
		for(iele=0; iele<numNeum; ++iele) {
			double kappa = neum[iele]->el->getProperty()->kappaHelm; // PJSA 1-15-2008
			neum[iele]->dofs(*this->dsa,edofs);
			neum[iele]->neumVector(this->nodes,elementNeumanScatterForce,
			                       kappa,direction[0],direction[1],direction[2],
			                       pointSourceFlag);
			for(i=0;i<neum[iele]->numDofs();i++) {
				int cn = c_dsa->getRCN(edofs[i]);
				if(cn >= 0) {
					ScalarTypes::addScalar(force[cn], elementNeumanScatterForce[i].real(),
					                       elementNeumanScatterForce[i].imag());
				}
			}
		}

		ComplexVector elementWetInterfaceScatterForce(this->maxNumDOFs,0.0);
		int numWetDof;
		for(iele=0; iele<numWet; ++iele) {
			HelmElement *he = dynamic_cast<HelmElement *>(wet[iele]->el);
			if (he==0 && wet[iele]->el2==0) {
				numWetDof = wet[iele]->numSolidDofs();
				wet[iele]->solidDofs(*this->dsa,edofs);
			} else if (he!=0 && wet[iele]->el2==0) {
				numWetDof = wet[iele]->numDofs();
				wet[iele]->dofs(*this->dsa,edofs);
			} else {
				numWetDof = wet[iele]->numWetDofs();
				wet[iele]->wetDofs(*this->dsa,edofs);
			}

			complex<double> kappaw =
				sqrt(geoSource->shiftVal())/wet[iele]->soundSpeed;
			wet[iele]->wetInterfaceVector(this->nodes,
			                              elementWetInterfaceScatterForce, real(kappaw),
			                              direction[0], direction[1], direction[2],0,
			                              pointSourceFlag);
			for(i=0;i<numWetDof;i++) {
				int cn = c_dsa->getRCN(edofs[i]);
				if(cn >= 0)
					ScalarTypes::addScalar(force[cn],
					                       elementWetInterfaceScatterForce[i].real(),
					                       elementWetInterfaceScatterForce[i].imag());
			}
		}
	}

	// ... ADD GRAVITY FORCES
	if(gravityFlag()) addGravityForce<Scalar>(force);

	// ... ADD THERMAL FORCES
	if(thermalFlag() && !sinfo.isNonLinExtF()) addThermalForce<Scalar>(force);

	// ... ADD PRESSURE LOAD
	if(!sinfo.isNonLinExtF()) addPressureForce<Scalar>(force);

	// ... ADD LMPC RHS
	if(!sinfo.isNonLin() || sinfo.getNLInfo().linearelastic) addMpcRhs<Scalar>(force);

	// scale RHS force for coupled domains
	if(sinfo.isCoupled) {
		int cdofs[6];  DofSet structdofs = DofSet::XYZdisp | DofSet::XYZrot;
		for(i=0; i<numnodes; ++i) {
			c_dsa->number(i, structdofs, cdofs);
			for(int j=0; j<6; ++j)
				if(cdofs[j] > -1) force[cdofs[j]] *= cscale_factor;
		}
	}

	// COMPUTE NON-HOMOGENEOUS FORCE CONTRIBUTION
	// IN THE CASE OF NONLINEAR, NON-HOMOGENEOUS (PRESCRIBED) FORCES
	// ARE TAKEN CARE OF USING THE GEOMSTATE CLASS, NOT BY
	// MODIFYING THE RHS VECTOR
	if(kuc && !sinfo.isNonLin()) {
		GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);

		// CONSTRUCT NON-HOMONGENOUS DIRICHLET BC VECTOR (PRESCRIBED)
		for(i=0; i<numDirichlet; ++i) {
			int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
			if(dof < 0) continue;
			dof = c_dsa->invRCN(dof);
			if(dof >= 0) {
				if(sinfo.isCoupled && dbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof], dbc[i].val/coupledScaling);
				else ScalarTypes::initScalar(Vc[dof], dbc[i].val);
			}
		}

		// CONSTRUCT NON-HOMONGENOUS COMPLEX DIRICHLET BC VECTOR
		ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
		for(i=0; i<numComplexDirichlet; ++i) {
			int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
			if(dof2 < 0) continue;
			dof2 = c_dsa->invRCN(dof2);
			if(dof2 >= 0) {
				if(sinfo.isCoupled && cdbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval/coupledScaling, cdbcMRHS[i].imval/coupledScaling);
				else ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
			}
		}

		kuc->multSubtract(Vc, force);
	}
}

template<class Scalar>
void
Domain::computeReactionForce(GenVector<Scalar> &fc, GenVector<Scalar> &Vu,
                             GenSparseMatrix<Scalar> *_kuc, GenSparseMatrix<Scalar> *_kcc)
{
	// TODO include external force on the constrained dofs
	GenCuCSparse<Scalar> *kuc = dynamic_cast<GenCuCSparse<Scalar> *>(_kuc);
	if(kuc) kuc->transposeMultNew(Vu.data(), fc.data()); // fc = kuc^T * Vu
	else fc.zero();

	GenCuCSparse<Scalar> *kcc = dynamic_cast<GenCuCSparse<Scalar> *>(_kcc);
	if(kcc) {
		GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);

		for(int i=0; i<numDirichlet; ++i) {
			int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
			if(dof < 0) continue;
			dof = c_dsa->invRCN(dof);
			if(dof >= 0) {
				if(sinfo.isCoupled && dbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof], dbc[i].val/coupledScaling); else
					ScalarTypes::initScalar(Vc[dof], dbc[i].val);
			}
		}

		ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
		for(int i=0; i<numComplexDirichlet; ++i) {
			int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
			if(dof2 < 0) continue;
			dof2 = c_dsa->invRCN(dof2);
			if(dof2 >= 0) {
				if(sinfo.isCoupled && cdbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval/coupledScaling, cdbcMRHS[i].imval/coupledScaling); else
					ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
			}
		}

		kcc->multAddNew(Vc.data(), fc.data()); // fc += Kcc * Vc
	}
}

template<class Scalar>
void
Domain::buildRHSForce(GenVector<Scalar> &force, GenVector<Scalar> &tmp,
                      GenSparseMatrix<Scalar> *kuc, GenSparseMatrix<Scalar> *muc,
                      GenSparseMatrix<Scalar> **cuc_deriv,
                      GenSparseMatrix<Scalar> **kuc_deriv,
                      GenSparseMatrix<Scalar> **kuc_arubber_l,
                      GenSparseMatrix<Scalar> **kuc_arubber_m,
                      double omega, double delta_omega, GeomState *gs)
{
	if(! dynamic_cast<GenSubDomain<Scalar>*> (this))
		checkSommerTypeBC(this);
	// ... ZERO FORCE VECTOR INITIALLY
	force.zero();

	// ... COMPUTE EXTERNAL FORCE FROM REAL NEUMAN BC
	int i;
	for(i=0; i < numNeuman; ++i) {
		int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
		if(dof < 0) continue;
		switch(nbc[i].type) {
			case(BCond::Forces) : case(BCond::Flux) : case(BCond::Convection) : case(BCond::Hneu) : {
				double loadFactor = domain->getLoadFactor(nbc[i].loadsetid);
				ScalarTypes::addScalar(force[dof], loadFactor*nbc[i].val);
			} break;
			default :
				ScalarTypes::addScalar(force[dof], nbc[i].val);
		}
	}

	// ... COMPUTE EXTERNAL FORCE FROM COMPLEX NEUMAN BC
	for(i=0; i < numComplexNeuman; ++i) {
		int dof  = c_dsa->locate(cnbc[i].nnum, (1 << cnbc[i].dofnum));
		if(dof < 0) continue;
/*
    switch(nbc[i].type) {
      case(BCond::Forces) : case(BCond::Flux) : case(BCond::Convection) : case(BCond::Hneu) : {
        double loadFactor = domain->getLoadFactor(cnbc[i].loadsetid);
        ScalarTypes::addScalar(force[dof], loadFactor*cnbc[i].reval, loadFactor*cnbc[i].imval);
      } break;
      default :*/
		ScalarTypes::addScalar(force[dof], cnbc[i].reval, cnbc[i].imval);
//    }
	}

	// PJSA: new FETI-H stuff
	if (implicitFlag) {
		int i, iele;
		double *direction = getWaveDirection();
// RT
		ComplexVector elementNeumanScatterForce(this->maxNumDOFs,0.0);
		int* edofs = (int*) dbg_alloca(this->maxNumDOFs*sizeof(int));
		for(iele=0; iele<numNeum; ++iele) {
			double kappa = omega/ScalarTypes::Real(neum[iele]->soundSpeed);
			//double kappa = neum[iele]->el->getProperty()->kappaHelm; // PJSA 1-15-2008
			neum[iele]->dofs(*this->dsa,edofs);
			neum[iele]->neumVector(this->nodes,elementNeumanScatterForce,
			                       kappa,direction[0],direction[1],direction[2],
			                       pointSourceFlag);
			for(i=0;i<neum[iele]->numDofs();i++) {
				int cn = c_dsa->getRCN(edofs[i]);
				if(cn >= 0) {
					ScalarTypes::addScalar(force[cn], elementNeumanScatterForce[i].real(),
					                       elementNeumanScatterForce[i].imag());
				}
			}
		}
	}

	if (implicitFlag) {

		double *direction = getWaveDirection();

		int iele;
		ComplexVector elementWetInterfaceScatterForce(this->maxNumDOFs,0.0);
		int* edofs = (int*) alloca(this->maxNumDOFs*4*sizeof(int));
		int numWetDof;
		for(iele=0; iele<numWet; ++iele) {
			HelmElement *he = dynamic_cast<HelmElement *>(wet[iele]->el);
			if (he==0 && wet[iele]->el2==0) {
				numWetDof = wet[iele]->numSolidDofs();
				wet[iele]->solidDofs(*this->dsa,edofs);
			} else if (he!=0 && wet[iele]->el2==0) {
				numWetDof = wet[iele]->numDofs();
				wet[iele]->dofs(*this->dsa,edofs);
			} else {
				numWetDof = wet[iele]->numWetDofs();
				wet[iele]->wetDofs(*this->dsa,edofs);
			}

			complex<double> kappaw =
				omega/wet[iele]->soundSpeed;
			wet[iele]->wetInterfaceVector(this->nodes,
			                              elementWetInterfaceScatterForce, real(kappaw),
			                              direction[0], direction[1], direction[2],0,
			                              pointSourceFlag);
			int i;
			for(i=0;i<numWetDof;i++) {
				int cn = c_dsa->getRCN(edofs[i]);
				if(cn >= 0)
					ScalarTypes::addScalar(force[cn],
					                       elementWetInterfaceScatterForce[i].real(),
					                       elementWetInterfaceScatterForce[i].imag());
			}
		}
	}

	// ... ADD GRAVITY FORCES
	if(gravityFlag()) addGravityForce<Scalar>(force);

	// ... ADD THERMAL FORCES
	if(thermalFlag() && !sinfo.isNonLinExtF()) addThermalForce<Scalar>(force);

	// ... ADD PRESSURE LOAD
	if(!sinfo.isNonLinExtF()) addPressureForce<Scalar>(force);

	GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);

	// CONSTRUCT NON-HOMONGENOUS DIRICHLET BC VECTOR (PRESCRIBED)
	for(i=0; i<numDirichlet; ++i) {
		int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
		if(dof < 0) continue;
		dof = c_dsa->invRCN(dof);
		if(dof >= 0) {
			if(sinfo.isCoupled && dbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof], dbc[i].val/coupledScaling);
			else ScalarTypes::initScalar(Vc[dof], dbc[i].val);
		}
	}

	// CONSTRUCT NON-HOMONGENOUS COMPLEX DIRICHLET BC VECTOR
	ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
	for(i=0; i<numComplexDirichlet; ++i) {
		int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
		if(dof2 < 0) continue;
		dof2 = c_dsa->invRCN(dof2);
		if(dof2 >= 0) {
			if(sinfo.isCoupled && cdbc[i].dofnum < 6) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval/coupledScaling, cdbcMRHS[i].imval/coupledScaling);
			else ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
		}
	}

	// scale RHS force for coupled domains
	if(sinfo.isCoupled) {
		int cdofs[6];  DofSet structdofs = DofSet::XYZdisp | DofSet::XYZrot;
		for(i=0; i<numnodes; ++i) {
			c_dsa->number(i, structdofs, cdofs);
			for(int j=0; j<6; ++j)
				if(cdofs[j] > -1) force[cdofs[j]] *= cscale_factor;
		}
	}

	// COMPUTE NON-HOMOGENEOUS FORCE CONTRIBUTION
	// IN THE CASE OF NONLINEAR, NON-HOMOGENEOUS (PRESCRIBED) FORCES
	// ARE TAKEN CARE OF USING THE GEOMSTATE CLASS, NOT BY
	// MODIFYING THE RHS VECTOR
	if(probType() != SolverInfo::NonLinStatic &&
	   probType() != SolverInfo::NonLinDynam  &&
	   probType() != SolverInfo::ArcLength) {
		if(kuc) {
			kuc->multSubtract(Vc, force);
		}
		complex<double> *lambda=0, *mu=0, *deltalambda=0, *deltamu=0;
		int num_arubber = geoSource->num_arubber;
		if (num_arubber>0)  {
			lambda = new complex<double>[num_arubber];
			mu = new complex<double>[num_arubber];
			deltalambda = new complex<double>[num_arubber];
			deltamu = new complex<double>[num_arubber];
		}
		geoSource->getARubberLambdaMu(omega,deltalambda,deltamu);
		geoSource->getARubberLambdaMu(omega-delta_omega,lambda,mu);
		for(int ir=0;ir<num_arubber;ir++) {
			deltalambda[ir] -= lambda[ir];
			deltamu[ir] -= mu[ir];
		}

// Acoustic rubber    
		for(int ir=0;ir<num_arubber;ir++) {
			tmp.zero();
			if(kuc_arubber_l) if (kuc_arubber_l[ir]) kuc_arubber_l[ir]->multSubtract(Vc,tmp);
			Scalar c = 0;
			ScalarTypes::addComplex(c,deltalambda[ir]);
			tmp *= c;
			force += tmp;
			tmp.zero();
			if (kuc_arubber_m) if (kuc_arubber_m[ir]) kuc_arubber_m[ir]->multSubtract(Vc,tmp);
			c = 0;
			ScalarTypes::addComplex(c,deltamu[ir]);
			tmp *= c;
			force += tmp;
		}

		if (muc) {
			tmp.zero();
			muc->multSubtract(Vc, tmp);
			tmp *= (omega-delta_omega)*(omega-delta_omega)-omega*omega;
			force += tmp;
		}
		if(cuc_deriv) if (cuc_deriv[0]) {
				tmp.zero();
				cuc_deriv[0]->multSubtract(Vc, tmp);
				Scalar c;
				ScalarTypes::initScalar(c,0.0,delta_omega);
				tmp *= c;
				force += tmp;
			}
	}
}

template<class Scalar>
void
Domain::buildFreqSweepRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *muc,
                               GenSparseMatrix<Scalar> **cuc_deriv,
                               GenSparseMatrix<Scalar> **kuc_deriv,
                               int iRHS, double omega)
{
	if(iRHS < 1) return; // this shouldn't happen

	// PJSA 5-15-06
	if (implicitFlag) {
		double *direction = getWaveDirection();
		ComplexVector elementNeumanScatterForce(this->maxNumDOFs,0.0);
		for(int iele=0; iele<numNeum; ++iele) {
			int *dofs = neum[iele]->dofs(*this->dsa);
			double kappa = neum[iele]->el->getProperty()->kappaHelm; // PJSA 1-15-2008
			neum[iele]->neumVectorDeriv(this->nodes,elementNeumanScatterForce,
			                            kappa,direction[0],direction[1],direction[2],iRHS, pointSourceFlag);
			for(int i=0;i<neum[iele]->numDofs();i++) {
				int cn = c_dsa->getRCN(dofs[i]);
				if(cn >= 0) {
					ScalarTypes::addScalar(force[cn],
					                       elementNeumanScatterForce[i].real() /
					                       pow(real(neum[iele]->soundSpeed),iRHS),
					                       elementNeumanScatterForce[i].imag() /
					                       pow(real(neum[iele]->soundSpeed),iRHS));
				}
			}
		}
	}

	if (implicitFlag) {

		double *direction = getWaveDirection();

		int iele;
		ComplexVector elementWetInterfaceScatterForce(this->maxNumDOFs,0.0);
		int* edofs = (int*) alloca(this->maxNumDOFs*4*sizeof(int));
		int numWetDof;
		for(iele=0; iele<numWet; ++iele) {
			HelmElement *he = dynamic_cast<HelmElement *>(wet[iele]->el);
			if (he==0 && wet[iele]->el2==0) {
				numWetDof = wet[iele]->numSolidDofs();
				wet[iele]->solidDofs(*this->dsa,edofs);
			} else if (he!=0 && wet[iele]->el2==0) {
				numWetDof = wet[iele]->numDofs();
				wet[iele]->dofs(*this->dsa,edofs);
			} else {
				numWetDof = wet[iele]->numWetDofs();
				wet[iele]->wetDofs(*this->dsa,edofs);
			}

			complex<double> kappaw =
				sqrt(geoSource->shiftVal())/wet[iele]->soundSpeed;
			wet[iele]->wetInterfaceVector(this->nodes,
			                              elementWetInterfaceScatterForce, real(kappaw),
			                              direction[0], direction[1], direction[2],iRHS,
			                              pointSourceFlag);
			if (he==0 && wet[iele]->el2==0) {
				for(int i=0;i<numWetDof;i++)
					elementWetInterfaceScatterForce[i] *= cscale_factor;
			} else if (he!=0 && wet[iele]->el2==0) {
			} else {
				for(int i=0;i<numWetDof/4;i++) {
					elementWetInterfaceScatterForce[4*i+0] *= cscale_factor;
					elementWetInterfaceScatterForce[4*i+1] *= cscale_factor;
					elementWetInterfaceScatterForce[4*i+2] *= cscale_factor;
				}
			}

			int i;
			for(i=0;i<numWetDof;i++) {
				int cn = c_dsa->getRCN(edofs[i]);
				if(cn >= 0)
					ScalarTypes::addScalar(force[cn],
					                       elementWetInterfaceScatterForce[i].real()/
					                       pow(real(wet[iele]->soundSpeed),iRHS),
					                       elementWetInterfaceScatterForce[i].imag()/
					                       pow(real(wet[iele]->soundSpeed),iRHS));
			}
		}
	}

	int i;
	GenVector<Scalar> Vc(numDirichlet+numComplexDirichlet, 0.0);
	GenVector<Scalar> Vc_tmp(numDirichlet+numComplexDirichlet);

	// CONSTRUCT NON-HOMONGENOUS DIRICHLET BC VECTOR (PRESCRIBED)
	for(i=0; i<numDirichlet; ++i) {
		int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
		if(dof < 0) continue;
		dof = c_dsa->invRCN(dof);
		if(dof >= 0) ScalarTypes::initScalar(Vc[dof], dbc[i].val);
	}

	// CONSTRUCT NON-HOMONGENOUS COMPLEX DIRICHLET BC VECTOR
	ComplexBCond *cdbcMRHS = cdbc + iWaveDir * numComplexDirichlet;
	for(i=0; i<numComplexDirichlet; ++i) {
		int dof2 = dsa->locate(cdbc[i].nnum,(1 << cdbc[i].dofnum));
		if(dof2 < 0) continue;
		dof2 = c_dsa->invRCN(dof2);
		if(dof2 >= 0) ScalarTypes::initScalar(Vc[dof2], cdbcMRHS[i].reval, cdbcMRHS[i].imval);
	}

	if(muc && iRHS <= 2) {
		if(iRHS==1) Vc_tmp.linC(Vc, -2.0*omega);
		else Vc_tmp.linC(Vc,-2.0);
		muc->multSubtract(Vc_tmp, force);
	}

	if(cuc_deriv && cuc_deriv[iRHS-1]) cuc_deriv[iRHS-1]->multSubtract(Vc, force);
	if(kuc_deriv && kuc_deriv[iRHS-1]) kuc_deriv[iRHS-1]->multSubtract(Vc, force);
}


template<class Scalar>
void
Domain::buildDeltaK(double w0, double w, GenSparseMatrix<Scalar> *deltaK,
                    GenSparseMatrix<Scalar> *deltaKuc) {
	// Acoustic rubber elements
	complex<double> *krarray = new complex<double>[3*maxNumDOFs*maxNumDOFs];
	complex<double> *krarray2 = new complex<double>[3*maxNumDOFs*maxNumDOFs];
	for(int iele = 0; iele < numele; ++iele) {
		StructProp *prop = packedEset[iele]->getProperty();
		if(!prop || (prop && (prop->E0==0.0 || prop->mu0==0.0)) ) continue;

		if(matrixTimers) matrixTimers->formTime -= getTime();
		int N = 0; // number of derivatives
		packedEset[iele]->aRubberStiffnessDerivs(nodes, krarray,N,w0);
		packedEset[iele]->aRubberStiffnessDerivs(nodes, krarray2,N,w);
		int ndof = packedEset[iele]->numDofs();
		for(int i=0;i<ndof*ndof;i++)
			krarray2[i+2*ndof*ndof] -= krarray[i+2*ndof*ndof];
		if(sinfo.isCoupled) if(isStructureElement(iele)) {  // coupled scaling
				for(int i=0;i<ndof*ndof;i++)  krarray2[i+2*ndof*ndof] *= cscale_factor2;
			}

		FullSquareMatrixC kdel(ndof,krarray2+2*ndof*ndof);
		deltaK->add(kdel,(*allDOFs)[iele]);
		if (deltaKuc) deltaKuc->add(kdel,(*allDOFs)[iele]);

	}
	delete[] krarray;
	delete[] krarray2;
	fprintf(stderr,"Domain::buildDeltaK exiting\n");
}

template<class Scalar>
void
Domain::assembleATDROB(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops,double Kcoef)
{
//  fprintf(stderr, "ATDROB ok\n");
//  checkSommerTypeBC(this);
	double *v = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
	if(numScatter > 0) {
		int i;
		double temp = sinfo.ATDROBbeta/sinfo.ATDROBalpha;
		for(i=0; i<numScatter; i++) {
			int *dofs = scatter[i]->dofs(*dsa);
			FullSquareMatrix ms = scatter[i]->sommerMatrix(nodes,v);//sommerMatrix is negative definite...
			FullSquareMatrix mm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
			ms.multiply(mm,temp);
			updateMatrices<Scalar>(ops,K,dofs,0,&mm,0,Kcoef); //RT : 02112013 - needs to be done
		}
	}
}

template<class Scalar>
void
Domain::assembleSommer(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops)
{
	checkSommerTypeBC(this); // TODO check

	if(numSommer > 0) {
		if(sommer[0]->dim() == 3) {
			if(sommerfeldType==1) sommerfeldType = 3;
			if(sommerfeldType==2) sommerfeldType = 4;
		}

		if(sommerfeldType == 1 )
			fprintf(stderr, " ... 1st order Bayliss-Turkel       ...\n");
		else if(sommerfeldType == 2 )
			fprintf(stderr, " ... 2nd order Bayliss-Turkel       ...\n");
		else if(sommerfeldType == 3 )
			fprintf(stderr, " ... 1st order 3D Bayliss-Turkel    ...\n");
		else if(sommerfeldType == 4 )
			fprintf(stderr, " ... 2nd order 3D Bayliss-Turkel    ...\n");

		if(sommerfeldType == 1 || sommerfeldType == 2)
			getCurvatures(this);
		if(curvatureFlag != 1) {
			if(sommerfeldType == 3 || sommerfeldType == 4)
				getCurvatures3D(this);
		}

		double *v = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
		double *vbt = (double *) dbg_alloca(maxNumDOFs*maxNumDOFs*sizeof(double));
		GenSubDomain<Scalar> *subCast = dynamic_cast<GenSubDomain<Scalar>*>(this);
		bool mdds_flag = (K && subCast && sinfo.solvercntl->type == SolverSelection::Direct); // multidomain direct solver

		// This loops adds the contribution of the terms emanating
		// from the non-reflecting boundary conditions
		//
		for(int i=0; i<numSommer; i++) {
			ComplexD *bt2Matrix = 0;
			ComplexD **bt2nMatrix = 0;
			int *dofs = sommer[i]->dofs(*dsa);
			int *dofs_mdds = dofs;
			if(mdds_flag) {
				// TODO Get rid of this horror!!! Why do we go to global????
				auto glnums = const_cast<int *>(subCast->getGlNodes().data());
				sommer[i]->renum(glnums);
				dofs_mdds = sommer[i]->dofs(*domain->getDSA());

				sommer[i]->renum(subCast->getGlobalToLocalNode());
			}
			FullSquareMatrix ms = sommer[i]->sommerMatrix(nodes,v);
			FullSquareMatrix mm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
			double kappa = sommer[i]->el->getProperty()->kappaHelm; // PJSA 1-15-2008
			double ss = real(sommer[i]->el->getProperty()->soundSpeed);

			ms.multiply(mm,kappa);

			// "Zero-order" term
			updateMatrices(ops,K,dofs,dofs_mdds,0,&mm);

			double psi; // curvature of the boundary
			double HH, KK;
			// 1st order Bayliss-Turkel boundary condition
			if(sommerfeldType == 1 ) {
				psi = curvatures[i];
				HH = psi/2.0;
				ms.multiply(mm,-HH);
				updateMatrices(ops,K,dofs,dofs_mdds,&mm,0);
			}

				// 2nd order Bayliss-Turkel boundary condition
			else if(sommerfeldType == 2) {
				FullSquareMatrix ks = sommer[i]->turkelMatrix(nodes,vbt);
				psi = curvatures[i];
				// fprintf(stderr,"curvature %f\n",psi);

				ComplexD cm = ComplexD(-0.5*psi,0.0) + psi*psi/8.0/ComplexD(psi,-kappa);
				ComplexD cs = -0.5/ComplexD(psi,-kappa);
				double c1,c2,c3,c4;
				c1 = imag(cm);
				c2 = real(cm);
				c3 = imag(cs);
				c4 = real(cs);

				FullSquareMatrix mm1(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
				// mm = ms*c2 + ks*c4;
				ms.multiply(mm,c2);
				ks.multiply(mm1,c4);
				mm += mm1;
				updateMatrices(ops,K,dofs,dofs_mdds,&mm,0);
				// mm = ms*c1 + ks*c3;
				ms.multiply(mm,c1);
				ks.multiply(mm1,c3);
				mm += mm1;
				updateMatrices(ops,K,dofs,dofs_mdds,0,&mm);
			}

				// 1st order 3D Bayliss-Turkel boundary condition
			else if(sommerfeldType == 3) {
				if(curvatureFlag != 1) {
					HH = 0.0;
					for(int iNode=0;iNode<somElemToNode->num(i);iNode++) {
						int iNodeNumber = (*somElemToNode)[i][iNode];
						int iSommerNode = nodeToSommerNodeMap[iNodeNumber];
						HH += curvaturesH[iSommerNode];
					}
					HH /= somElemToNode->num(i);
				}
				else {
					HH = 1.0/curvatureConst1;
				}
				ms.multiply(mm,-HH);
				updateMatrices(ops,K,dofs,dofs_mdds,&mm,0);
			}

				// 2nd order 3D Bayliss-Turkel boundary condition
			else if(sommerfeldType == 4) {
				HH = 0.0;
				KK = 0.0;
				int iNode;
				if(curvatureFlag != 1) {
					for(iNode=0;iNode<somElemToNode->num(i);iNode++) {
						int iNodeNumber = (*somElemToNode)[i][iNode];
						int iSommerNode = nodeToSommerNodeMap[iNodeNumber];
						HH += curvaturesH[iSommerNode];
						KK += curvaturesK[iSommerNode];
					}
					HH /= somElemToNode->num(i);
					KK /= somElemToNode->num(i);
				}
				else {
					HH = 0.5*(1.0/curvatureConst1+1.0/curvatureConst1);
					KK = 1.0/curvatureConst1*1.0/curvatureConst1;
				}
				ComplexD ii=ComplexD(0.0, 1.0);
				ComplexD cm = -ii/2.0/kappa/(1.0+ii*2.0*HH/kappa)*(KK-HH*HH) - ComplexD(HH,0.0);

				double c1,c2;
				c1 = imag(cm);
				c2 = real(cm);

				if(curvatureFlag != 2) {
					ms.multiply(mm,c2);
					updateMatrices(ops,K,dofs,dofs_mdds,&mm,0);
					ms.multiply(mm,c1);
					updateMatrices(ops,K,dofs,dofs_mdds,0,&mm);
				}

				bt2Matrix = new DComplex [ms.dim()*ms.dim()*sizeof(DComplex)];
				for(int j=0; j<ms.dim()*ms.dim(); ++j) bt2Matrix[j] = 0.0;  // PJSA
				if(solInfo().doFreqSweep) {
					int N = solInfo().getSweepParams()->nFreqSweepRHS - 1; // number of derivatives
					bt2nMatrix = new DComplex * [N];
					for(int n=1; n<=N; ++n) {
						bt2nMatrix[n-1] = new DComplex [ms.dim()*ms.dim()*sizeof(DComplex)];
						for(int j=0; j<ms.dim()*ms.dim(); ++j) bt2nMatrix[n-1][j] = 0.0;  // PJSA
					}
				}
				FullSquareMatrix ksRe(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
				FullSquareMatrix ksIm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));

				if(curvatureFlag == 2) {
					double HHH[3];
					double KKK[3];
					for(iNode=0;iNode<somElemToNode->num(i);iNode++) {
						int iNodeNumber = (*somElemToNode)[i][iNode];
						int iSommerNode = nodeToSommerNodeMap[iNodeNumber];
						HHH[iNode] = curvaturesH[iSommerNode];
						KKK[iNode] = curvaturesK[iSommerNode];
					}

					sommer[i]->sommerMatrixEllipsoid(nodes, kappa, HHH, KKK, bt2Matrix);

					int kDof=0;
					int iDof, jDof;
					for (iDof=0; iDof<ms.dim(); iDof++) {
						for (jDof=0; jDof<ms.dim(); jDof++) {
							ksRe[iDof][jDof] = real(bt2Matrix[kDof]);
							ksIm[iDof][jDof] = imag(bt2Matrix[kDof]);
							kDof++;
						}
					}

					updateMatrices(ops,K,dofs,dofs_mdds,&ksRe,&ksIm);
				}

				double curv_e[3];
				double curv_f[3];
				double curv_g[3];
				double tau1[3][3];
				double tau2[3][3];
				if (curvatureFlag != 1) {
					int n1 = (*somElemToNode)[i][0];
					int n2 = (*somElemToNode)[i][1];
					int n3 = (*somElemToNode)[i][2];
					int nn[3];
					nn[0] = nodeToSommerNodeMap[n1];
					nn[1] = nodeToSommerNodeMap[n2];
					nn[2] = nodeToSommerNodeMap[n3];
					double *curv_normal[3];
					int iiNode;
					for (iiNode=0; iiNode<3; iiNode++) {
						int curvNode = nn[iiNode];
						curv_e[iiNode]      = curvatures_e[curvNode];
						curv_f[iiNode]      = curvatures_f[curvNode];
						curv_g[iiNode]      = curvatures_g[curvNode];
						curv_normal[iiNode] = curvatures_normal[nn[iiNode]];
						getTau(curv_normal[iiNode], tau1[iiNode], tau2[iiNode]);

						// This is for experiments with exact curvatures for ellipsoid
						if (curvatureFlag == 2) {
							tau1[iiNode][0] = curvatures_tau1[curvNode][0];
							tau1[iiNode][1] = curvatures_tau1[curvNode][1];
							tau1[iiNode][2] = curvatures_tau1[curvNode][2];
							tau2[iiNode][0] = curvatures_tau2[curvNode][0];
							tau2[iiNode][1] = curvatures_tau2[curvNode][1];
							tau2[iiNode][2] = curvatures_tau2[curvNode][2];
						}
					}
				}
				else {
					curv_e[0] = curv_e[1] = curv_e[2] = 1.0/curvatureConst1;
					curv_g[0] = curv_g[1] = curv_g[2] = 1.0/curvatureConst1;
					curv_f[0] = curv_f[1] = curv_f[2] = 0.0;
					int nds[3],iiNode;
					sommer[i]->nodes(nds);
					for (iiNode=0; iiNode<3; iiNode++) {
						double nrmal[3];
						Node nd = nodes.getNode(nds[iiNode]);
						double l = sqrt(nd.x*nd.x+nd.y*nd.y+nd.z*nd.z);
						nrmal[0] = nd.x/l;
						nrmal[1] = nd.y/l;
						nrmal[2] = nd.z/l;
						getTau(nrmal, tau1[iiNode], tau2[iiNode]);
					}
				}
				if (curvatureFlag != 1) {
					sommer[i]->BT2(nodes, curv_e, curv_f, curv_g, tau1, tau2, kappa, bt2Matrix);
					if(solInfo().doFreqSweep) {
						int N = solInfo().getSweepParams()->nFreqSweepRHS - 1; // number of derivatives
						for(int n=1; n<=N; ++n)
							sommer[i]->BT2n(nodes, curv_e, curv_f, curv_g, tau1, tau2, kappa, bt2nMatrix[n-1], n);
					}
				}
				else {
					sommer[i]->sphereBT2(nodes, curvatureConst1 , kappa, bt2Matrix);
				}
				if (curvatureFlag == 2) {
					sommer[i]->ellipsoidBT2(nodes, curvatureConst2, curvatureConst1, kappa, bt2Matrix);
				}

				int kDof=0;
				int iDof, jDof;
				for (iDof=0; iDof<ms.dim(); iDof++) {
					for (jDof=0; jDof<ms.dim(); jDof++) {
						ksRe[iDof][jDof] = real(-1.0/(2.0*ii*kappa)*bt2Matrix[kDof]);
						ksIm[iDof][jDof] = imag(-1.0/(2.0*ii*kappa)*bt2Matrix[kDof]);
						kDof++;
					}
				}

				updateMatrices(ops,K,dofs,dofs_mdds,&ksRe,&ksIm);
			}

			if(solInfo().doFreqSweep) {
				switch(sommerfeldType) {
					case 0:
					case 1:
					case 3:
						updateDampingMatrices(ops, dofs, 0, &ms, ss, 1); // zero- and first-order sommerfeld
						break;
					case 4:
						computeSommerDerivatives(HH, KK, curvatureFlag, dofs, ms, bt2nMatrix, kappa, ss, ops); // 3D second-order sommerfeld
						break;
					case 2:
					default:
						std::cerr << " *** ERROR: Sommerfeld type " << sommerfeldType << " is not supported for frequency sweep \n";
						break;
				}
			}
			if(bt2Matrix) delete [] bt2Matrix;
			if(bt2nMatrix) {
				int N = solInfo().getSweepParams()->nFreqSweepRHS - 1; // number of derivatives
				for(int n=1; n<=N; ++n) if(bt2nMatrix[n-1]) delete [] bt2nMatrix[n-1];
				delete [] bt2nMatrix;
			}
			delete [] dofs;
		}
	}
}

template<class Scalar>
void
Domain::computeSommerDerivatives(double HH, double KK, int curvatureFlag, int *dofs, FullSquareMatrix &ms,
                                 DComplex **bt2nMatrix, double kappa, double ss, AllOps<Scalar> *ops)
{
	// PJSA 5-26-05
	// this function is for 3D second-order sommerfeld
	FullSquareMatrix mm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
	ComplexD ii=ComplexD(0.0, 1.0);
	int N = solInfo().getSweepParams()->nFreqSweepRHS - 1;  // number of derivatives
	for(int n=1; n<=N; ++n) {
		DComplex cm = (n==1) ? ii : 0;
		cm -= ((KK-HH*HH)*ii/2.0 * pow(-1.0,n)*double(DFactorial(n))/pow((kappa+ii*2.0*HH),n+1));
		if(curvatureFlag != 2) {
			ms.multiply(mm,real(cm)); // mm = real part of cm*ms
			updateDampingMatrices(ops,dofs,&mm,0,ss,n);
			ms.multiply(mm,imag(cm)); // mm = imaginary part of cm*ms
			updateDampingMatrices(ops,dofs,0,&mm,ss,n);
		}
		else std::cerr << " *** WARNING: 3D 2nd order Sommerfeld with curvatureFlag 2 is not supported for frequency sweep \n";
	}

	FullSquareMatrix ksRe(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
	FullSquareMatrix ksIm(ms.dim(),(double*)dbg_alloca(ms.dim()*ms.dim()*sizeof(double)));
	for(int n=1; n<=N; ++n) {
		DComplex cm = -pow(1.0/(2.0*ii*kappa),n+1)*pow(-2.0*ii,n)*double(DFactorial(n));
		int kDof = 0;
		for(int iDof=0; iDof<ms.dim(); iDof++) {
			for(int jDof=0; jDof<ms.dim(); jDof++) {
				DComplex cmbt = cm*bt2nMatrix[n-1][kDof];
				ksRe[iDof][jDof] = real(cmbt);
				ksIm[iDof][jDof] = imag(cmbt);
				kDof++;
			}
		}
		updateDampingMatrices(ops,dofs,&ksRe,&ksIm,ss,n);
	}
}

template<class Scalar>
void
Domain::updateMatrices(AllOps<Scalar> *ops, GenSparseMatrix<Scalar> *Z, int *dofs, int *dofs_mdds,
                       FullSquareMatrix *reEl, FullSquareMatrix *imEl, double Kcoef)
{
	if((sinfo.isATDARB()) || (sinfo.ATDROBalpha != 0.0)) {
		if(reEl) {
			FullSquareMatrix temp(reEl->dim(),(double*)dbg_alloca(reEl->dim()*reEl->dim()*sizeof(double)));
			reEl->multiply(temp, Kcoef);
#if defined(_OPENMP)
#pragma omp critical
#endif
			if(Z) Z->add(temp, dofs_mdds);
			if(ops && ops->spp) ops->spp->add(temp, dofs);
			if(ops && ops->Kuc) ops->Kuc->add(temp, dofs);
			if(ops && ops->Kcc) ops->Kcc->add(temp, dofs);
		}
		if(imEl) {
			FullSquareMatrix temp(imEl->dim(),(double*)dbg_alloca(imEl->dim()*imEl->dim()*sizeof(double)));
			imEl->multiply(temp, Kcoef);
#if defined(_OPENMP)
#pragma omp critical
#endif
			if(Z) Z->addImaginary(*imEl, dofs_mdds);
			if(ops && ops->spp) ops->spp->addImaginary(temp, dofs);
			if(ops && ops->Kuc) ops->Kuc->addImaginary(temp, dofs);
			if(ops && ops->Kcc) ops->Kcc->addImaginary(temp, dofs);
		}
	}
	else {
		if(reEl) {
#if defined(_OPENMP)
#pragma omp critical
#endif
			if(Z) Z->add(*reEl, dofs_mdds);
			if(ops && ops->spp) ops->spp->add(*reEl, dofs);
			if(ops && ops->Kuc) ops->Kuc->add(*reEl, dofs);
			if(ops && ops->Kcc) ops->Kcc->add(*reEl, dofs);
		}
		if(imEl) {
#if defined(_OPENMP)
#pragma omp critical
#endif
			if(Z) Z->addImaginary(*imEl, dofs_mdds);
			if(ops && ops->spp) ops->spp->addImaginary(*imEl, dofs);
			if(ops && ops->Kuc) ops->Kuc->addImaginary(*imEl, dofs);
			if(ops && ops->Kcc) ops->Kcc->addImaginary(*imEl, dofs);
		}
	}
}


template<class Scalar>
void
Domain::updateDampingMatrices(AllOps<Scalar> *ops, int *dofs, FullSquareMatrix *reEl,
                              FullSquareMatrix *imEl, double ss, int n)
{
	if(reEl) {
		*reEl /= pow(ss,n);
		if(ops && ops->C_deriv && ops->C_deriv[n-1]) ops->C_deriv[n-1]->add(*reEl, dofs);
		if(ops && ops->Cuc_deriv && ops->Cuc_deriv[n-1]) ops->Cuc_deriv[n-1]->add(*reEl, dofs);
	}
	if(imEl) {
		*imEl /= pow(ss,n);
		if(ops && ops->C_deriv && ops->C_deriv[n-1]) ops->C_deriv[n-1]->addImaginary(*imEl, dofs);
		if(ops && ops->Cuc_deriv && ops->Cuc_deriv[n-1]) ops->Cuc_deriv[n-1]->addImaginary(*imEl, dofs);
	}
}

//-------------------------------------------------------------------------------------
template<class Scalar>
int Domain::processDispTypeOutputs(OutputInfo &oinfo, Scalar (*glDisp)[11], int numNodes,
                                   int i, double time, double freq, int printFlag)  {

	int success = 0;
	double tag;
	switch (oinfo.type)  {

		case OutputInfo::Displacement:
			tag = time;
			success = 1;
		case OutputInfo::EigenPair:  {
			if (success == 0)
				tag = freq;
			if (oinfo.nodeNumber == -1) { // all nodes, or select group
				geoSource->outputNodeVectors(i, glDisp, numNodes, tag);
			}
			else { // one node
				int inode = (outFlag == 1) ? nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
				geoSource->outputNodeVectors(i, &(glDisp[inode]), 1, tag);
			}
			success = 1;
		}
			break;
		case OutputInfo::Disp6DOF:
			tag = time;
			success = 1;
		case OutputInfo::EigenPair6:  {
			if (success == 0)
				tag = freq;
			if (oinfo.nodeNumber == -1) { // all nodes, or select group
				geoSource->outputNodeVectors6(i, glDisp, numNodes, tag);
			}
			else { // one node
				int inode = (outFlag == 1) ? nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
				geoSource->outputNodeVectors6(i, &(glDisp[inode]), 1, tag);
			}
			success = 1;
		}
			break;
		default:
			break;
	}

	Scalar *globVal = 0;
	if (success == 0)  {
		globVal = new Scalar[numNodes];
		int dof = -1;

		switch (oinfo.type)  {

			case OutputInfo::DispX:
				if(dof==-1) dof = 0;
			case OutputInfo::DispY:
				if(dof==-1) dof = 1;
			case OutputInfo::DispZ:
				if(dof==-1) dof = 2;
			case OutputInfo::RotX:
				if(dof==-1) dof = 3;
			case OutputInfo::RotY:
				if(dof==-1) dof = 4;
			case OutputInfo::RotZ:
				if(dof==-1) dof = 5;

				for (int iNode=0; iNode<numNodes; ++iNode)
					globVal[iNode] = glDisp[iNode][dof];
				success = 1;
			case OutputInfo::DispMod:
				if (success == 0)  {
					for (int iNode = 0; iNode < numNodes; ++iNode)
						globVal[iNode] = ScalarTypes::sqrt(glDisp[iNode][0]*glDisp[iNode][0] +
						                                   glDisp[iNode][1]*glDisp[iNode][1] +
						                                   glDisp[iNode][2]*glDisp[iNode][2]);
					success = 1;
				}
			case OutputInfo::RotMod:
				if (success == 0)  {
					for (int iNode = 0; iNode < numNodes; ++iNode)
						globVal[iNode] = ScalarTypes::sqrt(glDisp[iNode][3]*glDisp[iNode][3] +
						                                   glDisp[iNode][4]*glDisp[iNode][4] +
						                                   glDisp[iNode][5]*glDisp[iNode][5]);
					success = 1;
				}
			case OutputInfo::TotMod:
				if (success == 0)  {
					for (int iNode = 0; iNode < numNodes; ++iNode)
						globVal[iNode] = ScalarTypes::sqrt(glDisp[iNode][0]*glDisp[iNode][0] +
						                                   glDisp[iNode][1]*glDisp[iNode][1] +
						                                   glDisp[iNode][2]*glDisp[iNode][2] +
						                                   glDisp[iNode][3]*glDisp[iNode][3] +
						                                   glDisp[iNode][4]*glDisp[iNode][4] +
						                                   glDisp[iNode][5]*glDisp[iNode][5]);
					success = 1;
				}
				if (oinfo.nodeNumber == -1)
					geoSource->outputNodeScalars(i, globVal, numNodes, time);
				else {
					int inode = (outFlag == 1) ? nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
					geoSource->outputNodeScalars(i, &(globVal[inode]), 1, time);
				}
				break;
			case OutputInfo::Temperature:
				if (dof==-1) { dof = 6; tag = time; }
			case OutputInfo::AcousticPressure:
				if (dof==-1) { dof = 7; tag = time; }
			case OutputInfo::EigenPressure:
			case OutputInfo::HelmholtzModes:
			case OutputInfo::Helmholtz:
				if (dof==-1) { dof = 7; tag = freq; }

				for (int iNode=0; iNode<numNodes; ++iNode)
					globVal[iNode] = glDisp[iNode][dof];

				success = 1;
			case OutputInfo::EigenSlosh:
				if (success == 0)  {
					if (dof == -1) dof = 10;
					tag = freq;
					for (int iNode=0; iNode<numNodes; ++iNode)
						globVal[iNode] = glDisp[iNode][dof];
					success = 1;
				}
				if(oinfo.nodeNumber == -1)
					geoSource->outputNodeScalars(i, globVal, numNodes, tag);
				else {
					int inode = (outFlag == 1) ? nodeTable[oinfo.nodeNumber]-1 : oinfo.nodeNumber;
					geoSource->outputNodeScalars(i, &(globVal[inode]), 1, tag);
				}
				break;
			default:
				break;
		}
		delete [] globVal;
	}

	return success;
}

template<class Scalar>
int Domain::processOutput(OutputInfo::Type &type, GenVector<Scalar> &d_n, Scalar *bcx, int i,
                          double time, double freq, int printFlag)  {

	int success = 1;
	switch(type) {
		case OutputInfo::SloshDisplacement:
			getSloshDispAll(d_n, bcx, i, freq);
			break;
		case OutputInfo::SloshDispX:
			getSloshDisp(d_n, bcx, i, SLDX, freq);
			break;
		case OutputInfo::SloshDispY:
			getSloshDisp(d_n, bcx, i, SLDY, freq);
			break;
		case OutputInfo::SloshDispZ:
			getSloshDisp(d_n, bcx, i, SLDZ, freq);
			break;
		case OutputInfo::YModulus:
			getElementAttr(i,YOUNG);
			break;
		case OutputInfo::MDensity:
			getElementAttr(i,MDENS);
			break;
		case OutputInfo::Thicknes:
			getElementAttr(i,THICK);
			break;
		case OutputInfo::StressXX:
			getStressStrain(d_n,bcx,i,SXX, time);
			break;
		case OutputInfo::StressYY:
			getStressStrain(d_n,bcx,i,SYY, time);
			break;
		case OutputInfo::StressZZ:
			getStressStrain(d_n,bcx,i,SZZ, time);
			break;
		case OutputInfo::StressXY:
			getStressStrain(d_n,bcx,i,SXY, time);
			break;
		case OutputInfo::StressYZ:
			getStressStrain(d_n,bcx,i,SYZ, time);
			break;
		case OutputInfo::StressXZ:
			getStressStrain(d_n,bcx,i,SXZ, time);
			break;
		case OutputInfo::StrainXX:
			getStressStrain(d_n,bcx,i,EXX, time);
			break;
		case OutputInfo::StrainYY:
			getStressStrain(d_n,bcx,i,EYY, time);
			break;
		case OutputInfo::StrainZZ:
			getStressStrain(d_n,bcx,i,EZZ, time);
			break;
		case OutputInfo::StrainXY:
			getStressStrain(d_n,bcx,i,EXY, time);
			break;
		case OutputInfo::StrainYZ:
			getStressStrain(d_n,bcx,i,EYZ, time);
			break;
		case OutputInfo::StrainXZ:
			getStressStrain(d_n,bcx,i,EXZ, time);
			break;
		case OutputInfo::StressVM:
			getStressStrain(d_n,bcx,i,VON, time);
			break;
		case OutputInfo::AggrStVM:
			getStressStrain(d_n,bcx,i,AGGREGATEDVON, time);
			break;
		case OutputInfo::Damage:
			getStressStrain(d_n,bcx,i,DAMAGE, time);
			break;
		case OutputInfo::EquivalentPlasticStrain:
			getStressStrain(d_n, bcx, i, EQPLSTRN, time);
			break;
		case OutputInfo::StressPR1:
			getPrincipalStress(d_n,bcx,i,PSTRESS1,time);
			break;
		case OutputInfo::StressPR2:
			getPrincipalStress(d_n,bcx,i,PSTRESS2,time);
			break;
		case OutputInfo::StressPR3:
			getPrincipalStress(d_n,bcx,i,PSTRESS3,time);
			break;
		case OutputInfo::StrainPR1:
			getPrincipalStress(d_n,bcx,i,PSTRAIN1,time);
			break;
		case OutputInfo::StrainPR2:
			getPrincipalStress(d_n,bcx,i,PSTRAIN2,time);
			break;
		case OutputInfo::StrainPR3:
			getPrincipalStress(d_n,bcx,i,PSTRAIN3,time);
			break;
		case OutputInfo::InXForce:
			getElementForces(d_n, bcx, i, INX, time);
			break;
		case OutputInfo::InYForce:
			getElementForces(d_n, bcx, i, INY, time);
			break;
		case OutputInfo::InZForce:
			getElementForces(d_n, bcx, i, INZ, time);
			break;
		case OutputInfo::AXMoment:
			getElementForces(d_n, bcx, i, AXM, time);
			break;
		case OutputInfo::AYMoment:
			getElementForces(d_n, bcx, i, AYM, time);
			break;
		case OutputInfo::AZMoment:
			getElementForces(d_n, bcx, i, AZM, time);
			break;
		case OutputInfo::StrainVM:
			getStressStrain(d_n,bcx,i,STRAINVON, time);
			break;
		case OutputInfo::HeatFlXX:
			getHeatFlux(d_n, bcx, i, HFLX);
			break;
		case OutputInfo::HeatFlXY:
			getHeatFlux(d_n, bcx, i, HFLY);
			break;
		case OutputInfo::HeatFlXZ:
			getHeatFlux(d_n, bcx, i, HFLZ);
			break;
		case OutputInfo::GrdTempX:
			getHeatFlux(d_n, bcx, i, GRTX);
			break;
		case OutputInfo::GrdTempY:
			getHeatFlux(d_n, bcx, i, GRTY);
			break;
		case OutputInfo::GrdTempZ:
			getHeatFlux(d_n, bcx, i, GRTZ);
			break;
		case OutputInfo::HeatFlX:
			getTrussHeatFlux(d_n, bcx, i, HFLX);
			break;
		case OutputInfo::GrdTemp:
			getTrussHeatFlux(d_n, bcx, i, GRTX);
			break;
		default:
			success = 0;
			break;
	}
	return success;
}

//-------------------------------------------------------------------------------------
template <>
void Domain::sensitivityPostProcessing(AllSensitivities<DComplex> &allSens, GenVector<DComplex> *sol,
                                       DComplex *bcx, GeomState *geomState, GeomState *refState, Corotator **allCorot) {
	filePrint(stderr, " ... WARNING : Domain::sensitivityPostProcessing is not implemented\n");
}

//-------------------------------------------------------------------------------------
template <>
void Domain::sensitivityPostProcessing(AllSensitivities<double> &allSens, GenVector<double> *sol, double *bcx,
                                       GeomState *geomState, GeomState *refState, Corotator **allCorot) {
#ifdef USE_EIGEN3
	OutputInfo *oinfo = geoSource->getOutputInfo();
	int numOutInfo = geoSource->getNumOutInfo();
	if(firstOutput) geoSource->openOutputFiles();
	for(int i = 0; i < numOutInfo; ++i)  {
		if(oinfo[i].type == OutputInfo::WeigThic) geoSource->outputSensitivityScalars(i, allSens.weightWRTthick, allSens.weight);
		if(oinfo[i].type == OutputInfo::WeigShap) geoSource->outputSensitivityScalars(i, allSens.weightWRTshape, allSens.weight);
		if(oinfo[i].type == OutputInfo::AGstShap) geoSource->outputSensitivityScalars(i, allSens.aggregatedVonMisesWRTshape,
		                                                                              *aggregatedStress,
		                                                                              allSens.dwrAggregatedStressVM);
		if(oinfo[i].type == OutputInfo::AGstThic) geoSource->outputSensitivityScalars(i, allSens.aggregatedVonMisesWRTthick,
		                                                                              *aggregatedStress,
		                                                                              allSens.dwrAggregatedStressVM);
		if(oinfo[i].type == OutputInfo::VMstThic) {
			if(solInfo().sensitivityMethod == SolverInfo::Direct)
				geoSource->outputSensitivityVectors(i, allSens.vonMisesWRTthick, 0.0, allSens.dwrStressVM);
			else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
				Vector stress(numNodes(),0.0);
				Vector stressWeight(numNodes(),0.0);
				if(sinfo.isNonLin()) {
					computeNormalizedNLVonMisesStress(*geomState,refState,allCorot,oinfo[i].surface,stress,stressWeight,false);
				} else {
					computeNormalizedVonMisesStress(*sol,bcx,oinfo[i].surface,stress,false);
				}
				int numThicknessGroup = getNumThicknessGroups();
				geoSource->outputSensitivityAdjointStressVectors(i, allSens.vonMisesWRTthick, stress.data(), 0, numThicknessGroup,
				                                                 stressNodes, allSens.dwrStressVM);
			}
		}
		if(oinfo[i].type == OutputInfo::VMstShap) {
			if(solInfo().sensitivityMethod == SolverInfo::Direct)
				geoSource->outputSensitivityVectors(i, allSens.vonMisesWRTshape, 0.0, allSens.dwrStressVM);
			else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
				Vector stress(numNodes(),0.0);
				Vector stressWeight(numNodes(),0.0);
				if(sinfo.isNonLin()) {
					computeNormalizedNLVonMisesStress(*geomState,refState,allCorot,oinfo[i].surface,stress,stressWeight,false);
				} else {
					computeNormalizedVonMisesStress(*sol,bcx,oinfo[i].surface,stress,false);
				}
				int numShapeVars = getNumShapeVars();
				geoSource->outputSensitivityAdjointStressVectors(i, allSens.vonMisesWRTshape, stress.data(), 0, numShapeVars,
				                                                 stressNodes, allSens.dwrStressVM);
			}
		}
		if(oinfo[i].type == OutputInfo::VMstMach) geoSource->outputSensitivityVectors(i, allSens.vonMisesWRTmach);
		if(oinfo[i].type == OutputInfo::VMstAlpha) geoSource->outputSensitivityVectors(i, allSens.vonMisesWRTalpha);
		if(oinfo[i].type == OutputInfo::VMstBeta) geoSource->outputSensitivityVectors(i, allSens.vonMisesWRTbeta);
		if(oinfo[i].type == OutputInfo::Statevector) {
			using namespace Rom;
			VecNodeDofConversion<6> vecNodeDof6Conversion(*getCDSA());
			if(allSens.lambdaDisp) {
				std::string fileName = std::string(oinfo[i].filename) + ".DispAdjoint";
				BasisOutputStream<6> adjointSnapFile(BasisFileId(fileName, BasisId::STATE, BasisId::SNAPSHOTS), vecNodeDof6Conversion,
				                                     (geoSource->getCheckFileInfo()->lastRestartFile != 0));
				for(int idof=0; idof<numTotalDispDofs; ++idof)
					adjointSnapFile << *allSens.lambdaDisp[idof];
			}
			if(allSens.lambdaStressVM) {
				std::string fileName = std::string(oinfo[i].filename) + ".StressVMAdjoint";
				BasisOutputStream<6> adjointSnapFile(BasisFileId(fileName, BasisId::STATE, BasisId::SNAPSHOTS), vecNodeDof6Conversion,
				                                     (geoSource->getCheckFileInfo()->lastRestartFile != 0));
				for(int inode=0; inode<numStressNodes; ++inode)
					adjointSnapFile << *allSens.lambdaStressVM[inode];
			}
			if(allSens.lambdaAggregatedStressVM) {
				std::string fileName = std::string(oinfo[i].filename) + ".AggregatedStressVMAdjoint";
				BasisOutputStream<6> adjointSnapFile(BasisFileId(fileName, BasisId::STATE, BasisId::SNAPSHOTS), vecNodeDof6Conversion,
				                                     (geoSource->getCheckFileInfo()->lastRestartFile != 0));
				adjointSnapFile << *allSens.lambdaAggregatedStressVM;
			}
		}
		if(oinfo[i].type == OutputInfo::DispThic) {
			int numThicknessGroup = getNumThicknessGroups();
			allSens.gdispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroup];
			if(solInfo().sensitivityMethod == SolverInfo::Direct) {
				for(int iparam=0; iparam<numThicknessGroup; ++iparam) {
					allSens.gdispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numnodes, 6);
					mergeDistributedDispSensitivity<double>(allSens.gdispWRTthick[iparam], allSens.dispWRTthick[iparam]);
				}
				geoSource->outputSensitivityDispVectors(i, allSens.gdispWRTthick, 0, numThicknessGroup, numnodes);
			} else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
				double *disp = new double[numTotalDispDofs];
				for(int iparam=0; iparam<numThicknessGroup; ++iparam) {
					allSens.gdispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numTotalDispDofs, 1);
					mergeAdjointDistributedDispSensitivity<double>(allSens.gdispWRTthick[iparam], allSens.dispWRTthick[iparam], disp, sol);
				}
				geoSource->outputSensitivityAdjointDispVectors(i, allSens.gdispWRTthick, disp, 0, numThicknessGroup, dispNodes,
				                                               allSens.dwrDisp);
			}
		}
		if(oinfo[i].type == OutputInfo::DispShap) {
			int numShapeVars = getNumShapeVars();
			allSens.gdispWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
			if(solInfo().sensitivityMethod == SolverInfo::Direct) {
				for(int iparam=0; iparam<numShapeVars; ++iparam) {
					allSens.gdispWRTshape[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numnodes, 6);
					mergeDistributedDispSensitivity<double>(allSens.gdispWRTshape[iparam], allSens.dispWRTshape[iparam]);
				}
				geoSource->outputSensitivityDispVectors(i, allSens.gdispWRTshape, 0, numShapeVars, numnodes);
			} else if(solInfo().sensitivityMethod == SolverInfo::Adjoint) {
				double *disp = new double[numTotalDispDofs];
				for(int iparam=0; iparam<numShapeVars; ++iparam) {
					allSens.gdispWRTshape[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numTotalDispDofs, 1);
					mergeAdjointDistributedDispSensitivity<double>(allSens.gdispWRTshape[iparam], allSens.dispWRTshape[iparam], disp, sol);
				}
				geoSource->outputSensitivityAdjointDispVectors(i, allSens.gdispWRTshape, disp, 0, numShapeVars, dispNodes, allSens.dwrDisp);
			}
		}
		if(oinfo[i].type == OutputInfo::DispMach) {
			allSens.gdispWRTmach = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numnodes, 6);
			mergeDistributedDispSensitivity<double>(allSens.gdispWRTmach, allSens.dispWRTmach);
			geoSource->outputSensitivityDispVectors(i, allSens.gdispWRTmach, 0, numnodes);
		}
		if(oinfo[i].type == OutputInfo::DispAlph) {
			allSens.gdispWRTalpha = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numnodes, 6);
			mergeDistributedDispSensitivity<double>(allSens.gdispWRTalpha, allSens.dispWRTalpha);
			geoSource->outputSensitivityDispVectors(i, allSens.gdispWRTalpha, 0, numnodes);
		}
		if(oinfo[i].type == OutputInfo::DispBeta) {
			allSens.gdispWRTmach = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numnodes, 6);
			mergeDistributedDispSensitivity<double>(allSens.gdispWRTbeta, allSens.dispWRTbeta);
			geoSource->outputSensitivityDispVectors(i, allSens.gdispWRTbeta, 0, numnodes);
		}
	}
	firstOutput = false;
#endif
}

//-------------------------------------------------------------------------------------
#ifdef USE_EIGEN3
template<class Scalar>
void Domain::mergeDistributedDispSensitivity(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *gdispSen,
                                             Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *dispSen)
{
	int inode, nodeI;
	int realNode = -1;

	for(inode = 0; inode < numnodes; ++inode) {

		if(nodes[inode] == 0) continue;
		realNode++;
		nodeI = (outFlag) ? realNode : inode;

		int xLoc  = c_dsa->locate(inode, DofSet::Xdisp);
		int xLoc1 =   dsa->locate(inode, DofSet::Xdisp);

		if (xLoc >= 0)
			(*gdispSen)(nodeI,0) = (*dispSen)(xLoc,0);          // free
		else
			(*gdispSen)(nodeI,0) = 0.0;

		int yLoc  = c_dsa->locate(inode, DofSet::Ydisp);
		int yLoc1 =   dsa->locate(inode, DofSet::Ydisp);

		if (yLoc >= 0)
			(*gdispSen)(nodeI,1) = (*dispSen)(yLoc,0);
		else
			(*gdispSen)(nodeI,1) = 0.0;

		int zLoc  = c_dsa->locate(inode, DofSet::Zdisp);
		int zLoc1 =   dsa->locate(inode, DofSet::Zdisp);

		if (zLoc >= 0)
			(*gdispSen)(nodeI,2) = (*dispSen)(zLoc,0);
		else
			(*gdispSen)(nodeI,2) = 0.0;

		int xRot  = c_dsa->locate(inode, DofSet::Xrot);
		int xRot1 =   dsa->locate(inode, DofSet::Xrot);

		if (xRot >= 0)
			(*gdispSen)(nodeI,3) = (*dispSen)(xRot,0);
		else
			(*gdispSen)(nodeI,3) = 0.0;

		int yRot  = c_dsa->locate(inode, DofSet::Yrot);
		int yRot1 =   dsa->locate(inode, DofSet::Yrot);

		if (yRot >= 0)
			(*gdispSen)(nodeI,4) = (*dispSen)(yRot,0);
		else
			(*gdispSen)(nodeI,4) = 0.0;

		int zRot  = c_dsa->locate(inode, DofSet::Zrot);
		int zRot1 =   dsa->locate(inode, DofSet::Zrot);

		if (zRot >= 0)
			(*gdispSen)(nodeI,5) = (*dispSen)(zRot,0);
		else
			(*gdispSen)(nodeI,5) = 0.0;
/*
    int xTemp  = c_dsa->locate(inode, DofSet::Temp);
    int xTemp1 =   dsa->locate(inode, DofSet::Temp);

    if (xTemp >= 0)
      xyz[nodeI][6] = u[xTemp];
    else
      xyz[nodeI][6] = 0.0;

    int xHelm  = c_dsa->locate(inode, DofSet::Helm);
    int xHelm1 =   dsa->locate(inode, DofSet::Helm);

    if (xHelm >= 0)
      xyz[nodeI][7] = u[xHelm];
    else
      xyz[nodeI][7] = 0.0;

    int xPot  = c_dsa->locate(inode, DofSet::Potential);
    int xPot1 =   dsa->locate(inode, DofSet::Potential);

    if (xPot >= 0)
      xyz[nodeI][10] = u[xPot];
    else
      xyz[nodeI][10] = 0.0;

    // transform displacements and rotations (if present) from DOF_FRM to basic coordinates
    // and keep a copy of the original in xyz_loc
    if(!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
      && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
      if(xyz_loc) for(int j=0; j<11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
      bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
      transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
    }
*/
	}

}
#endif

//-------------------------------------------------------------------------------------
#ifdef USE_EIGEN3
template<class Scalar>
void Domain::mergeAdjointDistributedDispSensitivity(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *gdispSen,
                                                    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *dispSen,
                                                    Scalar *disp, GenVector<Scalar> *sol)
{
	int inode;
	int realNode = -1;

	int dispDofIndex = 0;
	for(inode = 0; inode < numDispNodes; ++inode) {
		int node = dispNodes[inode].nodeID, loc;
		int numDispDofs = dispNodes[inode].numdofs;
		for(int idof=0; idof<numDispDofs; ++idof) {
			int dof = dispNodes[inode].dofs[idof];
			loc = returnLocalDofNum(node, dof);
			if (loc >= 0) {
				(*gdispSen)(dispDofIndex,0) = (*dispSen)(dispDofIndex,0);
				if(sol) disp[dispDofIndex] = (*sol)[loc];
				else disp[dispDofIndex] = 0.0;
			} else {
				(*gdispSen)(dispDofIndex,0) = 0.0;
				disp[dispDofIndex] = 0.0;
			}
			dispDofIndex++;
		}
	}

}
#endif

//-------------------------------------------------------------------------------------
// Templated Post-processing for direct solver statics, frequency response, helmholtz and eigen
template<class Scalar>
void Domain::postProcessing(GenVector<Scalar> &sol, Scalar *bcx, GenVector<Scalar> &force,
                            int ndflag, int index, double time, double eigV,
                            GenSparseMatrix<Scalar> *kuc, GenSparseMatrix<Scalar> *kcc) {

	if(outFlag && !nodeTable) makeNodeTable(outFlag);
	int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes
	double freq;
	if (domain->probType() == SolverInfo::Modal) freq = eigV;
	else freq = domain->getFrequencyOrWavenumber();

	if (geoSource->isShifted() || domain->probType() == SolverInfo::Modal) time = freq;
	if (domain->solInfo().loadcases.size() > 0 && !domain->solInfo().doFreqSweep) time = domain->solInfo().loadcases.front();

	Scalar *globVal = 0;
	int numOutInfo = geoSource->getNumOutInfo();
	OutputInfo *oinfo = geoSource->getOutputInfo();

	if(numOutInfo && (firstOutput || domain->solInfo().loadcases.size() > 0) && ndflag<=1)
		filePrint(stderr," ... Postprocessing                 ...\n");

	// organize displacements
	int numNodeLim = std::max(numNodes,numnodes);

	Scalar (*xyz)[11] = new Scalar[numNodeLim][11];
	Scalar (*xyz_loc)[11] = (domain->solInfo().basicDofCoords) ? 0 : new Scalar[numNodeLim][11];
	int i;
	for(i = 0; i < numNodeLim; ++i)
		for (int j = 0 ; j < 11 ; j++)
			xyz[i][j] = 0.0;
	mergeDistributedDisp<Scalar>(xyz, sol.data(), bcx, xyz_loc);
	int numNodesOut = (outFlag) ? exactNumNodes : numNodes;

	// Open files and write file headers in first time step
	if(firstOutput) geoSource->openOutputFiles();

	int dof;
	int iNode;
	for(i = 0; i < numOutInfo; ++i)  {
		if(oinfo[i].sentype > 0) continue;
		if(oinfo[i].ndtype != ndflag) continue;
		if(ndflag !=0 && oinfo[i].type != OutputInfo::Disp6DOF && oinfo[i].type !=  OutputInfo::Displacement) continue;
		// if non-deterministic and NOT displacement
		if(oinfo[i].interval == 1
		   || oinfo[i].type == OutputInfo::Farfield
		   || oinfo[i].type == OutputInfo::Kirchhoff) {
			dof = -1;
			int success;
			if(oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
				success = processDispTypeOutputs(oinfo[i], xyz, numNodesOut, i, time, freq);
			else
				success = processDispTypeOutputs(oinfo[i], xyz_loc, numNodesOut, i, time, freq);
			if (success) continue;
			success = processOutput(oinfo[i].type, sol, bcx, i, time, freq, 0);
			if (success) continue;
			success = 1;
			switch(oinfo[i].type)  {
				case OutputInfo::Energies: {
					double Wext = ScalarTypes::Real(force*sol); // Wext = external energy
					double Wela =   0.5 * Wext;                 // Wela = elastic energy
					geoSource->outputEnergies(i, freq, Wext, 0.0, Wela, 0.0, 0.0, 0.0);
				} break;
				case OutputInfo::Farfield: case OutputInfo::Kirchhoff:
					outputFFP(sol, i);
					break;
				case OutputInfo::ModeError:
					break;  // This is handled in Problems.d/DynamDescr.C
					// The following 3 cases are not officially supported in manual
				case OutputInfo::ElemToNode:
					if(elemToNode)
						elemToNode->print(oinfo[i].filptr, oinfo[i].nodeNumber);
					break;
				case OutputInfo::NodeToElem:
					if(nodeToElem)
						nodeToElem->print(oinfo[i].filptr, oinfo[i].nodeNumber);
					break;
				case OutputInfo::NodeToNode:
					if(nodeToNode)
						nodeToNode->print(oinfo[i].filptr, oinfo[i].nodeNumber);
					break;
				case OutputInfo::Reactions: {
					GenVector<Scalar> fc(numDirichlet+numComplexDirichlet);
					computeReactionForce(fc, sol, kuc, kcc);
					Scalar (*rxyz)[3] = new Scalar[numNodeLim][3];
					DofSet dofs[3] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp };
					for(int inode = 0, realNode = -1; inode < nodes.size(); ++inode) {
						if(nodes[inode] == 0) continue;
						realNode++;
						int nodeI = (outFlag) ? realNode : inode;
						for(int k = 0; k < 3; ++k) {
							int dof =   dsa->locate(inode, dofs[k].list());
							int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
							rxyz[nodeI][k] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
						}
						if(oinfo[i].oframe == OutputInfo::Global && !domain->solInfo().basicDofCoords) {
							transformVectorInv(&rxyz[nodeI][0],inode,false);
						}
					}
					geoSource->outputNodeVectors(i, rxyz, numNodesOut, time);
					delete [] rxyz;
				} break;
				case OutputInfo::Reactions6: {
					GenVector<Scalar> fc(numDirichlet+numComplexDirichlet);
					computeReactionForce(fc, sol, kuc, kcc);
					Scalar (*rxyz)[6] = new Scalar[numNodeLim][6];
					DofSet dofs[6] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
					                   DofSet::Xrot, DofSet::Yrot, DofSet::Zrot };
					for(int inode = 0, realNode = -1; inode < nodes.size(); ++inode) {
						if(nodes[inode] == 0) continue;
						realNode++;
						int nodeI = (outFlag) ? realNode : inode;
						for(int k = 0; k < 6; ++k) {
							int dof =   dsa->locate(inode, dofs[k].list());
							int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
							rxyz[nodeI][k] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
						}
						if(oinfo[i].oframe == OutputInfo::Global && !domain->solInfo().basicDofCoords) {
							transformVectorInv(&rxyz[nodeI][0],inode,true);
						}
					}
					geoSource->outputNodeVectors6(i, rxyz, numNodesOut, time);
					delete [] rxyz;
				} break;
				case OutputInfo::HeatReactions: {
					GenVector<Scalar> fc(numDirichlet+numComplexDirichlet);
					computeReactionForce(fc, sol, kuc, kcc);
					Scalar *rxyz = new Scalar[numNodeLim];
					DofSet dofs[1] = { DofSet::Temp };
					for(int inode = 0, realNode = -1; inode < nodes.size(); ++inode) {
						if(nodes[inode] == 0) continue;
						realNode++;
						int nodeI = (outFlag) ? realNode : inode;
						for(int k = 0; k < 1; ++k) {
							int dof =   dsa->locate(inode, dofs[k].list());
							int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
							rxyz[nodeI] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
						}
					}
					geoSource->outputNodeScalars(i, rxyz, numNodesOut, time);
					delete [] rxyz;
				} break;
				case OutputInfo::Statevector:
					break;
				case OutputInfo::Velocity:
				case OutputInfo::Acceleration:
				case OutputInfo::Velocity6:
				case OutputInfo::Accel6:
					break;
				default:
					success = 0;
					break;
			}
			if (success == 0)
				fprintf(stderr, " *** WARNING: output %d is not supported \n", i);
		}
		if(globVal)
		{ delete [] globVal; globVal = 0; }
	}

	// --- Print Problem statistics to the screen -------------------------------
	if(firstOutput) {
		if (!domain->solInfo().doEigSweep) {

			// ... CALCULATE STRUCTURE MASS IF REQUESTED
			if(sinfo.massFlag)  {
				double mass = computeStructureMass();
				filePrint(stderr," ... Total System Mass = %10.4f ...\n", mass);
				filePrint(stderr," --------------------------------------\n");
			}
		}

		firstOutput = false;
	}

	if (xyz) delete [] xyz;
	if (xyz_loc) delete [] xyz_loc;
}

template<class Scalar>
void
Domain::addConstantForceSensitivity(GenVector<Scalar>& cnst_fSen, GenSparseMatrix<Scalar>* kuc)
{
	if(domain->gravityFlag()) addGravityForceSensitivity(cnst_fSen);
}

template <class Scalar>
void
Domain::computeConstantForce(GenVector<Scalar>& cnst_f, GenSparseMatrix<Scalar>* kuc)
{
	// This is called for linear and nonlinear dynamics
	// cnst_f is independent of t
	double loadFactor; // load amplification factor used for combination load case

	if(!dynamic_cast<GenSubDomain<Scalar>*>(this) && !sommerChecked) checkSommerTypeBC(this);

	cnst_f.zero();

	// ... COMPUTE FORCE FROM GRAVITY
	if(domain->gravityFlag()) addGravityForce(cnst_f);

	// ... COMPUTE FORCE FROM DISCRETE NEUMANN BOUNDARY CONDITIONS
	// note #1 when MFTT is present then FORCES contribution is not constant
	// note #2 when HFTT is present the FLUX contribution is not constant
	// note #3 see getStiffAndForce/getInternalForce for treatment of non-axial forces and all nodal moments in nonlinear analyses
	for(int i = 0; i < numNeuman; ++i) {
		if(sinfo.isNonLinExtF() && nbc[i].type == BCond::Forces && !(nbc[i].mtype == BCond::Axial && nbc[i].dofnum < 3)) continue;
		int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
		if(dof < 0) continue;
		switch(nbc[i].type) {
			case(BCond::Forces) : {
				if(!sinfo.isDynam() || !domain->getMFTT(nbc[i].loadsetid)) {
					double loadFactor = domain->getLoadFactor(nbc[i].loadsetid);
					cnst_f[dof] += loadFactor*nbc[i].val;
				}
			} break;
			case(BCond::Flux) : {
				if(!sinfo.isDynam() || !domain->getHFTT(nbc[i].loadsetid)) {
					double loadFactor = domain->getLoadFactor(nbc[i].loadsetid);
					cnst_f[dof] += loadFactor*nbc[i].val;
				}
			} break;
			case(BCond::Convection) : {
				double loadFactor = domain->getLoadFactor(nbc[i].loadsetid);
				cnst_f[dof] += loadFactor*nbc[i].val;
			} break;
			case(BCond::Actuators) : case(BCond::Usdf) : break;
			default : cnst_f[dof] += nbc[i].val;
		}
	}

	// ... COMPUTE FORCE FROM ACOUSTIC DISTRIBUTED NEUMANN BOUNDARY CONDITIONS
	if(sinfo.ATDDNBVal != 0.0 && !domain->getDefaultMFTT()) addAtddnbForce(cnst_f, 0);

	// ... COMPUTE FORCE FROM ACOUSTIC ROBIN BOUNDARY CONDITIONS
	if(sinfo.ATDROBalpha != 0.0 && !domain->getDefaultMFTT()) addAtdrobForce(cnst_f, 0);

	// ... COMPUTE FORCE FROM PRESSURE
	// note #1: even when MFTTs/CONWEP are present this term may now be constant
	// note #2: for NONLINEAR problems this term is not constant (see getStiffAndForce/getInternalForce)
	if(!sinfo.isNonLinExtF()) addPressureForce(cnst_f, 0);

	// ... ADD RHS FROM LMPCs for linear statics
	if((!sinfo.isNonLin() || sinfo.getNLInfo().linearelastic) && !sinfo.isDynam()) addMpcRhs(cnst_f);

	// ... COMPUTE FORCE FROM TEMPERATURES
	// note #1: for THERMOE problems TEMPERATURES are ignored
	// note #2: for NONLINEAR problems this term is not constant (see getStiffAndForce/getInternalForce)
	if(sinfo.thermalLoadFlag && !(sinfo.thermoeFlag >= 0) && !sinfo.isNonLinExtF()) addThermalForce(cnst_f);

	// ... COMPUTE FORCE FROM NON-HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
	// note #1: when USDD is present this is term is not constant (see computeExtForce)
	// note #2  for nonlinear this term is not constant (see getStiffAndForce/getInternalForce)
	if(numDirichlet && !(claw && claw->numUserDisp) && (!sinfo.isNonLin() || (sinfo.galerkinPodRom && sinfo.getNLInfo().linearelastic)) && kuc) {
		GenVector<Scalar> Vc(numDirichlet, 0.0);
		// construct the non-homogeneous dirichlet bc vector
		for(int i = 0; i < numDirichlet; ++i) {
			int dof = dsa->locate(dbc[i].nnum, (1 << dbc[i].dofnum));
			if(dof >= 0) {
				int dof2 = c_dsa->invRCN(dof);
				if(dof2 >= 0) Vc[dof2] = dbc[i].val;
			}
		}
		// compute the non-homogeneous force
		kuc->multSubtract(Vc, cnst_f);
	}
}

template <class Scalar>
void
Domain::computeExtForce(GenVector<Scalar>& f, double t, GenSparseMatrix<Scalar>* kuc,
                        ControlInterface *userSupFunc, GenSparseMatrix<Scalar>* cuc,
                        double tm, GenSparseMatrix<Scalar> *muc)
{
	f.zero();

	// ... COMPUTE FORCE FROM DISCRETE NEUMANN BOUNDARY CONDITIONS
	// note #1 when MFTT is not assigned FORCES contribution is constant (see computeConstantForce)
	// note #2 when HFTT is not assigned FLUX contribution is constant (see computeConstantForce)
	// note #3 see Domain::getStiffAndForce for treatment of non-axial forces and all nodal moments in nonlinear analyses
	if(numNeuman && (domain->getNumMFTT() || domain->getNumHFTT() || (claw && (claw->numUserForce || claw->numActuator)))) {
		for(int i = 0; i < numNeuman; ++i) {
			if(sinfo.isNonLinExtF() && (nbc[i].type == BCond::Forces || nbc[i].type == BCond::Usdf
			                            || nbc[i].type == BCond::Actuators) && !((nbc[i].mtype == BCond::Axial && nbc[i].dofnum < 3) || nbc[i].dofnum == 6)) continue;
			int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
			if(dof < 0) continue;
			switch(nbc[i].type) {
				case(BCond::Forces) : if(MFTTData *mftt = domain->getMFTT(nbc[i].loadsetid)) f[dof] += mftt->getVal(t)*nbc[i].val; break;
				case(BCond::Flux)   : if(MFTTData *hftt = domain->getHFTT(nbc[i].loadsetid)) f[dof] += hftt->getVal(t)*nbc[i].val; break;
				case(BCond::Actuators) : case(BCond::Usdf) : f[dof] += nbc[i].val; break;
				default : /* all other cases are constant */ ;
			}
		}
	}

	// COMPUTE FORCE FROM ACOUSTIC DISTRIBUTED NEUMANN BOUNDARY CONDITIONS
	// note #1: when one or more MFTTs are present this term may not be constant
	if(sinfo.ATDDNBVal != 0.0 && domain->getDefaultMFTT()) addAtddnbForce(f, 1, t);

	// COMPUTE FORCE FROM ACOUSTIC ROBIN BOUNDARY CONDITIONS
	// note #1: when one or more MFTTs are present this term may not be constant
	if(sinfo.ATDROBalpha != 0.0 && domain->getDefaultMFTT()) addAtdrobForce(f, 1, t);

	// COMPUTE FORCE FROM PRESSURE
	// note #1: when MFTT/CONWEP are present this term may not be constant
	// note #2: for NONLINEAR problems this term is follower (see getStiffAndForce/getInternalForce)
	if((domain->getNumMFTT() > 0 || sinfo.ConwepOnOff) && !sinfo.isNonLinExtF()) addPressureForce(f, 1, t);

	// ... ADD RHS FROM LMPCs for linear dynamics
	if((!sinfo.isNonLin() || sinfo.getNLInfo().linearelastic) && sinfo.isDynam()) addMpcRhs(f, t);

	// COMPUTE FORCE FROM THERMOE
	// note #1: for NONLINEAR problems this term is follower (see getStiffAndForce/getInternalForce)
	if(sinfo.thermoeFlag >= 0 && !sinfo.isNonLinExtF()) addThermalForce(f);

	// COMPUTE FORCE FROM NON-HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
	// note #1: when USDD is not present this term is constant (see computeConstantForce)
	// note #2: for nonlinear the contribution due to Kuc is follower (see getStiffAndForce/getInternalForce)
	// note #3: for linear and nonlinear dynamics the contribution due to Cuc and Muc is now included
	if(numDirichlet && (claw && claw->numUserDisp)) {
		GenVector<Scalar> Vc(numDirichlet, 0.0);
		// construct the non-homogeneous dirichlet bc vector
		for(int i = 0; i < numDirichlet; ++i) {
			int dof = dsa->locate(dbc[i].nnum, (1 << dbc[i].dofnum));
			if(dof < 0) continue;
			int dof2 = c_dsa->invRCN(dof);
			if(dof2 >= 0) Vc[dof2] = dbc[i].val;
		}

		// compute the non-homogeneous force due to Kuc
		if((!sinfo.isNonLin() || (sinfo.galerkinPodRom && sinfo.getNLInfo().linearelastic)) && kuc) kuc->multSubtract(Vc, f);

		if(sinfo.isDynam() && userSupFunc && claw && claw->numUserDisp > 0) {

			GenSubDomain<Scalar> *subCast = dynamic_cast<GenSubDomain<Scalar>*>(this);

			int glNumUserDisp = domain->getClaw()->numUserDisp;
			double *userDefineDisp = (double *) dbg_alloca( sizeof(double)*glNumUserDisp );
			double *userDefineVel  = (double *) dbg_alloca( sizeof(double)*glNumUserDisp );
			double *userDefineAcc  = (double *) dbg_alloca( sizeof(double)*glNumUserDisp );

			for(int i = 0; i < glNumUserDisp; ++i) {
				userDefineVel[i] = 0;
				userDefineAcc[i] = 0;
			}
			userSupFunc->usd_disp( tm, userDefineDisp, userDefineVel, userDefineAcc );

			Vc.zero();
			for(int i=0; i<claw->numUserDisp; ++i) {
				int dof = getDSA()->locate( claw->userDisp[i].nnum,
				                            1 << claw->userDisp[i].dofnum );
				if(dof < 0) continue;
				int dof1 = getCDSA()->invRCN( dof );
				if(dof1 >= 0) {
					int j = (subCast) ? subCast->getUserDispDataMap()[i] : i;
					Vc[dof1] = userDefineAcc[j];
				}
			}

			if(muc) muc->multSubtract(Vc, f); // fu -= Muc * a_c^{n+1-alpha_m}

			if(cuc) {

				for(int i = 0; i < glNumUserDisp; ++i) {
					userDefineVel[i] = 0;
					userDefineAcc[i] = 0;
				}
				userSupFunc->usd_disp( t, userDefineDisp, userDefineVel, userDefineAcc );

				Vc.zero();
				for(int i=0; i<claw->numUserDisp; ++i) {
					int dof = getDSA()->locate( claw->userDisp[i].nnum,
					                            1 << claw->userDisp[i].dofnum );
					if(dof < 0) continue;
					int dof1 = getCDSA()->invRCN( dof );
					if(dof1 >= 0) {
						int j = (subCast) ? subCast->getUserDispDataMap()[i] : i;
						Vc[dof1] = userDefineVel[j];
					}
				}

				cuc->multSubtract(Vc, f); // fu -= Cuc * v_c^{n+1-alpha_f}
			}
		}
	}
}

template <class Scalar>
void
Domain::computeExtForce4(GenVector<Scalar>& f, const GenVector<Scalar>& constantForce,
                         double t, GenSparseMatrix<Scalar>* kuc, ControlInterface *userSupFunc,
                         GenSparseMatrix<Scalar>* cuc, double tm, GenSparseMatrix<Scalar> *muc)
{
	// This is called for linear and nonlinear dynamics
	// doesn't include follower forces

	computeExtForce(f, t, kuc, userSupFunc, cuc, tm, muc);

	// ADD CONSTANT FORCE
	f += constantForce;
}
#endif
