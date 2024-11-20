#include <algorithm>
#include <iostream>
#include <Driver.d/Domain.h>
#include <Paral.d/MDDynam.h>
#include <Driver.d/Dynam.h>
#include <Paral.d/MDOp.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/Vector.h>
#include <Math.d/CuCSparse.h>
#include <Timers.d/GetTime.h>
#include <Control.d/ControlInterface.h>
#include <Threads.d/PHelper.h>
#include <Paral.d/SubDOp.h>

#ifdef DISTRIBUTED

#include <Dist.d/DistDom.h>

#endif

#include <Driver.d/DecDomain.h>
#include <Hetero.d/DistFlExchange.h>
#include <Feti.d/Feti.h>
#include <Utils.d/ModeData.h>
#include <Driver.d/SysState.h>
#include <Solvers.d/MultiDomainRbm.h>
#include <Corotational.d/DistrGeomState.h>

extern ModeData modeData;

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **_sd,
                             DistrVector *_v1, DistrVector *_v2, double _c1,
                             SubDOp *_Kuc, ControlInterface *_userSupFunc,
                             SubDOp *_Cuc, double _c2, SubDOp *_Muc) {
	f = _f;
	sd = _sd;
	v1 = _v1;
	v2 = _v2;
	c1 = _c1;
	Kuc = _Kuc;
	userSupFunc = _userSupFunc;
	Cuc = _Cuc;
	c2 = _c2;
	Muc = _Muc;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **_sd,
                             DistrVector *_v1, SubDOp *_Kuc) {
	v1 = _v1;
	f = _f;
	sd = _sd;
	Kuc = _Kuc;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **_sd) {
	f = _f;
	sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **_sd,
                             DistrVector *_v1, DistrVector *_v2, DistrVector *_v3) {
	v1 = _v1;
	v2 = _v2;
	v3 = _v3;
	f = _f;
	sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **_sd,
                             DistrVector *_v1, DistrVector *_v2, DistrVector *_v3, DistrVector *_v4) {
	v1 = _v1;
	v2 = _v2;
	v3 = _v3;
	v4 = _v4;
	f = _f;
	sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **_sd,
                             DistrVector *_v1, DistrVector *_v2) {
	v1 = _v1;
	v2 = _v2;
	f = _f;
	sd = _sd;
}


void
MultiDomainOp::runFor(int isub) {
	(this->*f)(isub);
}

void
MultiDomainOp::computeExtForce(int isub) {
	// Get the pointer to the part of the vector f corresponding to subdomain isub
	StackVector localf(v1->subData(isub), v1->subLen(isub));

	// Get the pointer to the part of the vector cnst_f corresponding to subdomain isub
	StackVector localg(v2->subData(isub), v2->subLen(isub));

	SparseMatrix *localKuc = (Kuc) ? (*Kuc)[isub] : 0;
	SparseMatrix *localCuc = (Cuc) ? (*Cuc)[isub] : 0;
	SparseMatrix *localMuc = (Muc) ? (*Muc)[isub] : 0;
	sd[isub]->computeExtForce4(localf, localg, c1, localKuc, userSupFunc, localCuc, c2, localMuc);
}

void
MultiDomainOp::getConstForce(int isub) {
	// Get the pointer to the part of the vector f correspoding to subdomain sNum
	StackVector f(v1->subData(isub), v1->subLen(isub));
	sd[isub]->computeConstantForce(f, (Kuc) ? (*Kuc)[isub] : NULL);
}

void
MultiDomainOp::makeAllDOFs(int isub) {
	sd[isub]->makeAllDOFs();
}

void
MultiDomainOp::getInitState(int isub) {
	StackVector disp(v1->subData(isub), v1->subLen(isub));
	StackVector veloc(v2->subData(isub), v1->subLen(isub));
	StackVector accel(v3->subData(isub), v1->subLen(isub));
	StackVector v_p(v4->subData(isub), v1->subLen(isub));

	if (geoSource->getCheckFileInfo()->lastRestartFile) {
		int extlen = (int) std::log10((double) sd[isub]->subNum() + 1) + 1;
		char *ext = new char[extlen + 2];
		sprintf(ext, "_%d", sd[isub]->subNum() + 1);
		sd[isub]->initDispVeloc(disp, veloc, accel, v_p, ext);
		delete[] ext;
	} else
		sd[isub]->initDispVeloc(disp, veloc, accel, v_p);
}

void
MultiDomDynPostProcessor::setPostProcessor(DistFlExchanger *exchanger) {
	distFlExchanger = exchanger;
}

void
MultiDomDynPostProcessor::setUserDefs(double **disps, double **vels) {
	usrDefDisps = disps;
	usrDefVels = vels;
}

void
MultiDomDynPostProcessor::setNodalTemps(DistrVector *_nodalTemps) {
	nodalTemps = _nodalTemps;
}

void
MultiDomDynPostProcessor::dynamOutput(int tIndex, double t, MDDynamMat &dynOps, DistrVector &distForce,
                                      DistrVector *distAeroF, SysState<DistrVector> &distState, DistrVector *distResF) {
	if (!times) times = new StaticTimers;
	startTimerMemory(times->output, times->memoryOutput);

	if (domain->solInfo().nRestart > 0) {
		for (int i = 0; i < decDomain->getNumSub(); ++i) {
			SubDomain *sd = decDomain->getSubDomain(i);
			int extlen = (int) std::log10((double) sd->subNum() + 1) + 1;
			char *ext = new char[extlen + 2];
			sprintf(ext, "_%d", sd->subNum() + 1);
			if (domain->solInfo().isNonLin()) {
				StackVector vel_ni(distState.getVeloc().subData(i), distState.getVeloc().subLen(i));
				StackVector acc_ni(distState.getAccel().subData(i), distState.getAccel().subLen(i));
				sd->writeRestartFile(t, tIndex, vel_ni, acc_ni, (*geomState)[i], ext);
			} else {
				StackVector d_ni(distState.getDisp().subData(i), distState.getDisp().subLen(i));
				StackVector v_ni(distState.getVeloc().subData(i), distState.getVeloc().subLen(i));
				StackVector v_pi(distState.getPrevVeloc().subData(i), distState.getPrevVeloc().subLen(i));
				sd->writeRestartFile(t, tIndex, d_ni, v_ni, v_pi, domain->solInfo().initExtForceNorm, ext);
			}
			delete[] ext;
		}
	}

	// Update bcx for time dependent prescribed displacements and velocities
	ControlLawInfo *claw = geoSource->getControlLaw();
	ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
	if (claw && claw->numUserDisp) {
		double *userDefineDisp = new double[claw->numUserDisp];
		double *userDefineVel = new double[claw->numUserDisp];
		double *userDefineAcc = new double[claw->numUserDisp];
		for (int i = 0; i < claw->numUserDisp; ++i) {
			userDefineVel[i] = 0;
			userDefineAcc[i] = 0;
		}
		userSupFunc->usd_disp(t, userDefineDisp, userDefineVel, userDefineAcc);
		paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC,
		           userDefineDisp, userDefineVel, userDefineAcc, false);
		if (domain->solInfo().ROMPostProcess) {
			execParal(decDomain->getNumSub(), this, &MultiDomDynPostProcessor::subUpdateGeomStateUSDD, userDefineDisp,
			          geomState, userDefineVel, userDefineAcc);
		}
		delete[] userDefineDisp;
		delete[] userDefineVel;
		delete[] userDefineAcc;
	}

	// Send displacements to fluid code (except explicit C0)
	SolverInfo &sinfo = domain->solInfo();
	if (sinfo.aeroFlag >= 0 && !sinfo.lastIt && tIndex != sinfo.initialTimeIndex &&
	    !(sinfo.newmarkBeta == 0 && sinfo.aeroFlag == 20)) {
		domain->getTimers().sendFluidTime -= getTime();
		// Send u + IDISP6 to fluid code.
		// IDISP6 is used to compute pre-stress effects.
		DistrVector d_n_aero(distState.getDisp());

		if (domain->solInfo().gepsFlg == 1) {

			for (int i = 0; i < decDomain->getNumSub(); ++i) {
				SubDomain *sd = decDomain->getSubDomain(i);
				BCond *iDis6 = sd->getInitDisp6();
				for (int j = 0; j < sd->numInitDisp6(); ++j) {
					int dof = sd->getCDSA()->locate(iDis6[j].nnum, 1 << iDis6[j].dofnum);
					if (dof >= 0)
						d_n_aero.subData(i)[dof] += iDis6[j].val;
				}
			}
		}

		SysState<DistrVector> state(d_n_aero, distState.getVeloc(), distState.getAccel(), distState.getPrevVeloc());

		distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
		if (verboseFlag) filePrint(stderr, " ... [E] Sent displacements         ...\n");
		domain->getTimers().sendFluidTime += getTime();
	}

	if (sinfo.aeroheatFlag >= 0 && tIndex != 0) {
		SysState<DistrVector> tempState(distState.getDisp(), distState.getVeloc(), distState.getPrevVeloc());

		distFlExchanger->sendTemperature(tempState);
		if (verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
	}

	if (sinfo.thermohFlag >= 0 && tIndex != 0) {

		for (int i = 0; i < decDomain->getNumSub(); ++i) {
			SubDomain *sd = decDomain->getSubDomain(i);
			for (int j = 0; j < sd->numNodes(); ++j) {
				int tloc = sd->getCDSA()->locate(j, DofSet::Temp);
				int tloc1 = sd->getDSA()->locate(j, DofSet::Temp);
				double temp = (tloc >= 0) ? distState.getDisp().subData(i)[tloc] : sd->getBcx()[tloc1];
				if (tloc1 < 0) temp = 0.0;
				nodalTemps->subData(i)[j] = temp;
			}
		}

		distFlExchanger->sendStrucTemp(*nodalTemps);
		if (verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
	}

	if (sinfo.isNonLin())
		decDomain->postProcessing(geomState, distForce, allCorot, t, &distState, distAeroF, geomState, reactions,
		                          &dynOps, distResF);
	else
		decDomain->postProcessing(distState.getDisp(), distForce, t, distAeroF, tIndex, &dynOps, &distState);
	stopTimerMemory(times->output, times->memoryOutput);
}

void
MultiDomDynPostProcessor::subUpdateGeomStateUSDD(int isub, double *userDefineDisp, DistrGeomState *geomState,
                                                 double *userDefineVel, double *userDefineAcc) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	ControlLawInfo *subClaw = sd->getClaw();
	if (subClaw) {
		if (subClaw->numUserDisp) {
			double *subUserDefineDisp = new double[subClaw->numUserDisp];
			double *subUserDefineVel = new double[subClaw->numUserDisp];
			double *subUserDefineAcc = new double[subClaw->numUserDisp];
			for (int i = 0; i < subClaw->numUserDisp; ++i) {
				int globalIndex = sd->getUserDispDataMap()[i];
				subUserDefineDisp[i] = userDefineDisp[globalIndex];
				subUserDefineVel[i] = userDefineVel[globalIndex];
				subUserDefineAcc[i] = userDefineAcc[globalIndex];
			}
			(*geomState)[isub]->updatePrescribedDisplacement(subUserDefineDisp, subClaw, sd->getNodes(),
			                                                 subUserDefineVel, subUserDefineAcc);
			delete[] subUserDefineDisp;
			delete[] subUserDefineVel;
			delete[] subUserDefineAcc;
		}
	}
}

MultiDomainDynam::~MultiDomainDynam() {
	int nsub = decDomain->getNumSub();
	delete times;
	if (geomState) delete geomState;
	if (refState) delete refState;
	if (kelArray) {
		for (int i = 0; i < nsub; ++i) delete[] kelArray[i];
		delete[] kelArray;
	}
	if (melArray) {
		for (int i = 0; i < nsub; ++i) delete[] melArray[i];
		delete[] melArray;
	}
	if (allCorot) {
		for (int i = 0; i < nsub; ++i) {
			if (allCorot[i]) {
				for (int iElem = 0; iElem < decDomain->getSubDomain(i)->numElements(); ++iElem) {
					if (allCorot[i][iElem] && (allCorot[i][iElem] != dynamic_cast<Corotator *>(decDomain->getSubDomain(
						i)->getElementSet()[iElem])))
						delete allCorot[i][iElem];
				}
				delete[] allCorot[i];
			}
		}
		delete[] allCorot;
	}
	if (usrDefDisps) {
		for (int i = 0; i < decDomain->getNumSub(); ++i) delete[] usrDefDisps[i];
		delete[] usrDefDisps;
	}
	if (usrDefVels) {
		for (int i = 0; i < decDomain->getNumSub(); ++i) delete[] usrDefVels[i];
		delete[] usrDefVels;
	}
	delete decDomain;
	if (reactions) delete reactions;
	if (prevFrc) delete prevFrc;
	if (prevFrcBackup) delete prevFrcBackup;
	if (aeroForce) delete aeroForce;
	if (distFlExchanger) delete distFlExchanger;
}

MDDynamMat *
MultiDomainDynam::buildOps(double coeM, double coeC, double coeK) {
	// Have each subdomain create their operators, then put the
	// dynamic matrices in the Feti Solver
	dynMat = new MDDynamMat;

	times->getFetiSolverTime -= getTime();
	decDomain->buildOps(*dynMat, coeM, coeC, coeK, (Rbm **) 0, kelArray, true, melArray);

	if (domain->tdenforceFlag()) {
		domain->MakeNodalMass(dynMat->M, decDomain->getAllSubDomains());
	}

	int useRbmFilter = (domain->solInfo().isNonLin()) ? 0 : domain->solInfo().filterFlags;
	if (useRbmFilter || domain->solInfo().rbmflg) {
		MultiDomainRbm<double> *rigidBodyModes = decDomain->constructRbm();
		if (useRbmFilter) {
			filePrint(stderr, " ... RBM Filter Level %d Requested   ...\n", useRbmFilter);
			projector_prep(rigidBodyModes, dynMat->M);
		}
		delete rigidBodyModes;
	}

	times->getFetiSolverTime += getTime();
	return dynMat;
}

MultiDomainDynam::MultiDomainDynam(Domain *d)
	: MultiDomainBase(d->solInfo()) {
	domain = d;

#ifdef DISTRIBUTED
	decDomain = new GenDistrDomain<double>(domain);
#else
	decDomain = new GenDecDomain<double>(domain);
#endif
	times = new StaticTimers;

	claw = 0;
	userSupFunc = 0;
	kelArray = 0;
	melArray = 0;
	allCorot = 0;
	geomState = 0;
	refState = 0;
	dynMat = 0;
	reactions = 0;
	prevFrc = 0;
	prevFrcBackup = 0;
	aeroForce = 0;
	distFlExchanger = 0;
	dynMat = 0;
	usrDefDisps = 0;
	usrDefVels = 0;
}

const DistrInfo &
MultiDomainDynam::solVecInfo() const {
	return decDomain->solVecInfo();
}

const DistrInfo &
MultiDomainDynam::masterSolVecInfo() const {
	return decDomain->masterSolVecInfo();
}

DistrInfo &
MultiDomainDynam::bcInfo() {
	// prescribed boundary condition distributed vector information
	return *decDomain->pbcVectorInfo();
}

void
MultiDomainDynam::processLastOutput() {
	OutputInfo *oinfo = geoSource->getOutputInfo();
	domain->solInfo().lastIt = true;
	for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
		oinfo[iOut].interval = 1;
}

void
MultiDomainDynam::preProcess() {
	times->preProcess -= getTime();

	// Makes local renumbering, connectivities and dofsets
	decDomain->preProcess();

	// Make all element's dofs
	MultiDomainOp mdop(&MultiDomainOp::makeAllDOFs, decDomain->getAllSubDomains());
#ifdef DISTRIBUTED
	execParal(decDomain->getNumSub(), &mdop, &MultiDomainOp::runFor);
#else
	threadManager->execParal(decDomain->getNumSub(), &mdop);
#endif
	times->preProcess += getTime();

	// Check for user supplied routines (control, force or displacement)
	claw = geoSource->getControlLaw();
	userSupFunc = geoSource->getUserSuppliedFunction();

	// Make the geomState and refState (used for prestress, explicit nonlinear and contact)
	if ((domain->solInfo().gepsFlg == 1 && domain->numInitDisp6() > 0) || domain->solInfo().isNonLin() ||
	    domain->tdenforceFlag()) {
		times->timeGeom -= getTime();
		geomState = new DistrGeomState(decDomain);
		refState = new DistrGeomState(decDomain);
		times->timeGeom += getTime();
	}

	// Update geomState with prescribed dirichlet boundary conditions (explicit nonlinear and contact)
	if (domain->solInfo().isNonLin() || domain->tdenforceFlag())
		execParal(decDomain->getNumSub(), this, &MultiDomainDynam::initSubPrescribedDisplacement);

	// Make corotators and kelArray (used for prestress and explicit nonlinear)
	if ((domain->solInfo().gepsFlg == 1 && domain->numInitDisp6() > 0) || domain->solInfo().isNonLin()) {
		times->corotatorTime -= getTime();
		allCorot = new Corotator **[decDomain->getNumSub()];
		execParal(decDomain->getNumSub(), this, &MultiDomainDynam::makeSubCorotators);
		times->corotatorTime += getTime();

		times->kelArrayTime -= getTime();
		kelArray = new FullSquareMatrix *[decDomain->getNumSub()];
		if (domain->solInfo().isNonLin() && (domain->solInfo().newmarkBeta == 0 || domain->solInfo().samplingPodRom
		                                     || domain->solInfo().svdPodRom || domain->solInfo().ROMPostProcess))
			melArray = new FullSquareMatrix *[decDomain->getNumSub()];
		execParal(decDomain->getNumSub(), this, &MultiDomainDynam::makeSubElementArrays);
		times->kelArrayTime += getTime();
	}

	// Initialization for contact
	if (domain->tdenforceFlag())
		domain->InitializeDynamicContactSearch(decDomain->getNumSub(), decDomain->getAllSubDomains());
	else if (domain->solInfo().isNonLin() && domain->GetnContactSurfacePairs()) {
		filePrint(stderr,
		          " *** ERROR: \"tdenforce off\" is not supported for multi-domain nonlinear dynamics. Exiting...\n");
		exit(-1);
	}

	// Allocate vector to store reaction forces
	if (!reactions) reactions = new DistrVector(*decDomain->pbcVectorInfo());
}

void
MultiDomainDynam::makeSubCorotators(int isub) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	int numele = sd->numElements();
	allCorot[isub] = new Corotator *[numele];
	sd->createCorotators(allCorot[isub]);
}

void
MultiDomainDynam::makeSubElementArrays(int isub) {
	SubDomain *sd = decDomain->getSubDomain(isub);

	// allocate the element mass and/or stiffness array
	if (melArray) sd->createKelArray(kelArray[isub], melArray[isub]);
	else sd->createKelArray(kelArray[isub]);

	// update geomState with IDISP6 if GEPS is requested (geometric prestress / linear only)
	if ((sd->numInitDisp6() > 0) && (domain->solInfo().gepsFlg == 1)) // GEPS
		(*geomState)[isub]->updatePrescribedDisplacement(sd->getInitDisp6(), sd->numInitDisp6(), sd->getNodes());

	// build the element stiffness matrices.
	if (!domain->solInfo().ROMPostProcess && !domain->solInfo().galerkinPodRom && !domain->solInfo().samplingPodRom) {
		// Note: for explicit nonlinear ROMs the initial tangent stiffness is required for stability timestep computation
		//       but it is computed later after nodal inertia assembly (see DistrExplicitPodProjectionNonLinDynamicBase::preProcess)
		Vector elementInternalForce(sd->maxNumDOF(), 0.0);
		Vector residual(sd->numUncon(), 0.0);
		sd->getStiffAndForce(*(*geomState)[isub], elementInternalForce, allCorot[isub], kelArray[isub], residual,
		                     1.0, 0.0, (*geomState)[isub], (Vector *) NULL, ((melArray) ? melArray[isub] : NULL));
	}
}

void
MultiDomainDynam::initSubPrescribedDisplacement(int isub) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	if (sd->nDirichlet() > 0)
		(*geomState)[isub]->updatePrescribedDisplacement(sd->getDBC(), sd->nDirichlet(), sd->getNodes());
}

void
MultiDomainDynam::getTimes(double &dt, double &tmax) {
	dt = domain->solInfo().getTimeStep();
	tmax = domain->solInfo().tmax;
}

void
MultiDomainDynam::getInitialTime(int &initTimeIndex, double &initTime) {
	initTimeIndex = domain->solInfo().initialTimeIndex;
	initTime = domain->solInfo().initialTime;
}

double
MultiDomainDynam::getInitialForceNorm() {
	return domain->solInfo().initExtForceNorm;
}

int
MultiDomainDynam::getTimeIntegration() {
	return domain->solInfo().timeIntegration;
}

int
MultiDomainDynam::getFilterFlag() {
	int filterFlag = domain->solInfo().filterFlags;
	// note: only level 1 rbm filter is used for nonlinear
	return (domain->solInfo().isNonLin()) ? std::min(filterFlag, 1) : filterFlag;
}

int *
MultiDomainDynam::boundary() {
	filePrint(stderr, "Paral.d/MDDynam.C: boundary not implemented here\n");
	return 0;
}

double *
MultiDomainDynam::boundaryValue() {
	filePrint(stderr, "Paral.d/MDDynam.C: boundaryValue not implemented here\n");
	return 0;
}

Domain *
MultiDomainDynam::getDomain() {
	return domain;
}

void
MultiDomainDynam::getNewMarkParameters(double &beta, double &gamma,
                                       double &alphaf, double &alpham) {
	beta = domain->solInfo().newmarkBeta;
	gamma = domain->solInfo().newmarkGamma;
	alphaf = domain->solInfo().newmarkAlphaF;
	alpham = domain->solInfo().newmarkAlphaM;
}

void
MultiDomainDynam::getQuasiStaticParameters(double &maxVel, double &delta) {
	maxVel = domain->solInfo().qsMaxvel;
	delta = domain->solInfo().delta;
}

void
MultiDomainDynam::getSteadyStateParam(int &steadyFlag, int &steadyMin,
                                      int &steadyMax, double &steadyTol) {
	steadyFlag = domain->solInfo().steadyFlag;
	steadyMin = domain->solInfo().steadyMin;
	steadyMax = domain->solInfo().steadyMax;
	steadyTol = domain->solInfo().steadyTol;
}

void
MultiDomainDynam::getSensitivityStateParam(double &sensitivityTol, double &ratioSensitivityTol) {
	sensitivityTol = domain->solInfo().sensitivityTol;
	ratioSensitivityTol = domain->solInfo().ratioSensitivityTol;
}

void
MultiDomainDynam::getContactForce(DistrVector &d_n, DistrVector &dinc, DistrVector &ctc_f, double t_n_p, double dt,
                                  double dt_old) {
	ctc_f.zero();
	if(t_n_p < domain->solInfo().tdenforceInitia || t_n_p >= domain->solInfo().tdenforceFinal) return;
	if (domain->tdenforceFlag()) {
		times->tdenforceTime -= getTime();

		times->updateSurfsTime -= getTime();
		domain->UpdateSurfaceTopology(decDomain->getNumSub(), decDomain->getAllSubDomains()); // remove deleted elements
		domain->UpdateSurfaces(geomState, 1, decDomain->getAllSubDomains()); // update to current configuration

		// copy and update the current state (geomState) to the predicted state
		DistrGeomState *predictedState = new DistrGeomState(*geomState);
		if (domain->solInfo().isNonLin()) {
			predictedState->update(dinc, 1);
		} else {
			DistrVector d_n_p(decDomain->solVecInfo());
			d_n_p = d_n + dinc;
			execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subExplicitUpdate, d_n_p, predictedState);
		}
		times->updateSurfsTime += getTime();

		// update the prescribed displacements to their correct value at the time of the predictor
		if (claw && userSupFunc && claw->numUserDisp) {
			double *userDefineDisp = new double[claw->numUserDisp];
			double *userDefineVel = new double[claw->numUserDisp];
			double *userDefineAcc = new double[claw->numUserDisp];
			for (int i = 0; i < claw->numUserDisp; ++i) {
				userDefineVel[i] = 0;
				userDefineAcc[i] = 0;
			}
			userSupFunc->usd_disp(t_n_p, userDefineDisp, userDefineVel, userDefineAcc);
			execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subUpdateGeomStateUSDD, userDefineDisp,
			          predictedState,
			          userDefineVel, userDefineAcc);
			delete[] userDefineDisp;
			delete[] userDefineVel;
			delete[] userDefineAcc;
		}

		times->updateSurfsTime -= getTime();
		domain->UpdateSurfaces(predictedState, 2, decDomain->getAllSubDomains()); // update to predicted configuration
		times->updateSurfsTime += getTime();

		times->contactSearchTime -= getTime();
		domain->PerformDynamicContactSearch(dt_old, dt);
		times->contactSearchTime += getTime();

		times->contactForcesTime -= getTime();
		domain->AddContactForces(dt_old, dt, ctc_f);
		times->contactForcesTime += getTime();

		delete predictedState;
		times->tdenforceTime += getTime();
	}
}

void
MultiDomainDynam::updateState(double dt_n_h, DistrVector &v_n_h, DistrVector &d_n) {
	if (domain->solInfo().isNonLin()) {
		*refState = *geomState; // (AN) update refState values
		DistrVector dinc(solVecInfo());
		dinc = dt_n_h * v_n_h;
		geomState->update(dinc, 1);
		geomState->setVelocity(v_n_h);
		geomState->get_tot_displacement(d_n, false);
	}
}

void
MultiDomainDynam::subUpdateStates(int isub, DistrGeomState *refState, DistrGeomState *geomState, double time) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	GeomState *subRefState = (refState) ? (*refState)[isub] : 0;
	sd->updateStates(subRefState, *(*geomState)[isub], allCorot[isub], time);
}

void
MultiDomainDynam::computeExtForce2(SysState<DistrVector> &distState,
                                   DistrVector &f, DistrVector &cnst_f, int tIndex,
                                   double t, DistrVector *aero_f,
                                   double gamma, double alphaf) {
	times->formRhs -= getTime();
	SolverInfo &sinfo = domain->solInfo();

	// compute USDD prescribed displacements
	double *userDefineDisp = 0;
	double *userDefineVel = 0;
	double *userDefineAcc = 0;
	if (claw && userSupFunc) {
		if (claw->numUserDisp) {
			userDefineDisp = new double[claw->numUserDisp];
			userDefineVel = new double[claw->numUserDisp];
			userDefineAcc = new double[claw->numUserDisp];
			for (int i = 0; i < claw->numUserDisp; ++i) {
				userDefineVel[i] = 0;
				userDefineAcc[i] = 0;
			}
			userSupFunc->usd_disp(t, userDefineDisp, userDefineVel, userDefineAcc);
			paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC,
			           userDefineDisp, userDefineVel, userDefineAcc, false); // update bcx, vcx, acx
		}
	}

	// finish update of geomState. note that for nonlinear problems the positiion and rotation nodal variables
	// have already been updated in updateDisplacement
	if (sinfo.isNonLin() || domain->tdenforceFlag()) {
		if (!sinfo.isNonLin()) {
			execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subExplicitUpdate, distState.getDisp(),
			          geomState);
		}
		execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subUpdateGeomStateUSDD, userDefineDisp, geomState,
		          userDefineVel, userDefineAcc);
	}

	// update nodal temperatures for thermoe problem
	if (sinfo.thermoeFlag >= 0 && tIndex >= 0) {
		distFlExchanger->getStrucTemp(nodalTemps->data());
		if (verboseFlag) filePrint(stderr, " ... [E] Received temperatures     ...\n");
		if (geomState) geomState->setNodalTemperatures(*nodalTemps);
	}

	// add f(t) to cnst_f
	double dt = sinfo.getTimeStep();
	double alpham = sinfo.newmarkAlphaM;
	double t0 = sinfo.initialTime;
	double tm = (t == t0) ? t0 : t + dt * (alphaf - alpham);
	MultiDomainOp mdop(&MultiDomainOp::computeExtForce,
	                   decDomain->getAllSubDomains(), &f, &cnst_f, t, dynMat->Kuc, userSupFunc, dynMat->Cuc, tm,
	                   dynMat->Muc);
	threadManager->execParal(decDomain->getNumSub(), &mdop);
	if (userDefineDisp) delete[] userDefineDisp;
	if (userDefineVel) delete[] userDefineVel;
	if (userDefineAcc) delete[] userDefineAcc;

	// add USDF forces
	if (claw && userSupFunc) {
		if (claw->numUserForce) {
			double *userDefineForce = new double[claw->numUserForce];
			userSupFunc->usd_forc(t, userDefineForce);
			decDomain->addUserForce(f, userDefineForce);
			delete[] userDefineForce;
		}
	}

	// add ACTUATOR forces
	if (claw && userSupFunc) {
		if (claw->numActuator) {
			double *ctrdisp = new double[claw->numSensor];
			double *ctrvel = new double[claw->numSensor];
			double *ctracc = new double[claw->numSensor];
			double *ctrfrc = new double[claw->numActuator];
#ifdef DISTRIBUTED
			for (int i = 0; i < claw->numSensor; ++i)
				ctrdisp[i] = ctrvel[i] = ctracc[i] = std::numeric_limits<double>::min();
#endif
			DistrVector &disp = distState.getDisp();
			DistrVector &vel = distState.getVeloc();
			DistrVector &acc = distState.getAccel();
			decDomain->extractControlData(disp, vel, acc, ctrdisp, ctrvel, ctracc);
#ifdef DISTRIBUTED
			structCom->globalMax(claw->numSensor, ctrdisp);
			structCom->globalMax(claw->numSensor, ctrvel);
			structCom->globalMax(claw->numSensor, ctracc);
#endif
			userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t);
			decDomain->addCtrl(f, ctrfrc);
			delete[] ctrdisp;
			delete[] ctrvel;
			delete[] ctracc;
			delete[] ctrfrc;
		}
	}

	// add aeroelastic forces from fluid dynamics code
	if (sinfo.aeroFlag >= 0 && tIndex >= 0 &&
	      !(geoSource->getCheckFileInfo()->hotRestart() && sinfo.aeroFlag == 20 && !sinfo.dyna3d_compat &&
	      tIndex == sinfo.initialTimeIndex)) {
		if(tIndex % sinfo.subcycle == 0) {

		domain->getTimers().receiveFluidTime -= getTime();
		aeroForce->zero();
		int iscollocated;
		double tFluid = distFlExchanger->getFluidLoad(*aeroForce, tIndex, t,
		                                              alphaf, iscollocated);
		if (verboseFlag) filePrint(stderr, " ... [E] Received fluid load        ...\n");

		if (sinfo.aeroFlag == 20) {
			if (prevIndex >= 0)
				aero_f->linC(0.5, *aeroForce, 0.5, *prevFrc);
			else
				*aero_f = *aeroForce;
		} else {
			if (iscollocated == 0) {
				if (prevIndex >= 0) {
					*aeroForce *= (1 / gamma);
					aeroForce->linAdd(((gamma - 1.0) / gamma), *prevFrc);
				}
			}

			double alpha = 1.0 - alphaf;
			if (prevIndex < 0) alpha = 1.0;

			aero_f->linC(alpha, *aeroForce, (1.0 - alpha), *prevFrc);
		}

		*prevFrc = *aeroForce;
		prevTime = tFluid;
		prevIndex = tIndex;

		if (sinfo.aeroFlag == 20 && sinfo.dyna3d_compat) {
			if (sinfo.stop_AeroF) sinfo.stop_AeroS = true;
			double dt = sinfo.getTimeStep()*sinfo.subcycle;
			if (tIndex == sinfo.initialTimeIndex) sinfo.t_AeroF = sinfo.initialTime + 1.5 * dt;
			else sinfo.t_AeroF += dt;
			double maxTime_AeroF = sinfo.tmax - 0.5 * dt;
			double sendtim;
			if ((sinfo.t_AeroF < (maxTime_AeroF - 0.01 * dt)) || tIndex == sinfo.initialTimeIndex)
				sendtim = 1e100;
			else {
				sendtim = 0;
				sinfo.stop_AeroF = true;
			}
			int restartinc = std::max(sinfo.nRestart, 0);
			distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep()*sinfo.subcycle, sendtim, restartinc,
			                           sinfo.isCollocated, sinfo.alphas);
			if (tIndex == 0) // Send the parameter a second time for fluid iteration 1 to 2
				distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep()*sinfo.subcycle, sendtim, restartinc,
				                           sinfo.isCollocated, sinfo.alphas);
		}
		domain->getTimers().receiveFluidTime += getTime();
		}
		f += *aero_f;
	}

	// add aerothermal fluxes from fluid dynamics code
	if (sinfo.aeroheatFlag >= 0 && tIndex >= 0) {

		aeroForce->zero();
		double tFluid = distFlExchanger->getFluidFlux(*aeroForce, tIndex, t);
		if (verboseFlag) filePrint(stderr, " ... [T] Received fluid fluxes      ...\n");

		/*  Compute fluid flux at n+1/2, since we use midpoint rule in thermal */

		int useProjector = sinfo.filterFlags;

		if (tIndex == 0)
			f += *aeroForce;
		else {
			if (useProjector) f = *aeroForce;
			else
				f.linAdd(0.5, *aeroForce, 0.5, *prevFrc);
		}

		*prevFrc = *aeroForce;
	}

	// apply rbmfilter projection
	if (sinfo.filterFlags) trProject(f);

	if (tIndex == 1)
		sinfo.initExtForceNorm = f.norm();

	times->formRhs += getTime();

}

void
MultiDomainDynam::getConstForce(DistrVector &v) {
	times->formRhs -= getTime();
	MultiDomainOp mdop(&MultiDomainOp::getConstForce, decDomain->getAllSubDomains(), &v, (dynMat) ? dynMat->Kuc : NULL);
	threadManager->execParal(decDomain->getNumSub(), &mdop);
	times->formRhs += getTime();
}

void
MultiDomainDynam::addConstForceSensitivity(DistrVector &v) {
	filePrint(stderr, " ... MultiDomainDynam::addConstForceSensitivity is not implemented\n");
}

void
MultiDomainDynam::getInitState(SysState<DistrVector> &state) {
	// initialize state with IDISP/IDISP6/IVEL/IACC or RESTART
	MultiDomainOp mdop(&MultiDomainOp::getInitState, decDomain->getAllSubDomains(),
	                   &state.getDisp(), &state.getVeloc(), &state.getAccel(),
	                   &state.getPrevVeloc());
	threadManager->execParal(decDomain->getNumSub(), &mdop);

	if (sinfo.filterFlags) {
		project(state.getDisp());
		project(state.getVeloc());
		project(state.getAccel());
		project(state.getPrevVeloc());
	}

	if (geoSource->getCheckFileInfo()->lastRestartFile) {
		filePrint(stderr, " ... Restarting From a Previous Run ...\n");
		if (domain->solInfo().isNonLin()) {
			for (int i = 0; i < decDomain->getNumSub(); ++i) {
				SubDomain *sd = decDomain->getSubDomain(i);
				StackVector d_ni(state.getDisp().subData(i), state.getDisp().subLen(i));
				StackVector v_ni(state.getVeloc().subData(i), state.getVeloc().subLen(i));
				StackVector a_ni(state.getAccel().subData(i), state.getAccel().subLen(i));
				StackVector v_pi(state.getPrevVeloc().subData(i), state.getPrevVeloc().subLen(i));
				int extlen = (int) std::log10((double) sd->subNum() + 1) + 1;
				char *ext = new char[extlen + 2];
				sprintf(ext, "_%d", sd->subNum() + 1);
				sd->readRestartFile(d_ni, v_ni, a_ni, v_pi, sd->getBcx(), sd->getVcx(), *((*geomState)[i]), ext);
				delete[] ext;
				sd->updateStates((*geomState)[i], *((*geomState)[i]), allCorot[i], sd->solInfo().initialTime);
			}
		}
		domain->solInfo().initialTimeIndex = decDomain->getSubDomain(0)->solInfo().initialTimeIndex;
		domain->solInfo().initialTime = decDomain->getSubDomain(0)->solInfo().initialTime;
		domain->solInfo().initExtForceNorm = decDomain->getSubDomain(0)->solInfo().initExtForceNorm;
	} else if (geomState) {
		geomState->update(state.getDisp());
		geomState->setVelocityAndAcceleration(state.getVeloc(), state.getAccel());
	}

	// if we have a user supplied function, give it the initial state at the sensors
	// .. first update bcx, vcx in case any of the sensors have prescribed displacements
	if (claw && userSupFunc) {
		if (claw->numUserDisp) {
			double *userDefineDisp = new double[claw->numUserDisp];
			double *userDefineVel = new double[claw->numUserDisp];
			double *userDefineAcc = new double[claw->numUserDisp];
			for (int i = 0; i < claw->numUserDisp; ++i) {
				userDefineVel[i] = 0;
				userDefineAcc[i] = 0;
			}
			userSupFunc->usd_disp(domain->solInfo().initialTime, userDefineDisp, userDefineVel, userDefineAcc);
			paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC,
			           userDefineDisp, userDefineVel, userDefineAcc, false);
			delete[] userDefineDisp;
			delete[] userDefineVel;
			delete[] userDefineAcc;
		}
		if (claw->numSensor) {
			double *ctrdisp = new double[claw->numSensor];
			double *ctrvel = new double[claw->numSensor];
			double *ctracc = new double[claw->numSensor];
#ifdef DISTRIBUTED
			for (int i = 0; i < claw->numSensor; ++i)
				ctrdisp[i] = ctrvel[i] = ctracc[i] = std::numeric_limits<double>::min();
#endif
			DistrVector &disp = state.getDisp();
			DistrVector &vel = state.getVeloc();
			DistrVector &acc = state.getAccel();
			decDomain->extractControlData(disp, vel, acc, ctrdisp, ctrvel, ctracc);
#ifdef DISTRIBUTED
			structCom->globalMax(claw->numSensor, ctrdisp);
			structCom->globalMax(claw->numSensor, ctrvel);
			structCom->globalMax(claw->numSensor, ctracc);
#endif
			userSupFunc->init(ctrdisp, ctrvel, ctracc);
			delete[] ctrdisp;
			delete[] ctrvel;
			delete[] ctracc;
		}
	}
}

MultiDomDynPostProcessor *
MultiDomainDynam::getPostProcessor() {
	if (domain->solInfo().aeroFlag >= 0) {
		mddPostPro = new MultiDomDynPostProcessor(decDomain, distFlExchanger, times, geomState, allCorot, melArray,
		                                          reactions);
		return mddPostPro;
	} else {
		mddPostPro = new MultiDomDynPostProcessor(decDomain, times, geomState, allCorot, melArray, reactions);
	}
	return mddPostPro;
}

void
MultiDomainDynam::printTimers(MDDynamMat *dynOps, double timeLoop) {
	times->numSubdomain = decDomain->getNumSub();
	//filePrint(stderr," ... Print Timers                   ... \n");

	for (int i = 0; i < decDomain->getNumSub(); ++i) {
		domain->getTimers().formTime += decDomain->getSubDomain(i)->getTimers().formTime;
		domain->getTimers().assemble += decDomain->getSubDomain(i)->getTimers().assemble;
		domain->getTimers().formRhs += decDomain->getSubDomain(i)->getTimers().formRhs;
	}

	if (isFeti(domain->solInfo().solvercntl->type) && domain->solInfo().solvercntl->fetiInfo.version == 3) {
		times->printFetiDPtimers(domain->getTimers(),
		                         dynOps->dynMat->getSolutionTime(),
		                         domain->solInfo(),
		                         dynOps->dynMat->getTimers(),
		                         geoSource->getCheckFileInfo()[0],
		                         domain);
	} else {
		times->printStaticTimers(domain->getTimers(),
		                         dynOps->dynMat->getSolutionTime(),
		                         domain->solInfo(),
		                         dynOps->dynMat->getTimers(),
		                         geoSource->getCheckFileInfo()[0],
		                         domain);
	}

/*
   switch(domain->solInfo().solvercntl->fetiInfo.version) {
     default:
     case FetiInfo::feti1:
     case FetiInfo::feti2:
       times->printStaticTimers(domain->getTimers(),
                                dynOps->dynMat->getSolutionTime(),
                                domain->solInfo() ,
                                dynOps->dynMat->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;
                                                                                                 
     case FetiInfo::fetidp:
       times->printFetiDPtimers(domain->getTimers(),
                                dynOps->dynMat->getSolutionTime(),
                                domain->solInfo() ,
                                dynOps->dynMat->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;
   }
*/
}

void
MultiDomainDynam::getRayleighCoef(double &alpha) {
	alpha = domain->solInfo().alphaDamp;
}

SubDOp *
MultiDomainDynam::getpK(MDDynamMat *dynOps) {
	return dynOps->K;
}

SubDOp *
MultiDomainDynam::getpM(MDDynamMat *dynOps) {
	return dynOps->M;
}

SubDOp *
MultiDomainDynam::getpC(MDDynamMat *dynOps) {
	return dynOps->C;
}

void
MultiDomainDynam::computeStabilityTimeStep(double &dt, MDDynamMat &dMat) {
	double dt_c;
	int eid_c;
	if (domain->solInfo().isNonLin()) {
		dt_c = std::numeric_limits<double>::infinity();
		eid_c = -1;
		for (int i = 0; i < decDomain->getNumSub(); ++i) {
			int eid_ci;
			double dt_ci = decDomain->getSubDomain(i)->computeStabilityTimeStep(kelArray[i], melArray[i],
			                                                                    (*geomState)[i], eid_ci);
			if (dt_ci < dt_c) eid_c = eid_ci;
			dt_c = std::min(dt_c, dt_ci);
		}
#ifdef DISTRIBUTED
		double dt_cp = dt_c;
		dt_c = structCom->globalMin(dt_c);
		if (dt_cp != dt_c) eid_c = -1;
		eid_c = structCom->globalMax(eid_c);
#endif
	} else
		dt_c = decDomain->computeStabilityTimeStep(dMat);

	if (dt_c == std::numeric_limits<double>::infinity()) {
		filePrint(stderr, " **************************************\n");
		filePrint(stderr, " Stability max. timestep could not be  \n");
		filePrint(stderr, " determined for this model.            \n");
		if (domain->solInfo().isNonLin() && eid_c > -1) {
			filePrint(stderr, " Element with inf. time step = %7d\n", eid_c + 1);
		}
		filePrint(stderr, " Specified time step is selected\n");
		filePrint(stderr, " **************************************\n");
		domain->solInfo().stable = 0;
	} else {
		filePrint(stderr, " **************************************\n");
		if (domain->solInfo().modifiedWaveEquation) {
			dt_c = 1.73205 * dt_c;
			filePrint(stderr, " CONDITIONALLY STABLE MODIFIED WAVE EQUATION \n");
		} else
			filePrint(stderr, " CONDITIONALLY STABLE NEWMARK ALGORITHM \n");
		filePrint(stderr, " --------------------------------------\n");
		filePrint(stderr, " Specified time step      = %10.4e\n", dt);
		filePrint(stderr, " Stability max. time step = %10.4e\n", dt_c);
		if (domain->solInfo().isNonLin()) {
			filePrint(stderr, " Element with min. time step = %7d\n", eid_c + 1);
		}
		filePrint(stderr, " **************************************\n");
		if ((domain->solInfo().stable == 1 && dt_c < dt) || domain->solInfo().stable == 2) {
			dt = dt_c;
			filePrint(stderr, " Stability max. time step is selected\n");
		} else
			filePrint(stderr, " Specified time step is selected\n");
		filePrint(stderr, " **************************************\n");
	}

	for (int i = 0; i < decDomain->getNumSub(); ++i)
		decDomain->getSubDomain(i)->solInfo().setTimeStep(dt);
	domain->solInfo().setTimeStep(dt);
}

int
MultiDomainDynam::getModeDecompFlag() {
	return domain->solInfo().modeDecompFlag;
}

void
MultiDomainDynam::modeDecompPreProcess(SparseMatrix *M) {
	filePrint(stderr, "Paral.d/MDDynam.C: modeDecompPreProcess not implemented here\n");
}

void
MultiDomainDynam::modeDecomp(double t, int tIndex, DistrVector &d_n) {
	filePrint(stderr, "Paral.d/MDDynam.C: modeDecomp not implemented here\n");
}

void
MultiDomainDynam::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {
	if (domain->solInfo().isNonLin()) {
		execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subGetInternalForce, f, t, tIndex);
	} else {
		f.zero();
		execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subGetKtimesU, d, f);
	}

	if (domain->solInfo().timeIntegration == 1) decDomain->getSolVecAssembler()->assemble(f); // quasistatic only
}

void
MultiDomainDynam::getFollowerForce(DistrVector &f, double t, int tIndex) {
	execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subGetFollowerForce, f, t, tIndex);
}

void
MultiDomainDynam::pull_back(DistrVector &f) {
	if (domain->solInfo().isNonLin() && !domain->solInfo().galerkinPodRom &&
	    !domain->solInfo().getNLInfo().linearelastic) {
		// Transform both moments and forces to convected frame: f = [R^T  I ]*f
		//                                                           [ I  R^T]
		geomState->pull_back(f);
	}
}

void
MultiDomainDynam::push_forward(DistrVector &a) {
	if (domain->solInfo().isNonLin() && !domain->solInfo().galerkinPodRom &&
	    !domain->solInfo().getNLInfo().linearelastic) {
		// Transform 2nd time-derivative of displacement to spatial frame: a = [R I]*a
		//                                                                     [I I]
		// Note: the angular accelerations are deliberately not transformed.
		geomState->push_forward(a);
	}
}

void
MultiDomainDynam::subExplicitUpdate(int isub, DistrVector &d, DistrGeomState *geomState) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	StackVector subd(d.subData(isub), d.subLen(isub));
	(*geomState)[isub]->explicitUpdate(sd->getNodes(), subd);
}

void
MultiDomainDynam::subUpdateGeomStateUSDD(int isub, double *userDefineDisp, DistrGeomState *geomState,
                                         double *userDefineVel, double *userDefineAcc) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	ControlLawInfo *subClaw = sd->getClaw();
	if (subClaw) {
		if (subClaw->numUserDisp) {
			double *subUserDefineDisp = new double[subClaw->numUserDisp];
			double *subUserDefineVel = new double[subClaw->numUserDisp];
			double *subUserDefineAcc = new double[subClaw->numUserDisp];
			for (int i = 0; i < subClaw->numUserDisp; ++i) {
				int globalIndex = sd->getUserDispDataMap()[i];
				subUserDefineDisp[i] = userDefineDisp[globalIndex];
				subUserDefineVel[i] = userDefineVel[globalIndex];
				subUserDefineAcc[i] = userDefineAcc[globalIndex];
			}
			(*geomState)[isub]->updatePrescribedDisplacement(subUserDefineDisp, subClaw, sd->getNodes(),
			                                                 subUserDefineVel, subUserDefineAcc);
			delete[] subUserDefineDisp;
			delete[] subUserDefineVel;
			delete[] subUserDefineAcc;
		}
	}
}

void
MultiDomainDynam::subUpdateUsrDefDispsAndVels(int isub, double *userDefineDisp, double *userDefineVel) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	ControlLawInfo *claw = sd->getClaw();
	DofSetArray *dsa = sd->getDSA();
	int *locToGlUserDispMap = sd->getUserDispDataMap();
	double *bcx = usrDefDisps[isub];
	double *vcx = usrDefVels[isub];

	for (int i = 0; i < claw->numUserDisp; ++i) {
		int dof = dsa->locate(claw->userDisp[i].nnum, 1 << claw->userDisp[i].dofnum);
		if (dof >= 0) {
			bcx[dof] = userDefineDisp[locToGlUserDispMap[i]];
			vcx[dof] = userDefineVel[locToGlUserDispMap[i]];
		}
	}
}

void
MultiDomainDynam::subGetInternalForce(int isub, DistrVector &f, double &t, int &tIndex) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	Vector residual(f.subLen(isub), 0.0);
	Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)
	StackVector *subReactions = NULL;
	if (reactions) {
		subReactions = new StackVector(reactions->subData(isub), reactions->subLen(isub));
		subReactions->zero();
	}

	// NOTE #1: for explicit nonlinear dynamics, geomState and refState are the same object -- AN: no longer true
	// NOTE #2: by convention, the internal variables associated with a nonlinear constitutive relation are not updated
	//          when getStiffAndForce is called, so we have to call updateStates.
	if (domain->solInfo().newmarkBeta == 0 && domain->solInfo().stable && domain->solInfo().isNonLin() &&
	    tIndex % domain->solInfo().stable_freq == 0) {
		sd->getStiffAndForce(*(*geomState)[isub], eIF, allCorot[isub], kelArray[isub], residual, 1.0, t,
		                     (*refState)[isub],
		                     subReactions, melArray[isub]);
/* PJSA 10/12/2014 this is done in getStiffAndForce now because it needs to be done before handleElementDeletion.
    sd->updateStates((*geomState)[isub], *(*geomState)[isub], allCorot[isub]);
*/
	} else {
		sd->getInternalForce(*(*geomState)[isub], eIF, allCorot[isub], kelArray[isub], residual, 1.0, t,
		                     (*refState)[isub],
		                     subReactions, (melArray) ? melArray[isub] : NULL);
	}
	StackVector subf(f.subData(isub), f.subLen(isub));
	subf.linC(residual, -1.0); // f = -residual
	if (subReactions) delete subReactions;
}

void
MultiDomainDynam::subGetKtimesU(int isub, DistrVector &d, DistrVector &f) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	StackVector subf(f.subData(isub), f.subLen(isub));
	StackVector subd(d.subData(isub), d.subLen(isub));
	sd->getKtimesU(subd, (double *) 0, subf, 1.0, (kelArray) ? kelArray[isub] : (FullSquareMatrix *) 0);
}

void
MultiDomainDynam::subGetFollowerForce(int isub, DistrVector &f, double &t, int &tIndex) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	StackVector subf(f.subData(isub), f.subLen(isub));
	Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)
	sd->getFollowerForce(*(*geomState)[isub], eIF, allCorot[isub], kelArray[isub], subf, 1.0, t, (Vector *) NULL,
	                     false);
}

void
MultiDomainDynam::computeTimeInfo() {
	// Time integration information

	// Get total time and time step size and store them
	double totalTime = domain->solInfo().tmax;
	double dt = domain->solInfo().getTimeStep();

	// Compute maximum number of steps
	int maxStep = (int) ((totalTime + 0.49 * dt) / dt);

	// Compute time remainder
	double remainder = totalTime - maxStep * dt;
	if (std::abs(remainder) > 0.01 * dt) {
		domain->solInfo().tmax = maxStep * dt;
		filePrint(stderr, " Warning: Total time is being changed to : %e\n", domain->solInfo().tmax);
	}
}

int
MultiDomainDynam::aeroSensitivityPreProcess(DistrVector &disp, DistrVector &vel, DistrVector &accel,
                                            DistrVector &lastVel) {
	filePrint(stderr, " ... MultiDomainDynam::aeroSensitivityPreProcess has not been implemented\n");
	return 0;
}

int
MultiDomainDynam::sendDisplacements(DistrVector &disp, DistrVector &vel, DistrVector &accel, DistrVector &lastVel) {
	filePrint(stderr, " ... MultiDomainDynam::sendDisplacements has not been implemented\n");
	return 0;
}

int
MultiDomainDynam::aeroPreProcess(DistrVector &disp, DistrVector &vel,
                                 DistrVector &accel, DistrVector &lastVel) {
	// get solver info
	SolverInfo &sinfo = domain->solInfo();

	// Initialize previous force data
	// Reexamine for the case of restart
	prevFrc = new DistrVector(solVecInfo());
	prevFrc->zero();
	prevFrcBackup = new DistrVector(solVecInfo());
	prevIndex = -1;
	prevTime = 0;

	// Initialize the aeroforce vector
	aeroForce = new DistrVector(solVecInfo());

	if (sinfo.aeroFlag < 0)
		return 0;

	auto cpuToSub = geoSource->getCpuToSub();

	// get cpu id
#ifdef USE_MPI
	int myId = structCom->myID();
#else
	int myId = 0;
#endif

	int numLocSub = cpuToSub->num(myId);

	SubDomain **subdomain = decDomain->getAllSubDomains();

	// allocate for pointer arrays
	CoordSet **cs = new CoordSet *[numLocSub];
	Elemset **elemSet = new Elemset *[numLocSub];
	DofSetArray **cdsa = new DofSetArray *[numLocSub];
	DofSetArray **dsa = new DofSetArray *[numLocSub];
	usrDefDisps = new double *[numLocSub];
	usrDefVels = new double *[numLocSub];

	int iSub;
	for (iSub = 0; iSub < numLocSub; iSub++) {

		// assemble coordsets in this mpi
		cs[iSub] = &subdomain[iSub]->getNodes();

		// assemble element sets in this mpi
		elemSet[iSub] = &subdomain[iSub]->getElementSet();

		// assemble constrained and unconstrained dofset arrays in this mpi
		cdsa[iSub] = subdomain[iSub]->getCDSA();
		dsa[iSub] = subdomain[iSub]->getDSA();

		// allocate and initialize for the user defined disps and vels
		int numDofs = dsa[iSub]->size();
		usrDefDisps[iSub] = new double[numDofs];
		usrDefVels[iSub] = new double[numDofs];

		for (int iDof = 0; iDof < numDofs; iDof++) {
			usrDefDisps[iSub][iDof] = 0.0;
			usrDefVels[iSub][iDof] = 0.0;
		}
	}

	int numOutInfo = geoSource->getNumOutInfo();
	OutputInfo *oinfo = geoSource->getOutputInfo();

	int flag = 0;

	// Check if aero forces are requested for output
	int iInfo;
	for (iInfo = 0; iInfo < numOutInfo; ++iInfo) {
		if (oinfo[iInfo].type == OutputInfo::AeroForce) {
			flag = 1;
			break;
		}
	}

	// create distributed fluid exchanger
	OutputInfo *oinfo_aero = (flag) ? oinfo + iInfo : NULL;
	std::set<int> &aeroEmbeddedSurfaceId = domain->GetAeroEmbedSurfaceId();
	if (aeroEmbeddedSurfaceId.size() != 0) {
		int iSurf = -1;
		for (int i = 0; i < domain->getNumSurfs(); i++)
			if (aeroEmbeddedSurfaceId.find((*domain->viewSurfEntities())[i]->ID()) != aeroEmbeddedSurfaceId.end()) {
				iSurf = i;
				break; // only allows one surface.
			}
		if (iSurf < 0) {
			fprintf(stderr, "ERROR: Embedded wet surface not found! Aborting...\n");
			exit(-1);
		}
		distFlExchanger = new DistFlExchanger(cs, elemSet, (*domain->viewSurfEntities())[iSurf],
		                                      &domain->getNodes(), domain->getNodeToElem(),
		                                      decDomain->getElemToSub().get(), subdomain,
		                                      cdsa, dsa, oinfo_aero, domain->solInfo().elementDeletion);
	} else {
		if (domain->solInfo().elementDeletion) {
			filePrint(stderr, " *** WARNING: The C0 algorithm and an embedded surface id must be specified\n"
			                  "     under AERO for an aeroelastic analysis with element deletion, otherwise\n"
			                  "     Aero-F will not be notified of any topological changes in structure.\n");
		}
		distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa, oinfo_aero);
	}
	mddPostPro->setPostProcessor(distFlExchanger);
	mddPostPro->setUserDefs(usrDefDisps, usrDefVels);

	// negotiate with the fluid code
	distFlExchanger->negotiate();

	int restartinc = (sinfo.nRestart >= 0) ? (sinfo.nRestart) : 0;

	DistrVector dispAero(disp);

	if (sinfo.gepsFlg == 1) {
		// If we are in the first time step, and we initialized with
		// IDISP6, do not send IDISP6
		if (domain->numInitDisp() == 0 && sinfo.zeroInitialDisp != 1) {
			filePrint(stderr, " ... DO NOT SEND IDISP6             ...\n");
		} else {
			filePrint(stderr, " ... SENDING IDISP6                 ...\n");
			for (iSub = 0; iSub < numLocSub; iSub++) {
				BCond *iDis6 = subdomain[iSub]->getInitDisp6();
				for (int i = 0; i < subdomain[iSub]->numInitDisp6(); ++i) {
					int dof = cdsa[iSub]->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
					if (dof >= 0)
						dispAero[dof] += iDis6[i].val;
				}
			}
		}
	}

	SysState<DistrVector> state(dispAero, vel, accel, lastVel);

	if (sinfo.aeroFlag == 8) {
		distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.mppFactor,
		                           restartinc, sinfo.isCollocated, sinfo.alphas);
		distFlExchanger->sendModeFreq(modeData.frequencies, modeData.numModes);
		if (verboseFlag) filePrint(stderr, " ... [E] Sent parameters and mode frequencies ...\n");
		distFlExchanger->sendModeShapes(modeData.numModes, modeData.numNodes,
		                                modeData.modes, state, sinfo.mppFactor);
		if (verboseFlag) filePrint(stderr, " ... [E] Sent mode shapes           ...\n");
	} else {
		double aero_tmax = sinfo.tmax;
		if (sinfo.newmarkBeta == 0 && !sinfo.dyna3d_compat) aero_tmax += sinfo.getTimeStep();
		double aero_dt = (sinfo.dyna3d_compat) ? 0 : sinfo.getTimeStep();
		distFlExchanger->sendParam(sinfo.aeroFlag, aero_dt, aero_tmax, restartinc,
		                           sinfo.isCollocated, sinfo.alphas);
		if (verboseFlag) filePrint(stderr, " ... [E] Sent parameters            ...\n");

		// initialize the Parity
		if (sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
			distFlExchanger->initRcvParity(1);
			distFlExchanger->initSndParity(1);
		} else {
			distFlExchanger->initRcvParity(-1);
			distFlExchanger->initSndParity(-1);
		}

		if (sinfo.aeroFlag == 20 && sinfo.dyna3d_compat) {
			distFlExchanger->sendSubcyclingInfo(0);
			distFlExchanger->sendNoStructure();
		}

		// send initial displacements
    if(!(geoSource->getCheckFileInfo()->hotRestart() && sinfo.aeroFlag == 20)) {
			distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
			if (verboseFlag) filePrint(stderr, " ... [E] Sent initial displacements ...\n");
		}

		if (sinfo.aeroFlag == 1) { // Ping pong only
			filePrint(stderr, "Ping Pong Only requested. Structure code exiting\n");
		}
	}

	return sinfo.aeroFlag;
}

int
MultiDomainDynam::cmdCom(int cmdFlag) {
	return distFlExchanger->cmdCom(cmdFlag);
}

int
MultiDomainDynam::getAeroAlg() {
	return domain->solInfo().aeroFlag;
}

void
MultiDomainDynam::aeroSend(double time, DistrVector &d_n, DistrVector &v_n, DistrVector &a_n, DistrVector &v_p) {
	startTimerMemory(times->output, times->memoryOutput);

	domain->getTimers().sendFluidTime -= getTime();
	SysState<DistrVector> state(d_n, v_n, a_n, v_p);

	if (claw && userSupFunc) {
		if (claw->numUserDisp) { // USDD
			// Note: the approprate value of "time" passed into this function should be t^{n+½} for A6 and C0, and
			// t^{n+1} otherwise, where t^n denotes the time at the end of the current structure timestep. Note that
			// the predictor in FlExchanger::sendDisplacements is not applied to prescribed displacements; we directly
			// compute here the desired values of the prescribed displacements/velocities rather than predicting them.
			double *userDefineDisp = new double[claw->numUserDisp];
			double *userDefineVel = new double[claw->numUserDisp];
			double *userDefineAcc = new double[claw->numUserDisp];
			for (int i = 0; i < claw->numUserDisp; ++i) {
				userDefineVel[i] = 0;
				userDefineAcc[i] = 0;
			}
			userSupFunc->usd_disp(time, userDefineDisp, userDefineVel, userDefineAcc);
			// update usrDefDisps, usrDefVels
			execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subUpdateUsrDefDispsAndVels, userDefineDisp,
			          userDefineVel);
			delete[] userDefineDisp;
			delete[] userDefineVel;
			delete[] userDefineAcc;
		}
	}

	if (domain->solInfo().dyna3d_compat) {
		std::set<int> &aeroEmbeddedSurfaceId = domain->GetAeroEmbedSurfaceId();
		int numDeletedElements = 0;
		if (aeroEmbeddedSurfaceId.size() != 0 && domain->solInfo().elementDeletion) {
			for (int i = 0; i < decDomain->getNumSub(); ++i) {
				numDeletedElements += decDomain->getSubDomain(i)->getNewDeletedElements().size();
			}
#ifdef DISTRIBUTED
			numDeletedElements = structCom->globalSum(numDeletedElements);
#endif
		}
		if (numDeletedElements > 0) {
			distFlExchanger->sendNewStructure();
		} else {
			distFlExchanger->sendNoStructure();
		}
	}

	distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
	domain->getTimers().sendFluidTime += getTime();
	if (verboseFlag) filePrint(stderr, " ... [E] Sent displacements         ...\n");

	stopTimerMemory(times->output, times->memoryOutput);
}

void
MultiDomainDynam::a5TimeLoopCheck(int &parity, double &t, double dt) {
	if (domain->solInfo().aeroFlag == 5) {
		if (!parity) t -= dt;
		parity = (parity ? 0 : 1);
	}
}

void
MultiDomainDynam::a5StatusRevise(int parity, SysState<DistrVector> &curState, SysState<DistrVector> &bkState) {
	if (domain->solInfo().aeroFlag == 5) {
		if (parity) { // restore
			*prevFrc = *prevFrcBackup;
			prevIndex = prevIndexBackup;
			prevTime = prevTimeBackup;
			curState.getDisp() = bkState.getDisp();
			curState.getVeloc() = bkState.getVeloc();
			curState.getAccel() = bkState.getAccel();
			curState.getPrevVeloc() = bkState.getPrevVeloc();
		} else { // backup
			*prevFrcBackup = *prevFrc;
			prevIndexBackup = prevIndex;
			prevTimeBackup = prevTime;
			bkState.getDisp() = curState.getDisp();
			bkState.getVeloc() = curState.getVeloc();
			bkState.getAccel() = curState.getAccel();
			bkState.getPrevVeloc() = curState.getPrevVeloc();
		}
	}
}

void
MultiDomainDynam::thermoePreProcess(DistrVector &, DistrVector &, DistrVector &) {
	if (domain->solInfo().thermoeFlag >= 0) {

		auto cpuToSub = geoSource->getCpuToSub();
		int myId = structCom->myID();
		int numLocSub = cpuToSub->num(myId);
		SubDomain **subdomain = decDomain->getAllSubDomains();

		// if sinfo.aeroFlag >= 0, flExchanger has already been initialize before,
		// thus, only when sinfo.aeroFlag < 0 is necessary.
		if (domain->solInfo().aeroFlag < 0) {

			// allocate for pointer arrays
			CoordSet **cs = new CoordSet *[numLocSub];
			Elemset **elemSet = new Elemset *[numLocSub];
			DofSetArray **cdsa = new DofSetArray *[numLocSub];
			DofSetArray **dsa = new DofSetArray *[numLocSub];
			usrDefDisps = new double *[numLocSub];
			usrDefVels = new double *[numLocSub];

			int iSub;
			for (iSub = 0; iSub < numLocSub; iSub++) {

				// assemble coordsets in this mpi
				cs[iSub] = &subdomain[iSub]->getNodes();

				// assemble element sets in this mpi
				elemSet[iSub] = &subdomain[iSub]->getElementSet();

				// assemble constrained and unconstrained dofset arrays in this mpi
				cdsa[iSub] = subdomain[iSub]->getCDSA();
				dsa[iSub] = subdomain[iSub]->getDSA();

				// allocate and initialize for the user defined disps and vels
				int numDofs = dsa[iSub]->size();
				usrDefDisps[iSub] = new double[numDofs];
				usrDefVels[iSub] = new double[numDofs];

				for (int iDof = 0; iDof < numDofs; iDof++) {
					usrDefDisps[iSub][iDof] = 0.0;
					usrDefVels[iSub][iDof] = 0.0;
				}
			}

			distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);
			mddPostPro->setPostProcessor(distFlExchanger);
			mddPostPro->setUserDefs(usrDefDisps, usrDefVels);
		}

		nodalTemps = new DistrVector(decDomain->ndVecInfo());
		for (int iSub = 0; iSub < numLocSub; iSub++) subdomain[iSub]->temprcvd = nodalTemps->subData(iSub);
		int buffLen = nodalTemps->size();

		distFlExchanger->thermoread(buffLen);

		distFlExchanger->getStrucTemp(nodalTemps->data());
		if (verboseFlag) filePrint(stderr, " ... [E] Received initial temperatures ...\n");
		if (geomState) geomState->setNodalTemperatures(*nodalTemps);
	}
}

void
MultiDomainDynam::thermohPreProcess(DistrVector &d, DistrVector &, DistrVector &) {
	if (domain->solInfo().thermohFlag >= 0) {

		auto cpuToSub = geoSource->getCpuToSub();
		int myId = structCom->myID();
		int numLocSub = cpuToSub->num(myId);
		SubDomain **subdomain = decDomain->getAllSubDomains();

		// if sinfo.aeroheatFlag >= 0, flExchanger has already been initialize before,
		// thus, only when sinfo.aeroheatFlag < 0 is necessary.
		if (domain->solInfo().aeroheatFlag < 0) {

			// allocate for pointer arrays
			CoordSet **cs = new CoordSet *[numLocSub];
			Elemset **elemSet = new Elemset *[numLocSub];
			DofSetArray **cdsa = new DofSetArray *[numLocSub];
			DofSetArray **dsa = new DofSetArray *[numLocSub];
			usrDefDisps = new double *[numLocSub];
			usrDefVels = new double *[numLocSub];

			int iSub;
			for (iSub = 0; iSub < numLocSub; iSub++) {

				// assemble coordsets in this mpi
				cs[iSub] = &subdomain[iSub]->getNodes();

				// assemble element sets in this mpi
				elemSet[iSub] = &subdomain[iSub]->getElementSet();

				// assemble constrained and unconstrained dofset arrays in this mpi
				cdsa[iSub] = subdomain[iSub]->getCDSA();
				dsa[iSub] = subdomain[iSub]->getDSA();

				// allocate and initialize for the user defined disps and vels
				int numDofs = dsa[iSub]->size();
				usrDefDisps[iSub] = new double[numDofs];
				usrDefVels[iSub] = new double[numDofs];

				for (int iDof = 0; iDof < numDofs; iDof++) {
					usrDefDisps[iSub][iDof] = 0.0;
					usrDefVels[iSub][iDof] = 0.0;
				}
			}

			distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);
			mddPostPro->setPostProcessor(distFlExchanger);
			mddPostPro->setUserDefs(usrDefDisps, usrDefVels);
		}

		nodalTemps = new DistrVector(decDomain->ndVecInfo());
		int buffLen = nodalTemps->size();
		mddPostPro->setNodalTemps(nodalTemps);

		distFlExchanger->thermoread(buffLen);

		for (int i = 0; i < decDomain->getNumSub(); ++i) {
			SubDomain *sd = decDomain->getSubDomain(i);
			for (int j = 0; j < sd->numNodes(); ++j) {
				int tloc = sd->getCDSA()->locate(j, DofSet::Temp);
				int tloc1 = sd->getDSA()->locate(j, DofSet::Temp);
				double temp = (tloc >= 0) ? d.subData(i)[tloc] : sd->getBcx()[tloc1];
				if (tloc1 < 0) temp = 0.0;
				nodalTemps->subData(i)[j] = temp;
			}
		}

		distFlExchanger->sendStrucTemp(*nodalTemps);
		if (verboseFlag) filePrint(stderr, " ... [T] Sent initial temperatures  ...\n");
	}

}

int
MultiDomainDynam::getThermoeFlag() {
	return domain->solInfo().thermoeFlag;
}

int
MultiDomainDynam::getThermohFlag() {
	return domain->solInfo().thermohFlag;
}

void
MultiDomainDynam::aeroHeatPreProcess(DistrVector &disp, DistrVector &vel, DistrVector &lastVel) {
	// get solver info
	SolverInfo &sinfo = domain->solInfo();

	// Initialize previous force data
	// Reexamine for the case of restart
	prevFrc = new DistrVector(solVecInfo());
	prevFrc->zero();
	prevFrcBackup = new DistrVector(solVecInfo());
	prevIndex = -1;
	prevTime = 0;

	// Initialize the aeroforce vector
	aeroForce = new DistrVector(solVecInfo());

	if (sinfo.aeroheatFlag < 0)
		return;

	auto cpuToSub = geoSource->getCpuToSub();

	// get cpu id
#ifdef USE_MPI
	int myId = structCom->myID();
#else
	int myId = 0;
#endif

	int numLocSub = cpuToSub->num(myId);

	SubDomain **subdomain = decDomain->getAllSubDomains();

	// allocate for pointer arrays
	CoordSet **cs = new CoordSet *[numLocSub];
	Elemset **elemSet = new Elemset *[numLocSub];
	DofSetArray **cdsa = new DofSetArray *[numLocSub];
	DofSetArray **dsa = new DofSetArray *[numLocSub];
	usrDefDisps = new double *[numLocSub];
	usrDefVels = new double *[numLocSub];

	int iSub;
	for (iSub = 0; iSub < numLocSub; iSub++) {

		// assemble coordsets in this mpi
		cs[iSub] = &subdomain[iSub]->getNodes();

		// assemble element sets in this mpi
		elemSet[iSub] = &subdomain[iSub]->getElementSet();

		// assemble constrained and unconstrained dofset arrays in this mpi
		cdsa[iSub] = subdomain[iSub]->getCDSA();
		dsa[iSub] = subdomain[iSub]->getDSA();

		// allocate and initialize for the user defined disps
		int numDofs = dsa[iSub]->size();
		usrDefDisps[iSub] = new double[numDofs];

		for (int iDof = 0; iDof < numDofs; iDof++) {
			usrDefDisps[iSub][iDof] = 0.0;
		}
	}

	int numOutInfo = geoSource->getNumOutInfo();
	OutputInfo *oinfo = geoSource->getOutputInfo();

	int flag = 0;

	// Check if aero fluxes are requested for output
	int iInfo;
	for (iInfo = 0; iInfo < numOutInfo; ++iInfo) {
		if (oinfo[iInfo].type == OutputInfo::AeroForce) {
			flag = 1;
			break;
		}
	}

	// create distributed fluid exchanger
	if (flag)
		distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa, oinfo + iInfo);
	else
		distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);

	mddPostPro->setPostProcessor(distFlExchanger);
	mddPostPro->setUserDefs(usrDefDisps, usrDefVels);

	// negotiate with the fluid code
	distFlExchanger->negotiate();

	int restartinc = (sinfo.nRestart >= 0) ? (sinfo.nRestart) : 0;

	SysState<DistrVector> state(disp, vel, lastVel);

	distFlExchanger->sendTempParam(sinfo.aeroheatFlag, sinfo.getTimeStep(), sinfo.tmax, restartinc,
	                               sinfo.alphat);
	if (verboseFlag) filePrint(stderr, " ... [T] Sent parameters            ...\n");

	// send initial displacements
	distFlExchanger->sendTemperature(state);
	if (verboseFlag) filePrint(stderr, " ... [T] Sent initial temperatures  ...\n");
}

int
MultiDomainDynam::getAeroheatFlag() {
	return domain->solInfo().aeroheatFlag;
}

void
MultiDomainDynam::getGravityForce(DistrVector &f) {
	execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subGetGravityForce, f);
}

void
MultiDomainDynam::subGetGravityForce(int isub, DistrVector &f) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	StackVector subf(f.subData(isub), f.subLen(isub));
	subf.zero();
	sd->addGravityForce(subf);
}

void
MultiDomainDynam::getUnamplifiedExtForce(DistrVector &f, int loadcaseid) {
	execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subGetUnamplifiedExtForce, f, loadcaseid);
}

void
MultiDomainDynam::subGetUnamplifiedExtForce(int isub, DistrVector &f, int loadcaseid) {
	SubDomain *sd = decDomain->getSubDomain(isub);
	StackVector subf(f.subData(isub), f.subLen(isub));
	sd->computeUnamplifiedExtForce(subf, loadcaseid);
}

void
MultiDomainDynam::getAeroelasticForceSensitivity(int, double, DistrVector *, double, double) {
	filePrint(stderr, " ... MultiDomainDynam::getAeroelasticForceSensitivity is not implemented\n");
	exit(-1);
}

void
MultiDomainDynam::preProcessSA() {
	filePrint(stderr, " ... MultiDomainDynam::preProcessSA is not implemented\n");
	exit(-1);
}

void
MultiDomainDynam::postProcessSA(MDDynamMat *, DistrVector &) {
	filePrint(stderr, " ... MultiDomainDynam::postProcessSA is not implemented\n");
	exit(-1);
}

void
MultiDomainDynam::sensitivityPostProcessing(DistrVector *) {
	filePrint(stderr, " ... MultiDomainDynam::sensitivityPostProcessing is not implemented\n");
	exit(-1);
}

#include <Paral.d/MDNLQStatic.h>
#include <Driver.d/NLStaticProbType.h>

void
MultiDomainDynam::solveAndUpdate(DistrVector &force, DistrVector &dinc, DistrVector &d, double relaxFac, double time) {
	int numElemStates = geomState->getTotalNumElemStates();
	if (!refState) {
		// For the first coupling cycle refState is the initial state as defined by either IDISP or restart, if specified.
		refState = new DistrGeomState(*geomState);
	} else if (numElemStates == 0) {
		// In this case dlambda is only used for the first cycle.
		domain->solInfo().getNLInfo().dlambda = domain->solInfo().getNLInfo().maxLambda = 1.0;
		for (int i = 0; i < decDomain->getNumSub(); ++i)
			decDomain->getSubDomain(i)->solInfo().getNLInfo().dlambda = decDomain->getSubDomain(
				i)->solInfo().getNLInfo().maxLambda = 1.0;
	}

	MDNLQStatic nlstatic(domain, decDomain, force, refState);
	NLStaticSolver<ParallelSolver, DistrVector, MultiDomainPostProcessor, MDNLQStatic, DistrGeomState> nlsolver(
		&nlstatic);
	nlsolver.solve();

	nlsolver.getGeomState()->get_inc_displacement(dinc, *geomState, false);
	if (numElemStates == 0) {
		// In this case refState now stores the solution of the previous non-linear solve, and will be used as the initial
		// guess for the next coupling cycle's non-linear solve.
		*refState = *nlsolver.getGeomState();
	}

	dinc *= relaxFac;
	geomState->update(dinc);
	if (numElemStates != 0) {
		execParal(decDomain->getNumSub(), this, &MultiDomainDynam::subUpdateStates, refState, geomState, time);
	}

	geomState->get_tot_displacement(d, false);
}

