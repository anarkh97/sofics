#include <Paral.d/DomainGroupTask.h>
#include <iostream>
#include <Driver.d/SysState.h>
#include <Paral.d/MDDynam.h>
#include <Threads.d/Paral.h>
#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>
#include <Driver.d/DecDomain.h>
#include <Driver.d/SubDomain.h>
#include <Paral.d/MDOp.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/Vector.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/DiagMatrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Timers.d/GetTime.h>
#include <Control.d/ControlInterface.h>
#include <Threads.d/PHelper.h>
#include <Paral.d/GenMS.h>
#include <Solvers.d/SolverFactory.h>
#ifdef DISTRIBUTED
#include <Utils.d/DistHelper.h>
#endif

class IntFullM;

extern SolverInfo &solInfo;
extern GeoSource *geoSource;

template<class Scalar>
GenDomainGroupTask<Scalar>::GenDomainGroupTask(int _nsub, GenSubDomain<Scalar> **_sd, double _cm,
                                               double _cc, double _ck, Rbm **_rbms, FullSquareMatrix **_kelArray,
                                               double _alpha, double _beta, int _numSommer, int _solvertype,
                                               FSCommunicator *_com, FullSquareMatrix **_melArray,
                                               FullSquareMatrix **_celArray, MatrixTimers &_mt)
	: mt(_mt)
{
	nsub = _nsub;
	sd = _sd;
	dynMats = new GenSolver<Scalar> *[nsub];
	spMats  = new GenSparseMatrix<Scalar> *[nsub];
	M    = new GenSparseMatrix<Scalar> *[nsub];
	Muc  = new GenSparseMatrix<Scalar> *[nsub];
	Mcc  = new GenSparseMatrix<Scalar> *[nsub];
	C    = new GenSparseMatrix<Scalar> *[nsub];
	Cuc  = new GenSparseMatrix<Scalar> *[nsub];
	Ccc  = new GenSparseMatrix<Scalar> *[nsub];
	C_deriv    = new GenSparseMatrix<Scalar> **[nsub];
	Cuc_deriv    = new GenSparseMatrix<Scalar> **[nsub];
	K_deriv    = new GenSparseMatrix<Scalar> **[nsub];
	Kuc_deriv    = new GenSparseMatrix<Scalar> **[nsub];
	num_K_deriv = 0;
	K_arubber_l = new GenSparseMatrix<Scalar> **[nsub];
	K_arubber_m = new GenSparseMatrix<Scalar> **[nsub];
	Kuc_arubber_l = new GenSparseMatrix<Scalar> **[nsub];
	Kuc_arubber_m = new GenSparseMatrix<Scalar> **[nsub];
	num_K_arubber = 0;
	K    = new GenSparseMatrix<Scalar> *[nsub];
	if(domain->solInfo().solvercntl->precond) {
		spp = new GenSparseMatrix<Scalar> *[nsub];
		sps = new GenSolver<Scalar> *[nsub];
	}
	else {
		spp = 0;
		sps = 0;
	}
	rbms = _rbms;
	kelArray = _kelArray;
	melArray = _melArray;
	celArray = _celArray;
	Kuc  = new GenSparseMatrix<Scalar> *[nsub];
	coeM    = _cm;
	coeC    = _cc;
	coeK    = _ck;
	numSommer = _numSommer;
	alpha   = _alpha;
	beta    = _beta;
	solvertype = _solvertype;
	com = _com;
	makeC = (alpha != 0.0 || beta != 0.0 || (numSommer > 0) || domain->getElementSet().hasDamping());
// RT - 053013 - to enable multiple impedance section, build C_deriv whenever C
//  makeC_deriv = (makeC && solInfo.doFreqSweep && solInfo.getSweepParams()->nFreqSweepRHS > 1);
	makeC_deriv = (makeC && solInfo.doFreqSweep);
}

template<class Scalar>
GenDomainGroupTask<Scalar>::~GenDomainGroupTask()
{
	// delete [] dynMats;
	//delete [] spMats;
	//delete [] rbms;
	// don't delete K,Kuc,C,Cuc,M,Muc
}

template<class Scalar>
void
GenDomainGroupTask<Scalar>::runForWB(int isub, bool make_feti)
{
	mt.constructTime -= getTime();
	DofSetArray     *dsa = sd[isub]->getDSA();
	ConstrainedDSA *cdsa = sd[isub]->getCDSA();

	K[isub] = 0;
	Kuc[isub] = 0;
	M[isub] = 0;
	Muc[isub] = 0;
	Mcc[isub] = 0;
	C[isub] = 0;
	Cuc[isub] = 0;
	Ccc[isub] = 0;
	C_deriv[isub] = 0;
	Cuc_deriv[isub] = 0;
	K_deriv[isub] = 0;
	Kuc_deriv[isub] = 0;
	K_arubber_l[isub] = 0;
	Kuc_arubber_l[isub] = 0;
	K_arubber_m[isub] = 0;
	Kuc_arubber_m[isub] = 0;
	if(spp) spp[isub] = 0;
	if(sps) sps[isub] = 0;

	if((cdsa->size() - dsa->size()) != 0)
		Kuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();

	if(solInfo.isDynam() || solInfo.doFreqSweep || solInfo.probType == SolverInfo::Modal || solInfo.probType == SolverInfo::PodRomOffline) {

		// XML Need to introduce Mcc
		if(solInfo.isCoupled)
			K[isub] = sd[isub]->template constructNBSparseMatrix<Scalar>(); // unsymmetric
		else
			K[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();

		if(solInfo.newmarkBeta == 0.0) { // explict dynamics
			if(solvertype != 10) {
				int numN = sd[isub]->numNodes();
				// import a diagonal connectivity for the mass matrix
				Connectivity connForMass(numN);
				connForMass = connForMass.modify();
				M[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>(cdsa, &connForMass);
			}
			else {
				// if only spectral elements are used, the mass matrix is purely diagonal for implicit and explicit
				M[isub] = new GenDiagMatrix<Scalar>(cdsa);
			}
		}
		else {
			if(solInfo.isCoupled)
				M[isub] = sd[isub]->template constructNBSparseMatrix<Scalar>();
			else
				M[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();
		}

		if((cdsa->size() - dsa->size()) != 0) {
			Muc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
			Mcc[isub] = sd[isub]->template constructCCSparse<Scalar>();
		}

		if(makeC) {
			C[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();

			if((cdsa->size() - dsa->size()) != 0) {
				Cuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
				Ccc[isub] = sd[isub]->template constructCCSparse<Scalar>();
			}

			if(makeC_deriv) {
				int numC_deriv, numRHS;
				numRHS = solInfo.getSweepParams()->nFreqSweepRHS;
				if((numSommer > 0) && ((sd[isub]->sommerfeldType == 2) || (sd[isub]->sommerfeldType == 4)))
					numC_deriv = numRHS - 1;
				else
					numC_deriv = 1;
				C_deriv[isub] = new GenSparseMatrix<Scalar> * [numRHS - 1];
				for(int n = 0; n < numC_deriv; ++n) {
					C_deriv[isub][n] = sd[isub]->template constructDBSparseMatrix<Scalar>();
				}
				for(int n = numC_deriv; n < numRHS - 1; ++n)
					C_deriv[isub][n] = 0;
				if(cdsa->size() > 0 && (cdsa->size() - dsa->size()) != 0) {
					Cuc_deriv[isub] = new GenSparseMatrix<Scalar> * [numRHS - 1];
					for(int n = 0; n < numC_deriv; ++n) {
						Cuc_deriv[isub][n] = sd[isub]->template constructCuCSparse<Scalar>();
					}
					for(int n = numC_deriv; n < numRHS - 1; ++n)
						Cuc_deriv[isub][n] = 0;
				}
				num_K_deriv = solInfo.doFreqSweep ? solInfo.getSweepParams()->nFreqSweepRHS-1:0;
				K_deriv[isub] = new GenSparseMatrix<Scalar> * [num_K_deriv+1];
				for(int n = 0; n <= num_K_deriv; ++n) {
					K_deriv[isub][n] = sd[isub]->template constructDBSparseMatrix<Scalar>();
				}
				if(cdsa->size() > 0 && (cdsa->size() - dsa->size()) != 0) {
					Kuc_deriv[isub] = new GenSparseMatrix<Scalar> * [num_K_deriv+1];
					for(int n = 0; n <= num_K_deriv; ++n) {
						Kuc_deriv[isub][n] = sd[isub]->template constructCuCSparse<Scalar>();
					}
				}
				num_K_arubber = geoSource->num_arubber;
				K_arubber_l[isub] = new GenSparseMatrix<Scalar> * [num_K_arubber];
				K_arubber_m[isub] = new GenSparseMatrix<Scalar> * [num_K_arubber];
				for(int n = 0; n < num_K_arubber; ++n) {
					K_arubber_l[isub][n] =
						sd[isub]->template constructDBSparseMatrix<Scalar>();
					K_arubber_m[isub][n] =
						sd[isub]->template constructDBSparseMatrix<Scalar>();
				}
				if(cdsa->size() > 0 && (cdsa->size() - dsa->size()) != 0) {
					Kuc_arubber_l[isub] = new GenSparseMatrix<Scalar> * [num_K_arubber];
					Kuc_arubber_m[isub] = new GenSparseMatrix<Scalar> * [num_K_arubber];
					for(int n = 0; n < num_K_arubber; ++n) {
						Kuc_arubber_l[isub][n] = sd[isub]->template constructCuCSparse<Scalar>();
						Kuc_arubber_m[isub][n] = sd[isub]->template constructCuCSparse<Scalar>();
					}
				}
			}
		}
		else if(solInfo.ATDARBFlag != -2.0) { // for acoustic damping
			if(solvertype != 10) {
				// build the modified damping matrix for implicit and explicit
				C[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>(cdsa, sd[isub]->getNodeToNode_sommer());
			}
			else
				C[isub] = new GenDiagMatrix<Scalar>(cdsa);
			if((cdsa->size() - dsa->size()) != 0)
				Cuc[isub] = sd[isub]->template constructCuCSparse<Scalar>();
		}
	}
	else if(solInfo.filterQ == 0 && (solInfo.filterFlags || solInfo.hzemFilterFlag || solInfo.slzemFilterFlag)) {
		M[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();
	}

	// builds the datastructures for Kii, Kib, Kbb
	if(isFeti(domain->solInfo().solvercntl->type) && make_feti) { // FETI
		if(sd[isub]->numMPCs() > 0)
			sd[isub]->makeKbbMpc();
		else {
			sd[isub]->makeKbb(sd[isub]->getCCDSA());
		}
	}

	GenMultiSparse<Scalar> *allMats = 0;
	if(make_feti) {
		if(isFeti(domain->solInfo().solvercntl->type) || domain->solInfo().solvercntl->type == SolverSelection::BlockDiag) {
			// construct local solver for the subdomain
			SolverCntl &local_cntl = *domain->solInfo().solvercntl->fetiInfo.local_cntl;
			if(local_cntl.type == SolverSelection::Direct || local_cntl.type == SolverSelection::Iterative) {
				dynMats[isub] = GenSolverFactory<Scalar>::getFactory()
					->createSolver(sd[isub]->getNodeToNode(), sd[isub]->getDSA(), sd[isub]->getCCDSA(),
					               local_cntl, spMats[isub], (Rbm*) NULL, spp[isub], sps[isub]);
				// Wrap the sparse matrix.
				spMats[isub] = new SolverWrapper<Scalar>(sd[isub]->getNodeToNode(), sd[isub]->getDSA(), sd[isub]->getCCDSA(),
					spMats[isub]);
			}
			else if(isFeti(local_cntl.type)) { // local solver is feti
				std::cerr << "using FETI-DP solver for local problem\n";
				BCond *corner_dbc = new BCond[3*sd[isub]->numCorners()];
				for(int i = 0; i < sd[isub]->numCorners(); ++i)
					for(int j = 0; j < 3; ++j)
						corner_dbc[3*i+j].setData(sd[isub]->getLocalCornerNodes()[i], j, 0.0);
				sd[isub]->setDirichlet(3*sd[isub]->numCorners(), corner_dbc);
#ifdef USE_MPI
				GenDecDomain<Scalar> *decSubDomain =
					new GenDecDomain<Scalar>(sd[isub],
					                         new Communicator(CommunicatorHandle{(MPI_Comm)MPI_COMM_SELF}));
#else
				GenDecDomain<Scalar> *decSubDomain = new GenDecDomain<Scalar>(sd[isub], NULL);
#endif
				decSubDomain->preProcess();
				GenMDDynamMat<Scalar> ops;
				decSubDomain->buildOps(ops, 0.0, 0.0, 1.0);
				exit(-1); // XXX this is not finished yet
				//dynMats[isub] = ops.dynMat;
				//((GenFetiDPSolver<Scalar> *) Krr)->initL(cc_dsa);
			}
		}
    else if(domain->solInfo().solvercntl->type == SolverSelection::Iterative||
            (domain->solInfo().solvercntl->type == SolverSelection::Direct && domain->solInfo().solvercntl->subtype == 13)) {
			dynMats[isub] = 0;
			switch(domain->solInfo().solvercntl->iterSubtype) {
				case 2 :
					spMats[isub] = sd[isub]->template constructNBSparseMatrix<Scalar>();
					break;
				default:
				case 3 :
					spMats[isub] = sd[isub]->template constructDBSparseMatrix<Scalar>();
					break;
#ifdef USE_EIGEN3
				case 4:
					spMats[isub] = sd[isub]->template constructEiSparseMatrix<Scalar,Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >();
					break;
#endif
			}
			if(domain->solInfo().solvercntl->precond == 1) {
				GenDiagMatrix<Scalar> *dm = new GenDiagMatrix<Scalar>(sd[isub]->getCDSA());
				spp[isub] = (GenSparseMatrix<Scalar>*) dm;
				sps[isub] = (GenSolver<Scalar>*) dm;
			}
		}
		else if(domain->solInfo().solvercntl->type == static_cast<SolverSelection>(-1)) {
			dynMats[isub] = nullptr;
			spMats[isub] = nullptr;
		}

		if(isFeti(domain->solInfo().solvercntl->type) && domain->solInfo().getFetiInfo().version == FetiInfo::fetidp) {
			sd[isub]->constructKcc();
			sd[isub]->constructKrc();
		}

		if(isFeti(domain->solInfo().solvercntl->type)) {
			if(geoSource->isShifted() && domain->solInfo().getFetiInfo().prectype == FetiInfo::nonshifted)
				allMats = new GenMultiSparse<Scalar>(spMats[isub], sd[isub]->Krc.get(), sd[isub]->Kcc.get());
			else
				allMats = new GenMultiSparse<Scalar>(spMats[isub], sd[isub]->KiiSparse.get(), sd[isub]->Kbb.get(),
				                                     sd[isub]->Kib.get(), sd[isub]->Krc.get(), sd[isub]->Kcc.get());
		}
	}

	AllOps<Scalar> allOps;

	if(geoSource->isShifted() && solInfo.getFetiInfo().prectype == FetiInfo::nonshifted)
		allOps.K = new GenMultiSparse<Scalar>(K[isub], sd[isub]->KiiSparse.get(), sd[isub]->Kbb.get(), sd[isub]->Kib.get());
	else
		allOps.K = K[isub];
	allOps.C = C[isub];
	allOps.Cuc = Cuc[isub];
	allOps.Ccc = Ccc[isub];
	allOps.M = M[isub];
	allOps.Muc = Muc[isub];
	allOps.Mcc = Mcc[isub];
	allOps.Kuc = Kuc[isub];
	allOps.C_deriv = C_deriv[isub];
	allOps.Cuc_deriv = Cuc_deriv[isub];
	allOps.K_deriv = K_deriv[isub];
	allOps.Kuc_deriv = Kuc_deriv[isub];
	allOps.n_Kderiv = num_K_deriv;
	allOps.K_arubber_l = K_arubber_l[isub];
	allOps.K_arubber_m = K_arubber_m[isub];
	allOps.Kuc_arubber_l = Kuc_arubber_l[isub];
	allOps.Kuc_arubber_m = Kuc_arubber_m[isub];
	allOps.num_K_arubber = num_K_arubber;
	mt.constructTime += getTime();

	allOps.spp = (spp) ? spp[isub] : 0;
	FullSquareMatrix *subKelArray = (kelArray) ? kelArray[isub] : 0;
	FullSquareMatrix *subMelArray = (melArray) ? melArray[isub] : 0;
	FullSquareMatrix *subCelArray = (celArray) ? celArray[isub] : 0;
	if(isFeti(domain->solInfo().solvercntl->type))
		sd[isub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, allMats, subKelArray, subMelArray, subCelArray);
	else
		sd[isub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, spMats[isub], subKelArray, subMelArray, subCelArray);

	if(allMats) delete allMats;
}

template class GenDomainGroupTask<double>;
template class GenDomainGroupTask<std::complex<double>>;
