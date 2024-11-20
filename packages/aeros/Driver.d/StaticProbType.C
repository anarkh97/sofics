#include <Timers.d/StaticTimers.h>
#include <Utils.d/DistHelper.h>
#include <cassert>
#include <iostream>

extern int verboseFlag;

template < class Scalar,
		class OpSolver,
		class VecType,
		class PostProcessor,
		class ProblemDescriptor,
		class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
		PostProcessor, ProblemDescriptor, ComplexVecType >
::solve()
{
	probDesc->preProcess();
#ifdef USE_EIGEN3
	if(domain->solInfo().sensitivity) {
		probDesc->preProcessSA();
		if(!domain->runSAwAnalysis) {
			AllSensitivities<Scalar> *allSens = probDesc->getAllSensitivities();
			domain->sensitivityPostProcessing(*allSens);
			return;
		}
	}
#endif

	rhs = new VecType(probDesc->solVecInfo());
	rhs->zero();
	sol = new VecType(probDesc->solVecInfo());
	sol->zero();
	VecType *savedSol = new VecType(probDesc->solVecInfo());
	VecType *savedRhs = new VecType(probDesc->solVecInfo());


	if((domain->solInfo().inpc || domain->solInfo().noninpc) && !isFeti(domain->solInfo().solvercntl->type))
		sol->setn(domain->numUncon());

	//solver=0; postProcessor=0;
	allOps=0; postProcessor=0;
	allOps = probDesc->getAllOps();
	postProcessor = probDesc->getPostProcessor();

	//solver = probDesc->getSolver();

	if(domain->solInfo().doFreqSweep) { // multiple RHS solve
		for(int is=0;is<domain->set_of_frequencies.size();is++) {

			domain->setFrequencySet(is);
			domain->solInfo().curSweepParam = is;
			domain->solInfo().setDamping(domain->solInfo().getSweepParams()->betaD,
			                             domain->solInfo().getSweepParams()->alphaD);
			SPropContainer& sProps = geoSource->getStructProps();
			double min_w = *min_element(domain->coarse_frequencies->begin(),
			                            domain->coarse_frequencies->end() );
			double max_w = *max_element(domain->coarse_frequencies->begin(),
			                            domain->coarse_frequencies->end() );
			for(int iProp=0;iProp<geoSource->getNumProps();iProp++) {
				domain->updateSDETAF(&sProps[iProp],min_w+1e-12*(max_w-min_w));
				domain->updateRUBDAFT(&sProps[iProp],min_w+1e-12*(max_w-min_w));
				// Make sure you are to the right...
			}

			int padeN = domain->solInfo().getSweepParams()->padeN;  // 1 for 1-point pade or taylor, 2 for 2-point pade, etc
			filePrint(stderr, " ... Frequency Sweep Analysis       ... \n");
			filePrint(stderr, " ... Number of coarse freqencies = %3d     ... \n", domain->coarse_frequencies->size());
			if(padeN > 1)
				filePrint(stderr, " ... Number of interpolated freqencies = %3d     ... \n", domain->frequencies->size());
			else
				filePrint(stderr, " ... Number of extrapolated freqencies = %3d     ... \n", domain->frequencies->size());
			int nRHS = domain->solInfo().getSweepParams()->nFreqSweepRHS;
			if(nRHS==0) filePrint(stderr, " *** ERROR: nRHS = 0 \n");

			Scalar *VhKV = 0, *VhMV =0, *VhCV =0,
					**VhK_arubber_lV = 0, **VhK_arubber_mV = 0;

// RT 070513
			if (domain->solInfo().getSweepParams()->freqSweepMethod==SweepParams::Taylor &&
			    domain->solInfo().getSweepParams()->nFreqSweepRHS==1 &&
			    (domain->solInfo().loadcases.size()>1 || domain->numWaveDirections>0 ) )
			{
				bool first_time = true;
				while(domain->coarse_frequencies->size() > 0) {
					double wc = domain->coarse_frequencies->front();
					domain->setSavedFreq(wc); // Not needed ??
					domain->isCoarseGridSolve = true;
					geoSource->setOmega(wc);
					if(!first_time || is>0) rebuildSolver(wc);
					else first_time = false;

					if (domain->solInfo().loadcases.size()>1) {
						std::list<int> loadcases_backup;
						loadcases_backup = domain->solInfo().loadcases;
						if(domain->solInfo().loadcases.size() > 0) {
							while(domain->solInfo().loadcases.size() > 0) {
								probDesc->getRHS(*rhs);
								allOps->sysSolver->solve(*rhs,*sol);
								if(domain->solInfo().isCoupled) scaleDisp(*sol);
								bool printTimers =
										(domain->solInfo().loadcases.size() > 1 ||
										 domain->coarse_frequencies->size() > 1 ||
										 is < domain->set_of_frequencies.size()-1 ) ? false : true;
								postProcessor->staticOutput(*sol, *rhs, printTimers);
								domain->solInfo().loadcases.pop_front();
							}
						}
						domain->solInfo().loadcases = loadcases_backup;
					}
					else if (domain->numWaveDirections>0) {
						for(int iDir=0;iDir<domain->numWaveDirections;iDir++) {
							probDesc->setIWaveDir(iDir);
							probDesc->getRHS(*rhs);
							allOps->sysSolver->solve(*rhs,*sol);
							if(domain->solInfo().isCoupled) scaleDisp(*sol);
							bool printTimers =
									( domain->coarse_frequencies->size() > 1 ||
									  iDir  <domain->numWaveDirections-1 )  ? false : true;
							postProcessor->staticOutput(*sol, *rhs,printTimers);
						}
					}
					domain->coarse_frequencies->pop_front();
				}
			} else

// RTRT
			if (domain->solInfo().getSweepParams()->isAdaptSweep &&
			    domain->solInfo().getSweepParams()->adaptSweep.deltaRHS>=0 &&
			    domain->solInfo().loadcases.size()<=1 && domain->numWaveDirections<=1)
			{
				probDesc->setIWaveDir(0);

				int maxP,numS;
				double w1,w2,atol;
				int dgp_flag;
				int minRHS,maxRHS,deltaRHS;
				maxP = domain->solInfo().getSweepParams()->adaptSweep.maxP;
				numS = domain->solInfo().getSweepParams()->adaptSweep.numS;
				dgp_flag = domain->solInfo().getSweepParams()->adaptSweep.dgp_flag;
				w1 = domain->solInfo().getSweepParams()->adaptSweep.w1;
				w2 = domain->solInfo().getSweepParams()->adaptSweep.w2;
				atol = domain->solInfo().getSweepParams()->adaptSweep.atol;
				minRHS = domain->solInfo().getSweepParams()->adaptSweep.minRHS;
				maxRHS = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
				deltaRHS = domain->solInfo().getSweepParams()->adaptSweep.deltaRHS;
				filePrint(stderr,"Adaptive derivatives %d %d %d.\n",minRHS,maxRHS,deltaRHS);
				std::vector<double> w(maxP);
				double *newW = new double[maxP];
				int maxV = maxP*maxRHS;
				int numP = 0;

				VecType **sol_prev = new VecType * [maxV+maxRHS+1];
				VecType **GP_solprev = new VecType * [maxV+maxRHS+1];
				for(int i = 0; i < maxV+maxRHS+1; ++i) {
					sol_prev[i] = new VecType(probDesc->solVecInfo());
					GP_solprev[i] = sol_prev[i];
				}
				VecType **GP_orth_solprev = new VecType * [maxV];
				for(int i = 0; i < maxV; ++i)
					GP_orth_solprev[i] = new VecType(probDesc->solVecInfo());
//       VecType **GP_solprev = new VecType * [maxV];
//       for (int ii = 0; ii < maxP; ++ii)
//         for (int iRHS = 0; iRHS < nRHS; ++iRHS)
//           GP_solprev[iRHS + ii*nRHS] = sol_prev[iRHS+1 + ii*(nRHS+1)];
				VecType **aa = new VecType * [maxV];
				VecType **bb = new VecType * [maxV];
				VecType **cc = new VecType * [maxV];
				for(int i = 0; i < maxV; ++i) {
					aa[i] = new VecType(probDesc->solVecInfo());
					bb[i] = new VecType(probDesc->solVecInfo());
					cc[i] = new VecType(probDesc->solVecInfo());
				}
				int nOrtho = 0;
//       int ncheck = 6;
				int ncheck = 18;
				double *wcheck = new double[ncheck];
				ncheck = 9;
#define ADAPT2
#ifdef ADAPT2
				for(int ic=0;ic<ncheck;ic++)
					wcheck[ic] = w1 + double(ic+1)*(w2-w1)/double(ncheck+1);
				w[numP] = w1;
				geoSource->setOmega(w[numP]);
				if (is>0) rebuildSolver(w[numP]);
				adaptGP(dgp_flag,minRHS,maxRHS,deltaRHS,nOrtho,
				        sol, GP_solprev, GP_orth_solprev,
				        aa,bb,cc,
				        VhKV, VhMV, VhCV, VhK_arubber_lV, VhK_arubber_mV,
				        ncheck, wcheck, w[numP], ((w1+w2)/2.0)/w[numP], atol);
				numP++;

				w[numP] = w2;
				geoSource->setOmega(w[numP]);
				rebuildSolver(w[numP]);
				double oldw = w[numP];
				adaptGP(dgp_flag,minRHS,maxRHS,deltaRHS,nOrtho,
				        sol, GP_solprev, GP_orth_solprev,
				        aa,bb,cc,
				        VhKV, VhMV, VhCV, VhK_arubber_lV, VhK_arubber_mV,
				        ncheck, wcheck, oldw, ((w1+w2)/2.0)/w[numP], atol);
				numP++;
#else
				w[numP] = (w1+w2)/2.0;
     for(int ic=0;ic<ncheck/2;ic++)
       wcheck[ic] = w1 +
                     double(ic+1)*(w[numP]-w1)/double(ncheck/2+1);
     for(int ic=0;ic<ncheck/2;ic++)
       wcheck[ic+ncheck/2] = w2 -
                             double(ic+1)*(w2-w[numP])/double(ncheck/2+1);
     adaptGP(dgp_flag,minRHS,maxRHS,deltaRHS,nOrtho,
             sol, GP_solprev, GP_orth_solprev,
             aa,bb,cc,
             VhKV, VhMV, VhCV, VhK_arubber_lV, VhK_arubber_mV,
             ncheck, wcheck, w[numP], ((w1+w2)/2.0)/w[numP], atol);
     numP++;
#endif

				sort(w.begin(),w.begin()+numP-1);
				while (1) {
					int numNewP = 0;
#ifdef ADAPT2
					int extrap = 0;
#else
					int extrap = 2;
#endif
					for(int i=0;i<numP-1+extrap;i++) {
						double resmax = 0.0;
						double wmax = 0.0;
//         int numres = 19;
						int numres = 9;
						double sres[numres];
						for(int j=0;j<numres;j++) {
#ifdef ADAPT2
//             double wc = (w[i]+w[i+1])/2.0 +
//               (double(j)-double(numres-1)/2.0)/double(numres)*(w[i+1]-w[i]);
							double wc = w[i] + double(j+1)/double(numres+1)*(w[i+1]-w[i]);
#else
							double wl,wr;
           if (i==0) { wl = w1; wr = w[i]; }
           else if (i==numP) { wl = w[i-1]; wr = w2; }
           else { wl = w[i-1]; wr = w[i]; }
           double wc = wl + double(j+1)/double(numres+1)*(wr-wl);
#endif

							sres[j] = adaptGPSolRes(dgp_flag,nOrtho,sol,
							                        GP_solprev, GP_orth_solprev,
							                        aa,bb,cc,
							                        VhKV, VhMV, VhCV,
							                        VhK_arubber_lV, VhK_arubber_mV,
							                        wc, wc-oldw);
							if (resmax<sres[j]) { resmax = sres[j]; wmax = wc; }
						}
						if (resmax>atol) {
							newW[numNewP] = wmax;
							numNewP++;
						}
					}
					int oldNumP = numP;
					double timex = 0.0;
					if (numNewP>0) {
						for(int i=0;i<numNewP;i++) {
							w[numP] = newW[i];
							geoSource->setOmega(w[numP]);
							timex -= getTime();
							rebuildSolver(w[numP]);
							oldw = w[numP];
							timex += getTime();
#ifdef ADAPT2
							int ii;
							for(ii=0;ii<oldNumP;ii++) if (w[ii]>w[numP]) break;
							double wl = w[ii-1];
							double wr = w[ii];
#else
							double wl;
           double wr; 
           int ii;
           if (w[numP]<w[0]) { wl = w1; wr = w[0]; }
           else if (w[numP]>w[oldNumP-1]) { wl = w[oldNumP-1]; wr = w2; }
           else {
             for(ii=0;ii<oldNumP;ii++) if (w[ii]>w[numP]) break;
             wl = w[ii-1];
             wr = w[ii];
           }
#endif
							ncheck = 18;
							for(int ic=0;ic<ncheck/2;ic++)
								wcheck[ic] = wl +
								             double(ic+1)*(w[numP]-wl)/double(ncheck/2+1);
							for(int ic=0;ic<ncheck/2;ic++)
								wcheck[ic+ncheck/2] = wr -
								                      double(ic+1)*(wr-w[numP])/double(ncheck/2+1);

							adaptGP(dgp_flag,minRHS,maxRHS,deltaRHS,nOrtho,
							        sol, GP_solprev, GP_orth_solprev,
							        aa,bb,cc,
							        VhKV, VhMV, VhCV, VhK_arubber_lV, VhK_arubber_mV,
							        ncheck, wcheck, oldw, ((w1+w2)/2.0)/w[numP], atol);
							numP++;
							if (numP>=maxP) break;
						}
						sort(w.begin(),w.begin()+numP);
						if (numP>=maxP) break;
					}
					else break;
					if (verboseFlag) filePrint(stderr,"Rebuild time: %e\n",timex);
				}


				domain->isCoarseGridSolve = false;
				for(int i=0;i<=numS;i++) {
					double wc = w1 + (double(i))/double(numS)*(w2-w1);
					geoSource->setOmega(wc);
					double res;
					res = adaptGPSolRes(dgp_flag,nOrtho,sol,GP_solprev, GP_orth_solprev,
					                    aa,bb,cc,
					                    VhKV, VhMV, VhCV, VhK_arubber_lV, VhK_arubber_mV,
					                    wc, wc-oldw);
					domain->frequencies->push_front(wc);
					if(domain->solInfo().isAcousticHelm()) {
						SPropContainer& sProps = geoSource->getStructProps();
						for(int iProp=0;iProp<geoSource->getNumProps();iProp++) {
							if(sProps[iProp].kappaHelm!=0.0 || sProps[iProp].kappaHelmImag!=0.0) {
								complex<double> k1 = wc/sProps[iProp].soundSpeed;
								sProps[iProp].kappaHelm = real(k1);
								sProps[iProp].kappaHelmImag = imag(k1);
							}
						}
					}
					postProcessor->staticOutput(*sol, *rhs,  i==numS);
					domain->frequencies->pop_front();
				}
			} else if (domain->solInfo().getSweepParams()->isAdaptSweep &&
			           domain->solInfo().getSweepParams()->adaptSweep.deltaRHS<0)
			{
				adaptWindowSweep();
			} else if (domain->solInfo().getSweepParams()->isAdaptSweep &&
			           (domain->solInfo().loadcases.size()>1 ||
			            domain->numWaveDirections>1))
			{
				adaptSweep();
			} else {
// RTRT

				// some initialization
				probDesc->setIWaveDir(0);
				VecType **sol_prev = new VecType * [(nRHS+1)*padeN];
				for(int i = 0; i < (nRHS+1)*padeN; ++i)
					sol_prev[i] = new VecType(probDesc->solVecInfo());
				double *h = new double[padeN];

				//--- UH --- Data for Pade Lanczos
				double *PadeLanczos_VtKV = 0, *PadeLanczos_Vtb = 0;
				VecType **PadeLanczos_solprev = 0;
				if (domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::PadeLanczos) {
					//--- Allocate new set of pointers to have them without any gap
					PadeLanczos_solprev = new VecType * [nRHS*padeN];
					for (int ii = 0; ii < padeN; ++ii) {
						for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
							if (domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::PadeLanczos)
								PadeLanczos_solprev[iRHS + ii*nRHS] = sol_prev[iRHS + ii*(nRHS+1)];
							else
								PadeLanczos_solprev[iRHS + ii*nRHS] = sol_prev[iRHS+1 + ii*(nRHS+1)];
						}
					}
				} // if (domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos)
				//--- UH --- Data for Pade Lanczos
				VecType **GP_solprev = 0, **GP_orth_sol_prev = 0;

				if ( domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::GalProjection ||
				     domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::KrylovGalProjection ||
				     domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::WCAWEGalProjection) {
					GP_orth_sol_prev = new VecType * [nRHS*padeN];
					for(int i = 0; i < nRHS*padeN; ++i)
						GP_orth_sol_prev[i] = new VecType(probDesc->solVecInfo());
					//--- Allocate new set of pointers to have them without any gap
					GP_solprev = new VecType * [nRHS*padeN];
					for (int ii = 0; ii < padeN; ++ii)
						for (int iRHS = 0; iRHS < nRHS; ++iRHS)
							GP_solprev[iRHS + ii*nRHS] = sol_prev[iRHS+1 + ii*(nRHS+1)];
				}

				// loop over coarse grid
				bool first_time = true;
				bool savesol = false, savedsol = false;
				bool gpReorthoFlag = true;
				double w0;
				int count = 0;
				while(domain->coarse_frequencies->size() > 0) {
					int offset = count*(nRHS+1);
					sol_prev[offset]->zero();
					domain->isCoarseGridSolve = true;
					double wc = domain->coarse_frequencies->front();
					if(count == 0) w0 = wc;

					// if this isn't the first coarse solve then rebuild solver with new frequency
					if(!first_time || is>0) rebuildSolver(wc);
					else first_time = false;
					gpReorthoFlag = true;

					//----- UH ------
					if (domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::PadeLanczos) {
						filePrint(stderr,"\n ... Build Subspace M-orthonormal Basis      ... \n");
						PadeLanczos_BuildSubspace(nRHS, PadeLanczos_solprev + count * nRHS,
						                          PadeLanczos_solprev, count * nRHS);

					} else  if (domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::KrylovGalProjection) {
						krylovGalProjection(nRHS,count*nRHS,sol,
						                    GP_solprev, GP_orth_sol_prev,
						                    VhKV, VhMV, VhCV, wc, 0.0);
					} else {
						for(int iRHS = 0; iRHS < nRHS; iRHS++) {
							filePrint(stderr,"\n ... Solving RHS #%3d               ...\n",iRHS);
							if(iRHS > 0) {
								// PJSA: modified 5-29-04 now passing u and all derivatives to getFreqSweepRHS (needed for higher order sommerfeld)
								// note: sol_prev[0] = 0, sol_prev[1] = u, sol_prev[2] = u^(1), etc.
								*sol_prev[offset+iRHS] = *sol;
								probDesc->getFreqSweepRHS(rhs, sol_prev+offset, iRHS);
							}
							else probDesc->getRHS(*rhs);
							sol->zero();
							allOps->sysSolver->solve(*rhs, *sol);
//         forceContinuity(*sol);

/*
         filePrint(stderr,"\n ... Dumping RHS #%3d               ...\n",iRHS);
         *rhs = *sol;
          scaleDisp(*rhs);
         postProcessor->staticOutput(*rhs, *rhs,false); 
*/
							rhs->zero();
						} // for (int iRHS = 0; iRHS < nRHS; iRHS++)
						*sol_prev[offset+nRHS] = *sol;
					}
					//----- UH ------
					// PJSA 9-22-06: fix for coupled multi point pade
					if(domain->solInfo().isCoupled)
						if (domain->solInfo().getSweepParams()->freqSweepMethod != SweepParams::KrylovGalProjection && domain->solInfo().getSweepParams()->freqSweepMethod != SweepParams::WCAWEGalProjection)
							for(int i=1; i<(nRHS+1); ++i) scaleDisp(*sol_prev[offset+i]);
					bool printTimers = ((domain->coarse_frequencies->size()+domain->frequencies->size()) > 1 || is < domain->set_of_frequencies.size()-1  ) ? false : true;

					domain->setSavedFreq(domain->coarse_frequencies->front());
					//----- UH ------
					if (domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::PadeLanczos) {
						//--- We destroy the arrays PadeLanczos_VtKV and PadeLanczos_VtB when V is incomplete.
						//--- This is not optimal as we recompute several times some coefficients.
						//--- Outputting should be done only when the basis V is complete.
						if(PadeLanczos_VtKV) { delete[] PadeLanczos_VtKV; PadeLanczos_VtKV = 0; }
						if(PadeLanczos_Vtb) { delete[] PadeLanczos_Vtb; PadeLanczos_Vtb = 0; }

						//--- sol_prev[offset] stores the normalized solution vector !
						//--- So we need to recompute the solution.
						//--- Note that wc is referred only on the first call to this routine.
						//--- It is used to get K (as K is actually storing K - wc*wc*M)
						PadeLanczos_Evaluate((count+1)*nRHS, PadeLanczos_solprev,
						                     PadeLanczos_VtKV, PadeLanczos_Vtb, wc, sol);
						if (savesol)  {
							savedsol = true;
							*savedSol = *sol;
							*savedRhs = *rhs;
						}
						else {
							savesol = true;
							postProcessor->staticOutput(*sol, *rhs, printTimers);
						}
					} else if (domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::KrylovGalProjection || domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::WCAWEGalProjection) {
						if (savesol)  {
							savedsol = true;
							*savedSol = *sol;
							*savedRhs = *rhs;
						}
						else {
							savesol = true;
							postProcessor->staticOutput(*sol, *rhs, printTimers);
						}
					} else {
						if (savesol)  {
							savedsol = true;
							*savedSol = *sol_prev[offset+1];
							*savedRhs = *rhs;
						}
						else {
							savesol = true;
							postProcessor->staticOutput(*sol_prev[offset+1], *rhs, printTimers);
						}
					}
					//----- UH ------

					if(domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::Pade ||
					   domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::Fourier) {
						double hi = wc-w0;
						h[count] = hi;
					}
					count++;

					if(!(domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::PadeLanczos &&
					     padeN > 1 && count == padeN && domain->coarse_frequencies->size() > 1))
						domain->coarse_frequencies->pop_front();

					if(count == padeN) {  // this means there are enough coarse solves to solve for pade polynomial coefs & reconstruct
						filePrint(stderr, " --------------------------------------\n");
						count = 0;
						double max_freq_before_rebuild = wc;
						int ncoarse = domain->coarse_frequencies->size();
						if((ncoarse > 0) && (padeN==1))
							max_freq_before_rebuild = wc + (domain->coarse_frequencies->front() - wc)/2.0;

						// now loop over all the freqencies & reconstruct solutions for each
						// starting with 2nd (1st has already been solved & output)
						StaticTimers *times = probDesc->getStaticTimers();
						double xtime = 0.0;
						xtime -= getTime();

						while((domain->frequencies->size() > 0) && ((domain->frequencies->front() < max_freq_before_rebuild)
						                                            || (domain->coarse_frequencies->size() == 0))) {

							startTimerMemory(times->timeFreqSweep, times->memoryFreqSweep);
							domain->isCoarseGridSolve = false;
							double w = domain->frequencies->front();
							double deltaw;
							if (padeN>1) deltaw = w - w0;
							else deltaw = w - wc;
							filePrint(stderr, " ... Reconstructing solution for f = %f ...\n", w/(2.0*PI));
							//--------
							switch(domain->solInfo().getSweepParams()->freqSweepMethod) {
								case SweepParams::Taylor :
								{
									double _IsConverged = 1.0e-12;
									*sol = *sol_prev[1];
									double prevSolSqNorm, solSqNorm = sol->sqNorm();
									double relError = 0.0;
									//std::cerr << "  i = 0, sol->sqNorm() = " << sol->sqNorm() << std::endl;
									for(int i=1; i<nRHS; ++i) {
										prevSolSqNorm = solSqNorm;
										sol->linAdd(1.0/DFactorial(i)*pow(deltaw,i), *sol_prev[i+1]);
										solSqNorm = sol->sqNorm();
										relError = fabs((solSqNorm-prevSolSqNorm)/prevSolSqNorm);
										//std::cerr << "  i = " << i << ", solSqNorm = " << solSqNorm
										//     << ", factor = " << 1.0/double(DFactorial(i))*pow(deltaw,i)
										//     << ", sol_prev[i+1]->sqNorm() = " << sol_prev[i+1]->sqNorm()
										//     << ", relError = " << relError << std::endl;
										if(relError < _IsConverged) break;
									}
								} break;
								case SweepParams::Pade :
									pade(sol, sol_prev, h, deltaw);
									break;
								case SweepParams::Pade1 :
									pade1(sol, sol_prev, deltaw);
									break;
									//--- UH --- Evaluate the Pade approximation
								case SweepParams::PadeLanczos :
									PadeLanczos_Evaluate(padeN * nRHS, PadeLanczos_solprev,
									                     PadeLanczos_VtKV, PadeLanczos_Vtb, w, sol);
									break;
									//--- UH ---
								case SweepParams::GalProjection:
									galProjection(gpReorthoFlag,nRHS*padeN,sol,GP_orth_sol_prev,
									              GP_solprev, VhKV, VhMV, VhCV,
									              VhK_arubber_lV, VhK_arubber_mV,
									              w, w-wc);
									gpReorthoFlag = false;
									break;
								case SweepParams::KrylovGalProjection:
									krylovGalProjection(0,nRHS*padeN,sol,GP_solprev, GP_orth_sol_prev,
									                    VhKV, VhMV, VhCV, w, w-wc);
									break;
								case SweepParams::Fourier :
									fourier(sol, sol_prev, h, deltaw);
									break;
							} // switch (domain->solInfo().freqSweepMethod)
							//--------
							geoSource->setImpe(w); // for output file and/or restart
							bool printTimers = ((ncoarse+domain->frequencies->size()) > 1) ? false : true;
							stopTimerMemory(times->timeFreqSweep, times->memoryFreqSweep);

							postProcessor->staticOutput(*sol, *rhs, printTimers);
							domain->frequencies->pop_front();
							if(padeN == 1 && w < w0 && !domain->frequencies->empty() && domain->frequencies->front() > w0) {
								domain->isCoarseGridSolve = true;
								postProcessor->staticOutput(*savedSol, *savedRhs, printTimers);
								savedsol = false;
							}
						} // while((domain->frequencies->size() > 0) ... )
						xtime += getTime();
						if (verboseFlag) filePrint(stderr,"Projection  time: %e\n",xtime);

						// print out saved solution
						if (savedsol)  {
							domain->isCoarseGridSolve = true;
							postProcessor->staticOutput(*savedSol, *savedRhs, printTimers);
						}
						else {
							savesol = true;
						}

						if(domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::PadeLanczos) {
						}
						else {
							if((padeN > 1) && (ncoarse > 0)) {
								count = (padeN-ncoarse > 1) ? padeN-ncoarse : 1;
								int coffset = (padeN-count)*(nRHS+1);
								for(int i = 0; i < (nRHS+1)*count; ++i) *sol_prev[i] = *sol_prev[coffset+i];
								if (domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::KrylovGalProjection || domain->solInfo().getSweepParams()->freqSweepMethod == SweepParams::WCAWEGalProjection)
									for(int i = 0; i < nRHS*count; ++i)
										*GP_orth_sol_prev[i] = *GP_orth_sol_prev[(padeN-count)*nRHS+i];
								double deltah = h[padeN-count];
								w0 += deltah;
								for(int i=0; i<count; ++i) h[i] = h[padeN-count+i] - deltah;
							} // if((padeN > 1) && (ncoarse > 0))
						}

					} // if (count == padeN)

				} // while(domain->coarse_frequencies->size() > 0)
				for(int i = 0; i < (nRHS+1)*padeN; ++i)
					delete sol_prev[i];
				delete [] sol_prev;
				if (domain->solInfo().getSweepParams()->freqSweepMethod  == SweepParams::GalProjection) {
					for(int i = 0; i < nRHS*padeN; ++i)
						delete GP_orth_sol_prev[i];
					delete [] GP_orth_sol_prev;
				}
				delete [] h;
				//--- UH --- Data for Pade Lanczos
				if (PadeLanczos_VtKV)
					delete[] PadeLanczos_VtKV;
				if (PadeLanczos_Vtb)
					delete[] PadeLanczos_Vtb;
				if (PadeLanczos_solprev)
					delete[] PadeLanczos_solprev;
				//--- UH --- Data for Pade Lanczos
				if (GP_solprev) delete[] GP_solprev;
				if (VhKV) delete[] VhKV;
				if (VhMV) delete[] VhMV;
				if (VhCV) delete[] VhCV;
			}
		}
	} // if(domain->solInfo().doFreqSweep)
	else if(domain->solInfo().probType == SolverInfo::HelmholtzDirSweep) { // PJSA: never tested this in FEM/DPH
		int iRHS;
		int nDir = domain->getNumWaveDirections();
		if (nDir==0) nDir =1;
		for (iRHS = 0; iRHS < nDir; iRHS++) {
			probDesc->setIWaveDir(iRHS);
			if(nDir > 1) filePrint(stderr,"\n Running RHS #%d\n",iRHS);
			probDesc->getRHS(*rhs);
			allOps->sysSolver->solve(*rhs, *sol);
			postProcessor->staticOutput(*sol, *rhs);
			rhs->zero();
			sol->zero();
		}
	}
	else if(domain->solInfo().noninpc) {
		// initialization
		probDesc->getRHS(*rhs);
		VecType *psi_u = new VecType(probDesc->solVecInfo(sfem->getP()));
		if(!isFeti(domain->solInfo().solvercntl->type))
			psi_u->setn(domain->numUncon());

		psi_u->zero();
		for(int i=0; i<domain->solInfo().nsample; ++i) {
			std::cout << std::endl << "Realization number  " << i+1 << std::endl;
			if(i>0) {
				sfem->genXiPsi(i); // seed=i
				probDesc->assignRandMat();
				probDesc->rebuildSolver();
				postProcessor->setSolver(allOps->sysSolver);
			}
			allOps->sysSolver->solve(*rhs,*sol);
			sfem_noninpc->update_Psi_u(sol,psi_u);
			sol->zero();
		}
		delete sol;
		sol = 0;
		sfem_noninpc->compute_Coefs(psi_u);
		VecType *sol1 = new VecType(probDesc->solVecInfo(1));
		if(!isFeti(domain->solInfo().solvercntl->type))
			sol1->setn(domain->numUncon());

		sfem_noninpc->computeMean(psi_u,sol1);
		postProcessor->staticOutput(*sol1, *rhs, false, 1);
		delete sol1;
		VecType *sol2 = new VecType(probDesc->solVecInfo(1));
		if(!isFeti(domain->solInfo().solvercntl->type))
			sol2->setn(domain->numUncon());

		sfem_noninpc->computeStdDev(psi_u,sol2);
		postProcessor->staticOutput(*sol2, *rhs, false, 2);
		delete sol2;

		VecType *sol3 = new VecType(probDesc->solVecInfo(1));
		if(!isFeti(domain->solInfo().solvercntl->type))
			sol3->setn(domain->numUncon());
		int nsample_output = sfem->getnsamp_out();
		for(int i=0; i< nsample_output; ++i) {
			sfem_noninpc->computePdf(i,psi_u,sol3);
			postProcessor->staticOutput(*sol3, *rhs, true, 3);
		}
		delete sol3;
		stochStress(psi_u);
	}
 else if(domain->solInfo().inpc) {
   probDesc->getRHSinpc(*rhs);
   allOps->sysSolver->solve(*rhs,*sol);
//   postProcessor->staticOutput(*sol, *rhs, true, 0); /

   VecType *sol1 = new VecType(probDesc->solVecInfo(1));
   if(!isFeti(domain->solInfo().solvercntl->type))
   	sol1->setn(domain->numUncon());
   sfem_inpc->computeMean(sol,sol1);
   postProcessor->staticOutput(*sol1, *rhs, true, 1);
   delete sol1;
 
   VecType *sol2 = new VecType(probDesc->solVecInfo(1)); 
   if(!isFeti(domain->solInfo().solvercntl->type))
   	sol2->setn(domain->numUncon());
   sfem_inpc->computeStdDev(sol,sol2);
   postProcessor->staticOutput(*sol2, *rhs, false, 2);
   delete sol2;

   VecType *sol3 = new VecType(probDesc->solVecInfo(1));
   if(!isFeti(domain->solInfo().solvercntl->type))
   	sol3->setn(domain->numUncon());
   int nsample_output = sfem->getnsamp_out();
   for(int i=0; i<nsample_output; ++i) {
     sfem_inpc->computePdf(i,sol,sol3); 
     postProcessor->staticOutput(*sol3, *rhs, false, 3);
   }
   delete sol3;

   probDesc->retrieveElemset(); 
   delete allOps->sysSolver; 
   allOps->sysSolver = 0;
   stochStress(sol); // YYY Commented out temporarily
 }
 else if(domain->solInfo().loadcases.size() > 0) {
   while(domain->solInfo().loadcases.size() > 0) {
     probDesc->getRHS(*rhs);
     allOps->sysSolver->solve(*rhs,*sol);
     if(domain->solInfo().isCoupled) scaleDisp(*sol); // PJSA 9-22-06
     bool printTimers = (domain->solInfo().loadcases.size() > 1) ? false : true;
     postProcessor->staticOutput(*sol, *rhs, printTimers);
     domain->solInfo().loadcases.pop_front();
   }
 }
 else if (domain->numWaveDirections>0) {
     for(int iDir=0;iDir<domain->numWaveDirections;iDir++) {
     probDesc->setIWaveDir(iDir);
     probDesc->getRHS(*rhs);
     allOps->sysSolver->solve(*rhs,*sol);
     if(domain->solInfo().isCoupled) scaleDisp(*sol); // PJSA 9-22-06
     bool printTimers = (iDir  <domain->numWaveDirections-1) ? false : true;
     postProcessor->staticOutput(*sol, *rhs,printTimers);
     }
 }
 else {
   probDesc->getRHS(*rhs);
   try {
     allOps->sysSolver->solve(*rhs,*sol);
   }
   catch(std::runtime_error& e) {
     std::cerr << "exception: " << e.what() << std::endl;
   }
   if(domain->solInfo().isCoupled) scaleDisp(*sol); // PJSA 9-22-06
   postProcessor->staticOutput(*sol, *rhs);
 }

#ifdef USE_EIGEN3
 if(domain->solInfo().sensitivity) { 
   probDesc->postProcessSA(*sol);
 }
#endif

 geoSource->closeOutputFiles(); 
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::stochStress(VecType *sol)
{
 
//  int n = domain->numUncon();
//  int P = sfem->getP();
//  int n_node = geoSource->numNode();
//  int npsize=P*n_node;

  MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>* mltintg = new MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>();


//  int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes
  Scalar *globVal = 0;
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int qmd=4;
  int dof;
//  int iNode;
  for(int i = 0; i < numOutInfo; ++i)  { 
    if(oinfo[i].ndtype == 0) continue; // YYY fix ndflag issue
    if(oinfo[i].interval == 1) {
      dof = -1;

/*      if(domain->solInfo().noninpc)
        mltintg->computeStressStat(qmd,sol,i,oinfo[i].type,postProcessor,probDesc,oinfo[i].ndtype);
      else
        mltintg->computeStressStat(qmd,sol,i,oinfo[i].type,postProcessor,probDesc,oinfo[i].ndtype);
*/
      switch(oinfo[i].type)  {
     
        default:
        //  fprintf(stderr, " *** WARNING: output %d is not supported \n", i);
          break;

        case OutputInfo::StressXX:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SXX,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SXX,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressYY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SYY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SYY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressZZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SZZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SZZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressXY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SXY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SXY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressYZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SYZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SYZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressXZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SXZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SXZ,postProcessor,probDesc,oinfo[i].ndtype);
          break; 
        case OutputInfo::StressVM:
            if(domain->solInfo().noninpc)
             mltintg->computeStressStat(qmd,sol,i,VON,postProcessor,probDesc,oinfo[i].ndtype);
            else
             mltintg->computeStressStat(qmd,sol,i,VON,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainXX:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EXX,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EXX,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainYY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EYY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EYY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainZZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EZZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EZZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainXY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EXY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EXY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainYZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EYZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EYZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainXZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EXZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EXZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainVM:
            if(domain->solInfo().noninpc)
             mltintg->computeStressStat(qmd,sol,i,STRAINVON,postProcessor,probDesc,oinfo[i].ndtype);
            else
             mltintg->computeStressStat(qmd,sol,i,STRAINVON,postProcessor,probDesc,oinfo[i].ndtype);
          break;

      }
    } 
    if(globVal)
      { delete [] globVal; globVal = 0; }
  } 
  delete mltintg;
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::rebuildSolver(double w1)
{
  dbg_alloca(0);
  geoSource->setOmega(w1);
  filePrint(stderr,"\n ... Rebuilding LHS with frequency = %f ... \n", w1/(2.0*PI));
  if(domain->solInfo().isAcousticHelm()) {  // set new wave number for acoustic helmholtz
    SPropContainer& sProps = geoSource->getStructProps();
    for(int iProp=0;iProp<geoSource->getNumProps();iProp++) {
       if(sProps[iProp].kappaHelm!=0.0 || sProps[iProp].kappaHelmImag!=0.0) {
         complex<double> k1 = w1/sProps[iProp].soundSpeed;
         sProps[iProp].kappaHelm = real(k1);
         sProps[iProp].kappaHelmImag = imag(k1);
       } 
    }
  }
  probDesc->rebuildSolver();
  postProcessor->setSolver(allOps->sysSolver);
  //solver = probDesc->getSolver();

  if(a) {
    for(int i = 0; i < ia; ++i) delete a[i];
    delete [] a;
    a = 0;
  }
  if(b) {
    for(int i = 0; i < ib; ++i) delete b[i];
    delete [] b;
    b = 0;
  }
  if(P) { delete P; P = 0; }
  if(Q) { delete Q; Q = 0; }

  if(ca) {
    for(int i = 0; i < ia; ++i) delete ca[i];
    delete [] ca;
    ca = 0;
  }
  if(cb) {
    for(int i = 0; i < ib; ++i) delete cb[i];
    delete [] cb;
    cb = 0;
  }
  if(cP) { delete cP; cP = 0; }
  if(cQ) { delete cQ; cQ = 0; }

}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::scaleDisp(VecType &u)
{
  probDesc->scaleDisp(u);
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::scaleInvDisp(VecType &u)
{
  probDesc->scaleInvDisp(u);
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::scaleDisp(VecType &u, double alpha)
{
  probDesc->scaleDisp(u, alpha);
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::forceContinuity(VecType &u)
{
  probDesc->forceContinuity(u);
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::forceAssemble(VecType &u)
{
  probDesc->forceAssemble(u);
}


//------------------------------------------------------------------------------

