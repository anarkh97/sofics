#include <Timers.d/GetTime.h>
#include <unistd.h>
#include <Paral.d/SubDOp.h>
#include <Driver.d/SysState.h>
#include <Problems.d/DynamProbTraits.h>

#include <limits>

extern int contactPrintFlag;

//-------------------------------------------------------------------------------------------
template<class VecType>
SysState<VecType> & SysState<VecType>::operator=(const SysState<VecType> &v2)
{
   d_n = v2.getDispConst();
   v_n = v2.getVelocConst();
   a_n = v2.getAccelConst();
   v_n_p = v2.getPrevVelocConst();

   return *this;

}

//-------------------------------------------------------------------------------------------

template <class VecType,
          class ProblemDescriptor> 
NewmarkWorkVec<VecType,ProblemDescriptor>::NewmarkWorkVec(int _typ, ProblemDescriptor *probDesc)
{
   typ=_typ;

   switch (typ) {

     case -1:
       d_n_p  = new VecType( probDesc->solVecInfo() );
       o_n_p  = new VecType( probDesc->masterSolVecInfo() );
       ext_f  = new VecType( probDesc->solVecInfo() );
       rhs    = new VecType( probDesc->solVecInfo() );
       break;
     case 0:
       tmp1   = new VecType( probDesc->solVecInfo() );
       tmp2   = new VecType( probDesc->solVecInfo() );
       fint   = new VecType( probDesc->solVecInfo() );
       ext_f  = new VecType( probDesc->solVecInfo() );
       v_n_h  = new VecType( probDesc->solVecInfo() );
       break;
     case 1:
       d_n_p  = new VecType( probDesc->solVecInfo() );
       v_n_p  = new VecType( probDesc->solVecInfo() );
       a_n_p  = new VecType( probDesc->solVecInfo() );
       rhs    = new VecType( probDesc->solVecInfo() );
       ext_f  = new VecType( probDesc->solVecInfo() );
       d_n_h  = new VecType( probDesc->solVecInfo() );
       v_n_h  = new VecType( probDesc->solVecInfo() );
       Md_n_h = new VecType( probDesc->solVecInfo() );
       Cd_n_h = new VecType( probDesc->solVecInfo() );
       tmp1   = new VecType( probDesc->solVecInfo() );
       break;
   }
}

//-----------------------------------------------------------------------------
template <class VecType,
          class ProblemDescriptor>
NewmarkWorkVec<VecType,ProblemDescriptor> & NewmarkWorkVec<VecType,ProblemDescriptor>:: operator=(const NewmarkWorkVec<VecType,ProblemDescriptor> &v)
{

     if (typ!=v.typ) {
           
         switch (typ) {
                                                                                                  
	     case -1:
       		delete  d_n_p;
       		delete  o_n_p;
       		delete  ext_f;
       		delete  rhs;
       		break;
     	     case 0:
       		delete  tmp1;
       		delete  tmp2;
       		delete  fint;
                delete  ext_f;
                delete  v_n_h;
       		break;
     	     case 1:
       		delete  d_n_p;
       		delete  v_n_p;
       		delete  a_n_p;
       		delete  rhs;
       		delete  ext_f;
       		delete  d_n_h;
       		delete  v_n_h;
       		delete  Md_n_h;
       		delete  Cd_n_h;
       		delete  tmp1;
       		break;
      	 }
         typ=v.typ;
     }

     switch (typ) {
                                                                                                  
         case -1:
            d_n_p  = new VecType();
            *d_n_p = v.get_d_n_pConst();
            o_n_p  = new VecType();
            *o_n_p = v.get_o_n_pConst();
            ext_f  = new VecType();
            *ext_f = v.get_ext_fConst();
            rhs    = new VecType();
            *rhs   = v.get_rhsConst();
            break;
         case 0:
            tmp1  = new VecType();
            *tmp1 = v.get_tmp1Const();
            tmp2  = new VecType();
            *tmp2 = v.get_tmp2Const();
            fint  = new VecType();
            *fint = v.get_fintConst();
            ext_f  = new VecType();
            *ext_f = v.get_ext_fConst();
            v_n_h  = new VecType();
            *v_n_h = v.get_v_n_hConst();
            break;
         case 1:
            d_n_p  = new VecType();
            *d_n_p = v.get_d_n_pConst();
            v_n_p  = new VecType();
            *v_n_p = v.get_v_n_pConst();
            a_n_p  = new VecType();
            *a_n_p = v.get_a_n_pConst();
            rhs    = new VecType();
            *rhs   = v.get_rhsConst();
            ext_f  = new VecType();
            *ext_f = v.get_ext_fConst();
            d_n_h  = new VecType();
            *d_n_h = v.get_d_n_hConst();
            v_n_h  = new VecType();
            *v_n_h = v.get_v_n_hConst();
            Md_n_h = new VecType();
            *Md_n_h = v.get_Md_n_hConst();
            Cd_n_h = new VecType();
            *Cd_n_h = v.get_Cd_n_hConst();
            tmp1   = new VecType();
            *tmp1  = v.get_tmp1Const();
            break;
   }

}
//------------------------------------------------------------------------------

template <class VecType, 
          class ProblemDescriptor> 
NewmarkWorkVec<VecType,ProblemDescriptor>::~NewmarkWorkVec()
{
   switch (typ) {
   
     case -1:
       delete  d_n_p;
       delete  o_n_p;
       delete  ext_f;
       delete  rhs;
       break;
     case 0:
       delete  tmp1;
       delete  tmp2;
       delete  fint;
       delete  ext_f;
       delete  v_n_h;
       break;
     case 1:
       delete  d_n_p;  
       delete  v_n_p;  
       delete  a_n_p;  
       delete  rhs;    
       delete  ext_f; 
       delete  d_n_h;  
       delete  v_n_h;
       delete  Md_n_h;
       delete  Cd_n_h;
       delete  tmp1;
       break;
   }
}

//------------------------------------------------------------------------------

template < 
     class DynOps, 
     class VecType, 
     class PostProcessor, 
     class ProblemDescriptor,
     class Scalar>
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>::DynamicSolver(ProblemDescriptor *PrbD)
 : lambda_nSen(0), d_nSen(0), v_nSen(0), a_nSen(0), v_pSen(0), rhsSen(0), aeroForceSen(0), curSenState(0)
{
    probDesc = PrbD; 
    aeroAlg = -1;
    dynOps = 0;
    postProcessor = 0;
    aeroForce = 0;
}

//------------------------------------------------------------------------------

template < 
     class DynOps, 
     class VecType, 
     class PostProcessor, 
     class ProblemDescriptor,
     class Scalar>
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>::~DynamicSolver() {
  if(postProcessor) delete postProcessor;
  if(dynOps) delete dynOps;
  if(aeroForce) delete aeroForce;
}

//------------------------------------------------------------------------------

template <
     class DynOps,             // Data Structure for K, C, M and dynMat
     class VecType,            // Vector type used
     class PostProcessor,      
     class ProblemDescriptor,
     class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::solve()
{
   // Construct renumbering, connectivities, dofsets
   // and boundary conditions. Suppose to be common for
   // all the different dynamic schemes.

   probDesc->preProcess();

#ifdef USE_EIGEN3
   if(domain->solInfo().sensitivity) {
     probDesc->preProcessSA();
     if(!domain->runSAwAnalysis) {
       AllSensitivities<double> *allSens = probDesc->getAllSensitivities();
       domain->sensitivityPostProcessing(*allSens);
       return;
     }
   }
#endif

   postProcessor = probDesc->getPostProcessor();

   // Allocate vectors for displacement, velocity, 
   //                      acceleration and last velocity
   d_n = new VecType( probDesc->solVecInfo() );
   v_n = new VecType( probDesc->solVecInfo() );
   a_n = new VecType( probDesc->solVecInfo() );
   v_p = new VecType( probDesc->solVecInfo() );

   // Get time loop information 
   aeroAlg = probDesc->getAeroAlg();
   if(aeroAlg >= 0 || probDesc->getThermoeFlag() >= 0) {
     probDesc->computeTimeInfo(); // computes a new tmax if necessary
   }
   probDesc->getTimes( dt, tmax );

   // Set up initial conditions
   curState = new SysState<VecType>( *d_n, *v_n, *a_n, *v_p);
   probDesc->getInitState( *curState );
#ifdef USE_EIGEN3
   if(domain->solInfo().sensitivity) {
     if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
       d_nSen = new VecType( probDesc->solVecInfo() );
       v_nSen = new VecType( probDesc->solVecInfo() );
       a_nSen = new VecType( probDesc->solVecInfo() );
       v_pSen = new VecType( probDesc->solVecInfo() );
       *d_nSen = *v_nSen = *a_nSen = *v_pSen = 0.0;
       curSenState = new SysState<VecType>( *d_nSen, *v_nSen, *a_nSen, *v_pSen);
     }
     else if(domain->solInfo().sensitivityMethod == SolverInfo::Adjoint) {
       lambda_nSen = new VecType( probDesc->solVecInfo() );
       v_nSen = new VecType( probDesc->solVecInfo() );
       a_nSen = new VecType( probDesc->solVecInfo() );
       v_pSen = new VecType( probDesc->solVecInfo() );
       *lambda_nSen = *v_nSen = *a_nSen = *v_pSen = 0.0;
       curSenState = new SysState<VecType>( *lambda_nSen, *v_nSen, *a_nSen, *v_pSen);
     }
     probDesc->getSensitivityStateParam(sensitivityTol,ratioSensitivityTol);
   }
#endif

   // The aeroPreProcess is done later now, so that the correct
   // time-step is passed to the fluid in the explicit case
   // However, if we are doing a ping-pong, first we run aeroPreProcess and return 
   if(aeroAlg == 1 || aeroAlg == 8) {
     probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p ); // [S] sent initial displacements ...
     return;
   }

   aeroForce = (aeroAlg >= 0 || probDesc->getThermoeFlag() >= 0) ? new VecType(probDesc->solVecInfo()) : 0;
   if(aeroForce) aeroForce->zero();
   aeroForceSen = (domain->solInfo().sensitivity) ? new VecType(probDesc->solVecInfo()) : 0;
   if(aeroForceSen) aeroForceSen->zero();

   // Build time independent forces i.e. gravity force, pressure force
   constForce = new VecType( probDesc->solVecInfo() );

   // Get SteadyState Flag and Parameters
   probDesc->getSteadyStateParam(steadyFlag, steadyMin, steadyMax, steadyTol);

   // Get Time Integration Scheme
   algType = probDesc->getTimeIntegration();
   if(aeroAlg == 10 && algType != 1) {
      fprintf(stderr, "WARNING: B0 AERO type only valid for quasi-static.  Running QUASISTATICS for 1 iteration\n");
      domain->solInfo().timeIntegration = algType = 1;
   }
   
   double timeLoop = -getTime();
   switch ( algType )
   {
     // Newmark
     case 0:

       // Get Newmark Parameter  
       probDesc->getNewMarkParameters( beta, gamma, alphaf, alpham );
       
       // ... Newmark Beta == 0 -> CENTRAL DIFFERENCES ALGORITHM
       // ... Newmark Beta != 0 -> NEWMARK ALGORITHM
       //     This is now the Generalized Alpha Method, of which
       //     the NEWMARK algorithm is a subset.

       if(beta == 0.0) {

         // Defining Working Arrays   
         workVec = new NewmarkWorkVec<VecType,ProblemDescriptor>(0,probDesc);
 
         if(domain->solInfo().order != 1 && gamma != 0.5) {
           filePrint(stderr," ### WATCH gamma for explicit Newmark is set to 0.5 ###\n");
           gamma=0.5;
         }

         // Build Necessary Operators (K, M)
         bool fourthOrder = probDesc->getDomain()->solInfo().modifiedWaveEquation;
         if (fourthOrder) 
           dynOps = probDesc->buildOps(1.0, dt*gamma, dt*dt*beta);
         else
           dynOps = probDesc->buildOps(1.0, 0.0, 0.0);

         probDesc->getConstForce( *constForce );
 
         // Check stability time step
         if(domain->solInfo().stable && (aeroAlg < 0 || domain->solInfo().dyna3d_compat)) probDesc->computeStabilityTimeStep(dt, *dynOps);

         if(aeroAlg == 20) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_n ); // See Eq. 51 of C.Farhat et al. IJNME(2010) Robust and provably ... 
         else
           if(aeroAlg >= 0) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p );
         if(probDesc->getThermoeFlag() >= 0) probDesc->thermoePreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getAeroheatFlag() >= 0) probDesc->aeroHeatPreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getThermohFlag() >= 0) probDesc->thermohPreProcess(*d_n, *v_n, *v_p);
         
         // Time Integration Loop 
         explicitNewmarkLoop( *curState, *constForce, *dynOps, *workVec, dt, tmax);
       } 
       else {

         if(aeroAlg >= 0) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p );
         if(probDesc->getThermoeFlag() >= 0) probDesc->thermoePreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getAeroheatFlag() >= 0) probDesc->aeroHeatPreProcess(*d_n, *v_n, *v_p);
         if(probDesc->getThermohFlag() >= 0) probDesc->thermohPreProcess(*d_n, *v_n, *v_p);

         // Defining Working Arrays   
         workVec = new NewmarkWorkVec<VecType,ProblemDescriptor>(1,probDesc);

         if(domain->solInfo().order == 1) { // heat
           // build K, M and dynK = (M + (dt*gamma)*K
           dynOps = probDesc->buildOps(1, 0, dt*gamma);
         }
         else { // mech and acou
           // Build K, M, C and dynK = ((1-alpham)/(1-alphaf))*M + (dt*gamma)*C + (dt*dt*beta)*K
           dynOps = probDesc->buildOps(((1.0-alpham)/(1.0-alphaf)), dt*gamma, dt*dt*beta);
         }

         probDesc->getConstForce( *constForce );

         // Time Integration Loop 
         implicitNewmarkLoop( *curState, *constForce, *dynOps, *workVec, dt, tmax);
       }
       break;
       
     // Quasi-Static
     case 1:
       if(aeroAlg >= 0) probDesc->aeroPreProcess( *d_n, *v_n, *a_n, *v_p ); // [S] sent initial displacements ...
       if(probDesc->getThermoeFlag() >= 0) probDesc->thermoePreProcess(*d_n, *v_n, *v_p);
       if(probDesc->getAeroheatFlag() >= 0) probDesc->aeroHeatPreProcess(*d_n, *v_n, *v_p);
       if(probDesc->getThermohFlag() >= 0) probDesc->thermohPreProcess(*d_n, *v_n, *v_p);

       probDesc->getQuasiStaticParameters(maxVel, delta);
       
       // Defining Working Arrays   
       workVec = new NewmarkWorkVec<VecType,ProblemDescriptor>(-1,probDesc);
       if(domain->solInfo().sensitivity) { workSenVec = new NewmarkWorkVec<VecType,ProblemDescriptor>(-1,probDesc); }

       // Build Necessary Operators (only K!)
       dynOps = probDesc->buildOps(0,0,1); 
 
       probDesc->getConstForce( *constForce );
       
       // Quasi-Static Loop
       quasistaticLoop( *curState, *constForce, *dynOps, *workVec, dt, tmax, aeroAlg);   

       // Aeroelastic Sensitivity Quasi-Static 
#ifdef USE_EIGEN3
       SensitivityInfo *senInfo = probDesc->getSensitivityInfo();
       int numThicknessGroups = domain->getNumThicknessGroups();
       int numShapeVars = domain->getNumShapeVars();
       int numStructParamTypes = int(bool(numThicknessGroups))+int(bool(numShapeVars));
       int numStructQuantTypes = domain->getNumSensitivityQuantityTypes();
       int numFluidQuantTypes(0);
       if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
         probDesc->sendNumParam(numStructParamTypes, 0, ratioSensitivityTol*sensitivityTol);
       } else { // Adjoint method
         probDesc->sendNumParam(numStructQuantTypes, 0, ratioSensitivityTol*sensitivityTol);
         probDesc->getNumParam(numFluidQuantTypes);
         filePrint(stderr, " ... In structure,  numFluidQuantTypes is %d\n", numFluidQuantTypes);
       }
       
       if(domain->solInfo().sensitivity) { 
         double wall0 = -getTime(); 
         probDesc->postProcessSA(dynOps,*d_n);
         wall0 += getTime();
         std::cerr << " ****** wall clock time for sensitivity pre-computation is " << wall0/1000.0 << std::endl;
         AllSensitivities<double> *allSens = probDesc->getAllSensitivities();
         double walls = -getTime();
         for(int isen = 0; isen < probDesc->getNumSensitivities(); ++isen) {  // structure sensitivities
           switch (senInfo[isen].type) {

             case SensitivityInfo::DisplacementWRTshape:

               if( numShapeVars > 0 ) { 
                 if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
                   probDesc->sendNumParam(numShapeVars, 1, ratioSensitivityTol*sensitivityTol);
                   if(!allSens->dispWRTshape) { 
                     allSens->dispWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
                     for(int ishap=0; ishap< numShapeVars; ++ishap) {
                       allSens->dispWRTshape[ishap] = new Eigen::Matrix<double, Eigen::Dynamic,
                                                                                Eigen::Dynamic> (domain->numUncon(),1);
                       rhsSen = new VecType( probDesc->solVecInfo() );
                       rhsSen->copy(allSens->linearstaticWRTshape[ishap]->data());
                       (*rhsSen) *= -1; 
                       aeroForceSen->zero();
                       *d_nSen = 0.0;
                       aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
                       *allSens->dispWRTshape[ishap] = Eigen::Map<Eigen::Matrix<double,
                                                                  Eigen::Dynamic,Eigen::Dynamic>>(d_nSen->data(),domain->numUncon(),1);
                       delete rhsSen;
                     }
                   }
                 }
                 else if(domain->solInfo().sensitivityMethod == SolverInfo::Adjoint) {
                   
                   int numDispNodes = domain->getNumDispNodes();
                   int numTotalDispDofs = domain->getTotalNumDispDofs();
                   if(!allSens->lambdaDisp) {
                     probDesc->sendNumParam(numTotalDispDofs, 8,
                                            ratioSensitivityTol*sensitivityTol); // actvar = 8, structure displacement
                     computeLambdaDisp(numTotalDispDofs, numDispNodes);
                   }
                   allSens->dispWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
                   std::vector<DispNode> *dispNodes = domain->getDispNodes();
                   for(int ishap=0; ishap < numShapeVars; ++ishap) {
                     allSens->dispWRTshape[ishap] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numTotalDispDofs,1);
                     allSens->dispWRTshape[ishap]->setZero();
                     int dispDofIndex = 0;
                     for(int inode=0; inode < numDispNodes; ++inode) {
                       int numDispDofs = (*dispNodes)[inode].numdofs;
                       for(int idof=0; idof < numDispDofs; ++idof) {
                         (*allSens->dispWRTshape[ishap])(dispDofIndex,0) -=
                               allSens->lambdaDisp[dispDofIndex]->adjoint()*(allSens->linearstaticWRTshape[ishap]->col(0));
                         dispDofIndex++;
                       }
                     }
                   }
                 }
                  
               }
               break;

             case SensitivityInfo::DisplacementWRTthickness:

               if( numThicknessGroups > 0 ) {  
                 if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
                   probDesc->sendNumParam(numThicknessGroups, 5, ratioSensitivityTol*sensitivityTol);
                   if(!allSens->dispWRTthick) { 
                     allSens->dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
                     for(int iparam=0; iparam< numThicknessGroups; ++iparam) {
                       allSens->dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic,
                                                                                 Eigen::Dynamic>(domain->numUncon(),1);
                       rhsSen = new VecType( probDesc->solVecInfo() );
                       rhsSen->copy(allSens->linearstaticWRTthick[iparam]->data());
                       (*rhsSen) *= -1; 
                       aeroForceSen->zero();
                       *d_nSen = 0.0;
                       aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
                       *allSens->dispWRTthick[iparam] = Eigen::Map<Eigen::Matrix<double,
                                                                   Eigen::Dynamic,Eigen::Dynamic>>(d_nSen->data(),
                                                                                                   domain->numUncon(),1);
                       delete rhsSen;
                     }
                   }
                 }
                 else if(domain->solInfo().sensitivityMethod == SolverInfo::Adjoint) {
                   
                   int numDispNodes = domain->getNumDispNodes();
                   int numTotalDispDofs = domain->getTotalNumDispDofs();
                   if(!allSens->lambdaDisp) {
                     probDesc->sendNumParam(numTotalDispDofs,
                                            8, ratioSensitivityTol*sensitivityTol); // actvar = 8, structure displacement
                     computeLambdaDisp(numTotalDispDofs, numDispNodes);
                   }
                   allSens->dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
                   std::vector<DispNode> *dispNodes = domain->getDispNodes();
                   for(int iparam=0; iparam < numThicknessGroups; ++iparam) {
                     allSens->dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(numTotalDispDofs,1);
                     allSens->dispWRTthick[iparam]->setZero();
                     int dispDofIndex = 0;
                     for(int inode=0; inode < numDispNodes; ++inode) {
                       int numDispDofs = (*dispNodes)[inode].numdofs;
                       for(int idof=0; idof < numDispDofs; ++idof) {
                         (*allSens->dispWRTthick[iparam])(dispDofIndex,0) -=
                              allSens->lambdaDisp[dispDofIndex]->adjoint()*(allSens->linearstaticWRTthick[iparam]->col(0));
                         dispDofIndex++;
                       }
                     }
                   }
                 }
                  
               }
               break;

             case SensitivityInfo::AggregatedStressVMWRTthickness:

               if( numThicknessGroups > 0 ) {
                 if(domain->solInfo().sensitivityMethod != SolverInfo::Adjoint) { 
                   filePrint(stderr, " *** Error: aggregated stress sensitivity is not supported by direct sensitivity, exiting\n"); 
                   exit(-1); 
                 }
                 if(!allSens->lambdaAggregatedStressVM) computeLambdaAggregatedStressVM();
                 for(int iparam=0; iparam< numThicknessGroups; ++iparam)  
                   (*allSens->aggregatedVonMisesWRTthick)(iparam) -=
                            allSens->lambdaAggregatedStressVM->adjoint()*(allSens->linearstaticWRTthick[iparam]->col(0));
               }
               break;

             case SensitivityInfo::AggregatedStressVMWRTshape:

               if( numShapeVars > 0 ) {
                 if(domain->solInfo().sensitivityMethod != SolverInfo::Adjoint) { 
                   filePrint(stderr, " *** Error: aggregated stress sensitivity is not supported by direct sensitivity, exiting\n"); 
                   exit(-1); 
                 }
                 if(!allSens->lambdaAggregatedStressVM) computeLambdaAggregatedStressVM(); 
                 for(int ishap=0; ishap< numShapeVars; ++ishap)  
                   (*allSens->aggregatedVonMisesWRTshape)(ishap) -=
                          allSens->lambdaAggregatedStressVM->adjoint()*(allSens->linearstaticWRTshape[ishap]->col(0));
               }
               break;

             case SensitivityInfo::StressVMWRTshape:
              
               if( numShapeVars > 0 ) {  
                 if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
                   probDesc->sendNumParam(numShapeVars, 1, ratioSensitivityTol*sensitivityTol);
                   if(!allSens->dispWRTshape) { 
                     allSens->dispWRTshape = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numShapeVars];
                     for(int ishap=0; ishap< numShapeVars; ++ishap) {
                       allSens->dispWRTshape[ishap] = new Eigen::Matrix<double,
                                                                        Eigen::Dynamic, Eigen::Dynamic>(domain->numUncon(),1);
                       rhsSen = new VecType( probDesc->solVecInfo() );
                       rhsSen->copy(allSens->linearstaticWRTshape[ishap]->data());
                       (*rhsSen) *= -1; 
                       aeroForceSen->zero();
                       *d_nSen = 0.0;
                       aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
                       *allSens->dispWRTshape[ishap] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                                       Eigen::Dynamic >>(d_nSen->data(),
                                                                                        domain->numUncon(),1);
                       delete rhsSen;
                     }
                   }
                   for(int ishap=0; ishap< numShapeVars; ++ishap) 
                     allSens->vonMisesWRTshape->col(ishap) += *allSens->vonMisesWRTdisp * (*allSens->dispWRTshape[ishap]);
                 }
                 else if(domain->solInfo().sensitivityMethod == SolverInfo::Adjoint) {
                   int numStressNodes = domain->getNumStressNodes();
                   if(!allSens->lambdaStressVM) computeLambdaStressVM(numStressNodes);
                   for(int inode = 0; inode < numStressNodes; ++inode) 
                     for(int ishap=0; ishap< numShapeVars; ++ishap) { 
                       (*allSens->vonMisesWRTshape)(inode,ishap) -=
                             allSens->lambdaStressVM[inode]->adjoint()*(allSens->linearstaticWRTshape[ishap]->col(0));
                     }
                 }
               }
               break;

             case SensitivityInfo::StressVMWRTthickness:
             
               if( numThicknessGroups > 0 ) {
                 if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) {
                   probDesc->sendNumParam(numThicknessGroups, 5, ratioSensitivityTol*sensitivityTol);
                   allSens->dispWRTthick = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*[numThicknessGroups];
                   for(int iparam=0; iparam< numThicknessGroups; ++iparam) {
                     allSens->dispWRTthick[iparam] = new Eigen::Matrix<double, Eigen::Dynamic,
                                                                                Eigen::Dynamic>(domain->numUncon(),1);
                     rhsSen = new VecType( probDesc->solVecInfo() );
                     rhsSen->copy(allSens->linearstaticWRTthick[iparam]->data());
                     (*rhsSen) *= -1;
                     aeroForceSen->zero();
                     *d_nSen = 0.0;
                     aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
                     *allSens->dispWRTthick[iparam] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                                      Eigen::Dynamic > >(d_nSen->data(),
                                                                                                        domain->numUncon(),1);
                     allSens->vonMisesWRTthick->col(iparam) += *allSens->vonMisesWRTdisp * (*allSens->dispWRTthick[iparam]);
                     delete rhsSen;
                   }  
                 }
                 else if(domain->solInfo().sensitivityMethod == SolverInfo::Adjoint) {
                   int numStressNodes = domain->getNumStressNodes();
                   if(!allSens->lambdaStressVM) computeLambdaStressVM(numStressNodes);
                   for(int inode = 0; inode < numStressNodes; ++inode) 
                     for(int iparam=0; iparam < numThicknessGroups; ++iparam) 
                       (*allSens->vonMisesWRTthick)(inode,iparam) -=
                            allSens->lambdaStressVM[inode]->adjoint()*(allSens->linearstaticWRTthick[iparam]->col(0));
                 }
               }
               break;

             case SensitivityInfo::StressVMWRTmach:

               probDesc->sendNumParam(1, 2, ratioSensitivityTol*sensitivityTol);
               filePrint(stderr,"Sensitivity with respect to Mach number will be computed\n");
               allSens->dispWRTmach = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(domain->numUncon(),1);
               rhsSen = new VecType( probDesc->solVecInfo() );
               rhsSen->zero();
               aeroForceSen->copy(1e4);
               *d_nSen = 0.0;
               aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
               *allSens->dispWRTmach = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                       Eigen::Dynamic > >(d_nSen->data(),
                                                                                         domain->numUncon(),1);
               allSens->vonMisesWRTmach->col(0) += *allSens->vonMisesWRTdisp * (*allSens->dispWRTmach);
               delete rhsSen;
               break;

             case SensitivityInfo::StressVMWRTalpha:

               probDesc->sendNumParam(1, 3, ratioSensitivityTol*sensitivityTol);
               filePrint(stderr,"Sensitivity with respect to Angle of attack will be computed\n");
               allSens->dispWRTalpha = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(domain->numUncon(),1);
               rhsSen = new VecType( probDesc->solVecInfo() );
               rhsSen->zero();
               aeroForceSen->copy(1e4);
               *d_nSen = 0.0;
               aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
               *allSens->dispWRTalpha = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                        Eigen::Dynamic > >(d_nSen->data(),
                                                                                           domain->numUncon(),1);
               allSens->vonMisesWRTalpha->col(0) += *allSens->vonMisesWRTdisp * (*allSens->dispWRTalpha);
               delete rhsSen;
               break;

             case SensitivityInfo::StressVMWRTbeta:

               probDesc->sendNumParam(1, 4, ratioSensitivityTol*sensitivityTol);
               filePrint(stderr,"Sensitivity with respect to Yaw angle will be computed\n");
               allSens->dispWRTbeta = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(domain->numUncon(),1);
               rhsSen = new VecType( probDesc->solVecInfo() );
               rhsSen->zero();
               aeroForceSen->copy(1e4);
               *d_nSen = 0.0;
               aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
               *allSens->dispWRTbeta = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                        Eigen::Dynamic> >(d_nSen->data(),
                                                                                          domain->numUncon(),1);
               allSens->vonMisesWRTbeta->col(0) += *allSens->vonMisesWRTdisp * (*allSens->dispWRTbeta);
               delete rhsSen;
               break;

             default:
               break;
           }
         }
         probDesc->sendNumParam(numThicknessGroups+numShapeVars, 0, ratioSensitivityTol*sensitivityTol);
         filePrint(stderr, " ... numThicknessGroups = %d, numShapeVars = %d\n", numThicknessGroups, numShapeVars);
         for(int isen = 0; isen < numFluidQuantTypes; ++isen) { // fluid sensitivities
           computeLambdaFluidQuantity();   
           if(!allSens->linearstaticWRTshape && numShapeVars > 0) {
             filePrint(stderr, " *** Error! please comment out readsensitivity in your input file.\n"); exit(-1); 
           }          
           for(int ishap=0; ishap < numShapeVars; ++ishap) { //TODO: Now shape variable gets priority
             double fluidQuantity(0); 
             fluidQuantity -= (allSens->lambdaFluidQuantity->adjoint())*(allSens->linearstaticWRTshape[ishap]->col(0));
             probDesc->sendRelativeResidual(fluidQuantity);
           }
           if(!allSens->linearstaticWRTthick && numThicknessGroups > 0) {
             filePrint(stderr, " *** Error! please comment out thgrli in your input file.\n");  exit(-1);
           }
           for(int iparam=0; iparam < numThicknessGroups; ++iparam) { 
             double fluidQuantity(0); 
             fluidQuantity -= (allSens->lambdaFluidQuantity->adjoint())*(allSens->linearstaticWRTthick[iparam]->col(0));
             probDesc->sendRelativeResidual(fluidQuantity);
           }
         } 
         walls += getTime();
         std::cerr << " ****** wall clock time for sensitivity solve computation is " << walls/1000.0 << std::endl;
         probDesc->sensitivityPostProcessing(d_n);
       } 
#endif
       break;
   }
   timeLoop += getTime();
   probDesc->printTimers(dynOps, timeLoop);
   
   // Delete arrays
   delete curState;
   delete workVec;

   delete d_n;
   delete v_n;
   delete a_n;
   delete constForce;
   delete v_p;

   if(domain->solInfo().sensitivity) {
     if(aeroForceSen) delete aeroForceSen;
     if(curSenState) delete curSenState;
     if(workSenVec) delete workSenVec;
     if(lambda_nSen) delete lambda_nSen;
     if(d_nSen) delete d_nSen;
     if(v_nSen) delete v_nSen;
     if(a_nSen) delete a_nSen;
     if(v_pSen) delete v_pSen;
   }

}

template< class DynOps,        class VecType,
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::computeLambdaFluidQuantity()
{
  AllSensitivities<double> *allSens = probDesc->getAllSensitivities();
  if(!allSens->lambdaFluidQuantity) allSens->lambdaFluidQuantity = new Eigen::Matrix<double, Eigen::Dynamic, 1>(domain->numUncon());
  rhsSen = new VecType( probDesc->solVecInfo() );
  rhsSen->zero();
  aeroForceSen->zero();
  (*aeroForceSen)[0] = 1.0e2;
  *lambda_nSen = 0.0;
  aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
  *allSens->lambdaFluidQuantity = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                  Eigen::Dynamic> >(lambda_nSen->data(),
                                                                                    domain->numUncon(),1);
  delete rhsSen;
}

template< class DynOps,        class VecType,
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::computeLambdaDisp(int numTotalDispDofs, int numDispNodes)
{
  AllSensitivities<double> *allSens = probDesc->getAllSensitivities();
  std::vector<DispNode> *dispNodes = domain->getDispNodes();
  allSens->lambdaDisp = new Eigen::Matrix<double, Eigen::Dynamic, 1>*[numTotalDispDofs];
  int dispDofIndex = 0;
  for(int inode = 0; inode < numDispNodes; ++inode) {
    int node = (*dispNodes)[inode].nodeID, loc;
    int numDispDofs = (*dispNodes)[inode].numdofs;
    for(int idof = 0; idof<numDispDofs; ++idof) {
      allSens->lambdaDisp[dispDofIndex] = new Eigen::Matrix<double, Eigen::Dynamic, 1>(domain->numUncon());
      allSens->lambdaDisp[dispDofIndex]->setZero();
      rhsSen = new VecType( probDesc->solVecInfo() );
      Vector lambdadisp(domain->numUncon(),0.0);
      int dof = (*dispNodes)[inode].dofs[idof];
      loc = domain->returnLocalDofNum(node, dof);
      if (loc >= 0) { rhsSen->zero(); (*rhsSen)[loc] = 1.0; }
      else continue;
      aeroForceSen->zero();
      *lambda_nSen = 0.0;
      aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
      *allSens->lambdaDisp[dispDofIndex] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                           Eigen::Dynamic> >(lambda_nSen->data(),
                                                                                             domain->numUncon(),1);
      delete rhsSen;
      dispDofIndex++;
    }
  }
}

template< class DynOps,        class VecType,
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::computeLambdaAggregatedStressVM()
{
  AllSensitivities<double> *allSens = probDesc->getAllSensitivities();
  allSens->lambdaAggregatedStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>(domain->numUncon());
  probDesc->sendNumParam(1, 6, ratioSensitivityTol*sensitivityTol); // actvar = 6, aggregated von Mises stress
  rhsSen = new VecType( probDesc->solVecInfo() );
  for(int i=0; i<domain->numUncon(); ++i) (*rhsSen)[i] = (*allSens->aggregatedVonMisesWRTdisp)(i);
  aeroForceSen->zero();
  *lambda_nSen = 0.0;
  aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
  *allSens->lambdaAggregatedStressVM = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                       Eigen::Dynamic> >(lambda_nSen->data(),
                                                                                         domain->numUncon(),1);
  delete rhsSen;
}

template< class DynOps,        class VecType,
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::computeLambdaStressVM(int numStressNodes)
{
  AllSensitivities<double> *allSens = probDesc->getAllSensitivities();
  probDesc->sendNumParam(numStressNodes, 7, ratioSensitivityTol*sensitivityTol); // actvar = 7, von Mises stress
  std::vector<int> *stressNodes = domain->getStressNodes();
  allSens->lambdaStressVM = new Eigen::Matrix<double, Eigen::Dynamic, 1>*[numStressNodes];
  for(int inode = 0; inode < numStressNodes; ++inode) {
    allSens->lambdaStressVM[inode] = new Eigen::Matrix<double, Eigen::Dynamic, 1>(domain->numUncon());
    rhsSen = new VecType( probDesc->solVecInfo() );
    for(int i=0; i<domain->numUncon(); ++i) (*rhsSen)[i] = (*allSens->vonMisesWRTdisp)((*stressNodes)[inode],i);
    aeroForceSen->zero();
    *lambda_nSen = 0.0;
    aeroSensitivityQuasistaticLoop( *curSenState, *rhsSen, *dynOps, *workSenVec, dt, tmax, aeroAlg);
    *allSens->lambdaStressVM[inode] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                                                                      Eigen::Dynamic> >(lambda_nSen->data(),
                                                                                        domain->numUncon(),1);
    delete rhsSen;
  }
}

// -----------------------------------------------------------------------------//
//
// ... Quasi-Static Time Loop
//
// -----------------------------------------------------------------------------

template< class DynOps,        class VecType, 
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::quasistaticLoop(SysState<VecType>& curState, VecType& constForce,
                  DynOps& dynOps, NewmarkWorkVec<VecType,ProblemDescriptor>& workVec,
                  double dt, double tmax, int aeroFlg)
{ 
  filePrint(stderr, " ... Quasistatic loop               ... \n");
  if (aeroFlg == 10) steadyMax = 1;
  else if(steadyMax == 1) {
    filePrint(stderr, " ... Error! Maximum number of iterations must be greater than 1\n");
    exit(-1);
  }
   
  // get initial displacements
  VecType &d_n = curState.getDisp();
   
  // Initialize some vectors 
  VecType  &d_inc = workVec.get_d_n_p();
  VecType    &rhs = workVec.get_rhs();
  VecType  &ext_f = workVec.get_ext_f();
  VecType   &toto = workVec.get_o_n_p();

  // Initialize some parameters
  // Get initial Time
  int tIndex;
  int initIndex;
  double initialTime = 0.0;

  probDesc->getInitialTime(initIndex, initialTime);
  double initExtForceNorm = probDesc->getInitialForceNorm();
  tIndex = initIndex;

  int iSteady  = 0;

  double forceRef;

  double relaxFac = maxVel;
   
  // Output state of model

  postProcessor->dynamOutput( tIndex, initialTime, dynOps, ext_f, aeroForce, curState );

  //-----------------------------------------------------------------------
  // ... BEGIN MAIN TIME-LOOP
  //-----------------------------------------------------------------------

  double totalTime = -getTime();

  for (tIndex = tIndex+1; tIndex <= steadyMax; tIndex++) {

    // ... call projector for RBMs in case of rbmfilter level 2
    if (probDesc->getFilterFlag() == 2) probDesc->project( d_n );

    // ... compute external force (and receive fluid load)
    probDesc->computeExtForce2( curState, ext_f, constForce, tIndex, (double)tIndex*delta, aeroForce, 1.0, 0.0 );

    // ... build force reference norm 
    if (tIndex==initIndex+1) {
      if (initExtForceNorm == 0.0)
        forceRef=ext_f.norm();
      else 
        forceRef = initExtForceNorm;
      if(verboseFlag) filePrint(stderr, " ... Initial Force: %8.2e        ...\n", forceRef);
    }

    // ... build internal force 
    getInternalForce(dynOps, d_n, rhs, (double)tIndex*delta, tIndex);

    // ... check for convergence
    double relres = 0.0;
    toto = rhs-ext_f;
    if (forceRef != 0.0) relres = toto.norm()/forceRef;
    else {
      relres = toto.norm();
      filePrint(stdout, " ... WARNING: Reference External Force is zero, Relative residual is absolute error norm ...\n");
    }

    if(relres <= steadyTol && delta == 0) iSteady = 1;

    if(aeroAlg >= 0 || probDesc->getAeroheatFlag() >= 0) {
      filePrint(stderr," ... Pseudo-Step = %d  Rel. Res. = %10.4e ...\n",tIndex, relres);
  
      // command communication with fluid
      if(tIndex == steadyMax && !iSteady) {
        probDesc->processLastOutput();
        if(aeroFlg != 10) {
          postProcessor->dynamOutput( tIndex, (double)tIndex*delta, dynOps, ext_f, aeroForce, curState );
          probDesc->cmdCom(1);
          break;
        }
      }
      else {
        iSteady = probDesc->cmdCom(iSteady);
      }
    }

    // ... stop quasi-transient simulation if converged
    if(iSteady) {
      filePrint(stderr," ------------------------------------------------------\n");
      filePrint(stderr," ... Quasistatic Analysis Converged After %3d Steps ...\n",tIndex);
      filePrint(stderr," ------------------------------------------------------\n");
      probDesc->processLastOutput();
      postProcessor->dynamOutput( tIndex, (double)tIndex*delta, dynOps, ext_f, aeroForce, curState );
      break; 
    }
    else if (tIndex == steadyMax)
      probDesc->processLastOutput();

    // ... save load vector
    rhs = ext_f;

    if(domain->solInfo().isNonLin()) {

      // ... solve nonlinear system for current load, compute displacement increment
      //     and update solution, applying relaxation factor.
      probDesc->solveAndUpdate(rhs, d_inc, d_n, relaxFac, (double)tIndex*delta); // XXX consider delta != 0
    }
    else {

      // ... solve linear system for current load
      dynOps.dynMat->reSolve( rhs );

      // ... compute displacement increment
      d_inc.linC(rhs, -1.0, d_n);

      // ... update solution, applying relaxation factor
      d_n += relaxFac*d_inc;
    }

    // ... output current solution and send displacements to fluid
    postProcessor->dynamOutput( tIndex, (double)tIndex*delta, dynOps, ext_f, aeroForce, curState );
  }

  if (!iSteady && aeroAlg != 10) {
    filePrint(stderr," -------------------------------------------------------------\n");
    filePrint(stderr," ... Quasistatic Analysis Did Not Converge After %3d Steps ...\n",tIndex);
    filePrint(stderr," -------------------------------------------------------------\n");
  }

  // ... output CPU time spent in quasi-static loop
  totalTime += getTime();
#ifdef PRINT_TIMERS
  filePrint(stderr," ... Total Loop Time = %.2e s   ...\n",totalTime/1000.0);
#endif
}

// -----------------------------------------------------------------------------//
//
// ... Aeroelastic sensitivity Quasi-Static Time Loop
//
// -----------------------------------------------------------------------------

template< class DynOps,        class VecType, 
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::aeroSensitivityQuasistaticLoop(SysState<VecType>& curState, VecType& constForce,
                                 DynOps& dynOps, NewmarkWorkVec<VecType,ProblemDescriptor>& workVec,
                                 double dt, double tmax, int aeroFlg)
{ 
   filePrint(stderr, " ... Aeroelastic Sensitivity Quasistatic loop ... \n");
   if (aeroFlg == 10)  steadyMax = 1;
   else if(steadyMax == 1) {
     filePrint(stderr, " ... Error! Maximum number of iterations must be greater than 1\n");
     exit(-1);
   }
 
   // get initial displacements
   VecType &d_nSen = curState.getDisp();
   
   // Initialize some vectors 
   VecType  &d_n_pSen = workVec.get_d_n_p();
   VecType    &rhsSen = workVec.get_rhs();
   VecType  &ext_fSen = workVec.get_ext_f();

   // Initialize some parameters
   // Get initial Time
   int tIndex;
   int initIndex;
   double initialTime = 0.0;

   probDesc->getInitialTime(initIndex, initialTime);
   double initExtForceNorm = probDesc->getInitialForceNorm();
   tIndex = initIndex;

   int iSteady  = 0;

   double forceSenRef;

   double relaxFac = domain->solInfo().qsMaxvelSen;
   
   // Output state of model

   //-----------------------------------------------------------------------
   // ... BEGIN MAIN TIME-LOOP
   //-----------------------------------------------------------------------

   double totalTime = -getTime();

  for (tIndex = tIndex+1; tIndex <= steadyMax; tIndex++) {

    // ... call projector for RBMs in case of rbmfilter level 2
    if (probDesc->getFilterFlag() == 2) probDesc->project( d_nSen );

    // ... compute external force
    ext_fSen = constForce;
    if(domain->solInfo().sensitivityMethod == SolverInfo::Direct) ext_fSen += *aeroForceSen;
    else ext_fSen -= *aeroForceSen;

    // ... build force reference norm 
    if (tIndex==initIndex+1) {
      forceSenRef=ext_fSen.norm();
      filePrint(stderr, " ... Initial Force Norm: %10.6e     ...\n", forceSenRef);
    }

    // ... build internal force 
    probDesc->getInternalForce(d_nSen, rhsSen, (double)tIndex*delta, tIndex);

    // ... check for convergence
    double relres = 0.0;
    if (forceSenRef != 0.0) relres = norm(rhsSen-ext_fSen)/forceSenRef;
    else {
      relres = norm(rhsSen-ext_fSen);
      filePrint(stderr, " *** WARNING: Reference external force is zero, relative residual is absolute error norm\n");
    }

    if(relres <= sensitivityTol && delta == 0 && tIndex != 1) {
      // ... stop quasi-transient simulation if converged
      iSteady = 1;
      probDesc->cmdCom(iSteady);
      filePrint(stderr," ---------------------------------------------------------------------------\n");
      filePrint(stderr," ... Aeroelastic Sensitivity Quasistatic Analysis Converged After %d Steps ...\n",tIndex);
      filePrint(stderr," ---------------------------------------------------------------------------\n");
      break; 
    } else {
      probDesc->cmdCom(iSteady);
      // ... save load vector
      rhsSen=ext_fSen;

      // ... solve System for current load
      dynOps.dynMat->reSolve( rhsSen );

      // ... compute displacement increment;
      d_n_pSen.linC(rhsSen, -1.0, d_nSen);

      // ... apply relaxation factor
      d_n_pSen *= relaxFac;

      // ... update solution
      d_nSen   += d_n_pSen;

      probDesc->sendRelativeResidual(relres);
      probDesc->sendDisplacements(d_nSen, curState.getVeloc(), curState.getAccel(), curState.getPrevVeloc());
    }

    if(aeroAlg >= 0 || probDesc->getAeroheatFlag() >= 0) {
      filePrint(stderr," ... Pseudo-Step = %d  Rel. Res. = %10.4e ...\n",tIndex, relres);
 
      probDesc->getAeroelasticForceSensitivity(tIndex, (double)tIndex*delta, aeroForceSen, 1.0, 0.0); // [E] Received fluid load sensitivity...
      // command communication with fluid
      if(tIndex == steadyMax && !iSteady) { 
        if(aeroFlg != 10) {
          probDesc->cmdCom(1);
          break;
        }
      }
    }

  }

  if (!iSteady && aeroAlg != 10) {
    filePrint(stderr," ------------------------------------------------------------------------\n");
    filePrint(stderr," ... Aeroelastic Sensitivity Analysis Did Not Converge After %d Steps ...\n",tIndex);
    filePrint(stderr," ------------------------------------------------------------------------\n");
  }
//  probDesc->processLastOutput();
//  postProcessor->dynamOutput( tIndex, (double)tIndex*delta, dynOps, ext_fSen, aeroForceSen, curState );

  // ... output CPU time spent in quasi-static loop
  totalTime += getTime();
#ifdef PRINT_TIMERS
  filePrint(stderr," ... Total Loop Time = %.2e s   ...\n",totalTime/1000.0);
#endif
}

// -----------------------------------------------------------------------------//
// ... GENERAL IMPLICIT NEWMARK TIME LOOP
//
// DESCRIPTION: this routine implements the implicit Generalized Alpha method 
//              for time integration solution of the second-order differential equation:
//              M*(d^2u/dt^2) + C*(du/dt) + K*u = fext(u)
//
// WARNINGS   : 1) Second-order accuracy is obtained if and only if gamma = 1/2
//              2) Unconditionally stable when 2*beta >= gamma >= 1/2
//              3) This function does not work for beta == 0.0 ... we use the
//                 explicitNewmarkLoop function in this case
//
//  See J. Chung and G.M. Hulbert, "A time integration algorithm for structural
//  dynamics with improved numerical dissipation: the Generalized Alpha Method",
//  Journal of Applied Mechanics, no 60, 1993, pp 371-375.
//
//  This is a more general method that uses two parameters, alphaf and alpham to 
//  control numerical disspation.  In particular, it allows the user to limit
//  high frequency disspation while minimizing unwanted low frequency disspation.
//  The Newmark algorithm becomes a particular case of this method, achieved with
//  alphaf = alpham = 0.
// -----------------------------------------------------------------------------
template< class DynOps,        class VecType, 
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::implicitNewmarkLoop(SysState<VecType>& curState, VecType& constForce,
                      DynOps& dynOps, 
		      NewmarkWorkVec<VecType,ProblemDescriptor>& workVec,
                      double dt, double tmax)
{
   MatrixTimers &matrixTimers = domain->getTimers();
   if(domain->solInfo().order == 1) {
     if(gamma == 0.5)
       filePrint(stderr, " ... Implicit Midpoint Rule         ..."
                         " ... i.e. α = ½                     ...\n");
     else if(gamma == 1.0)
       filePrint(stderr, " ... Implicit Backward Euler Method ..."
                         " ... i.e. α = 1                     ...\n");
     else
       filePrint(stderr, " ... Imp. Generalized Midpoint Rule ...\n"
                         " ... with α = %5.3f                 ...\n", gamma);
   }
   else {
     if(beta == 0.25 && gamma == 0.5 && alphaf == 0.5 && alpham == 0.5)
       filePrint(stderr, " ... Implicit Midpoint Rule         ...\n"
                         " ... i.e. β = ¼, γ = ½, αf = αm = ½ ...\n");
     else if(beta == 0.25 && gamma == 0.5 && alphaf == 0 && alpham == 0)
       filePrint(stderr, " ... Implicit Newmark Method        ...\n"
                         " ... i.e. β = ¼, γ = ½, αf = αm = 0 ...\n");
     else
       filePrint(stderr, " ... Implicit Generalized-α Method  ...\n"
                         " ... with β  = %5.3f, γ  = %5.3f    ...\n"
                         " ...      αf = %5.3f, αm = %-6.3f   ...\n", beta, gamma, alphaf, alpham);
   }

   int parity = 0;
   SysState<VecType> *bkState = 0;
   // Allocate backup state for A5 algorithm
   VecType *d_bk, *v_bk, *a_bk, *v_p_bk;
   if(aeroAlg == 5) {
     d_bk = new VecType(probDesc->solVecInfo());
     v_bk = new VecType(probDesc->solVecInfo());
     a_bk = new VecType(probDesc->solVecInfo());
     v_p_bk = new VecType(probDesc->solVecInfo());
     bkState = new SysState<VecType>(*d_bk, *v_bk, *a_bk, *v_p_bk);
   }

   VecType &d_n = curState.getDisp();
   VecType &v_n = curState.getVeloc();
   VecType &a_n = curState.getAccel();
   VecType &v_p = curState.getPrevVeloc();

   // project initial displacements in case of rbmfilter
   if(probDesc->getFilterFlag() > 0) {
      probDesc->project( d_n );
      probDesc->project( v_n );
   }

   // Some vectors for implicit Newmark time loop
   VecType  &d_n_p = workVec.get_d_n_p();
   VecType  &v_n_p = workVec.get_v_n_p();
   VecType  &a_n_p = workVec.get_a_n_p();
   VecType   & rhs = workVec.get_rhs();
   VecType  &ext_f = workVec.get_ext_f();
   VecType  &d_n_h = workVec.get_d_n_h();
   VecType  &v_n_h = workVec.get_v_n_h();
   VecType &Md_n_h = workVec.get_Md_n_h();
   VecType &Cd_n_h = workVec.get_Cd_n_h();
   VecType   &tmp1 = workVec.get_tmp1();

   // Get initial time and time index
   int n = 0;
   double t = 0.0;
   probDesc->getInitialTime(n, t);

   // Compute the external force at the initial time: fext^0
   probDesc->computeExtForce2(curState, ext_f, constForce, -1, t, aeroForce, gamma, alphaf);

   // Compute the initial acceleration: a^0 = M^{-1}(fext^0 - Ku^0 - Cv^0)
   if(domain->solInfo().iacc_switch && dynOps.Msolver) {
     if(domain->solInfo().order == 1) {
       if(verboseFlag) filePrint(stderr," ... Computing initial first time derivative of temperature ...\n");
       matrixTimers.formRhs -= getTime();
       dynOps.K->mult(d_n, tmp1);
       v_n = ext_f - tmp1;
       matrixTimers.formRhs += getTime();
       dynOps.Msolver->reSolve(v_n);
     }
     else {
       if(verboseFlag) filePrint(stderr," ... Computing initial acceleration ...\n");
       matrixTimers.formRhs -= getTime();
       dynOps.K->mult(d_n, tmp1);
       a_n = ext_f - tmp1;
       if(dynOps.C) {
         dynOps.C->mult(v_n, tmp1);
         a_n -= tmp1;
       }
       matrixTimers.formRhs += getTime();
       dynOps.Msolver->reSolve(a_n);
       if(probDesc->getFilterFlag() == 2) probDesc->project(a_n);
     }
   }

   // Output initial state of model
   postProcessor->dynamOutput(n, t, dynOps, ext_f, aeroForce, curState);

   // ... BEGIN MAIN TIME-LOOP
   double s0 = -getTime(), s1 = -51, s2 = 0;
   char ch[4] = { '|', '/', '-', '\\' };
   int printNumber = (solInfo.printNumber > 0) ? solInfo.printNumber : std::numeric_limits<int>::max();
   if(aeroAlg < 0 && printNumber < std::numeric_limits<int>::max()) filePrint(stderr, " ⌈\x1B[33m Time Integration Loop In Progress: \x1B[0m⌉\n");

   for( ; t < tmax-0.01*dt; t += dt, s2 = s0+getTime()) {
     if(aeroAlg < 0 && (s2-s1 > printNumber)) {
       s1 = s2;
       filePrint(stderr, "\r ⌊\x1B[33m %c t = %9.3e Δt = %8.2e %3d%% \x1B[0m⌋",
                 ch[int(s1/250)%4], t+dt, dt, int((t+dt)/(tmax-0.01*dt)*100));
     }

     // ... For Aeroelastic A5 Algorithm, Do restore and backup here
     if(aeroAlg == 5) probDesc->a5StatusRevise(parity, curState, *bkState);

     // Mode decomposition of displacement
     if(probDesc->getModeDecompFlag()) probDesc->modeDecomp(t, n, d_n);

     // ... Construct force vector, includes time-independent constant force

     // ... Compute external force at time t+dt*(1-alphaf)
     // (if alphaf=0.5, external force is at time t+0.5*dt)
     probDesc->computeExtForce2(curState, ext_f, constForce, n, 
                                t+(dt*(1-alphaf)), aeroForce, gamma, alphaf);

     if(domain->solInfo().order == 1) { // heat (CURRENTLY ONLY IMPLEMENTED FOR alpham = alphaf = 1/2)
       // Solve for temperature: d^{n+1/2} = (M + gamma*dt*K)^{-1}(gamma*dt*f^{n+1/2} + M*(d^n+dt/2*(1-2*gamma)*v^n))
       matrixTimers.formRhs -= getTime();
       d_n_h = 1.0*d_n + (dt/2.0*(1.0-2.0*gamma))*v_n;
       dynOps.M->mult( d_n_h, Md_n_h );
       rhs = 1.0*Md_n_h + gamma*dt*ext_f;
       matrixTimers.formRhs += getTime();
       dynOps.dynMat->reSolve( rhs );

       // Extrapolate temperature solution to t^{n+1} : d^{n+1} = 2*d^{n+1/2} - d^n
       matrixTimers.updateState -= getTime();
       d_n = 2.0*rhs - 1.0*d_n;

       // Compute the first time derivative of temperature at t^{n+1}: v^{n+1} = 2/(gamma*dt)*(d^{n+1/2} - d^n) - (1-gamma)/(gamma)*v^n
       v_n_p = 2/(gamma*dt)*(d_n - rhs);
       if(gamma != 1.0) v_n_p -= (1.0-gamma)/gamma*v_n;
     }
     else { // mech, acou
       // ... Construct R.H.S. vector
       matrixTimers.formRhs -= getTime();
       // First: d_n_h = ((1-alpham)/(1-alphaf))*d_n + dt*(1-alpham)*v_n + dt*dt*((1-alpham)/2 - beta)*a_n
       d_n_h = ((1-alpham)/(1-alphaf))*d_n + dt*(1-alpham)*v_n + dt*dt*((1-alpham)/2 - beta)*a_n;

       // Second: Multiply by Mass Matrix M
       dynOps.M->mult( d_n_h, Md_n_h );

       // Third: Accumulate in rhs vector: rhs = Md_n_h + beta*dt*dt*ext_f
       rhs = 1.0*Md_n_h + beta*dt*dt*ext_f;

       if(dynOps.C) {
         // Fourth: d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma(1-alphaf))*v_n - dt*dt*dt*0.5(1-alphaf)*(2*beta - gamma)*a_n
         d_n_h = dt*gamma*d_n - dt*dt*(beta-gamma*(1-alphaf))*v_n - dt*dt*dt*0.5*(1-alphaf)*(2*beta - gamma)*a_n;

         // Fifth: Multiply by Damping Matrix C
         dynOps.C->mult( d_n_h, Cd_n_h );
  
         // Sixth: Accumulate in rhs vector 
         rhs += Cd_n_h;
       }
       matrixTimers.formRhs += getTime();

       dynOps.dynMat->reSolve( rhs ); // Now rhs contains d_(n+1-alphaf)

       // call projector for RBMs in case of rbmfilter level 2
       if (probDesc->getFilterFlag() == 2) probDesc->project( rhs );

       // one time step forward
       matrixTimers.updateState -= getTime();
       // d_n_p = 1/(1-alphaf)*[d_(n+1-alphaf)-alphaf*d_n] = d_(n+1)
       d_n_p = 1/(1-alphaf)*(1.0*rhs-alphaf*d_n);
   
       // a_n_p = 1/(dt^2*beta)*[d_(n+1)-d_n] - 1/(dt*beta)*v_n + (1-1/(2*beta))*a_n = a_(n+1)
       a_n_p = 1/(dt*dt*beta)*(d_n_p-d_n) -1/(dt*beta)*v_n + (1-1/(2*beta))*a_n;

       // v_n_h = gamma/(beta*dt)*[d_(n+1-alphaf) - d_n] + (1.0-(1.0-alphaf)*gamma/beta)*v_n + dt*(1.0-alphaf)*(2.0*beta-gamma)/(2*beta)*a_n
       v_n_h = gamma/(beta*dt)*(rhs-d_n) + (1-(1-alphaf)*gamma/beta)*v_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*a_n;

       // v_n_p = 1/(1-alphaf)*(v_n_h - alphaf*v_n) = v_(n+1)
       v_n_p = 1/(1-alphaf)*(1.0*v_n_h - alphaf*v_n);
     
       // Now swap v_n_p -> v_n and d_n_p -> d_n
       v_p = v_n;
       d_n.swap( d_n_p );
       a_n.swap( a_n_p );
     }
     v_n.swap( v_n_p );

     // Increment time index
     n++;
     matrixTimers.updateState += getTime();

     // ... current state is replaced by predicted value 
     // NOTE: current state is modified here only for output
     // it will be restored to the backup state at the
     // beginning of the next iteration
     if(!parity && aeroAlg == 5) {
       curState.getDisp().linC(0.5, curState.getDisp(), 0.5, bkState->getDisp());
       curState.getVeloc().linC(0.5, curState.getVeloc(), 0.5, bkState->getVeloc());
     }
     else {
       // FORCE PRINTING AT LAST ITERATION
       if(t+1.01*dt > tmax)  probDesc->processLastOutput();
     }

     postProcessor->dynamOutput(n, t+dt, dynOps, ext_f, aeroForce, curState);

     // ... For A5 Algorithm, do one time step back if necessary
     // add n so that time index is wound back as well as t
     if(aeroAlg == 5) probDesc->a5TimeLoopCheck( parity, t, dt );

   }
   if(aeroAlg < 0 && printNumber < std::numeric_limits<int>::max())
     filePrint(stderr, "\r ⌊\x1B[33m   t = %9.3e Δt = %8.2e 100%% \x1B[0m⌋\n", t, dt);

   if(aeroAlg == 5) {
     delete d_bk;
     delete v_bk;
     delete a_bk;
     delete v_p_bk;
     delete bkState;
   }

#ifdef PRINT_TIMERS
   if(verboseFlag) filePrint(stderr, " ... Total Loop Time = %.2e s   ...\n", s2/1000.0);
#endif
}

// -----------------------------------------------------------------------------
//
// ... CENTRAL DIFFERENCE TIME LOOP
//
// DESCRIPTION: this routine implements the explicit central difference time-integrator
//              for the solution of the second-order differential equation:
//              (a) LINEAR mech or acou: M*(d^2u/dt^2) + C*(du/dt) + K*u = fext(u)
//              (b) NONLINEAR mech: M*(d^2u/dt^2) + C*(du/dt) + fint(u) = fext(u)
//
// WARNINGS:    1. The mass matrix is automatically LUMPED (in main.C) unless MRATIO is set to 1
//              2. Viscous damping is supported, but to keep the scheme explicit the equilibrium 
//                 condition is expressed as M*a^{n+1} + C*v^{n+1/2} + fint(u^{n+1}) = fext^{n+1}
//                 where v^{n+1/2} = v^n + dt/2*a^n
//              3. Constraints must be enforced with the penalty method, except for contact/tied
//                 surfaces with multipliers which are supported using ACME
//
// -----------------------------------------------------------------------------
 
template< class DynOps,        class VecType, 
          class PostProcessor, class ProblemDescriptor,
          class Scalar>
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar>
::explicitNewmarkLoop(SysState<VecType>& curState, VecType& constForce, 
                      DynOps& dynOps, 
                      NewmarkWorkVec<VecType,ProblemDescriptor>& workVec,
                      double dt0, double tmax)
{
  MatrixTimers &matrixTimers = domain->getTimers();
  filePrint(stderr, " ... Exp. Central Difference Method ...\n"
                    " ... i.e. β = 0, γ = ½, αf = αm = 0 ...\n");

  int parity = 0;
  SysState<VecType> *bkState = 0;
  // Allocate backup state for A5 algorithm
  if(aeroAlg == 5) {
    VecType d_bk(probDesc->solVecInfo());
    VecType v_bk(probDesc->solVecInfo());
    VecType a_bk(probDesc->solVecInfo());
    VecType v_p_bk(probDesc->solVecInfo());
    bkState = new SysState<VecType>(d_bk, v_bk, a_bk, v_p_bk);
  }

  // Vectors we will need to use
  VecType &v_n = curState.getVeloc();
  VecType &v_p = curState.getPrevVeloc();
  VecType &d_n = curState.getDisp();
  VecType &a_n = curState.getAccel();
  VecType &v_n_h = workVec.get_v_n_h();
  VecType &fext = workVec.get_ext_f();
  VecType &fint = workVec.get_fint();
  VecType &tmp1 = workVec.get_tmp1();
  VecType &tmp2 = workVec.get_tmp2();
  VecType v_h_p(probDesc->solVecInfo());
  VecType *fext_p, *fint_p;

  if(domain->solInfo().check_energy_balance) {
    fext_p = new VecType(probDesc->solVecInfo());
    fint_p = new VecType(probDesc->solVecInfo());
  }
  double Wext = 0, Wint = 0, Wkin = 0;

  // We consider here an algorithm with a variable timestep
  double dt_n_h, dt_old; // deltat^{n+1/2} = t^{n+1} - t^n
  dt_old = dt_n_h = dt0;

  // project initial displacements in case of rbmfilter
  if(probDesc->getFilterFlag() > 0) {
    probDesc->project(d_n);
    probDesc->project(v_n);
  }
  v_n_h = v_n;
      
  handleDisplacement(*probDesc, d_n);
  handleVelocity(*probDesc, d_n);

  // Initialize time index (n) and time (t^n)
  int n = 0;
  double t_n = 0.0; // t^n
  double t_n_h; // t^{n+1/2} = 1/2*(t^{n+1}+t^n) = t^n + 1/2*deltat^{n+1/2}
  probDesc->getInitialTime(n, t_n);
  double t0 = t_n;

  // Get initial external force vector fext^0
  probDesc->computeExtForce2(curState, fext, constForce, n, t_n, aeroForce, 0.5, 0.0);

  // Compute the initial internal forces fint^0
  getInternalForce(dynOps, d_n, fint, t_n, n);

  // Compute the initial acceleration a^0 = M^{-1}(fext^0 - fint^0 - C*v^0)
  // note: for restarted nonlinear, the initial acceleration is read from the restart file
  if(solInfo.iacc_switch && !(geoSource->getCheckFileInfo()->lastRestartFile && domain->solInfo().isNonLin())) {
    if(verboseFlag) filePrint(stderr," ... Computing initial acceleration ...\n");
    domain->getTimers().formRhs -= getTime();
    if(dynOps.C) {
      dynOps.C->mult(v_n,tmp2);
      fint += tmp2;
    }
    a_n = fext - fint;
    domain->getTimers().formRhs += getTime();
    handleForce(*probDesc, fint);
    dynOps.dynMat->reSolve(a_n);

    if(domain->tdenforceFlag()) { // Contact corrector step: a^0 += M^{-1}*Fctc
      tmp1.linC(dt_n_h, v_n, 0.5*dt_n_h*dt_n_h, a_n); // predicted displacement d^1 = d^0 + dt^{1/2}*v^0 + dt^{1/2}*dt^{1/2}/2*a^0
      probDesc->getContactForce(d_n, tmp1, tmp2, t_n+dt_n_h, dt_n_h, dt_old);
      dynOps.dynMat->reSolve(tmp2);
      a_n += tmp2;
    }
    if(probDesc->getFilterFlag() == 2) probDesc->project(a_n);
  }
  handleAcceleration(*probDesc, a_n);

  // Output the state at t^0: d^0, v^0, a^0, fext^0
  postProcessor->dynamOutput(n, t_n, dynOps, fext, aeroForce, curState);

  // for using the modified wave equation in the acoustic problemType
  double coeff = 0.0;
  double coeff2 = -dt*dt/(probDesc->getDomain()->solInfo().modifiedWaveEquationCoef);
  bool fourthOrder = probDesc->getDomain()->solInfo().modifiedWaveEquation;

  // ... BEGIN MAIN TIME-LOOP
  double s0 = -getTime(), s1 = -51, s2 = 0;
  char ch[4] = { '|', '/', '-', '\\' };

  double eps1 = domain->solInfo().epsilon1;
  double eps2 = domain->solInfo().epsilon2;
#ifdef PRINT_ENERGIES
  std::ofstream energies("energies");
  if(domain->solInfo().isNonLin() && domain->solInfo().check_energy_balance)
    energies << "n " << "time " << "       Wkin           " << "   Wext          " << "    Wint         " << "   Wdis         " << "    Sum(W_i)     " 
             << " abs(Wkin+Wint-Wext) " << " eps1*max(We,Wi,Wk) " << "  dt " << std::endl;
#endif
  int printNumber = (solInfo.printNumber > 0) ? solInfo.printNumber : std::numeric_limits<int>::max();
  if(aeroAlg < 0 && printNumber < std::numeric_limits<int>::max()) filePrint(stderr, " ⌈\x1B[33m Time Integration Loop In Progress: \x1B[0m⌉\n");

  for( ; t_n < tmax-0.01*dt_n_h && !domain->solInfo().stop_AeroS; s2 = s0+getTime()) {

    // Time update:
    t_n_h = t_n + dt_n_h/2; // t^{n+1/2} = t^n + 1/2*deltat^{n+1/2}

    if(aeroAlg < 0 && (s2-s1 > printNumber)) {
      s1 = s2;
      filePrint(stderr, "\r ⌊\x1B[33m %c t = %9.3e Δt = %8.2e %3d%% \x1B[0m⌋",
                ch[int(s1/250)%4], t_n, dt_n_h, int((t_n-t0)/((tmax-t0)-0.01*dt_n_h)*100));
    }

    // First partial update nodal velocities:
    matrixTimers.updateState -= getTime();
    v_h_p = v_n_h;
    v_n_h.linC(v_n, t_n_h-t_n, a_n); // v^{n+1/2} = v^n + (t^{n+1/2} - t^n)*a^n
    matrixTimers.updateState += getTime();

    if (fourthOrder) { // TODO check this for case of variable timestep
      // this is as in the previous release of the FEM code (before august 28th 2008)
      //  the variables have been change to match the one declared above.
      // the alpha rayleigh damping is probably neglected
      //
      // External force
      probDesc->computeExtForce2(curState, fext, constForce, n+1, t_n+dt_n_h, aeroForce, 0.5, 0.0);

      d_n.linAdd(dt_n_h,v_n_h);
      handleDisplacement(*probDesc, d_n);

      // Internal force
      getInternalForce(dynOps, d_n, fint, t_n+dt_n_h, n+1);

      // 4th order loop see Cohen et al, Finite Elements in Analysis and Design 16(1994) pp 329-336
      // "Higher-order finite elements with mass lumping for the 1D wave equation"
      dynOps.dynMat->reSolve(fint);
      tmp1.linC(1.0,d_n,coeff2,fint);
      dynOps.K->mult(tmp1,fint);
      if (dynOps.C) {
        tmp1.linC(1.0,fext,-1.0,fint);
        dynOps.C->mult(v_n_h,tmp2);
        tmp1.linC(2.0,tmp1,-1.0,tmp2);
        dynOps.M->mult(v_n_h,tmp2);
        coeff=0.5*dt_n_h;
        v_n_h.linC(1.0,tmp2,coeff,tmp1);
        dynOps.dynMat->reSolve(v_n_h);
      } else {
        coeff=dt_n_h;
        a_n.linC(1.0,fext,-1.0,fint);

        dynOps.dynMat->reSolve(a_n);

        v_n_h.linC(1.0,v_n_h,coeff,a_n);
      }
      
      // Update the time index
      n += 1;

      // Output the state at t^n: d^n, v^n, a^n, fext^n
      postProcessor->dynamOutput(n, t_n+dt_n_h, dynOps, fext, aeroForce, curState);

    }
    else {
      if(aeroAlg == 5) probDesc->a5StatusRevise(parity, curState, *bkState);

      // Mode decomposition of displacement
      if(probDesc->getModeDecompFlag()) probDesc->modeDecomp(t_n, n, d_n);

      // Update the displacement at t^(n+1): d^{n+1} = d^n + dt^{n+1/2}*v^{n+1/2}
      matrixTimers.updateState -= getTime();
      if(domain->solInfo().isNonLin()) {
        probDesc->updateState(dt_n_h, v_n_h, d_n);
      }
      else d_n.linAdd(dt_n_h, v_n_h);
      matrixTimers.updateState += getTime();

      // C0: Send predicted displacement at t^{n+1.5} to fluid
      if(aeroAlg == 20 && !domain->solInfo().stop_AeroS) {
        int subcycle = domain->solInfo().subcycle;
        if((n + 1)%subcycle == 0) {
          probDesc->aeroSend(t_n_h+dt_n_h, d_n, v_n_h, a_n, v_h_p);
        }
      }

      // Compute the external force at t^{n+1}
      if(domain->solInfo().check_energy_balance) *fext_p = fext;
      
      probDesc->computeExtForce2(curState, fext, constForce, n+1, t_n+dt_n_h, aeroForce, 0.5, 0.0);
      handleDisplacement(*probDesc, d_n);

      // Compute the internal force at t^{n+1}
      if(domain->solInfo().check_energy_balance) *fint_p = fint;

      getInternalForce(dynOps, d_n, fint, t_n+dt_n_h, n+1);

      // Compute the acceleration at t^{n+1}: a^{n+1} = M^{-1}(fext^{n+1}-fint^{n+1}-C*v^{n+1/2})
      domain->getTimers().formRhs -= getTime();
      if(dynOps.C) {
         dynOps.C->mult(v_n_h, tmp1);
         fint.linAdd(1.0, tmp1);
      }
      a_n.linC(1.0, fext, -1.0, fint);
      domain->getTimers().formRhs += getTime();
      handleForce(*probDesc, fint);

      if(domain->solInfo().isNonLin()) probDesc->pull_back(a_n);
      dynOps.dynMat->reSolve(a_n);
      if(domain->solInfo().isNonLin()) probDesc->push_forward(a_n);

      if(domain->tdenforceFlag()) { // Contact corrector step
        for(int j = 0; j < domain->solInfo().tdenforceMaxItr; ++j) {
          tmp1.linC(dt_n_h, v_n_h, dt_n_h*dt_n_h, a_n); // predicted displacement d^{n+2} = d^{n+1} + dt^{n+1/2}*(v^{n+1/2} + dt^{n+1/2}*a^{n+1})
          probDesc->getContactForce(d_n, tmp1, tmp2, t_n+2*dt_n_h, dt_n_h, dt_old);
          double tmp2norm = tmp2.norm();
          if(contactPrintFlag && tmp2norm > 0) filePrint(stderr, "\n TDEnforcement: it = %d, ctc force = %e", j, tmp2.norm());
          if(tmp2norm < domain->solInfo().tdenforceTolAbs) break;
          dynOps.dynMat->reSolve(tmp2);
          a_n += tmp2;
        }
      }
      if(probDesc->getFilterFlag() == 2) probDesc->project(a_n);
      handleAcceleration(*probDesc, a_n);

      // Update the velocity at t^{n+1}: v^{n+1} = v^{n+1/2}+dt^{n+1/2}/2*a^{n+1}
      matrixTimers.updateState -= getTime();
      v_p = v_n;
      v_n.linC(1.0, v_n_h, 0.5*dt_n_h, a_n);
      matrixTimers.updateState += getTime();
      handleVelocity(*probDesc, v_n);

      // Energy balance check
      if(domain->solInfo().isNonLin() && domain->solInfo().check_energy_balance) {

        // 1. compute the kinetic energy
        if(domain->solInfo().activatePodRom)
          Wkin = 0.5*(v_n*v_n);
        else {
          dynOps.M->mult(v_n,tmp1); // tmp1 = M*v_n
          Wkin = 0.5*(v_n*tmp1);
        }

        // 2. compute the external energy. 
        tmp2 = *fext_p + fext;
        Wext += (0.5*dt_n_h)*(v_n_h*tmp2);

        // 3. compute the internal energy
        tmp2 = *fint_p + fint;
        Wint += (0.5*dt_n_h)*(v_n_h*tmp2);

#ifdef PRINT_ENERGIES
        // 4. compute the dissipation
        double Wdis;
        if(domain->solInfo().activatePodRom)
          Wdis = -0.125*dt_n_h*dt_n_h*(a_n*a_n);
        else {
          dynOps.M->mult(a_n,tmp1);
          Wdis = -0.125*dt_n_h*dt_n_h*(a_n*tmp1);
        }

        energies << setprecision(16) << n+1 << " " << t_n+dt_n_h << " "
                 << Wkin << " " << Wext << " " << Wint << " " << Wdis << " "
                 << Wkin+Wint-Wext+Wdis << " " << std::abs(Wkin+Wint-Wext) << " "
                 << eps1*std::max(Wext,std::max(Wint,Wkin)) << " " 
                 << dt_n_h << std::endl;
#endif
        // check |Wkin + Wint - Wext| <= 0.01*max(Wext,Wint,Wkin)
        if(std::abs(Wkin+Wint-Wext) <= std::max(eps1*std::max(Wext,std::max(Wint,Wkin)),eps2)) {
          std::cerr << "Warning: energy balance check failed at t = " << t_n << std::endl;
          //std::cerr << "exiting...\n";
          //exit(-1);
        }
      }

      // Update the time index
      n += 1;

      // ... current state is replaced by predicted value
      // NOTE: current state is modified here only for output
      // it will be restored to the backup state at the
      // beginning of the next iteration
      if(!parity && probDesc->getAeroAlg() == 5) {
        curState.getDisp().linC(0.5, curState.getDisp(), 0.5, bkState->getDisp());
        curState.getVeloc().linC(0.5, curState.getVeloc(), 0.5, bkState->getVeloc());
      }
      else {
        if(t_n+1.01*dt_n_h > tmax) probDesc->processLastOutput(); // force printing at last iteration
      }

      // Output the state at t^n: d^n, v^n, a^n, fext^n
      postProcessor->dynamOutput(n, t_n+dt_n_h, dynOps, fext, aeroForce, curState);

      // ... For A5 Algorithm, do one time step back if necessary
      // add n so that time index is wound back as well as t
      if(aeroAlg == 5) probDesc->a5TimeLoopCheck(parity, t_n, dt_n_h);

      // Time update:
      t_n += dt_n_h; // t^n = t^{n-1} + deltat^{n-1/2}

      // Choose a new time step deltat^{n+1/2}
      if(domain->solInfo().stable && (aeroAlg < 0 || domain->solInfo().dyna3d_compat) && domain->solInfo().isNonLin()
         && n%(domain->solInfo().subcycle*domain->solInfo().stable_freq) == 0) {
        filePrint(stderr,"\n");
        dt_old = dt_n_h;
        probDesc->computeStabilityTimeStep(dt_n_h, dynOps);
      }

    } 
  }
  if(aeroAlg < 0 && printNumber < std::numeric_limits<int>::max())
    filePrint(stderr, "\r ⌊\x1B[33m   t = %9.3e Δt = %8.2e 100%% \x1B[0m⌋\n", t_n, dt_n_h);

#ifdef PRINT_TIMERS
  if(verboseFlag) filePrint(stderr, " ... Total Loop Time = %.2e s   ...\n", s2/1000.0);
#endif

  if(domain->solInfo().check_energy_balance) {
    delete fext_p;
    delete fint_p;
  }

}

//------------------------------------------------------------------------------

template<
     class DynOps,             // Data Structure for K, C, M and dynMat
     class VecType,            // Vector type used
     class PostProcessor,      
     class ProblemDescriptor,
     class Scalar> 

int
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar >
::checkSteadyState(double time, double step, double criteria)
{
  int nstep = (int)(time / step);
  
  filePrint(stderr," ... Pseudo-Step = %d  Rel. Res. = %10.4e\n",nstep,criteria);

  if(nstep > steadyMax) return 2;

  return (criteria < steadyTol) ? 1 : 0;
}
  
//------------------------------------------------------------------------------

template<
     class DynOps,             // Data Structure for K, C, M and dynMat
     class VecType,            // Vector type used
     class PostProcessor,      
     class ProblemDescriptor,
     class Scalar> 
void
DynamicSolver< DynOps, VecType, PostProcessor, ProblemDescriptor, Scalar >
::getInternalForce(const DynOps &dynamOps, const VecType &disp, VecType &result, double time, int tIndex) {
  if (domain->solInfo().isNonLin()) {
    probDesc->getInternalForce(const_cast<VecType &>(disp), result, time, tIndex);
  } else {
    domain->getTimers().formRhs -= getTime();
    const_cast<DynOps &>(dynamOps).K->mult(const_cast<VecType &>(disp), result);
    domain->getTimers().formRhs += getTime();
  }
}

