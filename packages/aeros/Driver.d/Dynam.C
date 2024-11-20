#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <algorithm>

#ifndef TFLOP
#ifndef WINDOWS 
#include <dlfcn.h>
#endif
#endif

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <Math.d/DBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Driver.d/Dynam.h>
#include <Driver.d/Domain.h>
#include <Hetero.d/FlExchange.h>
#include <Utils.d/ModeData.h>
#include <Element.d/State.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Utils.d/DistHelper.h>
#include <Problems.d/ModalBase.h>

extern ModeData modeData; 
extern ModeData modeDataIDis;
extern ModeData modeDataIVel;
extern ModeData modeDataMode;
extern int verboseFlag;

void
Domain::initDispVeloc(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, const char* ext)
{
 // ... INITIALIZE DISPLACEMENT, VELOCITY, AND ACCELERATION 
 // ... VECTORS TO ZERO
 d_n = v_n = a_n = v_p = 0.0;

 // ... SET INITIAL VELOCITY
 if(numIVelModal) {
   filePrint(stderr, " ... Compute initial velocity from given modal basis ... \n");
   ModeData &modeData = (sinfo.ivel_modal_id == -1) ? ::modeData : modeDataIVel;
   modeData.addMultY(numIVelModal, iVelModal, v_n, c_dsa);
 }
 for(int i = 0; i < numIVel; ++i) {
   int dof = c_dsa->locate(iVel[i].nnum, 1 << iVel[i].dofnum);
   if(dof >= 0)
     v_n[dof] += iVel[i].val;
 }

 // If zeroInitialDisp is set, then return as the displacement is
 // already set to zero above.
 // note: geps is always zero for nonlinear
 if(sinfo.zeroInitialDisp == 0) {
   // ... SET INITIAL DISPLACEMENT FROM IDISP IF IDISP6 DOES NOT EXIST
   // ... OR IF WE ARE USING GEOMETRIC PRE-STRESS (GEPS)
   if(domain->numInitDisp6() == 0 || sinfo.gepsFlg == 1) { // note: always use global num to do this check
     if(numIDisModal) {
       filePrint(stderr, " ... Compute initial displacement from given modal basis ... \n");
       ModeData &modeData = (sinfo.idis_modal_id == -1) ? ::modeData : modeDataIDis;
       modeData.addMultY(numIDisModal, iDisModal, d_n, c_dsa);
     }
     for(int i = 0; i < numIDis; ++i) {
       int dof = c_dsa->locate(iDis[i].nnum, 1 << iDis[i].dofnum);
       if(dof >= 0) d_n[dof] += iDis[i].val;
     }
   }

   // ... SET INITIAL DISPLACEMENT FROM IDISP6
   // ... ONLY IF WE ARE NOT USING GEOMETRIC PRE-STRESS
   if(sinfo.gepsFlg == 0) {
     for(int i = 0; i < numIDis6; ++i) {
       int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
       if(dof >= 0)
         d_n[dof] = iDis6[i].val;
     }   
   }
 }

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(probType() != SolverInfo::NonLinDynam) {
   if(cinfo->lastRestartFile) {
     int fn;
     if(strlen(ext) != 0) {
       char *lastRestartFile = new char[strlen(cinfo->lastRestartFile)+strlen(ext)+1];
       strcpy(lastRestartFile, cinfo->lastRestartFile);
       strcat(lastRestartFile, ext);
       fn = open(lastRestartFile, O_RDONLY);
       delete [] lastRestartFile;
     } 
     else {
       fn = open(cinfo->lastRestartFile, O_RDONLY);
     }
     if(fn >= 0) {
       int vsize, restartTIndex;
       double restartT;
       int readSize = read(fn, &vsize, sizeof(int));
       if(vsize != d_n.size() || readSize != (int) sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       readSize = read(fn, &restartTIndex, sizeof(int));
       if(readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0) 
         sinfo.initialTimeIndex = 0;
       else
         sinfo.initialTimeIndex = restartTIndex;
       readSize = read(fn, &restartT, sizeof(double));
       if(readSize != sizeof(double))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTime = 0.0;
       else
         sinfo.initialTime = restartT;
       readSize = read(fn, d_n.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
     
       readSize = read(fn, v_n.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, v_p.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double))) {
         fprintf(stderr," *** WARNING: Older version of restart file"
                        " -- Missing velocity field is set to zero\n");
         v_p.zero();
       }
       readSize = read(fn, &(sinfo.initExtForceNorm), sizeof(double));
       if(readSize != sizeof(double))
         fprintf(stderr," *** ERROR: Inconsistent restart file: External Force Norm\n");

       close(fn);
     } else {
       perror(" *** ERROR: Restart file could not be opened: ");
     }
   
   } 
 }
}

void
Domain::writeRestartFile(double time, int timeIndex, Vector &d_n, 
                         Vector &v_n, Vector &v_p, double Fref, const char *ext)
{
// either test for pointer or frequency > 0
 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(timeIndex % sinfo.nRestart == 0 || time >= sinfo.tmax-0.1*sinfo.getTimeStep() || domain->solInfo().stop_AeroS) {
   int fn;
   if(strlen(ext) != 0) {
     char *currentRestartFile = new char[strlen(cinfo->currentRestartFile)+strlen(ext)+1];
     strcpy(currentRestartFile, cinfo->currentRestartFile);
     strcat(currentRestartFile, ext);
     fn = open(currentRestartFile, O_WRONLY | O_CREAT, 0666);
     delete [] currentRestartFile;
   } else
   fn = open(cinfo->currentRestartFile, O_WRONLY | O_CREAT, 0666);
   if(fn >= 0) {
     int vsize = d_n.size();
     int writeSize = write(fn, &vsize, sizeof(int));
     if(writeSize != sizeof(int))
       fprintf(stderr," *** ERROR: Writing restart file vector size\n");

     writeSize = write(fn, &timeIndex, sizeof(int));
     if(writeSize != sizeof(int))
       fprintf(stderr," *** ERROR: Writing restart file time index\n");

     writeSize = write(fn, &time, sizeof(double));
     if(writeSize != sizeof(double))
       fprintf(stderr," *** ERROR: Writing restart file time\n");

     writeSize = write(fn, d_n.data(), d_n.size()*sizeof(double));
     if(int(writeSize) != int(d_n.size()*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file displacement \n");

     writeSize = write(fn, v_n.data(), v_n.size()*sizeof(double));
     if(int(writeSize) != int(v_n.size()*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file velocity\n");

     writeSize = write(fn, v_p.data(), v_p.size()*sizeof(double));
     if(int(writeSize) != int(v_p.size()*sizeof(double)))
       fprintf(stderr," *** ERROR: Writing restart file prev. velocity\n");

     writeSize = write(fn, &Fref, sizeof(double));
     if(writeSize != sizeof(double))
       fprintf(stderr," *** ERROR: Writing restart file external force norm\n");
     close(fn);
   } else {
      perror(" *** ERROR: Restart file could not be opened: ");
      exit(-1);
   }
 }
}
//---------------------------------------------------------------------------------------------

void
Domain::getOrAddDofForPrint(bool ad, Vector& d_n, double* bcx, int iNode, double *xdata, 
                        int *dofx, double *ydata, int *dofy, double *zdata, int *dofz)
{
    // ad==true -> add :   d_n[xloc]+=data  
    // ad==false -> get:  data=d_n[xloc]
    
    //beginning of dynamOutput with: xdata=vx=variation, dn=dn_aero, bcx=bcx_aero 

    // c_dsa mean: dofsetarray without constrained dof ie free dof
    //   dsa mean: constrained dofsetarray

    // xloc>=0  : dof free
    // xloc1>=0 : dof exist but not free
    // xloc<0 && xloc1<0 : dof doesn't exist
  
    if (dofx){
      int xloc  = c_dsa->locate( iNode, *dofx);
      int xloc1 =   dsa->locate( iNode, *dofx);
      if(xloc >= 0)       { if(!ad) *xdata=d_n[xloc];  else d_n[xloc]+=*xdata; }
      else if(xloc1 >= 0 && bcx) { if(!ad) *xdata=bcx[xloc1]; else bcx[xloc1]+=*xdata;}
      else                { if(!ad) *xdata=0.0;}
    }
    
    if(dofy){
     int yloc  = c_dsa->locate( iNode, *dofy);
     int yloc1 =   dsa->locate( iNode, *dofy);
     if(yloc >= 0)       { if(!ad) *ydata=d_n[yloc];  else d_n[yloc]+=*ydata; }
     else if(yloc1 >= 0 && bcx) { if(!ad) *ydata=bcx[yloc1]; else bcx[yloc1]+=*ydata;}
     else                { if(!ad) *ydata=0.0;}
    }
   
    if(dofz){
      int zloc  = c_dsa->locate( iNode, *dofz);
      int zloc1 =   dsa->locate( iNode, *dofz);
      if(zloc >= 0)       { if(!ad) *zdata=d_n[zloc];  else d_n[zloc]+=*zdata; }
      else if(zloc1 >= 0 && bcx) { if(!ad) *zdata=bcx[zloc1]; else bcx[zloc1]+=*zdata;}
      else                { if(!ad) *zdata=0.0;}
    }
}
//-----------------------------------------------------------------------------------------------

void
Domain::addVariationOfShape_StructOpt(int iNode, CoordSet *nodescopy, double &x, double &y, double &z)
{
  // ... considering the variation of shape
  Node oldnod;      
  oldnod=nodescopy->getNode(iNode);
        
  double vx = nodes[iNode]->x - oldnod.x ;
  double vy = nodes[iNode]->y - oldnod.y ;
  double vz = nodes[iNode]->z - oldnod.z ;
       
  x=x+vx;   y=y+vy;    z=z+vz;
}
//------------------------------------------------------------------------------------------------

void
Domain::aeroSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState)
{
  // Send u + IDISP6 to fluid code.
  // IDISP6 is used to compute pre-stress effects.
  getTimers().sendFluidTime -= getTime();
  Vector d_n_aero(d_n);

  if(sinfo.gepsFlg == 1) {

    for(int i = 0; i < numIDis6; ++i) {
      int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
      if(dof >= 0)
        d_n_aero[dof] += iDis6[i].val;
    }
  }

  State state(c_dsa, dsa, bcx, vcx, d_n_aero, v_n, a_n, v_p);

  if(sinfo.dyna3d_compat) {
    if(aeroEmbeddedSurfaceId.size() != 0 && sinfo.elementDeletion && !newDeletedElements.empty()) {
      flExchanger->sendNewStructure(newDeletedElements);
    }
    else flExchanger->sendNoStructure();
  }
  flExchanger->sendDisplacements(state, -1, geomState);
  if(verboseFlag) filePrint(stderr, " ... [E] Sent displacements         ...\n");

  getTimers().sendFluidTime += getTime();
}

void
Domain::aeroheatSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState)
{
  State tempState(c_dsa, dsa, bcx, d_n, v_n, v_p);

  flExchanger->sendTemperature(tempState);
  if(verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
}

void
Domain::thermohSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState)
{
  int iNode;
  Vector tempsent(numnodes);

  for(iNode=0; iNode<numnodes; ++iNode) {
    int tloc  = c_dsa->locate( iNode, DofSet::Temp);
    int tloc1 =   dsa->locate( iNode, DofSet::Temp);
    double temp  = (tloc >= 0) ? d_n[tloc] : bcx[tloc1];
    if(tloc1 < 0) temp = 0.0;
    tempsent[iNode] = temp;
  }

  flExchanger->sendStrucTemp(tempsent);
  if(verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
}

void
Domain::buildAeroelasticForce(Vector& aero_f, PrevFrc& prevFrc, int tIndex, double t, double gamma, double alphaf, GeomState* geomState)
{
  // ... COMPUTE AEROELASTIC FORCE 
  getTimers().receiveFluidTime -= getTime();

  // ... Temporary variable for inter(extra)polated force
  double *tmpFmem = new double[numUncon()];
  StackVector tmpF(tmpFmem, numUncon());
  tmpF.zero();

  int iscollocated;
  double tFluid = flExchanger->getFluidLoad(tmpF, tIndex, t,
                                            alphaf, iscollocated, geomState);
  if(verboseFlag) filePrint(stderr, " ... [E] Received fluid load norm is %e ...\n", tmpF.norm());

  if(sinfo.aeroFlag == 20) {
    if(prevFrc.lastTIndex >= 0)
      aero_f.linC(0.5,tmpF,0.5,prevFrc.lastFluidLoad);
    else
      aero_f = tmpF;
  }
  else {
    if(iscollocated == 0) {
      if(prevFrc.lastTIndex >= 0) {
        tmpF *= (1/gamma);
        tmpF.linAdd(((gamma-1.0)/gamma), prevFrc.lastFluidLoad);
      }
    }

    double alpha = (prevFrc.lastTIndex < 0) ? 1.0 : 1.0-alphaf;
    aero_f.linC(alpha, tmpF, (1.0-alpha), prevFrc.lastFluidLoad);
  }
  prevFrc.lastFluidLoad = tmpF;
  prevFrc.lastFluidTime = tFluid;
  prevFrc.lastTIndex = tIndex;

  delete [] tmpFmem;

  if(sinfo.aeroFlag == 20 && sinfo.dyna3d_compat) {
    if(sinfo.stop_AeroF) sinfo.stop_AeroS = true;
    double dt = sinfo.getTimeStep()*sinfo.subcycle;
    if(tIndex == sinfo.initialTimeIndex) sinfo.t_AeroF = sinfo.initialTime + 1.5*dt;
    else sinfo.t_AeroF += dt;
    double maxTime_AeroF = sinfo.tmax-0.5*dt;
    double sendtim;
    if((sinfo.t_AeroF < (maxTime_AeroF-0.01*dt)) || tIndex == sinfo.initialTimeIndex)
      sendtim = 1e100;
    else {
      sendtim = 0;
      sinfo.stop_AeroF = true;
    }
    int restartinc = std::max(sinfo.nRestart, 0);
    flExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep()*sinfo.subcycle, sendtim, restartinc,
                           sinfo.isCollocated, sinfo.alphas,  sinfo.alphasv);
    if(tIndex == 0) // Send the parameter a second time for fluid iteration 1 to 2
      flExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep()*sinfo.subcycle, sendtim, restartinc,
                             sinfo.isCollocated, sinfo.alphas, sinfo.alphasv);
  }

  getTimers().receiveFluidTime += getTime();
}

void
Domain::buildAeroheatFlux(Vector &f, Vector &prev_f, int tIndex, double t)
{
  // ... ADD FLUID FLUX
  getTimers().receiveFluidTime -= getTime();

  // ... Temporary variable for inter(extra)polated force
  double *tmpFmem = (double *) dbg_alloca(sizeof(double)*numUncon());
  StackVector tmpF(tmpFmem, numUncon());

  double sflux = 0;
  double bfl ;

  tmpF.zero();
  flExchanger->getFluidFlux(tmpF, t, bfl);

  if(verboseFlag) filePrint(stderr, " ... [T] Received fluid fluxes      ...\n");

  int vectlen = tmpF.size();

/*  Compute fluid flux at n+1/2, since we use midpoint rule in thermal */

  int useProjector = domain->solInfo().filterFlags;

  if(tIndex == 0)
    f += tmpF;
  else {
    if(useProjector) f = tmpF;
    else
      f.linAdd(0.5, tmpF, 0.5, prev_f);
  }

  prev_f = tmpF;

/*  Print out sum of fluxes received at n+1 and at n+1/2 in a file
    specified by raerotfl in entry */

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  for(int i=0; i < numOutInfo; i++) {
    if(oinfo[i].interval != 0 && tIndex % oinfo[i].interval == 0) {
      switch(oinfo[i].type) {
        case OutputInfo::AeroFlux:
          fprintf(oinfo[i].filptr,"%e   ",t);
          fprintf(oinfo[i].filptr,"%e   ",bfl);
          for(int j=0; j<vectlen; j++)  {
            sflux += 0.5*(tmpF[j]+prev_f[j]);
          }
          fprintf(oinfo[i].filptr,"%e\n",sflux);
          fflush(oinfo[i].filptr);
        break;
        default: /* do nothing */
        break;
      }
    }
  }
}

void
Domain::thermoeComm()
{
  flExchanger->getStrucTemp(temprcvd);
  if(verboseFlag) filePrint(stderr, " ... [E] Received temperatures      ...\n");

  computeTDProps();
}

void
Domain::dynamOutput(int tIndex, double time, double *bcx, DynamMat& dMat, Vector& ext_f, Vector &aeroForce, 
                    Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* vcx, double* acx)
{
  if(sinfo.nRestart > 0 && !sinfo.modal && probType() != SolverInfo::NonLinDynam) {
    writeRestartFile(time, tIndex, d_n, v_n, v_p, sinfo.initExtForceNorm);
  }
  
  // Send to fluid code (except explicit C0)
  if(sinfo.aeroFlag >= 0 && !sinfo.lastIt && tIndex != sinfo.initialTimeIndex && !(sinfo.newmarkBeta == 0 && sinfo.aeroFlag == 20))
    aeroSend(d_n, v_n, a_n, v_p, bcx, vcx);

  if(sinfo.aeroheatFlag >= 0 && tIndex != 0 )
    aeroheatSend(d_n, v_n, a_n, v_p, bcx, vcx);

  if(sinfo.thermohFlag >=0 && tIndex != sinfo.initialTimeIndex) 
    thermohSend(d_n, v_n, a_n, v_p, bcx, vcx);

  int numOutInfo = geoSource->getNumOutInfo();

  // Open the file and write the header in the first time step
  if (tIndex == sinfo.initialTimeIndex) {
    if(numOutInfo > 0)
      geoSource->openOutputFiles();
  }

  // Check if there are any output files which need to be processed at this time
  if(geoSource->noOutput(tIndex)) return;

  this->dynamOutputImpl(tIndex, bcx, dMat, ext_f, aeroForce, d_n, v_n, a_n, v_p, vcx, acx, time, 0, numOutInfo);
} 

// Output for time-parallel linear dynamics
// There is a full set of files for each time-slice on the CPU
// Preconditions:
// 1) The function GeoSource::duplicateFilesForPita must have been called during the initialization
// 2) The corresponding output files must already be open
void
Domain::pitaDynamOutput(int tIndex, double *bcx, DynamMat& dMat, Vector& ext_f, Vector &aeroForce, 
                        Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* vcx, double* acx,
                        int sliceRank, double time)
{
  std::pair<int, int> requestIndices = geoSource->getTimeSliceOutputFileIndices(sliceRank);
  this->dynamOutputImpl(tIndex, bcx, dMat, ext_f, aeroForce, d_n, v_n, a_n, v_p, vcx, acx, time, requestIndices.first, requestIndices.second);
} 

// Perform output for files with indices in [firstRequest, lastRequest)
// Precondition: The corresponding output files must already be open
void
Domain::dynamOutputImpl(int tIndex, double *bcx, DynamMat& dMat, Vector& ext_f, Vector &aeroForce,
                        Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* vcx, double* acx,
                        double time, int firstRequest, int lastRequest)
{
  // Print out the displacement info
  if(outFlag && !nodeTable) makeNodeTable(outFlag);
  int numNodes = geoSource->numNode();  // don't print displacements for internal nodes
  int numNodeLim = std::max(numNodes,numnodes);
  double (*glDisp)[11] = new double[numNodeLim][11];
  double (*locDisp)[11] = (domain->solInfo().basicDofCoords) ? 0 : new double[numNodeLim][11];
  for (int i = 0; i < numNodeLim; ++i)
    for (int j = 0 ; j < 11 ; j++)
      glDisp[i][j] = 0.0;
  mergeDistributedDisp(glDisp, d_n.data(), bcx, locDisp);
  int numNodesOut = (outFlag) ? exactNumNodes : numNodes;

  for (int i = firstRequest; i < lastRequest; ++i) {
    enum {YOUNG,MDENS,THICK};
    int iNode, nodeI, realNode;
    int first_node, last_node, last_node_out;
    
    OutputInfo *oinfo = geoSource->getOutputInfo();
    if(sinfo.isNonLin() && (oinfo[i].isStressOrStrain() || oinfo[i].isRotation())) continue; // see Domain::postProcessing in NLStatic.C
    if(oinfo[i].sentype) continue;
    
    //CD: ad and get will be used in  addVariationOfShape_StructOpt and getOrAddDofForPrint which were 
    //    added to "clean" dynamOutput 
    bool get=false;
     
    if (oinfo[i].nodeNumber == -1) { first_node=0; last_node=numNodes; last_node_out=numNodesOut; }           
    else { first_node=oinfo[i].nodeNumber; last_node=last_node_out=first_node+1; }

    if ((oinfo[i].interval != 0) && (tIndex % oinfo[i].interval == 0)) {
      int dof=-1;
      int w = oinfo[i].width;
      int p = oinfo[i].precision;

      int success;
      if(oinfo[i].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
        success = processDispTypeOutputs(oinfo[i], glDisp, numNodesOut, i, time);
      else
        success = processDispTypeOutputs(oinfo[i], locDisp, numNodesOut, i, time);
      if (success) continue;
      success = processOutput(oinfo[i].type, d_n, bcx, i, time);
      if (success) continue;
      success = 1;

      int nNodes = last_node-first_node;
      int nNodesOut = last_node_out-first_node;
      switch(oinfo[i].type) {

        case OutputInfo::Velocity6: {
          double (*data)[6] = new double[nNodesOut][6];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data[nodeI], &DofSet::Xdisp,
                                data[nodeI]+1, &DofSet::Ydisp, data[nodeI]+2, &DofSet::Zdisp);

            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data[nodeI]+3, &DofSet::Xrot,
                                data[nodeI]+4, &DofSet::Yrot, data[nodeI]+5, &DofSet::Zrot);

            // transform velocity from DOF_FRM to basic coordinate frame
            if(oinfo[i].oframe == OutputInfo::Global) transformVectorInv(&(data[nodeI][0]), first_node+iNode, true);
          }
          geoSource->outputNodeVectors6(i, data, nNodesOut, time);
          delete [] data;
        }
          break;
        case OutputInfo::Velocity:  {
          double (*data)[3] = new double[nNodesOut][3];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data[nodeI], &DofSet::Xdisp,
                                data[nodeI]+1, &DofSet::Ydisp, data[nodeI]+2, &DofSet::Zdisp);

            // transform velocity from DOF_FRM to basic coordinate frame
            if(oinfo[i].oframe == OutputInfo::Global) transformVectorInv(&(data[nodeI][0]), first_node+iNode, false);
          }
          geoSource->outputNodeVectors(i, data, nNodesOut, time);
          delete [] data;
        }
          break;
        case OutputInfo::PressureFirstTimeDerivative: {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode) {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data+nodeI, &DofSet::Helm, 0, 0, 0, 0);
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;

        } 
          break;
        case OutputInfo::TemperatureFirstTimeDerivative: {
          double *data = new double[nNodesOut]; 
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode) {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            getOrAddDofForPrint(get, v_n, vcx, first_node+iNode, data+nodeI, &DofSet::Temp, 0,0,0,0);
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        }
          break;
        case OutputInfo::Accel6:  {
          double (*data)[6] = new double[nNodesOut][6];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            getOrAddDofForPrint(get, a_n, acx, first_node+iNode, data[nodeI],
                                &DofSet::Xdisp, data[nodeI]+1, &DofSet::Ydisp,
                                data[nodeI]+2, &DofSet::Zdisp);
            getOrAddDofForPrint(get, a_n, acx, first_node+iNode, data[nodeI]+3, 
                                &DofSet::Xrot, data[nodeI]+4, &DofSet::Yrot, data[nodeI]+5, 
                                &DofSet::Zrot);

             // transform acceleration from DOF_FRM to basic coordinate frame
            if(oinfo[i].oframe == OutputInfo::Global) transformVectorInv(&(data[nodeI][0]), first_node+iNode, true);
          }
          geoSource->outputNodeVectors6(i, data, nNodesOut, time);
          delete [] data;
        }
          break;
        case OutputInfo::Acceleration:  {
          double (*data)[3] = new double[nNodesOut][3];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            getOrAddDofForPrint(get, a_n, acx, first_node+iNode, data[nodeI], 
                                &DofSet::Xdisp, data[nodeI]+1, &DofSet::Ydisp, data[nodeI]+2, &DofSet::Zdisp);

            // transform acceleration from DOF_FRM to basic coordinate frame
            if(oinfo[i].oframe == OutputInfo::Global) transformVectorInv(&(data[nodeI][0]), first_node+iNode, false);
          }
          geoSource->outputNodeVectors(i, data, nNodesOut, time);
          delete [] data;
        }
          break;
        case OutputInfo::PressureSecondTimeDerivative: {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode) {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            getOrAddDofForPrint(get, a_n, (double *) 0, first_node+iNode, data+nodeI, &DofSet::Helm, 
                                0, 0, 0, 0);
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        }
          break;
        case OutputInfo::Energies: {
          double Wela, Wkin, error; 
          computeEnergies(d_n, ext_f, time, &aeroForce, &v_n, dMat.K, dMat.M, dMat.C, Wela, Wkin, error);
          geoSource->outputEnergies(i, time, Wext, Waero, Wela, Wkin, Wdmp, error);
        }
          break;
        case OutputInfo::AeroForce: break; // this is done in FlExchange.C
        case OutputInfo::AeroXForce:  {
          if(sinfo.aeroFlag >= 0) {
            double *data = new double[nNodesOut];
            realNode = -1;
            for (iNode = 0; iNode < nNodes; ++iNode)  {
              if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
              int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
              data[nodeI]  = (xloc >= 0) ? aeroForce[xloc] : 0.0;
            }
            geoSource->outputNodeScalars(i, data, nNodesOut, time);
            delete [] data;
          }
        } break;
        case OutputInfo::AeroYForce:  {
          if(sinfo.aeroFlag >= 0) {
            double *data = new double[nNodesOut];
            realNode = -1;
            for (iNode = 0; iNode < nNodes; ++iNode)  {
              if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
              int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
              data[nodeI]  = (yloc >= 0) ? aeroForce[yloc] : 0.0;
            }
            geoSource->outputNodeScalars(i, data, nNodesOut, time);
            delete [] data;
          }
        } break;
        case OutputInfo::AeroZForce:  {
          if(sinfo.aeroFlag >= 0) {
            double *data = new double[nNodesOut];
            realNode = -1;
            for (iNode = 0; iNode < nNodes; ++iNode)  {
              if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
              int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
              data[nodeI] = (zloc >= 0) ? aeroForce[zloc] : 0.0;
            }
            geoSource->outputNodeScalars(i, data, nNodesOut, time);
            delete [] data;
          }
        } break;
        case OutputInfo::AeroXMom:  {
          if(sinfo.aeroFlag >= 0) {
            double *data = new double[nNodesOut];
            realNode = -1;
            for (iNode = 0; iNode < nNodes; ++iNode)  {
              if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
              int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
              data[nodeI] = (xrot >= 0) ? aeroForce[xrot] : 0.0;
            }
            geoSource->outputNodeScalars(i, data, nNodesOut, time);
            delete [] data;
          }
        } break;
        case OutputInfo::AeroYMom:  {
          if(sinfo.aeroFlag >= 0) {
            double *data = new double[nNodesOut];
            realNode = -1;
            for (iNode = 0; iNode < nNodes; ++iNode)  {
              if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
              int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
              data[nodeI] = (yrot >= 0) ? aeroForce[yrot] : 0.0;
            }
            geoSource->outputNodeScalars(i, data, nNodesOut, time);
            delete [] data;
          }
        } break;
        case OutputInfo::AeroZMom:  {
          if(sinfo.aeroFlag >= 0) {
            double *data = new double[nNodesOut];
            realNode = -1;
            for (iNode = 0; iNode < nNodes; ++iNode)  {
              if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
              int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
              data[nodeI] = (zrot >= 0) ? aeroForce[zrot] : 0.0;
            }
            geoSource->outputNodeScalars(i, data, nNodesOut, time);
            delete [] data;
          }
        } break;
        case OutputInfo::ExternalXForce:  {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            int xloc  = c_dsa->locate(first_node+iNode, DofSet::Xdisp);
            data[nodeI]  = (xloc >= 0) ? ext_f[xloc] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        } break;
        case OutputInfo::ExternalYForce:  {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            int yloc  = c_dsa->locate(first_node+iNode, DofSet::Ydisp);
            data[nodeI]  = (yloc >= 0) ? ext_f[yloc] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        } break;
        case OutputInfo::ExternalZForce:  {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            int zloc  = c_dsa->locate(first_node+iNode, DofSet::Zdisp);
            data[nodeI] = (zloc >= 0) ? ext_f[zloc] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        } break;
        case OutputInfo::ExternalXMom:  {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            int xrot  = c_dsa->locate(first_node+iNode, DofSet::Xrot);
            data[nodeI] = (xrot >= 0) ? ext_f[xrot] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        } break;
        case OutputInfo::ExternalYMom:  {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            int yrot  = c_dsa->locate(first_node+iNode, DofSet::Yrot);
            data[nodeI] = (yrot >= 0) ? ext_f[yrot] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        } break;
        case OutputInfo::ExternalZMom:  {
          double *data = new double[nNodesOut];
          realNode = -1;
          for (iNode = 0; iNode < nNodes; ++iNode)  {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            int zrot  = c_dsa->locate(first_node+iNode, DofSet::Zrot);
            data[nodeI] = (zrot >= 0) ? ext_f[zrot] : 0.0;
          }
          geoSource->outputNodeScalars(i, data, nNodesOut, time);
          delete [] data;
        } break;
        case OutputInfo::Composit:
          getCompositeData(i,time);
          break;          
        case OutputInfo::TDEnforcement: {
          if(tdenforceFlag()) {
            // TODO outFlag == 1
            double *plot_data = new double[numNodes];
            if(oinfo[i].tdenforc_var == 1) // CONFACE
              for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.5;
            else
              for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.0;
            for(int iMortar=0; iMortar<nMortarCond; iMortar++) {
              MortarConds[iMortar]->get_plot_variable(oinfo[i].tdenforc_var,plot_data);
            }
            if(oinfo[i].nodeNumber == -1)
              geoSource->outputNodeScalars(i, plot_data, numNodes, time);
            else
              geoSource->outputNodeScalars(i, &plot_data[oinfo[i].nodeNumber], 1, time);
            delete [] plot_data;
          }
          else fprintf(stderr, " *** WARNING: output %d is not supported \n", i);
        } break;
        case OutputInfo::Reactions: {
          Vector fc(numDirichlet);
          computeReactionForce(fc, d_n, v_n, a_n, bcx, vcx, acx,
                               dMat.Kuc, dMat.Kcc, dMat.Cuc, dMat.Ccc, dMat.Muc, dMat.Mcc);
          double (*rxyz)[3] = new double[nNodesOut][3];
          DofSet dofs[3] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp };
          realNode = -1;
          for(int iNode = 0; iNode < nNodes; ++iNode) {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            for(int k = 0; k < 3; ++k) {
              int dof =   dsa->locate(first_node+iNode, dofs[k].list());
              int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
              rxyz[nodeI][k] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
            }
            if(oinfo[i].oframe == OutputInfo::Global && !domain->solInfo().basicDofCoords) {
              transformVectorInv(&rxyz[nodeI][0],iNode,false);
            }
          }
          geoSource->outputNodeVectors(i, rxyz, nNodesOut, time);
          delete [] rxyz;
          } break;
        case OutputInfo::Reactions6: {
          Vector fc(numDirichlet);
          computeReactionForce(fc, d_n, v_n, a_n, bcx, vcx, acx,
                               dMat.Kuc, dMat.Kcc, dMat.Cuc, dMat.Ccc, dMat.Muc, dMat.Mcc);
          double (*rxyz)[6] = new double[nNodesOut][6];
          DofSet dofs[6] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
                             DofSet::Xrot, DofSet::Yrot, DofSet::Zrot };
          realNode = -1;
          for(int iNode = 0; iNode < nNodes; ++iNode) {
            if(outFlag) { if(nodes[first_node+iNode] == 0) continue; nodeI = ++realNode; } else nodeI = iNode;
            for(int k = 0; k < 6; ++k) {
              int dof =   dsa->locate(first_node+iNode, dofs[k].list());
              int cdof = (dof >= 0) ? c_dsa->invRCN(dof) : -1;
              rxyz[nodeI][k] = (cdof >= 0) ? fc[cdof] : 0;     // constrained
            }
            if(oinfo[i].oframe == OutputInfo::Global && !domain->solInfo().basicDofCoords) {
              transformVectorInv(&rxyz[nodeI][0],iNode,true);
            }
          }
          geoSource->outputNodeVectors6(i, rxyz, nNodesOut, time);
          delete [] rxyz;
          } break;
        case OutputInfo::DeletedElements: {
          for(std::vector<std::pair<double,int> >::iterator it = outDeletedElements.begin(); it != outDeletedElements.end(); ++it) {
            filePrint(oinfo[i].filptr, " %12.6e  %9d          Undetermined\n", it->first, it->second+1);
            fflush(oinfo[i].filptr);
          }
          outDeletedElements.clear();
        } break;
        case OutputInfo::ModeError: // don't print warning message since these are
        case OutputInfo::ModeAlpha: // output in SingleDomainDynamic::modeDecomp
          break;
        case OutputInfo::Statevector: // don't print warning message for these either
        case OutputInfo::Velocvector:
        case OutputInfo::Accelvector:
        case OutputInfo::Residual:
        case OutputInfo::Jacobian:
        case OutputInfo::RobData:
        case OutputInfo::SampleMesh:
        case OutputInfo::ModalDsp:
        case OutputInfo::ModalExF:
        case OutputInfo::ModalMass:
        case OutputInfo::ModalStiffness:
        case OutputInfo::ModalDamping:
        case OutputInfo::ModalDynamicMatrix:
        case OutputInfo::ModalMatrices:
        case OutputInfo::ErrorIndicator:
          break;

        default:
          success = 0;
          break;
        
      }
      if (success == 0)
        fprintf(stderr, " *** WARNING: output %d is not supported \n", i);
    }
  }

  if (glDisp) delete [] glDisp;
  if (locDisp) delete [] locDisp;
  firstOutput = false;
}

//----------------------------------------------------------------------------------------------

ControlInterface* Domain::getUserSuppliedFunction() {

  ControlInterface *userSupFunc = 0;

#ifndef TFLOP
#ifndef WINDOWS
  if( claw ) {
    void *handle;


    dlerror(); // forget about the last error
    handle = dlopen(claw->fileName, RTLD_NOW);
    const char *errorMsg;
    if((errorMsg = dlerror() ) != 0) {
      fprintf(stderr," *** ERROR: in dynamic loading of %s: %s\n",
             claw->fileName,errorMsg);
      exit(-1);
    }

    ControlInterface ** fcp =
        (ControlInterface **) dlsym(handle, claw->routineName);

    if(fcp ==0) {
      fprintf(stderr," *** ERROR: in dynamic loading of %s: "
                     "control function not found\n",
                     claw->routineName);
      exit(-1);
    }

    userSupFunc = *fcp;
 }
#endif
#else
// Should have something here
#endif

 return userSupFunc;

}

void
Domain::aeroPreProcess(Vector& d_n, Vector& v_n, Vector& a_n,
                       Vector& v_p, double *bcx, double *vcx)
{
  if(sinfo.aeroFlag >= 0) {

    int numOutInfo = geoSource->getNumOutInfo();
    OutputInfo *oinfo = geoSource->getOutputInfo();

    int flag = 0;

    // Check if aero forces are requested for output
    int iInfo;
    for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
      if(oinfo[iInfo].type == OutputInfo::AeroForce) { 
        flag = 1;
        break;
      }
    }

    if(sinfo.aeroFlag == 20 && sinfo.newmarkBeta != 0) {
      filePrint(stderr, " *** ERROR: Requested AERO Algorithm is not available with implicit time-integrator. Aborting...\n");
      exit(-1);
    }
    if(!(sinfo.aeroFlag == 20 || sinfo.aeroFlag == 1) && sinfo.newmarkBeta == 0) {
      filePrint(stderr, " *** ERROR: Requested AERO Algorithm is not available with explicit time-integrator. Aborting...\n");
      exit(-1);
    }

    OutputInfo *oinfo_aero = (flag) ? oinfo+iInfo : NULL;
    if(aeroEmbeddedSurfaceId.size() != 0) {
      int iSurf = -1;
      for(int i=0; i<nSurfEntity; i++)
        if(aeroEmbeddedSurfaceId.find(SurfEntities[i]->ID()) != aeroEmbeddedSurfaceId.end()) {
          iSurf = i; 
          break; // only allows one surface.
        }
      if(iSurf<0) {
        filePrint(stderr, " *** ERROR: Embedded wet surface not found! Aborting...\n");
        exit(-1);
      }
      flExchanger = new FlExchanger(nodes, packedEset, SurfEntities[iSurf], c_dsa, oinfo_aero, sinfo.elementDeletion);
    }
    else {
      if(sinfo.elementDeletion) {
        filePrint(stderr," *** WARNING: The C0 algorithm and an embedded surface id must be specified\n"
                         "     under AERO for an aeroelastic analysis with element deletion, otherwise\n"
                         "     Aero-F will not be notified of any topological changes in structure.\n");
      }
      flExchanger = new FlExchanger(nodes, packedEset, c_dsa, oinfo_aero);
    }

    const char *matchFile = geoSource->getMatchFileName();
    if(matchFile == 0)
      matchFile = (char*) "MATCHER";

    flExchanger->read(0, matchFile);

    //KW: send the embedded wet surface to fluid 
    if(aeroEmbeddedSurfaceId.size() != 0) {
      flExchanger->sendEmbeddedWetSurface();
      if(verboseFlag) filePrint(stderr, " ... [E] Sent embedded wet surface  ...\n");
    }

    //XML New step of negotiation with fluid code
    flExchanger->negotiate();

    int restartinc = (solInfo().nRestart >= 0) ? (solInfo().nRestart) : 0;

    Vector d_n_aero(d_n);

    // Send u + IDISP6 to fluid code.
    // IDISP6 is used to compute pre-stress effects.
    if(sinfo.gepsFlg == 1 && sinfo.aeroFlag != 8) {
      // If we are in the first time step, and we initialized with
      // IDISP6, do not send IDISP6
      if(numIDis == 0 && sinfo.zeroInitialDisp != 1) {
        fprintf(stderr," ... DO NOT SEND IDISP6             ...\n");
      } else {
        fprintf(stderr," ... SENDING IDISP6                 ...\n");
        int i;
        for(i = 0; i < numIDis6; ++i) {
          int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
          if(dof >= 0)
            d_n_aero[dof] += iDis6[i].val;
        }
      }
    }
    
    State curState(c_dsa, dsa, bcx, vcx, d_n_aero, v_n, a_n, v_p);

    if(sinfo.aeroFlag == 8) {
      flExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.mppFactor,
                             restartinc, sinfo.isCollocated, sinfo.alphas, sinfo.alphasv);
      flExchanger->sendModeFreq(modeData.frequencies, modeData.numModes);
      if(verboseFlag) filePrint(stderr, " ... [E] Sent parameters and mode frequencies ...\n");
      if(sinfo.gepsFlg == 1) {
        flExchanger->sendModeShapes(modeData.numModes, modeData.numNodes,
                     modeData.modes, curState, sinfo.mppFactor, numIDis6, iDis6);
      } else flExchanger->sendModeShapes(modeData.numModes, modeData.numNodes, 
                       modeData.modes, curState, sinfo.mppFactor);
      if(verboseFlag) filePrint(stderr, " ... [E] Sent mode shapes           ...\n");
    }
    else {
      double aero_tmax = sinfo.tmax;
      if(sinfo.newmarkBeta == 0 && !sinfo.dyna3d_compat) aero_tmax += sinfo.getTimeStep();
      double aero_dt = (sinfo.dyna3d_compat) ? 0 : sinfo.getTimeStep();
      flExchanger->sendParam(sinfo.aeroFlag, aero_dt, aero_tmax, restartinc,
                             sinfo.isCollocated, sinfo.alphas, sinfo.alphasv);
      if(verboseFlag) filePrint(stderr, " ... [E] Sent parameters            ...\n");

      if(sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
        flExchanger->initRcvParity(1);
        flExchanger->initSndParity(1);
      } else {
        flExchanger->initRcvParity(-1);
        flExchanger->initSndParity(-1);
      }

      if(sinfo.aeroFlag == 20 && sinfo.dyna3d_compat) {
        flExchanger->sendSubcyclingInfo(0);
        flExchanger->sendNoStructure();
      }

      if(!(geoSource->getCheckFileInfo()->hotRestart() && (sinfo.aeroFlag == 20 ||
          (sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0)))) {
        flExchanger->sendDisplacements(curState);
        if(verboseFlag) filePrint(stderr, " ... [E] Sent initial displacements ...\n");
      }

      if(sinfo.aeroFlag == 1) { // Ping pong only
        fprintf(stderr, "Ping Pong Only requested. Structure code exiting\n");
      }
    }
  }
}

void
Domain::sendDisplacements(Vector& d_n, Vector& v_n, Vector& a_n,
                          Vector& v_p, double *bcx, double *vcx)
{
  if(sinfo.aeroFlag >= 0) {
    Vector d_n_aero(d_n);
    State curState(c_dsa, dsa, bcx, vcx, d_n_aero, v_n, a_n, v_p);
    flExchanger->sendDisplacements(curState);
  }
}

void
Domain::aeroSensitivityPreProcess(Vector& d_n, Vector& v_n, Vector& a_n,
                                  Vector& v_p, double *bcx, double *vcx)
{
  if(sinfo.aeroFlag >= 0) {

    int numOutInfo = geoSource->getNumOutInfo();
    OutputInfo *oinfo = geoSource->getOutputInfo();

    int flag = 0;

    // Check if aero forces are requested for output
    int iInfo;
    for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
      if(oinfo[iInfo].type == OutputInfo::AeroForce) { 
        flag = 1;
        break;
      }
    }

    if(sinfo.aeroFlag == 20 && sinfo.newmarkBeta != 0) {
      filePrint(stderr, " *** ERROR: Requested AERO Algorithm is not available with implicit time-integrator. Aborting...\n");
      exit(-1);
    }
    if(!(sinfo.aeroFlag == 20 || sinfo.aeroFlag == 1) && sinfo.newmarkBeta == 0) {
      filePrint(stderr, " *** ERROR: Requested AERO Algorithm is not available with explicit time-integrator. Aborting...\n");
      exit(-1);
    }

    OutputInfo *oinfo_aero = (flag) ? oinfo+iInfo : NULL;
    if(aeroEmbeddedSurfaceId.size() != 0) {
      int iSurf = -1;
      for(int i=0; i<nSurfEntity; i++)
        if(aeroEmbeddedSurfaceId.find(SurfEntities[i]->ID()) != aeroEmbeddedSurfaceId.end()) {
          iSurf = i; 
          break; // only allows one surface.
        }
      if(iSurf<0) {
        filePrint(stderr, " *** ERROR: Embedded wet surface not found! Aborting...\n");
        exit(-1);
      }
      flExchanger = new FlExchanger(nodes, packedEset, SurfEntities[iSurf], c_dsa, oinfo_aero, sinfo.elementDeletion);
    }
    else {
      if(sinfo.elementDeletion) {
        filePrint(stderr," *** WARNING: The C0 algorithm and an embedded surface id must be specified\n"
                         "     under AERO for an aeroelastic analysis with element deletion, otherwise\n"
                         "     Aero-F will not be notified of any topological changes in structure.\n");
      }
      flExchanger = new FlExchanger(nodes, packedEset, c_dsa, oinfo_aero);
    }

    const char *matchFile = geoSource->getMatchFileName();
    if(matchFile == 0)
      matchFile = (char*) "MATCHER";

    flExchanger->read(0, matchFile);

    //KW: send the embedded wet surface to fluid 
    if(aeroEmbeddedSurfaceId.size() != 0) {
      flExchanger->sendEmbeddedWetSurface();
      if(verboseFlag) filePrint(stderr, " ... [E] Sent embedded wet surface  ...\n");
    }

    //XML New step of negotiation with fluid code
    flExchanger->negotiate();

    int restartinc = (solInfo().nRestart >= 0) ? (solInfo().nRestart) : 0;

    Vector d_n_aero(d_n);

    // Send u + IDISP6 to fluid code.
    // IDISP6 is used to compute pre-stress effects.
    if(sinfo.gepsFlg == 1) {
      // If we are in the first time step, and we initialized with
      // IDISP6, do not send IDISP6
      if(numIDis == 0 && sinfo.zeroInitialDisp != 1) {
        fprintf(stderr," ... DO NOT SEND IDISP6             ...\n");
      } else {
        fprintf(stderr," ... SENDING IDISP6                 ...\n");
        int i;
        for(i = 0; i < numIDis6; ++i) {
          int dof = c_dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
          if(dof >= 0)
            d_n_aero[dof] += iDis6[i].val;
        }
      }
    }
    
    State curState(c_dsa, dsa, bcx, vcx, d_n_aero, v_n, a_n, v_p);

    if(sinfo.aeroFlag == 8) {
      flExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.mppFactor,
                             restartinc, sinfo.isCollocated, sinfo.alphas, sinfo.alphasv);
      flExchanger->sendModeFreq(modeData.frequencies, modeData.numModes);
      if(verboseFlag) filePrint(stderr, " ... [E] Sent parameters and mode frequencies ...\n");
      flExchanger->sendModeShapes(modeData.numModes, modeData.numNodes,
                   modeData.modes, curState, sinfo.mppFactor);
      if(verboseFlag) filePrint(stderr, " ... [E] Sent mode shapes           ...\n");
    }
    else {
//      double aero_tmax = sinfo.tmax;
//      if(sinfo.newmarkBeta == 0) aero_tmax += sinfo.getTimeStep();
//      flExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), aero_tmax, restartinc,
//                             sinfo.isCollocated, sinfo.alphas);
//      if(verboseFlag) fprintf(stderr," ... [E] Sent parameters ...\n");

      if(sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
        flExchanger->initRcvParity(1);
        flExchanger->initSndParity(1);
      } else {
        flExchanger->initRcvParity(-1);
        flExchanger->initSndParity(-1);
      }

      if(sinfo.aeroFlag == 20 && sinfo.dyna3d_compat) {
        flExchanger->sendSubcyclingInfo(0);
        flExchanger->sendNoStructure();
      }

      if(!(geoSource->getCheckFileInfo()->hotRestart() && (sinfo.aeroFlag == 20 ||
          (sinfo.isNonLin() && sinfo.isDynam() && sinfo.newmarkBeta != 0)))) {
        flExchanger->sendDisplacements(curState);
        if(verboseFlag) filePrint(stderr, " ... [E] Sent initial displacements ...\n");
      }

      if(sinfo.aeroFlag == 1) { // Ping pong only
        fprintf(stderr, "Ping Pong Only requested. Structure code exiting\n");
      }
    }
  }
}

void
Domain::thermoePreProcess()
{
  if(sinfo.thermoeFlag >=0) {

    int buffLen = numnodes;

    temprcvd = new double[numnodes]; // Initialize received temperature

    // if sinfo.aeroFlag >= 0, flExchanger has already been initialize before,
    // thus, only when sinfo.aeroFlag < 0 is necessary.
    if(sinfo.aeroFlag < 0)
      flExchanger = new FlExchanger(nodes, packedEset, c_dsa);

    flExchanger->thermoread(buffLen);
 
    flExchanger->getStrucTemp(temprcvd) ;
    if(verboseFlag) filePrint(stderr," ... [E] Received initial temperatures ...\n");

    computeTDProps();
  }
}

/*
** This subroutine computes the maximum stability time step,
** which is equal to 2.0 divided by the maximum natural frequency
** of the structure.  The maxium natural frequency is computed
** by power method.
*/

template<typename DynamMatType>
double
Domain::computeStabilityTimeStep(DynamMatType& dMat)
{
      using std::abs;
      using std::sqrt;

      double eigmax;
      double relTol    = sinfo.stable_tol; // stable_tol default is 1.0e-3
      double preeigmax = 0.0;

      int numdofs = dMat.K->dim();
      if(numdofs == 0) return std::numeric_limits<double>::infinity();

      int maxIte  = sinfo.stable_maxit; // stable_maxit default is 100

      Vector v(numdofs);
      Vector z(numdofs);

// Starts from an arbitrary array.
      int i,j;
      for (i=0; i<numdofs; ++i)
        v[i] = (double) (i+1) / (double) numdofs;

// Power iteration loop

      for (i=0; i<maxIte; ++i) {
        dMat.K->mult(v,z);

        for (j=0; j< numdofs; ++j)
          z[j] /= dMat.M->diag(j);

// Normalize

        double zmax = z[0];
        for (j=1; j< numdofs; ++j)
          if (abs(z[j])>zmax) zmax = abs(z[j]);

        eigmax = zmax;

        v = (1.0/zmax)*z;

        if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;

        preeigmax = eigmax;
      }

      // compute stability maximum time step
      double sdt = 2.0 / sqrt(eigmax);

      return sinfo.stable_cfl*sdt;
}

template double Domain::computeStabilityTimeStep(DynamMat& dMat);
template double Domain::computeStabilityTimeStep(ModalOps& dMat);

double Domain::computeStabilityTimeStepROM(GenFullSquareMatrix<double>& K_red)
{
      using std::abs;
      using std::sqrt;

      double eigmax;
      double relTol    = sinfo.stable_tol; // stable_tol default is 1.0e-3
      double preeigmax = 0.0;

      int numdofs = K_red.dim();
      int maxIte  = sinfo.stable_maxit; // stable_maxit default is 100

      Vector v(numdofs);
      Vector z(numdofs);

// Starts from an arbitrary array.
      int i,j;
      for (i=0; i<numdofs; ++i)
        v[i] = (double) (i+1) / (double) numdofs;

// Power iteration loop

      for (i=0; i<maxIte; ++i) {
        K_red.mult(v,z);

// Normalize

        double zmax = z[0];
        for (j=1; j< numdofs; ++j)
          if (abs(z[j])>zmax) zmax = abs(z[j]);

        eigmax = zmax;

        v = (1.0/zmax)*z;

        if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;

        preeigmax = eigmax;
      }

      // compute stability maximum time step
      double sdt = 2.0 / sqrt(eigmax);

      return sinfo.stable_cfl*sdt; 
}


double
Domain::computeStabilityTimeStep(FullSquareMatrix *kelArray, FullSquareMatrix *melArray, GeomState *geomState, int &eid)
{
  double sdt = std::numeric_limits<double>::infinity();
  eid = -1;
  for(int iele = 0; iele < numele; ++iele) {
    double sdt_iele = packedEset[iele]->computeStabilityTimeStep(kelArray[iele], melArray[iele], nodes, geomState,
                                                                 sinfo.stable_tol, sinfo.stable_maxit);
    if(sdt_iele < sdt) eid = packedEset[iele]->getGlNum();
    sdt = std::min(sdt, sdt_iele);
  }
  return sinfo.stable_cfl*sdt;
}

//------------------------------------------------------------------------------


double Domain::getKineticEnergy( Vector& vel, SparseMatrix * gMass )
{
  double energy=0.0;
  
  Vector tmpVec(c_dsa->size());
    
  gMass->mult(vel,tmpVec);

  energy = 0.5 * (vel * tmpVec);
  
  return energy;
}

//------------------------------------------------------------------------------

void
Domain::updateUsddInDbc(double* userDefineDisp, int* map)
{
  int j = 0;
  for(int i = 0; i < numDirichlet; ++i)
    if(dbc[i].type == BCond::Usdd) {
      int k = (map) ? map[j] : j;
      dbc[i].val = userDefineDisp[k];
      j++;
    }
}

void
Domain::updateUsdfInNbc(double* userDefineForce, int* map, double* weight)
{
  int j = 0;
  for(int i = 0; i < numNeuman; ++i) 
    if(nbc[i].type == BCond::Usdf) {
      int k = (map) ? map[j] : j;
      nbc[i].val = userDefineForce[k];
      if(weight) {
        int dof = c_dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
        nbc[i].val *= weight[dof];
      }
      j++;
    }
}

void
Domain::updateActuatorsInNbc(double* actuatorsForce, int* map, double* weight)
{
  int j = 0;
  for(int i = 0; i < numNeuman; ++i)   
    if(nbc[i].type == BCond::Actuators) {
      int k = (map) ? map[j] : j;
      nbc[i].val = actuatorsForce[k];
      if(weight) {
        int dof = c_dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
        nbc[i].val *= weight[dof];
      }
      j++;
    }
}

void
Domain::computeReactionForce(Vector &fc, Vector &Du, Vector &Vu, Vector &Au,
                             double *bcx, double *vcx, double *acx,
                             SparseMatrix *_kuc, SparseMatrix *_kcc,
                             SparseMatrix *_cuc, SparseMatrix *_ccc,
                             SparseMatrix *_muc, SparseMatrix *_mcc)
{
  // TODO include external force on the constrained dofs
  CuCSparse *kuc = dynamic_cast<CuCSparse *>(_kuc);

  if(kuc) kuc->transposeMultNew(Du.data(), fc.data()); // fc = Kuc^T * Du
  else fc.zero();

  CuCSparse *cuc = dynamic_cast<CuCSparse *>(_cuc);
  if(cuc) cuc->transposeMultAddNew(Vu.data(), fc.data()); // fc += Cuc^T * Vu

  CuCSparse *muc = dynamic_cast<CuCSparse *>(_muc);
  if(muc) muc->transposeMultAddNew(Au.data(), fc.data()); // fc += Muc^T * Vu

  CuCSparse *kcc = dynamic_cast<CuCSparse *>(_kcc);
  CuCSparse *ccc = dynamic_cast<CuCSparse *>(_ccc);
  CuCSparse *mcc = dynamic_cast<CuCSparse *>(_mcc);
  Vector Dc(numDirichlet, 0.0);
  Vector Vc(numDirichlet, 0.0);
  Vector Ac(numDirichlet, 0.0);

  for(int i=0; i<numDirichlet; ++i) {
    int dof = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
    if(dof < 0) continue;
    int cdof = c_dsa->invRCN(dof);
    if(cdof >= 0) {
      Dc[cdof] = bcx[dof];
      Vc[cdof] = vcx[dof];
      Ac[cdof] = (acx) ? acx[dof] : 0; 
    }
  }

  if(kcc) kcc->multAddNew(Dc.data(), fc.data()); // fc += Kcc * Dc
  if(ccc) ccc->multAddNew(Vc.data(), fc.data()); // fc += Ccc * Vc
  if(mcc) mcc->multAddNew(Ac.data(), fc.data()); // fc += Mcc * Ac
}

void
Domain::setModalEnergies(double _modalWela, double _modalWkin, double _modalWdmp)
{
  modalWela = _modalWela;
  modalWkin = _modalWkin;
  Wdmp = _modalWdmp;
}

void
Domain::computeEnergies(Vector &disp, Vector &force, double time, Vector *aeroForce, Vector *vel, SparseMatrix *K,
                        SparseMatrix *M, SparseMatrix *C, double &Wela, double &Wkin, double &error)
{
  double pWext = Wext, pWdmp = Wdmp;
  computeExtAndDmpEnergies(disp, force, time, aeroForce, vel, C);

  Vector tmpVec(numUncon());
  if(sinfo.modal) {
    Wkin = modalWkin;
    Wela = modalWela;
  }
  else {
    if(M) {
      M->mult(*vel, tmpVec);
      Wkin = 0.5 * ((*vel) * tmpVec);
    }
    if(K) {
      K->mult(disp, tmpVec);
      Wela = 0.5 * (disp * tmpVec);
    }
  }

  // XXX consider sign of Wdmp in this equation:
  error = (time == sinfo.initialTime) ? 0.0 : (Wela+Wkin+Wdmp-Wext)-(pWela+pWkin+pWdmp-pWext);

  pWela = Wela;
  pWkin = Wkin;
}

void
Domain::computeExtAndDmpEnergies(Vector &disp, Vector &force, double time, Vector *aeroForce,
                                 Vector *vel, SparseMatrix *C, Vector *folForce)
{
  // compute work done by external forces and dissipation due to viscous damping
  Vector tmpVec(numUncon());

  if(time == sinfo.initialTime) {
    Wext  = 0.0;
    Waero = 0.0;
    Wdmp  = 0.0;
    previousExtForce = new Vector(force);
    if(folForce) (*previousExtForce) += (*folForce);
    if(sinfo.aeroFlag >= 0) {
      previousAeroForce = new Vector(*aeroForce);
    }
    if(C) {
      C->mult((*vel), tmpVec);
      previousCq = new Vector(tmpVec);
    }
    previousDisp = new Vector(disp);
  }
  else {
    double c = solInfo().newmarkGamma;
    // NOTE: The formula for dW should be (1-c)*f^n + c*f^{n+1}. However, force 
    //       stores f^{n+1-alpha_f} and previousExtForce stores f^n.
    //       For this reason, first f^{n+1} is approximately computed as follows
    //       f_{n+1} ≈ 1/(1-alpha_f)*(f_{n+1-alpha_f} - alpha_f*f_n)
    double alphaf = (solInfo().newmarkBeta == 0) ? 0 : solInfo().newmarkAlphaF;
    tmpVec = 1/(1-alphaf)*(force - alphaf*(*previousExtForce));
    if(folForce) tmpVec += (*folForce);
    Wext += (c*tmpVec + (1.0-c)*(*previousExtForce)) * (disp - (*previousDisp));
    (*previousExtForce) = tmpVec;

    if(sinfo.aeroFlag >= 0) {
      tmpVec = 1/(1-alphaf)*((*aeroForce) - alphaf*(*previousAeroForce));
      Waero += (c*tmpVec + (1.0-c)*(*previousAeroForce)) *(disp - (*previousDisp));
      (*previousAeroForce) = tmpVec;
    }

    if(C) {
      C->mult((*vel), tmpVec);
      Wdmp += (c*tmpVec + (1.0-c)*(*previousCq))*(disp - (*previousDisp));
      (*previousCq) = tmpVec;
    }
    (*previousDisp) = disp;
  }
}

void
Domain::computeUnamplifiedExtForce(GenVector<double>& fcon, int loadsetid)
{
  // Compute constant part (f^con) of a time-dependent external force of the form
  // f^ext(t) = lambda(t)*f^con belonging to the specified loadset, where lambda(t) is the scalar amplification
  // factor obtained from an MFTT/HFTT.

  fcon.zero();

  for(int i = 0; i < numNeuman; ++i) {
    if(nbc[i].loadsetid != loadsetid) continue; // skip forces not belonging to the specified loadset
    if(sinfo.isNonLin() && nbc[i].type == BCond::Forces && 
       !(nbc[i].mtype == BCond::Axial && nbc[i].dofnum < 3)) continue; // skip configuration-dependent forces
    int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
    if(dof < 0) continue;
    switch(nbc[i].type) {
      case(BCond::Forces) : if(MFTTData *mftt = domain->getMFTT(nbc[i].loadsetid)) fcon[dof] += nbc[i].val; break;
      case(BCond::Flux)   : if(MFTTData *hftt = domain->getHFTT(nbc[i].loadsetid)) fcon[dof] += nbc[i].val; break;
      default : /* all other cases are not included, e.g. Usdf and Actuators */ ;
    }
  }
}
