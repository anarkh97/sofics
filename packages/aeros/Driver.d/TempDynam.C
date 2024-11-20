#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <iostream>
#include <algorithm>

#include <Driver.d/Domain.h>
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#include <Utils.d/linkfc.h>
#include <Driver.d/Dynam.h>
#include <Hetero.d/FlExchange.h>
#include <Corotational.d/GeomState.h>
#include <Element.d/State.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/GeoSource.h>

extern int verboseFlag;

void
Domain::getHeatFlux(Vector &tsol, double *bcx, int fileNumber, int hgIndex,
                    double time)
{
  // Postprocessing: Computes the NODAL heat fluxes or the temperature 
  //                 gradients of the QUAD element
  if(outFlag && !nodeTable) makeNodeTable(outFlag);
  int numNodes = (outFlag) ? exactNumNodes : geoSource->numNode();

  OutputInfo *oinfo = geoSource->getOutputInfo();
  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  int k;

  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;

  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;

  // ... ALLOCATE VECTORS AND INITIALIZE TO ZERO
  if(heatflux == 0) heatflux = new Vector(numNodes,0.0);
  if(elTemp == 0) elTemp = new Vector(maxNumDOFs,0.0);

  int iele;
  if(elheatflux == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    elheatflux = new Vector(maxNodesPerElement, 0.0);
  }

  // zero the vectors
  heatflux->zero();
  elTemp->zero();
  elheatflux->zero();

  for(iele=0; iele<numele; ++iele) {
     int NodesPerElement = elemToNode->num(iele);
     packedEset[iele]->nodes(nodeNumbers);

  // ... DETERMINE ELEMENT TEMPERATURE VECTOR
     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elTemp)[k] = tsol[cn];
        else
          (*elTemp)[k] = bcx[(*allDOFs)[iele][k]];
     }

// ... CALCULATE HEAT FLUX or TEMPERATURE GRADIENT VALUE FOR EACH NODE 
//     OF THE ELEMENT

      packedEset[iele]->computeHeatFluxes(*elheatflux, nodes, *elTemp, hgIndex);

// ... ASSEMBLE ELEMENT'S NODAL HEAT FLUXES OR TEMPERATURE GRADIENTS

     for(k=0; k<NodesPerElement; ++k) {
       int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
       (*heatflux)[node] += (*elheatflux)[k];
     }

    }

// ... PRINT HEAT FLUXES OR TEMPERATURE GRADIENTS DEPENDING ON hgIndex

    if (oinfo[fileNumber].nodeNumber == -1)
      geoSource->outputNodeScalars(fileNumber, heatflux->data(), numNodes, time);
    else
      geoSource->outputNodeScalars(fileNumber, heatflux->data()+oinfo[fileNumber].nodeNumber, 1, time);


}

void
Domain::getTrussHeatFlux(Vector &tsol, double *bcx, int fileNumber, int hgIndex,
                            double time)
{
 // Postprocessing: Computes the ELEMENT heat fluxes or the temperature 
 //                 gradients of the BAR element
 
  OutputInfo *oinfo = geoSource->getOutputInfo();
  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];
 
  int k;
 
  // ... OUTPUT FILE field width
  int w = oinfo[fileNumber].width;
 
  // ... OUTPUT FILE precision
  int p = oinfo[fileNumber].precision;
 
  // ... ALLOCATE VECTOR Element temperature  AND INITIALIZE TO ZERO
  if(elTemp == 0) elTemp = new Vector(maxNumDOFs,0.0);
 
  int iele;
  double eltrussflux = 0.;
 
  // zero the vectors
  elTemp->zero();
 
  // ... WRITE CURRENT TIME VALUE
  if(oinfo[fileNumber].nodeNumber == -1)
    fprintf(oinfo[fileNumber].filptr,"%*.*e\n",w,p,time);
 
  for(iele=0; iele<numele; ++iele) {
     //int NodesPerElement = elemToNode->num(iele);
     packedEset[iele]->nodes(nodeNumbers);
 
  // ... DETERMINE ELEMENT TEMPERATURE VECTOR
 
     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elTemp)[k] = tsol[cn];
        else
          (*elTemp)[k] = bcx[(*allDOFs)[iele][k]];
     }
 
// ... CALCULATE HEAT FLUX or TEMPERATURE GRADIENT  FOR EACH  ELEMENT
 
      packedEset[iele]->trussHeatFluxes(eltrussflux, nodes, *elTemp, hgIndex); 

// ... PRINT ELEMENT HEAT FLUXES OR TEMPERATURE GRADIENTS
 
       fprintf(oinfo[fileNumber].filptr,"% *.*E\n",w,p,eltrussflux);
  }
    fflush(oinfo[fileNumber].filptr);
} 

void
Domain::initTempVector(Vector& d_n, Vector& v_n, Vector& v_p)
{

 int i;

// ... INITIALIZE TEMPERATURE AND ITS FIRST DERIVATIVE
// ... VECTORS TO ZERO
 d_n = v_n = v_p = 0.0;

// ... INITIALIZE TEMPERATURE
 for(i = 0; i < numIDis; ++i) {
   int dof = c_dsa->locate(iDis[i].nnum, 1 << iDis[i].dofnum);
   if(dof >= 0)
     d_n[dof] = iDis[i].val;
 }

 ControlInfo *cinfo = geoSource->getCheckFileInfo();
 if(probType() != SolverInfo::NonLinDynam) {
   if(cinfo->lastRestartFile) {
     fprintf(stderr, " ... Restarting From a Previous Run ...\n");
     int fn = open(cinfo->lastRestartFile,O_RDONLY );
     if(fn >= 0) {
       int vsize, restartTIndex;
       double restartT;
       int readSize = read(fn, &vsize, sizeof(int));
       if(vsize != d_n.size() || readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       readSize = read(fn, &restartTIndex, sizeof(int));
       if(readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTimeIndex = 0;
       else
         sinfo.initialTimeIndex = restartTIndex;
       readSize = read(fn, &restartT, sizeof(double));
       if(readSize != int(sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTime = 0.0;
       else
         sinfo.initialTime = restartT;
       readSize = read(fn, d_n.data(), vsize*sizeof(double));
       if(readSize != vsize*int(sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, v_n.data(), vsize*sizeof(double));
       if(readSize != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, v_p.data(), vsize*sizeof(double));
       if(readSize != int(vsize*sizeof(double))) {
         fprintf(stderr," *** WARNING: Older version of restart file"
                        " -- Missing velocity field is set to zero\n");
         v_p.zero();
       }
       readSize = read(fn, &(sinfo.initExtForceNorm), sizeof(double));
       if(readSize != int(sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file: External Force Norm\n");

       close(fn);
     } else {
        perror(" *** ERROR: Restart file could not be opened: ");
        //exit(-1);
     }

   }
   //else {
   //    fprintf(stderr, " ... No restart                     ...\n");
   //}
 }


}


void
Domain::tempdynamOutput(int tIndex, double *bcx, DynamMat& dMat,
                        Vector& ext_f, Vector& d_n, Vector& v_n, Vector&v_p)
{
// The condition below (tIndex != 0) has to be maintained, because
//in TempProbType.C, we already sent the initial temperature during
//the pre-process. As this function is called after the latter for the
//sole purpose of initializing the printing output, the mentioned
//condition is neede.  

  int i;

  enum {HFLX=0, HFLY=1, HFLZ=2, GRTX=3, GRTY=4, GRTZ=5};

  double time = tIndex*sinfo.getTimeStep();

  if (sinfo.nRestart > 0 && !sinfo.modal) {
#ifndef SALINAS
    writeRestartFile(time, tIndex, d_n, v_n, v_p, sinfo.initExtForceNorm);
#endif
  }

  State tempState(c_dsa, dsa, bcx, d_n, v_n, v_p);

  if(sinfo.aeroheatFlag >= 0)
    if (tIndex != 0) {
      flExchanger->sendTemperature(tempState);
      if(verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
    }

  if(sinfo.thermohFlag >=0) 
    if (tIndex != 0) {
      int iNode;
      Vector tempsent(numnodes);

      for(iNode=0; iNode<numnodes; ++iNode) {
        int tloc  = c_dsa->locate( iNode, DofSet::Temp);
        int tloc1 =   dsa->locate( iNode, DofSet::Temp);
        double temp = (tloc >= 0) ? d_n[tloc] : bcx[tloc1];
        if(tloc1 < 0) temp = 0.0;
        tempsent[iNode] = temp; 
      }

      flExchanger->sendStrucTemp(tempsent);
      if(verboseFlag) filePrint(stderr, " ... [T] Sent temperatures          ...\n");
    }

// Open the file and write the header in the first time step

  if (tIndex==0) {
    geoSource->openOutputFiles();
  }
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

//   Print out the temperature info
  int numNodes = geoSource->numNode();

  for (i=0; i < numOutInfo; i++) {
    int iNode;
  if(oinfo[i].interval != 0 && tIndex % oinfo[i].interval == 0) {
    int w = oinfo[i].width;
    int p = oinfo[i].precision;

    switch(oinfo[i].type) {
       case OutputInfo::Temperature:
         if(oinfo[i].nodeNumber == -1) {
           fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
           for(iNode=0; iNode<numNodes; ++iNode) {
              int tloc  = c_dsa->locate(iNode, DofSet::Temp);
              int tloc1 =   dsa->locate(iNode, DofSet::Temp);
              double temp = (tloc >= 0) ? d_n[tloc] : bcx[tloc1];
              if(tloc1 < 0) temp = 0.0;
              fprintf(oinfo[i].filptr," % *.*E\n",w,p,temp);
            }
         } else {
           // if only one node is requested as output
           iNode = oinfo[i].nodeNumber;
           int tloc  = c_dsa->locate(iNode, DofSet::Temp);
           int tloc1 =   dsa->locate(iNode, DofSet::Temp);
           double temp = (tloc >= 0) ? d_n[tloc] : bcx[tloc1];
           if(tloc1 < 0) temp = 0.0;
           fprintf(oinfo[i].filptr,"  %e % *.*E\n", time, w, p, temp);
           }
         fflush(oinfo[i].filptr);
         break;
       case OutputInfo::TemperatureFirstTimeDerivative: // this is the first time derivative of the temperature
         if(oinfo[i].nodeNumber == -1) {
           fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
           for(iNode=0; iNode<numNodes; ++iNode) {
             int tloc  = c_dsa->locate(iNode, DofSet::Temp);
             int tloc1 =   dsa->locate(iNode, DofSet::Temp);
             double val = (tloc >= 0) ? v_n[tloc] : 0.0;
             if(tloc1 < 0) val = 0.0;
             fprintf(oinfo[i].filptr," % *.*E\n",w,p,val);
           }
         } else {
           // if only one node is requested as output
           iNode = oinfo[i].nodeNumber;
           int tloc  = c_dsa->locate(iNode, DofSet::Temp);
           int tloc1 =   dsa->locate(iNode, DofSet::Temp);
           double val = (tloc >= 0) ? v_n[tloc] : 0.0;
           if(tloc1 < 0) val = 0.0;
           fprintf(oinfo[i].filptr," %e % *.*E\n", time, w, p, val);
         }
         fflush(oinfo[i].filptr);
         break;
       case OutputInfo::HeatFlXX:
         getHeatFlux(d_n, bcx, i, HFLX, time);
         break;
       case OutputInfo::HeatFlXY:
         getHeatFlux(d_n, bcx, i, HFLY, time);
         break;
       case OutputInfo::HeatFlXZ:
         getHeatFlux(d_n, bcx, i, HFLZ, time);
         break;
       case OutputInfo::GrdTempX:
         getHeatFlux(d_n, bcx, i, GRTX, time);
         break;
       case OutputInfo::GrdTempY:
         getHeatFlux(d_n, bcx, i, GRTY, time);
         break;
       case OutputInfo::GrdTempZ:
         getHeatFlux(d_n, bcx, i, GRTZ, time);
         break;
       case OutputInfo::HeatFlX:
         getTrussHeatFlux(d_n, bcx, i, HFLX, time);
         break;
       case OutputInfo::GrdTemp:
         getTrussHeatFlux(d_n, bcx, i, GRTX, time);
         break;
       default:
         fprintf(stderr, " *** WARNING: output %d is not supported \n", i);
         break;
       }
     }

  }
}


void
Domain::computeExtForce(Vector &f, double t, int tIndex, SparseMatrix *kuc, Vector &prev_f)
{
  // ... ZERO THE FORCE VECTOR
  f.zero();

  // ... Temporary variable for inter(extra)polated force
  double *tmpFmem = (double *) dbg_alloca(sizeof(double)*numUncon());
  StackVector tmpF(tmpFmem, numUncon());

  // ... CONSTRUCT Vc
  Vector Vc(numDirichlet);

  // CONSTRUCT THE NON-HOMONGENOUS DIRICHLET BC VECTOR
  int i;
  for(i=0; i<numDirichlet; ++i) {
    int dof2 = dsa->locate(dbc[i].nnum,(1 << dbc[i].dofnum));
    dof2 = c_dsa->invRCN(dof2);
    double scaleFactor=1.0;
    if(dof2 >= 0) Vc[dof2] = scaleFactor*dbc[i].val;
  }

  // ... COMPUTE NON-HOMOGENEOUS FORCE ( f = f - kuc*Vc)
  if(kuc)
    kuc->multSubtract(Vc, f);

  // ... COMPUTE TIME DEPENDENT EXTERNAL FORCE
  for(i=0; i < numNeuman; ++i) {
    int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
    if(dof < 0) continue;
    switch(nbc[i].type) {
      case(BCond::Flux) : {
        MFTTData *mftt = domain->getHFTT(nbc[i].loadsetid);
        double loadFactor = (mftt && domain->solInfo().isDynam()) ? mftt->getVal(t) : domain->getLoadFactor(nbc[i].loadsetid);
        f[dof] += loadFactor*nbc[i].val;
      } break;
      case(BCond::Convection) : {
        double loadFactor = domain->getLoadFactor(nbc[i].loadsetid);
        f[dof] += loadFactor*nbc[i].val;
      } break;
      default : {
        f[dof] += nbc[i].val; 
      } break;
    }
  }

  // ... ADD FLUID FLUX
  if(sinfo.aeroheatFlag >= 0 && tIndex >= 0) {
    int j;
    double sflux = 0;
    double bfl ;

    tmpF.zero();
    flExchanger->getFluidFlux(tmpF, t, bfl);

    if(verboseFlag) filePrint(stderr, " ... [T] Received fluid fluxes      ...\n");

    int vectlen = tmpF.size();

//  Compute fluid flux at n+1/2, since we use midpoint rule in thermal

    int useProjector = domain->solInfo().filterFlags;

    if(tIndex == 0) 
      f += tmpF;
    else {
      if(useProjector) f = tmpF;
      else
        f.linAdd(0.5, tmpF, 0.5, prev_f);
    }

    prev_f = tmpF; 

//  Print out sum of fluxes received at n+1 and at n+1/2 in a file
//  specified by raerotfl in entry

    int numOutInfo = geoSource->getNumOutInfo();
    OutputInfo *oinfo = geoSource->getOutputInfo();

    for(i=0; i < numOutInfo; i++) {
      if(oinfo[i].interval != 0 && tIndex % oinfo[i].interval == 0) {
        switch(oinfo[i].type) {
          case OutputInfo::AeroFlux:
            fprintf(oinfo[i].filptr,"%e   ",t);
            fprintf(oinfo[i].filptr,"%e   ",bfl);
            for(j=0; j<vectlen; j++)  {
              sflux += 0.5*(tmpF[j]+prev_f[j]);
            } 
            fprintf(oinfo[i].filptr,"%e\n",sflux); 
            fflush(oinfo[i].filptr);
          break;
          default: // do nothing
          break;
        }
      }
    }
  }
}

void
Domain::aeroHeatPreProcess(Vector& d_n, Vector& v_n, Vector& v_p, double* bcx)
{
 if(sinfo.aeroheatFlag >= 0) {

   State tempState(c_dsa, dsa, bcx, d_n, v_n, v_p);

   int flag = 0;

   // Check if aero fluxes are requested for output
    
   int iInfo;
   int numOutInfo = geoSource->getNumOutInfo();
   OutputInfo *oinfo = geoSource->getOutputInfo();

   for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
     if(oinfo[iInfo].type == OutputInfo::AeroForce) {
       flag = 1;
       break;
     }
   }

   if(flag)
     flExchanger = new FlExchanger(nodes, packedEset, c_dsa, oinfo+iInfo);
   else
     flExchanger = new FlExchanger(nodes, packedEset, c_dsa );

   const char *matchFile = geoSource->getMatchFileName();
   if (matchFile == 0)
     matchFile = (char*) "MATCHER";

   flExchanger->read(0, matchFile);

   flExchanger->negotiate();

   int restartinc = (solInfo().nRestart >= 0) ? (solInfo().nRestart) : 0;

   flExchanger->sendTempParam(sinfo.aeroheatFlag, sinfo.getTimeStep(), 
                              sinfo.tmax, restartinc, sinfo.alphat);

   if(verboseFlag) filePrint(stderr, " ... [T] Sent temperature parameters...\n");

   flExchanger->sendTemperature(tempState);
   if(verboseFlag) filePrint(stderr, " ... [T] Sent initial temperatures  ...\n");
 }
}

void
Domain::thermohPreProcess(Vector& d_n, Vector& v_n, Vector& v_p, double* bcx)
{
  if(sinfo.thermohFlag >=0) {

     int tempbuffLen = numnodes;

     if (sinfo.aeroheatFlag < 0)
     flExchanger = new FlExchanger(nodes, packedEset, c_dsa );

     flExchanger->thermoread(tempbuffLen);

/* Send initial Temperature */

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
     if(verboseFlag) filePrint(stderr, " ... [T] Sent initial temperatures  ...\n");
  }
}

