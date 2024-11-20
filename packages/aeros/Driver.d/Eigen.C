#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <cmath>
#include <fstream>
#include <algorithm>

#include <Utils.d/dofset.h>
#include <Driver.d/Domain.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Math.d/VectorSet.h>
#include <Utils.d/resize_array.h>
#include <Utils.d/BinFileHandler.h>

#include <Driver.d/GeoSource.h>

void
Domain::getSloshDisp(Vector &sloshPotSol, double *bcx, int fileNumber, int hgIndex, double time)
{
  // Postprocessing: Computes the NODAL displacements for sloshing problem
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
  if(fluidDispSlosh == 0) fluidDispSlosh = new Vector(numNodes,0.0);
  if(elPotSlosh == 0) elPotSlosh = new Vector(maxNumDOFs,0.0);

  int iele;
  if(elFluidDispSlosh == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    elFluidDispSlosh = new Vector(maxNodesPerElement, 0.0);
  }

  // zero the vectors
  fluidDispSlosh->zero();
  elPotSlosh->zero();
  elFluidDispSlosh->zero();

  for(iele=0; iele<numele; ++iele) {
     int NodesPerElement = elemToNode->num(iele);
     packedEset[iele]->nodes(nodeNumbers);

  // ... DETERMINE ELEMENT POTENTIAL VECTOR
     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elPotSlosh)[k] = sloshPotSol[cn];
        else
          (*elPotSlosh)[k] = bcx[(*allDOFs)[iele][k]];
     }

// ... CALCULATE FLUID DISPLACEMENT VALUE FOR EACH NODE 
//     OF THE ELEMENT

      packedEset[iele]->computeSloshDisp(*elFluidDispSlosh, nodes, *elPotSlosh, hgIndex);

// ... ASSEMBLE ELEMENT'S NODAL FLUID DISPLACEMENTS

     for(k=0; k<NodesPerElement; ++k) {
       int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
       (*fluidDispSlosh)[node] += (*elFluidDispSlosh)[k];
     }
    }

// ... PRINT FLUID DISPLACEMENTS DEPENDING ON hgIndex
      
  if (oinfo[fileNumber].nodeNumber == -1)
    geoSource->outputNodeScalars(fileNumber, fluidDispSlosh->data(), numNodes, time);
  else  {
    geoSource->outputNodeScalars(fileNumber, fluidDispSlosh->data()+oinfo[fileNumber].nodeNumber, 1, time);
  }

}

void
Domain::getSloshDispAll(Vector &sloshPotSol, double *bcx, int fileNumber, double time)
{
  // Postprocessing: Computes the NODAL displacements for sloshing problem
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
  double (*fluidDispSloshAll)[3] = new double[numNodes][3];

  if(elPotSlosh == 0) elPotSlosh = new Vector(maxNumDOFs,0.0);

  int iele;
  if(elFluidDispSloshAll == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    elFluidDispSloshAll = new Vector(maxNodesPerElement*3, 0.0);
  }

  // zero the vectors
  elPotSlosh->zero();
  elFluidDispSloshAll->zero();

  for(iele=0; iele<numele; ++iele) {
     int NodesPerElement = elemToNode->num(iele);
     packedEset[iele]->nodes(nodeNumbers);

  // ... DETERMINE ELEMENT POTENTIAL VECTOR
     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elPotSlosh)[k] = sloshPotSol[cn];
        else
          (*elPotSlosh)[k] = bcx[(*allDOFs)[iele][k]];
     }

// ... CALCULATE FLUID DISPLACEMENT VALUE FOR EACH NODE OF THE ELEMENT
     packedEset[iele]->computeSloshDispAll(*elFluidDispSloshAll, nodes, *elPotSlosh);

// ... ASSEMBLE ELEMENT'S NODAL FLUID DISPLACEMENTS
     for(k=0; k<NodesPerElement; ++k) {
        int node = (outFlag) ? nodeTable[(*elemToNode)[iele][k]]-1 : (*elemToNode)[iele][k];
        for(int j = 0; j<3; ++j)
          fluidDispSloshAll[node][j] += (*elFluidDispSloshAll)[3*k+j];
     }
  }

// ... PRINT FLUID DISPLACEMENTS DEPENDING ON hgIndex
  if (oinfo[fileNumber].nodeNumber == -1)
    geoSource->outputNodeVectors(fileNumber, fluidDispSloshAll, numNodes, time);
  else
    geoSource->outputNodeVectors(fileNumber, fluidDispSloshAll+oinfo[fileNumber].nodeNumber, 1, time);

}

void
Domain::eigenOutput(Vector& eigenValues, VectorSet& eigenVectors, double* bcx, int convEig)
{
  const double pi = 3.141592653589793;
  
  int maxmode = (convEig && convEig < eigenValues.size()) ? convEig : eigenValues.size(); 

  if (sinfo.modeDecompFlag) {
    fprintf(stderr," ... Outputting EigenVectors in Binary file requested ...\n");
    fflush(stderr);
    BinFileHandler modefile("EIGENMODES" ,"w");
    modefile.write(&maxmode, 1);
    int veclength = eigenVectors[1].size();
    modefile.write(&veclength, 1);

    // Write eigenmodes in file EIGENMODES
    for (int imode=0; imode < maxmode; ++imode) {
      double *modes = eigenVectors[imode].data();
      modefile.write(modes, veclength);
    }
  }

  // Create a dummy Vector
  int veclength = eigenVectors[0].size();
  GenVector<double> dummyVector(veclength,0.0);

  for(int imode=0; imode < maxmode; ++imode) {
    if(domain->solInfo().buckling || solInfo().soltyp == 2)
      domain->postProcessing<double>(eigenVectors[imode],bcx,dummyVector,0, 0, 0.0, eigenValues[imode]);
    else if(domain->solInfo().sloshing)
      domain->postProcessing<double>(eigenVectors[imode],bcx,dummyVector,0, 0, 0.0, sqrt(eigenValues[imode])/(2.0*pi)*sqrt(gravitySloshing));
    else {
      double freq = (eigenValues[imode] < 0.0) ? 0.0 : sqrt(eigenValues[imode])/(2.0*pi);
      domain->postProcessing<double>(eigenVectors[imode],bcx,dummyVector,0, 0, 0.0, freq);
    }
  }

 // --- Print Problem statistics to screen ------------------------------
 if(!domain->solInfo().doEigSweep) {

   fprintf(stderr," --------------------------------------\n");

   // Print omega and frequency values to screen
   if(solInfo().buckling || solInfo().soltyp == 2) {
     fprintf(stderr," Mode\tLambda\n");
     fprintf(stderr," --------------------------------------\n");
     int imode;
     for(imode=0; imode<eigenValues.size(); ++imode)
       fprintf(stderr," %d\t%e\n",imode+1,eigenValues[imode]);
     fprintf(stderr," --------------------------------------\n");
   } else {
     if(domain->solInfo().sloshing) {
       fprintf(stderr," Mode\tLambda\t\tSloshing Frequency\n");
       fprintf(stderr," \t(=Omega^2/g)\t\t(=sqrt(Lambda*g)/(2*pi))\n");
     }
     else 
       fprintf(stderr," Mode\tOmega^2\t\tFrequency\n");
     fprintf(stderr," --------------------------------------\n");
     int imode;
     for(imode=0; imode<eigenValues.size(); ++imode) {
       if(domain->solInfo().sloshing)
         fprintf(stderr," %d\t%e\t%e\n",imode+1,eigenValues[imode], sqrt(eigenValues[imode])/(2.0*pi)*sqrt(gravitySloshing));
       else
         fprintf(stderr," %d\t%e\t%e\n",imode+1,eigenValues[imode], sqrt(eigenValues[imode])/(2.0*pi));
     }
     fprintf(stderr," --------------------------------------\n");
   }

 }
}

#ifdef USE_EIGEN3
void
Domain::eigenQROutput(Eigen::MatrixXd& Xmatrix, Eigen::MatrixXd& Qmatrix, Eigen::MatrixXd& Rmatrix)
{
  const char* Xoutput = domain->solInfo().xmatrixname;
  const char* Qoutput = domain->solInfo().qmatrixname;
  const char* Routput = domain->solInfo().rmatrixname;
  std::ofstream xout(Xoutput, std::ios::out);
  std::ofstream qout(Qoutput, std::ios::out);
  std::ofstream rout(Routput, std::ios::out);

  if(!xout) {
    std::cerr << "Error: cannot open file " << Xoutput << std::endl;
    exit(-1);
  }
  if(!qout) {
    std::cerr << "Error: cannot open file " << Qoutput << std::endl;
    exit(-1);
  }
  if(!rout) {
    std::cerr << "Error: cannot open file " << Routput << std::endl;
    exit(-1);
  }

  // write X, Q and R matrix
  Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, " ");
  xout << Xmatrix.rows() << std::endl << Xmatrix.cols() << std::endl;
  xout << Xmatrix.transpose().format(HeavyFmt) << std::endl;
  qout << Qmatrix.rows() << std::endl << Qmatrix.cols() << std::endl;
  qout << Qmatrix.transpose().format(HeavyFmt) << std::endl;
  rout << Rmatrix.rows() << std::endl << Rmatrix.cols() << std::endl;
  rout << Rmatrix.format(HeavyFmt) << std::endl;

  xout.close();
  qout.close();
  rout.close();
}
#endif
