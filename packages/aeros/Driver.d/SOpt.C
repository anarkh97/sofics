#include <cstdlib>
#include <cstdio>

#include <Driver.d/Domain.h>
#include <Hetero.d/FlExchange.h>

#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optsol.h>
#include <Structopt.d/Structopt.h>

extern Domain * domain;
extern int yyoptparse(void);

#ifdef AEROELASTIC
  #include <Utils.d/linkfc.h>
  extern "C"      {
    void _FORTRAN(inicom)(int &);
    void _FORTRAN(endcom)();
  }
#endif


//------------------------------------------------------------------------------

void Domain::structoptSolve() {

  optpro    = new Optpro;                    //Define Optimization Problem

  optsol    = new Optsol;                    //Define Solution Strategy

  structopt = new Structopt();               //Define Structural Opt.Prob.
                                             //alternatively: Fluid Opt.Prob.

  structoptInput();                          // Reading Optimiztion Input File

#ifdef AEROELASTIC

  int structID = 2;  

  _FORTRAN(inicom)(structID);

#endif
  
  structopt->build(domain,optpro,optsol);    //Build Structural Optimization    
                                             
  optsol->solve(optpro,structopt);           //Solve Optimization problem

  fclose(optunitout);                        //Closing Optimization Output
  fclose(optprotout);

#ifdef AEROELASTIC

  _FORTRAN(endcom)();

#endif
}

//------------------------------------------------------------------------------

void Domain::structoptInput() {

  //Reading Optimization Input-File

  FILE *optin = freopen(optinputfile,"r",stdin);
  if (optin == 0) {
    fprintf(stderr,"\n  *** ERROR: Could not open input file: %s\n",
            optinputfile);
    exit(-1);
  }

  //Checking for Errors during Reading
    
  int opterror = yyoptparse();

  if (opterror) {
  fprintf(stderr,
  "\n *** ERROR: Optimization-Input file contained errors. Aborting ... \n");
  exit(opterror);
  }
  
  //Open Outputfiles for Optimization
  
  char * basename = getBasename(optinputfile); 
  int fnamesize   = strlen(basename)+4;
  char * optfile  = (char*) malloc(sizeof(char*)*(fnamesize));
  strcpy(optfile,basename);
  strcat(optfile,".opt"); 
  optunitout = fopen(optfile,"w"); 

  optsol->fsize       = fnamesize;
  optsol->optprotfile = (char*) malloc(sizeof(char*)*(fnamesize));
  strcpy(optsol->optprotfile,basename);
  strcat(optsol->optprotfile,".nlp");   
  optprotout = fopen(optsol->optprotfile,"w"); 

  //Printing Problem and Solution Strategy Data

  optpro->print();                          
  optsol->print();                           
}

//------------------------------------------------------------------------------

void
Domain::postProcessing(Vector &sol, double *bcx, Vector& force, double time)
{

   double oldtime=sinfo.dt;
   sinfo.dt=time;

   enum {SXX=0,SYY=1,SZZ=2,SXY= 3,SYZ= 4,SXZ= 5,VON=6,
         EXX=7,EYY=8,EZZ=9,EXY=10,EYZ=11,EXZ=12,STRAINVON=13,
         VONTOP=14,VONBOT=15};

   enum {INX,INY,INZ,AXM,AYM,AZM};
   
   enum {YOUNG,MDENS,THICK};

   enum {PSTRESS1=0,PSTRESS2=1,PSTRESS3=2,
         PSTRAIN1=3,PSTRAIN2=4,PSTRAIN3=5};

   int iNode;
   
   Node oldnod;

 // Wext = external energy
 // Wela = elastic energy
 // Wkin = kinetic energy
 // Wdmp = damping energy
 // Total Energy = Wext+Wela+Wkin+Wdmp 

 double Wext=0.0,Wela=0.0,Wkin=0.0,Wdmp=0.0;

 int i;
 for(i=0; i<numOutInfo; ++i) 
 {
   if(oinfo[i].interval == 1) {
     int w = oinfo[i].width;
     int p = oinfo[i].precision;
     switch(oinfo[i].type)
     {
       default:
         break;
       case OutputInfo::Displacement:
         fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
         for(iNode=0; iNode<numnodes; ++iNode) {

            // If you want to skip nodes that do not exist.
            // i.e. node numbering gaps
            // if(dsa->firstdof(iNode) == -1) continue;

            int xloc  = c_dsa->locate( iNode, DofSet::Xdisp);
            int xloc1 =   dsa->locate( iNode, DofSet::Xdisp);

            double x,y,z;

	    if(xloc >= 0) 
              x = sol[xloc];
 	    else if(xloc1 >= 0)
              x = bcx[xloc1];
            else
              x = 0.0;

            int yloc  = c_dsa->locate( iNode, DofSet::Ydisp);
            int yloc1 =   dsa->locate( iNode, DofSet::Ydisp);

	    if(yloc >= 0) 
              y = sol[yloc];
 	    else if(xloc1 >= 0)
              y = bcx[yloc1];
            else
              y = 0.0;

            int zloc  = c_dsa->locate( iNode, DofSet::Zdisp);
            int zloc1 =   dsa->locate( iNode, DofSet::Zdisp);

	    if(zloc >= 0) 
              z = sol[zloc];
 	    else if(zloc1 >= 0)
              z = bcx[zloc1];
            else
              z = 0.0;

            // ... considering the variation of shape
   	    oldnod=nodescopy->getNode(iNode);
            double vx = nodes[iNode]->x - oldnod.x ;
            double vy = nodes[iNode]->y - oldnod.y ;
            double vz = nodes[iNode]->z - oldnod.z ;
	    
	    x=x+vx;   y=y+vy;    z=z+vz;

            fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                                   ,w,p,x,w,p,y,w,p,z);
         }
         fflush(oinfo[i].filptr);
         break;
       case OutputInfo::Disp6DOF:
         fprintf(oinfo[i].filptr,"  % *.*E\n",w,p,time);
         for(iNode=0; iNode<numnodes; ++iNode) {
            int xloc  = c_dsa->locate( iNode, DofSet::Xdisp);
            int xloc1 =   dsa->locate( iNode, DofSet::Xdisp);

            double x,y,z,xr,yr,zr;

            if(xloc >= 0)	// dof exists and is free
              x = sol[xloc];
            else if(xloc1 >= 0)	// dof exists and is constrained
              x = bcx[xloc1];
            else	        // dof does not exist
              x = 0.0;

            int yloc  = c_dsa->locate( iNode, DofSet::Ydisp);
            int yloc1 =   dsa->locate( iNode, DofSet::Ydisp);

            if(yloc >= 0)
              y = sol[yloc];
            else if(xloc1 >= 0)
              y = bcx[yloc1];
            else
              y = 0.0;

            int zloc  = c_dsa->locate( iNode, DofSet::Zdisp);
            int zloc1 =   dsa->locate( iNode, DofSet::Zdisp);

	    if(zloc >= 0) 
              z = sol[zloc];
 	    else if(zloc1 >= 0)
              z = bcx[zloc1];
            else
              z = 0.0;

            int xrot  = c_dsa->locate( iNode, DofSet::Xrot);
            int xrot1 =   dsa->locate( iNode, DofSet::Xrot);

            if(xrot >= 0)
              xr = sol[xrot];
            else if(xrot1 >= 0)
              xr = bcx[xrot1];
            else
              xr = 0.0;

            int yrot  = c_dsa->locate( iNode, DofSet::Yrot);
            int yrot1 =   dsa->locate( iNode, DofSet::Yrot);

            if(yrot >= 0)
              yr = sol[yrot];
            else if(yrot1 >= 0)
              yr = bcx[yrot1];
            else
              yr = 0.0;

            int zrot  = c_dsa->locate( iNode, DofSet::Zrot);
            int zrot1 =   dsa->locate( iNode, DofSet::Zrot);

            if(zrot >= 0)
              zr = sol[zrot];
            else if(zrot1 >= 0)
              zr = bcx[zrot1];
            else
              zr = 0.0;

            // ... considering the variation of shape
   	    oldnod=nodescopy->getNode(iNode);
            double vx = nodes[iNode]->x - oldnod.x ;
            double vy = nodes[iNode]->y - oldnod.y ;
            double vz = nodes[iNode]->z - oldnod.z ;
	    
	    x=x+vx;   y=y+vy;    z=z+vz;

            fprintf(oinfo[i].filptr,
              "%d % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n"
              ,iNode+1, w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
         }
         fflush(oinfo[i].filptr);
	 break;
       case OutputInfo::Temperature:
         fprintf(oinfo[i].filptr,"  %f\n",time);
         for(iNode=0; iNode<numnodes; ++iNode) {
            int tloc  = c_dsa->locate( iNode, DofSet::Temp);
            int tloc1 =   dsa->locate( iNode, DofSet::Temp);

	    double temp;

	    if(tloc >= 0)	// dof exists and is free
              temp = sol[tloc];
            else if(tloc1 >= 0)	// dof exists and is constrained
              temp = bcx[tloc1];
            else		// dof does not exist
              temp = 0.0;

            fprintf(oinfo[i].filptr,"% *.*E\t\n",w,p,temp);
          }
         fflush(oinfo[i].filptr);
         break;
       case OutputInfo::StressXX:
         getStressStrain(sol,bcx,i,SXX);
         break;
       case OutputInfo::StressYY:
         getStressStrain(sol,bcx,i,SYY);
         break;
       case OutputInfo::StressZZ:
         getStressStrain(sol,bcx,i,SZZ);
         break;
       case OutputInfo::StressXY:
         getStressStrain(sol,bcx,i,SXY);
         break;
       case OutputInfo::StressYZ:
         getStressStrain(sol,bcx,i,SYZ);
         break;
       case OutputInfo::StressXZ:
         getStressStrain(sol,bcx,i,SXZ);
         break;
       case OutputInfo::StrainXX:
         getStressStrain(sol,bcx,i,EXX);
         break;
       case OutputInfo::StrainYY:
         getStressStrain(sol,bcx,i,EYY);
         break;
       case OutputInfo::StrainZZ:
         getStressStrain(sol,bcx,i,EZZ);
         break;
       case OutputInfo::StrainXY:
         getStressStrain(sol,bcx,i,EXY);
         break;
       case OutputInfo::StrainYZ:
         getStressStrain(sol,bcx,i,EYZ);
         break;
       case OutputInfo::StrainXZ:
         getStressStrain(sol,bcx,i,EXZ);
         break;
       case OutputInfo::StressVM:
         getStressStrain(sol,bcx,i,VON);
         break;
       case OutputInfo::StressPR1:
         getPrincipalStress(sol,bcx,i,PSTRESS1);
         break;
       case OutputInfo::StressPR2:
         getPrincipalStress(sol,bcx,i,PSTRESS2);
         break;
       case OutputInfo::StressPR3:
         getPrincipalStress(sol,bcx,i,PSTRESS3);
         break;
       case OutputInfo::StrainPR1:
         getPrincipalStress(sol,bcx,i,PSTRAIN1);
         break;
       case OutputInfo::StrainPR2:
         getPrincipalStress(sol,bcx,i,PSTRAIN2);
         break;
       case OutputInfo::StrainPR3:
         getPrincipalStress(sol,bcx,i,PSTRAIN3);
         break;
       case OutputInfo::InXForce:
         getElementForces(sol, bcx, i, INX);
         break;
       case OutputInfo::InYForce:
         getElementForces(sol, bcx, i, INY);
         break;
       case OutputInfo::InZForce:
         getElementForces(sol, bcx, i, INZ);
         break;
       case OutputInfo::AXMoment:
         getElementForces(sol, bcx, i, AXM);
         break;
       case OutputInfo::AYMoment:
         getElementForces(sol, bcx, i, AYM);
         break;
       case OutputInfo::AZMoment:
         getElementForces(sol, bcx, i, AZM);
         break;
       case OutputInfo::Energies:
         // Wext = external energy
         Wext = force *  sol;
	 // Wela = elastic energy 
         Wela =   0.5 * Wext;

         fprintf(oinfo[i].filptr,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
           0.0,Wext,Wela,Wkin,Wdmp,Wext+Wela+Wkin+Wdmp);

         break;
       case OutputInfo::StrainVM:
	 getStressStrain(sol,bcx,i,STRAINVON);
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
       case OutputInfo::Composit:
         getCompositeData(i,time);
         break;
       case OutputInfo::Helmholtz:
         // ... PRINT (REAL) HELMHOLTZ SOLUTION
         for(iNode=0; iNode<numnodes; ++iNode) {
              int loc  = c_dsa->locate( iNode, DofSet::Helm);
              int loc1 =   dsa->locate( iNode, DofSet::Helm);

              double xHelm;
              if(loc >= 0)        // dof exists and is free
                xHelm = sol[loc];
              else if(loc1 >= 0)  // dof exists and is constrained
                xHelm = bcx[loc1];
              else                // dof does not exist
                xHelm = 0.0;

              fprintf(oinfo[i].filptr,"% *.*E\n",w,p,xHelm);
         }
         break;
       case OutputInfo::ShapeAtt:
         fprintf(oinfo[i].filptr,"  %f\n",time);
	 
         for(iNode=0; iNode<numnodes; ++iNode) {
	   oldnod=nodescopy->getNode(iNode);
           double x = nodes[iNode]->x - oldnod.x ;
           double y = nodes[iNode]->y - oldnod.y ;
           double z = nodes[iNode]->z - oldnod.z ;

           fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	 }
	 break;
       case OutputInfo::ShapeStc:
         fprintf(oinfo[i].filptr,"  %f\n",time);
	 
         for(iNode=0; iNode<numnodes; ++iNode) {
	   oldnod=nodescopy->getNode(iNode);
           double x = nodes[iNode]->x - oldnod.x ;
           double y = nodes[iNode]->y - oldnod.y ;
           double z = nodes[iNode]->z - oldnod.z ;

           fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	 }
	 break;
      }
    }
  } 

 sinfo.dt=oldtime;
}

//------------------------------------------------------------------------------

void
Domain::eigenOutput(Vector& eigenValues, VectorSet& eigenVectors,
                    double time)
{
  int i,iNode;
  const double pi = 3.141592653589793;

  enum {YOUNG,MDENS,THICK};

  Node oldnod;

  int iInfo;
  for(iInfo = 0; iInfo < numOutInfo; ++iInfo)
  {
    if(oinfo[iInfo].interval == 1) {
      int w = oinfo[iInfo].width;
      int p = oinfo[iInfo].precision;

      switch(oinfo[iInfo].type) {
      
      case OutputInfo::EigenPair:
      case OutputInfo::Disp6DOF:
 
        // ... Print eigenvalues and eigenvectors
        int imode,firstmode,maxmode;
      
        if (oinfo[iInfo].nodeNumber == -1 ) {
          firstmode=0;  
	  maxmode =sinfo.nEig; }
        else {
          firstmode=oinfo[iInfo].nodeNumber;
	  maxmode=firstmode+1; }
	  
        for(imode=firstmode; imode < maxmode; ++imode) {

           fprintf(oinfo[iInfo].filptr,"  %f\n",sqrt(eigenValues[imode])/(2.0*pi));

           for(i=0; i<numnodes; ++i) {

             // If you want to skip nodes that are not defined.
             if( dsa->firstdof(i) == -1 ) continue;

               int xloc = c_dsa->locate( i, DofSet::Xdisp);
               double x = (xloc >= 0) ? eigenVectors[imode][xloc] : 0;
 
               int yloc = c_dsa->locate( i, DofSet::Ydisp);
               double y = (yloc >= 0) ? eigenVectors[imode][yloc] : 0;

               int zloc = c_dsa->locate( i, DofSet::Zdisp);
               double z = (zloc >= 0) ? eigenVectors[imode][zloc] : 0;

               // ... considering the variation of shape
   	       oldnod=nodescopy->getNode(i);
               double vx = nodes[i]->x - oldnod.x ;
               double vy = nodes[i]->y - oldnod.y ;
               double vz = nodes[i]->z - oldnod.z ;
	   
	       x=x+vx;   y=y+vy;    z=z+vz;

               if(oinfo[iInfo].type == OutputInfo::Disp6DOF) {
                  xloc  = c_dsa->locate( i, DofSet::Xrot);
                  double xrot = (xloc >= 0) ? eigenVectors[imode][xloc] : 0;

                  yloc  = c_dsa->locate( i, DofSet::Yrot);
                  double yrot = (yloc >= 0) ? eigenVectors[imode][yloc] : 0;

                  zloc  = c_dsa->locate( i, DofSet::Zrot);
                  double zrot = (zloc >= 0) ? eigenVectors[imode][zloc] : 0;

                  fprintf(oinfo[iInfo].filptr,
                         "%d % *.*E\t% *.*E\t% *.*E\t% *.*E\t% *.*E\t% *.*E\n",
                         i+1,w,p,x,w,p,y,w,p,z,w,p,xrot,w,p,yrot,w,p,zrot);
               } else
                  fprintf(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",
                                         w,p,x,w,p,y,w,p,z);
             }
          }
          fflush(oinfo[iInfo].filptr);
          break;
        case OutputInfo::ShapeAtt:
          fprintf(oinfo[iInfo].filptr,"  %f\n",time);

          for(iNode=0; iNode<numnodes; ++iNode) {
	    oldnod=nodescopy->getNode(iNode);
            double x = nodes[iNode]->x - oldnod.x ;
            double y = nodes[iNode]->y - oldnod.y ;
            double z = nodes[iNode]->z - oldnod.z ;

            fprintf(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                                ,w,p,x,w,p,y,w,p,z);
	  }
	  break;
        case OutputInfo::ShapeStc:
          fprintf(oinfo[iInfo].filptr,"  %f\n",time);

          for(iNode=0; iNode<numnodes; ++iNode) {
	    oldnod=nodescopy->getNode(iNode);
            double x = nodes[iNode]->x - oldnod.x ;
            double y = nodes[iNode]->y - oldnod.y ;
            double z = nodes[iNode]->z - oldnod.z ;

            fprintf(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	  }
	  break;
        case OutputInfo::YModulus:
          getElementAttr(iInfo,YOUNG);
          break;
        case OutputInfo::MDensity:
          getElementAttr(iInfo,MDENS);
          break;
        case OutputInfo::Thicknes:
          getElementAttr(iInfo,THICK);
          break;
       }
     }  
   }

 // Print omega and frequency values to screen
 fprintf(stderr," Mode\tOmega^2\t\tFrequency\n");
 fprintf(stderr," --------------------------------------\n");
 int imode;
 for(imode=0; imode<sinfo.nEig; ++imode)
   fprintf(stderr," %d\t%e\t%e\n",imode+1,eigenValues[imode],
           sqrt(eigenValues[imode])/(2.0*pi));

}

//------------------------------------------------------------------------------

void
Domain::buildOptGrad() 
{

   // .... variable position of fe-nodes

   gradnodes = new CoordSet(numnodes);

   double *xyz = new double[3];
   
   xyz[0]=0.0; xyz[1]=0.0; xyz[2]=0.0; 
   
   int i;
   for (i=0;i<numnodes;i++) { gradnodes->nodeadd(i,xyz) ; }

   // .... variable material properties
   
   sgradprops = new StructProp[numprops];
   
   for (i=0;i<numprops;i++) {
      sgradprops[i].E         = 0.0; 
      sgradprops[i].A         = 0.0;
      sgradprops[i].nu        = 0.0;
      sgradprops[i].rho       = 0.0;
      sgradprops[i].eh        = 0.0;
      sgradprops[i].Ixx       = 0.0;
      sgradprops[i].c         = 0.0;
      sgradprops[i].alphaY    = 0.0;
      sgradprops[i].k         = 0.0;
      sgradprops[i].Q         = 0.0;
      sgradprops[i].W         = 0.0;
      sgradprops[i].P         = 0.0;
      sgradprops[i].Ta        = 0.0;
      sgradprops[i].kappaHelm = 0.0;
   }

   for (i=0; i<na; ++i) {
      Element *ele = eset[ attrib[i].nele ];
      if (attrib[i].attr < 0)
        ele->setProp(0);
      else
        ele->setGradProp(sgradprops+attrib[i].attr);
   }

   // .... variable composite properties

   for (i=0;i<numlayInfo;i++) {
     layInfo[i]->setGrad();
     layInfo[i]->zeroGrad();
   }

   for (i=0; i < na; ++i) {
      Element *ele = eset[ attrib[i].nele ];
      if (attrib[i].cmp_attr >= 0) {
        if(layInfo[attrib[i].cmp_attr] == 0) {
           continue;
	}
        ele->setCompositeGrad(layInfo[attrib[i].cmp_attr]->gradval());
      }
   }

   // Set up a derivative for each existing element frame

   EFrameData *cEFData = firstEFData;
   
   // Initialize pointer on list for derivative frames
   
   fstdEFData = 0;

   while(cEFData) {
     Element *ele = eset[cEFData->elnum];
     if(ele == 0) {
       fprintf(stderr,"\n *** WARNING: Frame was found for non existent"
                      " element %d\n", cEFData->elnum+1);
       cEFData = cEFData->next;
       continue;
     }

     // creating a linked list like for the EFrameData

     EFrameData *efd = new EFrameData;

     efd->next  = fstdEFData;
     fstdEFData = efd;
     efd->elnum = cEFData->elnum;
 
     // initializing with zero
  
     int i,j;
     for(i =0; i < 3; ++i)
       for(j =0; j < 3; ++j)
         efd->frame[i][j] = 0;
 
     // set up derivative in related element
  
     ele->setdFrame(&(efd->frame));
     
     cEFData = cEFData->next;
   }

}

//------------------------------------------------------------------------------

void
Domain::buildOptInf(int ivar, int numvar)
{
  if (ivar==0) optInf = new int*[numvar];
 
  actInf = new int[numele];
  
  int iele;
  for(iele=0; iele<numele; ++iele)
     actInf[iele]=packedEset[iele]->chkOptInf(*gradnodes);
     
  optInf[ivar]=actInf;   
}

//------------------------------------------------------------------------------

void
Domain::zeroGrad()
{
   Node gnode;
   
   int i;
   for (i=0;i<numnodes;i++) {
      (*gradnodes)[i]->x = 0.0 ; 
      (*gradnodes)[i]->y = 0.0 ;     
      (*gradnodes)[i]->z = 0.0 ;
   } 

   for (i=0;i<numprops;i++) {
      sgradprops[i].E         = 0.0; 
      sgradprops[i].A         = 0.0;
      sgradprops[i].nu        = 0.0;
      sgradprops[i].rho       = 0.0;
      sgradprops[i].eh        = 0.0;
      sgradprops[i].Ixx       = 0.0;
      sgradprops[i].Iyy       = 0.0;
      sgradprops[i].Izz       = 0.0;
      sgradprops[i].c         = 0.0;
      sgradprops[i].k         = 0.0;
      sgradprops[i].Q         = 0.0;
      sgradprops[i].W         = 0.0;
      sgradprops[i].P         = 0.0;
      sgradprops[i].Ta        = 0.0;
      sgradprops[i].kappaHelm = 0.0;
   }

   for (i=0;i<numlayInfo;i++) {
     layInfo[i]->zeroGrad();
   }

#ifdef AEROELASTIC
   flExchanger->zerogradFluidAttr();
#endif
}   

//------------------------------------------------------------------------------

void
Domain::computeTimePseudoLoad(Vector & transPseudo, int tIndex, double t,
                              Vector * gradAf)
{

// ... 1. add derivative of fluid pressure load

#ifdef AEROELASTIC
 if(sinfo.aeroFlag >= 0) { 

   double tFluid = flExchanger->getGradFluidLoad(nodes,transPseudo,tIndex,t);

   // ... save derivative of fluid load
   if ( gradAf ) (*gradAf)=transPseudo;
 }
#endif
}
