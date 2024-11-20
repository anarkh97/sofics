#include <cstdio>
#include <algorithm>
#include <Driver.d/DecDomain.h>
#include <Feti.d/Feti.h>
#include <Paral.d/MDDynam.h>
#include <Utils.d/BlockAlloc.h>
#include <Paral.d/DomainGroupTask.h>
#include <Paral.d/Assembler.h>
#include <Driver.d/DecDomainImpl.h>

template<>
double
GenDecDomain<double>::computeStabilityTimeStep(GenMDDynamMat<double>& dMat)
{
  using std::abs;
  double eigmax;
  double relTol    = domain->solInfo().stable_tol; // stable_tol default is 1.0e-3
  double preeigmax = 0.0;
  int numdofs = internalInfo->len;
  int maxIte  = domain->solInfo().stable_maxit; // stable_maxit default is 100
  GenDistrVector<double> v(*internalInfo);
  GenDistrVector<double> z(*internalInfo);
  // Starts from an arbitrary array.
  int i;
  for (i=0; i<numdofs; ++i)
    v.data()[i] = (double) (i+1) / (double) numdofs;

  Assembler *assembler = getSolVecAssembler();
  dMat.K->setAssembler(assembler); 
  dMat.M->setAssembler(assembler);
  // Power iteration loop
  for (i=0; i<maxIte; ++i) {

    dMat.K->mult(v,z);
    dMat.M->multInvertDiag(z);

    double zmax = z.infNorm();
    eigmax = zmax;
    v.linC(z,1.0/zmax);
    if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;
    preeigmax = eigmax;
  }
  // compute stability maximum time step
  double sdt = 2.0 / sqrt(eigmax);
  
  dMat.K->setAssembler((BasicAssembler *)0);
  dMat.M->setAssembler((BasicAssembler *)0);
  return domain->solInfo().stable_cfl*sdt;
}

template<>
void
GenDecDomain<double>::buildLocalFFP(int iSub, GenDistrVector<double> *u,
                                    double **ffp, int *numSample, double (*dir)[3], bool direction)
{
  fprintf(stderr, "WARNING: GenDecDomain<double>::buildLocalFFP not implemented \n");
}

template<>
void
GenDecDomain<DComplex>::buildLocalFFP(int iSub, GenDistrVector<DComplex> *u,
                                      DComplex **ffp, int *numSample, double (*dir)[3], bool direction)
{
 subDomain[iSub]->ffp(subDomain[iSub], *numSample, ffp[iSub], dir, u->subData(iSub), direction);
}

template<>
void
GenDecDomain<double>::buildFFP(GenDistrVector<double> &u, FILE *fffp, bool direction)
{
 fprintf(stderr, "WARNING: GenDecDomain<double>::buildFFP not implemented \n"); 
}

template<>
void
GenDecDomain<DComplex>::buildFFP(GenDistrVector<DComplex> &u, FILE *fffp, bool direction)
{
 if(direction && domain->numFFPDirections == 0) {
   int i,j;
   int nsint = std::max(2, domain->nffp);
   // Dimension of the problem
   int dim= domain->scatter[0]->dim();
   DComplex ffpCoef;
   if(dim!=3) ffpCoef = exp(DComplex(0.0,M_PI/4.0))/sqrt(8.0*M_PI*geoSource->kappa())*geoSource->global_average_rhof;
   else ffpCoef = DComplex(0.25/M_PI, 0.0)*geoSource->global_average_rhof;

   double (*ffpDir)[3];
   int numSamples;
   if(dim==2) {
     ffpDir = new double[nsint][3];
     for(i=0;i<nsint;i++) {
       ffpDir[i][0] = cos(2*M_PI*double(i)/double(nsint));
       ffpDir[i][1] = sin(2*M_PI*double(i)/double(nsint));
       ffpDir[i][2] = 0.0;
     }
     numSamples = nsint;
   } 
   else {
     int numTheta = nsint/2+1;
     ffpDir = new double[nsint*numTheta][3];
     for(i=0; i<numTheta; ++i) {
       double theta = M_PI*(-0.5+((double) i)/(numTheta-1.0));
       for(j=0;j<nsint;j++) {
         ffpDir[i*nsint+j][0] = cos(theta)*cos(2*j*M_PI/double(nsint));
         ffpDir[i*nsint+j][1] = cos(theta)*sin(2*j*M_PI/double(nsint));
         ffpDir[i*nsint+j][2] = sin(theta);
       }
     }
     numSamples = nsint*numTheta;
   }

   DComplex **localFFP = new DComplex * [numSub];
   for(i=0; i<numSub; ++i) {
     localFFP[i] = new DComplex[numSamples];
     for(j=0; j<numSamples; ++j) {
       localFFP[i][j] = DComplex(0.0,0.0);
     }
   }

   execParal(numSub, this, &GenDecDomain<DComplex>::buildLocalFFP, &u, localFFP, &numSamples, ffpDir, direction);

   for(i=1; i<numSub; ++i) 
     for(j=0; j<numSamples; ++j)
       localFFP[0][j] += localFFP[i][j];
   for(j=0; j<numSamples; ++j)
     localFFP[0][j] *= ffpCoef;
#ifdef DISTRIBUTED
   communicator->globalSum(numSamples, localFFP[0]);
   if(communicator->cpuNum() == 0) {
#endif
   if(dim==2) {
     for(j=0; j<nsint; ++j) {
       fprintf(fffp,"%e  %e  %.10e  %.10e\n",
               2*j*M_PI/double(nsint),0.0,
               ScalarTypes::Real(localFFP[0][j]), ScalarTypes::Imag(localFFP[0][j]));
     }
   }
   else {
     int numTheta = nsint/2+1;
     for(i=0; i<numTheta; ++i) {
       for(j=0; j<nsint; ++j) {
//         double y = 10.0 * log(2.0*M_PI*abs(localFFP[0][numSample])*
//                    abs(localFFP[0][numSample]))/log(10.0);
         fprintf(fffp,"%e  %e  %.10e  %.10e\n",
                 2*j*M_PI/double(nsint),M_PI*(-0.5+((double)i)/(numTheta-1.0)),
                 ScalarTypes::Real(localFFP[0][i*nsint+j]),
                 ScalarTypes::Imag(localFFP[0][i*nsint+j]));
       }
     }
   }
#ifdef DISTRIBUTED
   }
#endif
   for(i=0; i<numSub; ++i) delete [] localFFP[i];
   delete [] localFFP;
 } 
 else {
   // RT: new style input/output
   int i,j;
   int numSamples;
   double *evalDirLoc;

   if (direction) {
     numSamples = domain->numFFPDirections;
     evalDirLoc = domain->ffpDirections;
   } else {
     numSamples = domain->numKirchhoffLocations;
     evalDirLoc = domain->kirchhoffLocations;
   }

   DComplex **localFFP = new DComplex * [numSub];
   for(i=0; i<numSub; ++i) {
     localFFP[i] = new DComplex[numSamples];
     for(j=0; j<numSamples; ++j) {
       localFFP[i][j] = DComplex(0.0,0.0);
     }
   }

   double (*ffpDirLoc)[3] = new double[numSamples][3];
   for(i=0;i<numSamples;i++) {
     ffpDirLoc[i][0] = evalDirLoc[i*3+0];
     ffpDirLoc[i][1] = evalDirLoc[i*3+1];
     ffpDirLoc[i][2] = evalDirLoc[i*3+2];
   }

   execParal(numSub, this, &GenDecDomain<DComplex>::buildLocalFFP, &u, localFFP,
           &numSamples, ffpDirLoc, direction);

   for(i=1; i<numSub; ++i)
     for(j=0; j<numSamples; ++j)
       localFFP[0][j] += localFFP[i][j];

#ifdef DISTRIBUTED
   communicator->globalSum(numSamples, localFFP[0]);
   if(communicator->cpuNum() == 0) {
#endif

   for(i=0;i<numSamples;i++) {
     fprintf(fffp,"%e %e %e   %e %e\n",
             ffpDirLoc[i][0],
             ffpDirLoc[i][1],
             ffpDirLoc[i][2],
             real(localFFP[0][i]),imag(localFFP[0][i]));
   }
   fflush(fffp);
#ifdef DISTRIBUTED
   }
#endif

   for(i=0; i<numSub; ++i) delete [] localFFP[i];
   delete [] localFFP;
   delete [] ffpDirLoc;
 }
}

template class GenDecDomain<double>;
template class GenDecDomain<std::complex<double>>;