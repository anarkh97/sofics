#include <cstdio>
#include <Math.d/Skyline.d/SkyMatrix.C>
#include <Utils.d/SolverInfo.h>

extern SolverInfo &solInfo;

template<>
void
GenSkyMatrix<double>::Factor(Rbm *rigid)
{
   //filePrint(stderr," ... Skyline factor: neq %d size %d avg. band. %d \n",
   //         neqs(), size(), size() / neqs());

   int i,j;
   double *w = new double[7*numUncon];

   double *dummyZEM = w + numUncon;
   double *w1 = dummyZEM + numUncon;
   double *w2 = w1 + numUncon;

   double *rhs = 0;

   // ... INITIALIZE NUMBER OF OPERATIONS (NOPS)
   double nops = 0.0;

   // ... INITIALIZE # of GEOMETRIC RIGID BODY MODES
   int ngrbm = rigid->numRBM();

   // ... DECLARE A POINTER TO ZEM
   double *zem = NULL;

   // ... GET THE NUMBER OF COMPONENTS
   int numComp = rigid->numComponents();

   // ... INITIALIZE # OF ZERO ENERGY MODES (ZEM)
   int nTotZem = 0;
   int *nZemPerComp = (int *) dbg_alloca(numComp * sizeof(int));

   // ... ALLOCATE MEMORY FOR S1 AND S2
   // int sizeS = (numUncon - kstop)*(numUncon - kstop);
   int defblk = solInfo.solvercntl->sparse_defblk;
   double *s1 = (double *) dbg_alloca((defblk+4)*(defblk+4)*sizeof(double));
   double *s2 = (double *) dbg_alloca((defblk+4)*(defblk+4)*sizeof(double));

   // ... LOOP OVER COMPONENTS
   int n;
   for(n=0; n<numComp; ++n) {

     // ... CALL WITH FLAG = 1 FOR PARTIAL FACTORIZATION
     int flag = 1;

     // Get number of rbm per component
     int locNgrbm = rigid->numRBM(n);

     // Get number of dof per component
     int numDofPerComp = (numComp > 1) ? rigid->numDof(n) : numUncon;

     // Get first dof of the component
     int firstDofOfComp = rigid->firstDof(n);

     int kstop = 4*((numDofPerComp - defblk)/4);
     if(kstop < 0) kstop = 0;

     int k;
     if(firstDofOfComp != 0)
       for(k = 0; k < numDofPerComp; ++k)
          lacol[firstDofOfComp+k] -= firstDofOfComp;

     _FORTRAN(svbu4gmb)(skyA,dlp+firstDofOfComp,lacol+firstDofOfComp,rhs,
      w+firstDofOfComp,w1+firstDofOfComp,w2+firstDofOfComp,
      pivot+firstDofOfComp,TOLERANCE,numDofPerComp,flag,nops,nzem,
      dummyZEM,kstop,NULL,NULL,locNgrbm,seqid+firstDofOfComp);

     nTotZem += nzem;
     nZemPerComp[n] = nzem;
   }

   Vector *allrbms;
   if(nTotZem > 0) {
     allrbms = new Vector[nTotZem+ngrbm];
     Vector v(numUncon,0.0);
     int i;
     for(i=0; i<ngrbm; ++i)
        allrbms[i] = rbm->getGrbm(i);

     for(i=0; i<nzem; ++i)
       allrbms[i+ngrbm] = v;
   }

   nTotZem = 0;
   for(n=0; n<numComp; ++n) {
     // ... CALL WITH FLAG = 11 TO COMPLETE FACTORIZATION
     int flag = 11;

     // Get number of rbm per component
     int locNgrbm = rigid->numRBM(n);

     // Get number of dof per component
     int numDofPerComp = (numComp > 1) ? rigid->numDof(n) : numUncon;

     // Get first dof of the component
     int firstDofOfComp = rigid->firstDof(n);

     int kstop = 4*((numDofPerComp - defblk)/4);
     if(kstop < 0) kstop = 0;

     // ... UPDATE nzem
     int ZemPlusGrbm = nZemPerComp[n] + locNgrbm;

     // ... ALLOCATE MEMORY FOR nzem RIGID BODY MODES (even if we don't care)
     zem = (ZemPlusGrbm) ? (double *) dbg_alloca(numDofPerComp*ZemPlusGrbm*sizeof(double)) : NULL;

    _FORTRAN(svbu4gmb)(skyA,dlp+firstDofOfComp,lacol+firstDofOfComp,
        rhs,w+firstDofOfComp, w1+firstDofOfComp,w2+firstDofOfComp,
        pivot+firstDofOfComp,TOLERANCE,numDofPerComp,
                    flag,nops,ZemPlusGrbm,zem,kstop,s1,s2,locNgrbm,
                    seqid+firstDofOfComp);

     // Copy mechanisms into list of rbms
     for(i = 0; i < nZemPerComp[n]; ++i) {
       for(j = 0 ; j < numDofPerComp; ++j)
         allrbms[i+ngrbm+nTotZem][j+firstDofOfComp] = zem[i*numDofPerComp+j];
     }
    nTotZem += nZemPerComp[n];
    int k;
    if(firstDofOfComp != 0)
       for(k = 0; k < numDofPerComp; ++k)
          lacol[firstDofOfComp+k] += firstDofOfComp;
   }

// ... NOW COPY MECHANISM MODES TO GEOMETRIC RIGID BODY MODES IF NECESSARY
   nzem = nTotZem + ngrbm;
   if(print_nullity && nzem > 0)
     std::cerr << " ... Matrix is singular: size = " << numUncon << ", rank = " << numUncon-nzem << ", nullity = " << nzem
               << " (" << ngrbm << " grbm/hzem + " << nzem-ngrbm << " other) ...\n";

   if(nTotZem > 0) {
     rbm = new Rbm(allrbms,nzem,numUncon);
     myRbm = 1;
   }
   delete [] w;
}

template<>
void
GenSkyMatrix<DComplex>::Factor(Rbm *rigid)
{
 // GRBM factor not implemented for complex, so just use TRBM
 Factor();
}

template<>
double
GenSkyMatrix<double>::diag(int dof) const
{
  if(skyA[dlp[dof] - 1] <= 0)  // <= because taking square root elsewhere
    return (1.0);
  else
    return skyA[dlp[dof] - 1];
}

template<>
DComplex
GenSkyMatrix<DComplex>::diag(int dof) const
{
  return skyA[dlp[dof] - 1];
}

template<>
void
GenSkyMatrix<DComplex>::addImaginary(const FullSquareMatrix &ks, const int *dofs)
{
// Construct stiffness matrix K (skyA)
 int i, j, k;
 int kndof = ks.dim();                 // Element stiffness dimension
 for( i = 0; i < kndof; ++i ) {          // Loop over rows.
    if( rowColNum[dofs[i]] == -1 ) continue;    // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {              // Loop over columns.
       if( dofs[i] > dofs[j] ) continue;      // Work with upper symmetric half.
       if( rowColNum[dofs[j]] == -1 ) continue; // Skip constrained dofs
       k = dlp[rowColNum[dofs[j]]] - (rowColNum[dofs[j]]-rowColNum[dofs[i]])-1;
       skyA[k] += DComplex(0.0,ks[i][j]);
    }
 }
}

template<>
void
GenSkyMatrix<double>::addImaginary(const FullSquareMatrix &ks, const int *dofs)
{
  fprintf(stderr, "GenSkyMatrix<double>::addImaginary(...) is not implemented \n");
}

template<>
void
GenSkyMatrix<DComplex>::add(const FullSquareMatrixC &kel, const int *dofs)
{
// Construct stiffness matrix K (skyA)
 int i, j, ri, rj;
 int kndof = kel.dim();                	// Element stiffness dimension

 for( i = 0; i < kndof; ++i ) {           // Loop over rows.
    if( (ri = rowColNum[dofs[i]]) == -1 ) continue; // Skip constrained dofs
    for( j = 0; j < kndof; ++j ) {          // Loop over columns.
       if( dofs[i] > dofs[j] ) continue;    // Work with upper symmetric half.
       if( (rj = rowColNum[dofs[j]]) == -1 ) continue; // Skip constrained dofs
       skyA[dlp[rj] - rj + ri - 1] += kel[i][j];
    }
 }
}

template<>
void
GenSkyMatrix<double>::add(const FullSquareMatrixC &kel, const int *dofs)
{
  fprintf(stderr, "GenSkyMatrix<double>::add(const FullSquareMatrixC &kel, const int *dofs) is not implemented \n");
}

template<>
void
GenSkyMatrix<double>::print(FILE *f)
{
 int i;
 fprintf(f,"K(%d,%d) = %e;\n",1,1,skyA[0]);
 for(i=1; i<numUncon; ++i) {
   int numEntries = dlp[i] - dlp[i-1];
   int j;
   for(j=0; j<numEntries; ++j)
     fprintf(f,"K(%d,%d) = %e;\n",i-j+1,i+1,skyA[dlp[i]-1-j]);
 }
 fprintf(f,"Kt = transpose(K);\n");
 fprintf(f,"S = K+Kt;\n");
 fprintf(f,"D = diag(diag(K));\n");
 fprintf(f,"S = S - D;\n");
}

template<>
void
GenSkyMatrix<DComplex>::print(FILE *fid)
{
 int iCol;
 int iDof = 0;

 for(iCol=0; iCol<numUncon; ++iCol) {
   int iRow = iCol+iDof-dlp[iCol];
   for(; iDof < dlp[iCol]; ++iDof, iRow++) {
     fprintf(fid, "K(%d,%d) = %e ",iCol+1, iRow+2,real(skyA[iDof]));
     if (imag(skyA[iDof]) < 0.0)
        fprintf(fid," %e*i\n", imag(skyA[iDof]));
     else
        fprintf(fid,"+ %e*i\n", imag(skyA[iDof]));
   }
 }
}

// print the skyline matrix diagonal entries
template<>
void
GenSkyMatrix<double>::printDiagonals()
{
 fprintf(stderr,"------ Skyline Diagonals -----\n");
 int dof;
 for(dof=0; dof<numUncon; ++dof)
   fprintf(stderr,"K(%d,%d) = %e\n", dof,dof,skyA[dlp[dof] - 1]);
}

template<>
void
GenSkyMatrix<DComplex>::printDiagonals()
{
  fprintf(stderr, "GenSkyMatrix<DComplex>::printDiagonals() is not implemented \n");
}

template<>
void
GenSkyMatrix<double>::print1(int dof)
{
 int i;
 for(i = dlp[dof-1]; i < dlp[dof]; ++i)
   fprintf(stderr,"%d  %d  %e\n",dof+i-dlp[dof]+1, dof, skyA[i]);
}

template<>
void
GenSkyMatrix<DComplex>::print1(int dof)
{
  fprintf(stderr, "GenSkyMatrix<DComplex>::print1(int dof) is not implemented \n");
}

template<>
void
GenSkyMatrix<double>::printMatlab(int subNum)
{
  int i;
  char checkfile1[80];

  // File names
  sprintf(checkfile1, "K%d.FEM", subNum+1);

  // Open files
  FILE *fileReal;
  fileReal=fopen(checkfile1,"w");

  // DLP for the real part
  fprintf(fileReal,"%d \n", numUncon);
  for(i = 0; i < numUncon; i++)
    fprintf(fileReal," %d \n", dlp[i]);

  // Real part of the skyline
  fprintf (fileReal,"%d \n", dlp[numUncon-1]);
  for(i = 0; i < dlp[numUncon-1]; ++i)
    fprintf(fileReal, " %1.16e \n", skyA[i]);

  fclose(fileReal);
}

template<>
void
GenSkyMatrix<DComplex>::printMatlab(int subNum)
{
  fprintf(stderr, "GenSkyMatrix<DComplex>::printMatlab(int subNum) is not implemented \n");
}

template<>
void
GenSkyMatrix<double>::printMatlab(char *fileName)
{
 // Open files
 FILE *file;
 if(fileName)
   file=fopen(fileName,"w");
 else
   file=stderr;

 int i;
 /*
 for(i = 0; i < dlp[numUncon-1]; ++i)
    fprintf(file, " %1.16e \n", skyA[i]);
 */

 fprintf(file,"%% row col val\n");
 fprintf(file," %6d %6d %1.16e\n",1,1,skyA[0]);
 for(i=1; i<numUncon; ++i) {
   int numEntries = dlp[i] - dlp[i-1];
   int j;
   for(j=0; j<numEntries; ++j)
     fprintf(file," %6d %6d %1.16e\n",i-j+1,i+1,skyA[dlp[i]-1-j]);
 }
 if(fileName)
   fclose(file);
}

template<>
void
GenSkyMatrix<DComplex>::printMatlab(char *fileName)
{
  fprintf(stderr, "GenSkyMatrix<DComplex>::printMatlab(char *fileName) is not implemented \n");
}



template class GenSkyMatrix<double>;
template class GenSkyMatrix<complex<double>>;

