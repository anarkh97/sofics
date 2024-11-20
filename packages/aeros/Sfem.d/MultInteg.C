#include <Sfem.d/MultInteg.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/SolverInfo.h>
#include <iostream>

extern SolverInfo &solInfo;
extern GeoSource *geoSource;

template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor>
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::assignxw()
{
  // Gaussian quadrature points and weights
  rowx=9;
  colx=9;
  x=new double*[rowx];
  w=new double*[rowx];
  for (int i=0;i<rowx;i++) {
    x[i]=new double[colx];
    w[i]=new double[colx];
    for (int j=0;j<colx;j++) { 
      x[i][j]=0;
      w[i][j]=0;
    }
  }
  
  m=new int[rowx];
  m[0]=1;
  m[1]=2;
  m[2]=3;
  m[3]=4;
  m[4]=5;
  m[5]=6;
  m[6]=7;
  m[7]=8;
  m[8]=9;

  int i;
  i=0;
  x[i][0]=0;
  w[i][0]=1/0.56419;  

  i=1;
  x[i][0]=-1/sqrt(2.0);
  x[i][1]=1/sqrt(2.0);
  w[i][0]=0.88622693;
  w[i][1]=0.88622693;
 
  i=2;
  x[i][0]=-1.22474487;
  x[i][1]=0;
  x[i][2]=1.22474487;
  w[i][0]=0.29540898;
  w[i][1]=1.1816359;
  w[i][2]=0.29540898;

  i=3;
  x[i][0]=-1.65068012;
  x[i][1]=-0.52464762;
  x[i][2]=0.52464762;
  x[i][3]=1.65068012;
  w[i][0]=0.08131287;
  w[i][1]=0.80491409;
  w[i][2]=0.80491409;
  w[i][3]=0.08131287;

  i=4;
  x[i][0]=-2.02018287;
  x[i][1]=-0.95857246;
  x[i][2]=0;
  x[i][3]=0.95857246;
  x[i][4]=2.02018287;
  w[i][0]=0.01995324;
  w[i][1]=0.39361932;
  w[i][2]=0.94530872;
  w[i][3]=0.39361932;
  w[i][4]=0.01995324;
 
  i=5;
  x[i][0]=-2.35060497;
  x[i][1]=-1.33584907;
  x[i][2]=-0.43607741;
  x[i][3]=0.43607741;
  x[i][4]=1.33584907;
  x[i][5]=2.35060497;
  w[i][0]=0.00453;
  w[i][1]=0.15706732;
  w[i][2]=0.7246296;
  w[i][3]=0.7246296;
  w[i][4]=0.15706732;
  w[i][5]=0.00453;

  i=6;
  x[i][0]=-2.65196136;
  x[i][1]=-1.67355163;
  x[i][2]=-0.81628788;
  x[i][3]=0;
  x[i][4]=0.81628788;
  x[i][5]=1.67355163;
  x[i][6]=2.65196136;
  w[i][0]=0.0009718;
  w[i][1]=0.0545156;
  w[i][2]=0.4256073;
  w[i][3]=0.8102646;
  w[i][4]=0.4256073;
  w[i][5]=0.0545156;
  w[i][6]=0.0009718;

  i=7;
  x[i][0]=-2.93063742;
  x[i][1]=-1.98165676;
  x[i][2]=-1.15719371;
  x[i][3]=-0.381187;
  x[i][4]=0.381187;
  x[i][5]=1.15719371;
  x[i][6]=1.98165676;
  x[i][7]=2.93063742;
  w[i][0]=0.000199604;
  w[i][1]=0.01707798;
  w[i][2]=0.20780233;
  w[i][3]=0.66114701;
  w[i][4]=0.66114701;
  w[i][5]=0.20780233;
  w[i][6]=0.01707798;
  w[i][7]=0.000199604;

  i=8;
  x[i][0]=-3.1909932;
  x[i][1]=-2.26658058;
  x[i][2]=-1.46855329;
  x[i][3]=-0.72355102;
  x[i][4]=0;
  x[i][5]=0.72355102;
  x[i][6]=1.46855329;
  x[i][7]=2.26658058;
  x[i][8]=3.1909932;
  w[i][0]=0.00003961;
  w[i][1]=0.00494362;
  w[i][2]=0.08847453;
  w[i][3]=0.43265156;
  w[i][4]=0.72023522;
  w[i][5]=0.43265156;
  w[i][6]=0.08847453;
  w[i][7]=0.00494362;
  w[i][8]=0.00003961;


  for (int i=0;i<rowx;i++) {
    for (int j=0;j<colx;j++) {
      x[i][j]=x[i][j]*sqrt(2.0);
      w[i][j]=w[i][j]*0.56419;
    }
  }

}

template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::gensplit(int imd,int ttlp,int* grp,int cnt,int jtmp,int d1,double wt_smol)
{
// Function to split an integer n into d integers (permutation problem)
  cnt=cnt+1;
  int qq=jtmp-ttlp;
  ttlp=ttlp-1;

  for (int i=1;i<=qq;++i) {
    jtmp=jtmp-1;
    grp[cnt]=i;
    if (ttlp==0) {
     grp[cnt+1]=qq+1-i;
     // distribute d1 number of integers into a d-dimensional array
     int* targarray = new int[d];
     for (int ii=0;ii<d;ii++) targarray[ii]=0;
     int ilvl=1;  int icur=0;
     indist(grp,d1,targarray,icur,ilvl,wt_smol); // output is ivec
     delete [] targarray;
     targarray = 0;
    }
    else {
     gensplit(imd,ttlp,grp,cnt,jtmp,d1,wt_smol);
    }
  }
}


template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::indist(int* orgf,int d1f,int* trgf,int i_cur,int i_level,double wt_smol)
{
//Function for distribute d1 number of integers into a d-dimensional array, without disturbing the sequence

  for (int i=i_cur+1;i<=d-(d1f-i_level);i++) {
    trgf[i-1]=orgf[i_level-1];
    if (d1f==i_level) {
      double* tenx2 = new double[d]; // Apparently seems like defined twice, but actually this is local to this fn
      for (int ii=0;ii<d;ii++) tenx2[ii]=0;
      double* tenw2 = new double[d]; // Apparently seems like defined twice, but actually this is local to this fn
      for (int ii=0;ii<d;ii++) tenw2[ii]=0;
      //double wt_gauss=1;
      int count=-1; // changed from 0
      genten(count,tenx2,tenw2,trgf,wt_smol); // generates tensor product, evaluates function value, 
                                               // multiply by weight and update the integration result
      delete [] tenx2;
      tenx2 = 0;
      delete [] tenw2;
      tenw2 = 0;
    }
    else {
      indist(orgf,d1f,trgf,i,i_level+1,wt_smol);
    }
    trgf[i-1]=0;
  }
}


template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::genten(int cnt,double* tenxx,double* tenww,int* ivecc,double wt_extn)
{
// Function for generating tensor product and compute stresses at each point; needed for Smolyak and Kronecker
	cnt=cnt+1;
	int printflag = 1; // means don't print
	if (cnt<d) {
		for (int i=0;i<m[ivecc[cnt]];++i) {
			tenxx[cnt]=x[ivecc[cnt]][i];
			tenww[cnt]=w[ivecc[cnt]][i];
			genten(cnt,tenxx,tenww,ivecc,wt_extn);
			if (cnt==d-1) {
				for (int i=0;i<ndim;i++)  xi[i]=tenxx[i];
				build_psi();
				double prodw=1;
				for (int ii=0;ii<=cnt;++ii) prodw=prodw*tenww[ii];
				VecType* urealz = new VecType(probDesc_cur->solVecInfo(1));
				if(!isFeti(domain->solInfo().solvercntl->type))
					urealz->setn(domain->numUncon());
				urealz->zero();
				for(int ii=0;ii<P;ii++) {
					urealz->computeRealz(ii,psi[ii],(*mysol));
				}
				sfem->copyXi(xi); // xi in the argument is the local xi
//          if (stressIndex_cur < 7) sfem->assignRandMat();   //  only for stresses; YYY add principal stresses
				if (stressIndex_cur < 7) probDesc_cur->assignRandMat();   //  only for stresses; YYY add principal stresses
				postProcessor_cur->getStressStrain(urealz[0],fileNumber_cur,stressIndex_cur,time,printflag);
				integnodes++;
				res_cur = postProcessor_cur->getSfemStress(fileNumber_cur); // returning Scalar * :  current result
				if(ndflag_cur==1) { // mean
					for (int ii=0;ii<size_res;ii++) res[ii] = res[ii] + res_cur[ii]*prodw*wt_extn;
					delete urealz;
					urealz=0;
				}
				else if(ndflag_cur==2)
				{
					for (int ii=0;ii<size_res;ii++) res[ii] = res[ii] + pow(res_cur[ii],2)*prodw*wt_extn;
					for (int ii=0;ii<size_res;ii++) res[ii+size_res] = res[ii+size_res] + res_cur[ii]*prodw*wt_extn;
					delete urealz;
					urealz=0;
				}
				else std::cerr << "ndtype = " << ndflag_cur << " not supported by this routine" << std::endl;
			}
		}
	}
}


template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::computeStressStat(int qmd, VecType* sol, int fileNumber, int stressIndex, PostProcessor *postProcessor, ProblemDescriptor* probDesc, int ndflag)
{
  // The main function
  oinfo = geoSource->getOutputInfo();
  avgnum = oinfo[fileNumber].averageFlg;

  ndflag_cur=ndflag;
  fileNumber_cur=fileNumber;
  stressIndex_cur=stressIndex;
  postProcessor_cur = postProcessor;
  probDesc_cur = probDesc;
  postProcessor_cur->setsizeSfemStress(fileNumber); // Only for element-based output
  size_res = postProcessor_cur->getsizeSfemStress();

  makealpha();
  xi=new double[ndim]; // d=ndim
  psi=new double[P];
  mysol = sol;
  integnodes = 0;
  time = 0;

 if (ndflag_cur!=3) {
  if ((stressIndex_cur > 6) && (stressIndex_cur < 13))  directcomp();  // strains sxx syy szz sxy syz sxz
//  else  simulcomp();
//  else kroneckercomp();
//  else integmain(qmd);
  else  {
   int meth_int, param_int;
/*   ifstream readintegfile("integparamfile",ios::in);
   readintegfile >> meth_int; // 1 : Simulation, 2 : Kronecker, 3 : Smolyak
   readintegfile >> param_int; // Integration parameter, nosamp for simul, points in 1-D for kronecker, qmd for Smolyak
   readintegfile.close();*/
   meth_int=1;
   param_int=20;
   if(meth_int==1) simulcomp(param_int);
   else if(meth_int==2) kroneckercomp(param_int);
   else if(meth_int==3) integsmol(param_int);
   else std::cerr << "Invalid Integration method" << std::endl;
  }

 }
 else compdf();
// std::cerr << "Number of integration nodes = " << integnodes << std::endl;
// delete [] x; // YYY DG, should delete, but wrong syntax
// delete [] w;
// delete [] m;
}


template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::integsmol(int qmd)
{
  filePrint(stderr,"Integration using Smolyak cubature\n");

  int printflag = 1; // means don't print

  q=d+qmd;
  int lowlim;
  if (d>=q-d+1) lowlim = d;
  else lowlim = q-d+1;
/*  std::cerr << "d = " << d << std::endl;
  std::cerr << "qmd = " << qmd << std::endl;
  std::cerr << "lowlim = " << lowlim << std::endl;*/
  int uplim=q;

  if(ndflag_cur==1) //  mean
   {
    res = new Scalar[size_res];  
    for (int i=0;i<size_res;i++)  res[i]=0;
   }
  else if(ndflag_cur==2)
   {
    res = new Scalar[2*size_res];  
    for (int i=0;i<2*size_res;i++)  res[i]=0;
   }
  else
   std::cerr << "ndtype = " << ndflag_cur << " not supported currently" << std::endl;

  double weight;
  int d1, ilim, ilvl, icur,count,jtemp;
  int* grp1;
  int* targarray = new int[d];
  int* trgf = new int[d];
  double* tenx = new double[d];
  double* tenw = new double[d];

  for (int modi=lowlim; modi<=uplim;modi++) { // modi=|i|
    int  imd=modi-d;
    weight=pow((-1.0),(q-modi))*nchooser(d-1,q-modi);
    if(imd==0) { // (x1,x1, ... ,x1)
      for (int ii=0;ii<d;ii++) trgf[ii]=0; // replacement of indist()
      for (int ii=0;ii<d;ii++) tenx[ii]=0;
      count=-1; 
      for (int ii=0;ii<d;ii++) tenw[ii]=0; // tensor product of Gaussian quadrature weights
     genten(count,tenx,tenw,trgf,weight);
    }
    else {
      weight=pow((-1.0),(q-modi))*nchooser(d-1,q-modi);
      if(imd > d) {
	ilim=d;
      }
      else {
        ilim=imd;
      }
      for (int i=1;i<=ilim;++i) { // the min() is introduced to handle the cases q >= 2*d
        if(i==1) { 
          //if(!grp1) int* grp1;
          grp1=new int[i];
          grp1[0]=imd;
          d1=i;  
          for (int ii=0;ii<d;ii++) targarray[ii]=0;
          ilvl=1;  icur=0;
          indist(grp1,d1,targarray,icur,ilvl,weight); 
	  //          if(grp1) delete grp1;
          delete [] grp1;
          grp1=0;
        }
        else {
	  int totloop=i-2; // number of for loops
          for (int j=1;j<=imd-i+1;++j) {
            //if(!grp1) int* grp1;
            grp1=new int[i];
            for (int ii=0;ii<i;ii++) grp1[ii]=0;
            grp1[0]=j;
	    if(i==2) {
              grp1[1]=imd-j;
              // distribute d1 number of integers into a d-dimensional array
	      d1=i; 
              for (int ii=0;ii<d;ii++) targarray[ii]=0;
              ilvl=1;  icur=0;
              indist(grp1,d1,targarray,icur,ilvl,weight); 
            }
            else {
              count=0;
              jtemp=imd-j;
              gensplit(imd,totloop,grp1,count,jtemp,i,weight); // here d1=i, i.e. grp1 has 'i' elements
	      // Note : here indist(..) is called within gensplit(..)
            }
            delete [] grp1;
            grp1=0;
          }
        }
      }
    }
  }
  
  // compute std-dev
//  if(ndflag_cur==2)   for (int ii=0;ii<size_res;ii++) (*res)[ii]=ScalarTypes::sqrt((*res)[ii]-pow((*res)[size_res+ii],2));
  if(ndflag_cur==2)   for (int ii=0;ii<size_res;ii++) res[ii]=ScalarTypes::sqrt(std::abs(res[ii]-pow(res[size_res+ii],2)));
                                          //  abs() is used to avoid numerical error near zero-mean material
  postProcessor_cur->updateSfemStress(res,fileNumber_cur); // using postProcessor_cur
  printflag = 2; // means don't compute, only print
  postProcessor_cur->getStressStrain(mysol[0],fileNumber_cur,stressIndex_cur,time,printflag); // print result 

}

template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::simulcomp(int nsample_integ)
{
  //filePrint(stderr,"Integration by simulation\n");
  int printflag = 1; // means don't print
//  int nsample_integ = 5; //sfem->getnsamp_out(); // YYY take this data from oinfo, at least for stress and strain
  if(ndflag_cur==1) //  mean
   {
    res = new Scalar[size_res]; 
    for (int i=0;i<size_res;i++) res[i]=0;
   }
  else if(ndflag_cur==2)
   {
    res = new Scalar[2*size_res];
    for (int i=0;i<2*size_res;i++)  res[i]=0;
   }
  else {
   std::cerr << "ndtype = " << ndflag_cur << " not supported currently" << std::endl;
   return;
  }

  VecType* urealz = new VecType(probDesc_cur->solVecInfo(1));
  if(!isFeti(domain->solInfo().solvercntl->type)) urealz->setn(domain->numUncon());

  for(int i=0; i< nsample_integ; ++i) {
    genXiPsi(i);
    sfem->copyXi(xi);
//    if (stressIndex_cur < 7)  sfem->assignRandMat();   // only for stresses; YYY add principal stresses
    if (stressIndex_cur < 7)  probDesc_cur->assignRandMat();   // only for stresses; YYY add principal stresses
    urealz->zero();
    for(int ii=0;ii<P;ii++) {
      urealz->computeRealz(ii,psi[ii],(*mysol));
    }
    postProcessor_cur->getStressStrain(urealz[0],fileNumber_cur,stressIndex_cur,time,printflag);
//    std::cerr.flush();

    res_cur = postProcessor_cur->getSfemStress(fileNumber_cur);
     if(ndflag_cur==1) { // mean
       for (int ii=0;ii<size_res;ii++) res[ii]=res[ii]+res_cur[ii];
     }
     else if(ndflag_cur==2)
      {
       for (int ii=0;ii<size_res;ii++) res[ii]=res[ii]+pow(res_cur[ii],2);
       for (int ii=0;ii<size_res;ii++) res[ii+size_res]=res[ii+size_res] + res_cur[ii];
      }
     else std::cerr << "ndtype = " << ndflag_cur << " not supported by this routine" << std::endl;

  } // end of nsample_output loop

  delete urealz;
  urealz=0;
   

  if(ndflag_cur==1)  for (int i=0;i<size_res;i++)  res[i] = res[i]/Scalar(nsample_integ);
  else if(ndflag_cur==2)  for (int i=0;i<2*size_res;i++)  res[i] = res[i]/Scalar(nsample_integ);
  else std::cerr << "ndtype = " << ndflag_cur << " not supported currently" << std::endl;

  // compute std-dev
  if(ndflag_cur==2)   for (int ii=0;ii<size_res;ii++) res[ii] = ScalarTypes::sqrt(std::abs(res[ii]-pow(res[size_res+ii],2)));

// YYY delete res_cur
  postProcessor_cur->updateSfemStress(res,fileNumber_cur); // using postProcessor_cur
  printflag = 2; // means don't compute, only print
  postProcessor_cur->getStressStrain(mysol[0],fileNumber_cur,stressIndex_cur,time,printflag); // print result
  //delete []res; 
}


template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::kroneckercomp(int quadorder)
{
  filePrint(stderr,"Integration on the Kronecker tensor product grid\n");
  int printflag = 1; // means don't print
  if(ndflag_cur==1) //  mean
   {
    res = new Scalar[size_res]; 
    for (int i=0;i<size_res;i++)  res[i]=0;
   }
  else if(ndflag_cur==2)
   {
    res = new Scalar[2*size_res]; 
    for (int i=0;i<2*size_res;i++)  res[i]=0;
   }
  else
   std::cerr << "ndtype = " << ndflag_cur << " not supported currently" << std::endl;

  int* trgf = new int[d];
  for (int i=0;i<d;i++) trgf[i]=quadorder; // Define the 1-D rule
  double* tenx2 = new double[d]; // Apparently seems like defined twice, but actually this is local to this fn
  for (int ii=0;ii<d;ii++) tenx2[ii]=0;
  double* tenw2 = new double[d]; // Apparently seems like defined twice, but actually this is local to this fn
  for (int ii=0;ii<d;ii++) tenw2[ii]=0;
  //double wt_gauss=1;
  int count=-1; 
  double wt_sm=1;
  genten(count,tenx2,tenw2,trgf,wt_sm); // generates tensor product, evaluates function value,
                                        // multiply by weight and update the integration result
  delete [] tenx2;
  tenx2 = 0;
  delete [] tenw2;
  tenw2 = 0;
  delete [] trgf;
  trgf = 0;

  // compute std-dev
  if(ndflag_cur==2)   for (int ii=0;ii<size_res;ii++) res[ii]=ScalarTypes::sqrt(std::abs(res[ii]-pow(res[size_res+ii],2)));
// YYY delete res_cur
  postProcessor_cur->updateSfemStress(res,fileNumber_cur); // using postProcessor_cur
  printflag = 2;
  postProcessor_cur->getStressStrain(mysol[0],fileNumber_cur,stressIndex_cur,time,printflag); // print result
}


template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::directcomp()
{
  filePrint(stderr,"Using orthogonality of the chaos polynomials\n");
  res = new Scalar[size_res]; 
  for (int i=0;i<size_res;i++)  res[i]=0;
  int printflag = 0; // means compute and print, the usual one

  if (ndflag_cur==1) { // mean
   postProcessor_cur->getStressStrain(mysol->getBlock(0),fileNumber_cur,stressIndex_cur,time,printflag); // printflag is default 0
  }
  else if(ndflag_cur==2) {
   printflag = 1;
   psisqr=new double[P];
   build_psisqr();
   for (int i=1;i<P;i++) {
     postProcessor_cur->getStressStrain(mysol->getBlock(i),fileNumber_cur,stressIndex_cur,time,printflag);
//     postProcessor_cur->getStressStrain(mysol->getBlock(i),fileNumber_cur,stressIndex_cur,time,1);
     res_cur = postProcessor_cur->getSfemStress(fileNumber_cur);
     for (int j=0;j<size_res;j++) res[j]=res[j]+psisqr[i]*pow(res_cur[j],2);
   }
   for (int i=0;i<size_res;i++) res[i]=ScalarTypes::sqrt(res[i]);


// YYY delete res_cur
   postProcessor_cur->updateSfemStress(res,fileNumber_cur); // using postProcessor_cur
   printflag = 2; // means don't compute, only print

   postProcessor_cur->getStressStrain(mysol[0],fileNumber_cur,stressIndex_cur,time,printflag); // print result 
  }
  else std::cerr << "ndtype = " << ndflag_cur << " not supported currently" << std::endl;

}

template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
void MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::compdf()
{
  int printflag = 1;
  int nsample_output = sfem->getnsamp_out(); // YYY take this data from oinfo, at least for stress and strain

  Scalar* strchaos;
  if ((stressIndex_cur > 6) && (stressIndex_cur < 13)) {  // strains sxx syy szz sxy syz sxz
//   compute_and_save_coeficients
   strchaos = new Scalar(P*size_res); // chaos expansion of stress tensor components : get size_res through postProcessor_cur
   for (int i=0;i<P;i++) {
     postProcessor_cur->getStressStrain(mysol->getBlock(i),fileNumber_cur,stressIndex_cur,time,printflag); 
     res_cur = postProcessor_cur->getSfemStress(fileNumber_cur);
     for (int j=0;j<size_res;j++) strchaos[i*size_res+j] = res_cur[j];
   }
   delete mysol;
   mysol=0;
  }

  res = new Scalar[size_res]; // YYY check assignment

  for(int i=0; i< nsample_output; ++i) {
    for (int ii=0;ii<size_res;ii++)  res[ii]=0;
    genXiPsi(i);

    if ((stressIndex_cur > 6) && (stressIndex_cur < 13))  { // strains sxx syy szz sxy syz sxz  
     for (int k=0;k<P;k++)  for (int j=0;j<size_res;j++) res[j]=res[j]+psi[k]*strchaos[size_res*k+j]; // directly_from_expansion
    }
    else {
     sfem->copyXi(xi);
//     if (stressIndex_cur < 7)  sfem->assignRandMat();   // only for stresses; YYY add principal stresses
     if (stressIndex_cur < 7)  probDesc_cur->assignRandMat();   // only for stresses; YYY add principal stresses
     VecType* urealz = new VecType(probDesc_cur->solVecInfo(1));
     if(!isFeti(domain->solInfo().solvercntl->type)) urealz->setn(domain->numUncon());
     urealz->zero();
     for(int ii=0;ii<P;ii++) {
       urealz->computeRealz(ii,psi[ii],(*mysol));
     }
     postProcessor_cur->getStressStrain(urealz[0],fileNumber_cur,stressIndex_cur,time,printflag);
     res_cur = postProcessor_cur->getSfemStress(fileNumber_cur);
     delete urealz;
     urealz=0;
    } // end of if-else
    postProcessor_cur->updateSfemStress(res,fileNumber_cur); // YYY DG special care needed for printing pdf
    printflag = 2; // means don't compute, only print
    postProcessor_cur->getStressStrain(mysol[0],fileNumber_cur,stressIndex_cur,time,printflag); // YYY DG tricky for pdf, make sure it overwrites
    printflag=1;
  } // end of nsample_output loop

// YYY delete res_cur
}



template <class Scalar, class VecType, class PostProcessor, class ProblemDescriptor >
int MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>::nchooser(int n, int r) 
{
  int p;
  int q;
  if (r<n-r) {
    p=r;
    q=n-r;
  }
  else {
    p=n-r;
    q=r;
  }

  int num=1;
  for (int i=q+1;i<=n;i++)  num=num*i;
  int den = 1;
  for (int i=2; i<=p;i++)  den = den* i;
  int nd =num/den;
  return nd;
}

