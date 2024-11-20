#ifndef _MULTINTEG_H
#define _MULTINTEG_H

#include <Sfem.d/Sfem.h>

class OutputInfo;

extern Sfem *sfem;

template <class Scalar, 
          class VecType, 
          class PostProcessor,
          class ProblemDescriptor>

class MultInteg : public Sfem {
  double** x;
  double** w;
  int* m;
  int rowx;
  int colx;
  int size_res;
  int d; // dimension
  int q;
  int integnodes;
  int fileNumber_cur, stressIndex_cur;
  PostProcessor *postProcessor_cur;
  ProblemDescriptor* probDesc_cur;
  VecType* mysol;
  int ndflag_cur;
  OutputInfo *oinfo;
  int avgnum;
  Scalar* res;
  Scalar* res_cur;
  double time;
 public:
  MultInteg() {assignxw(); P=sfem->getP(); L= sfem->getL(); ndim=sfem->getndim(); d=ndim; output_order = sfem->getoutput_order();}
  ~MultInteg() {};//  delete res; res = 0; };
//  void setq(int qmd) {q=d+qmd;}
  int getd() {return d;}
  int getq() {return q;}
  void assignxw();
  void gensplit(int imd, int ttlp, int* grp, int cnt, int jtmp, int d1, double wt_smol);
  void indist(int* orgf, int d1f, int* trgf, int i_cur, int i_level, double wt_smol);
  void genten(int cnt, double* tenxx, double* tenww, int* ivecc, double wt_extn);
  Scalar* getres() {return res;}
  int nchooser(int n, int r);
  void computeStressStat(int qmd, VecType* sol, int fileNumber, int stressIndex, PostProcessor *postProcessor,
                         ProblemDescriptor* probDesc, int ndflag);
  void integsmol(int qmd); // Integration using Smolyak cubature
  void simulcomp(int nsample_integ);  // Integration by simulation
  void kroneckercomp(int quadorder); // Integration using Kronecker (i.e. full) tensor product rule
  void directcomp();
  void compdf();
};

#ifdef _TEMPLATE_FIX_
#include <Sfem.d/MultInteg.C>
#endif

#endif
