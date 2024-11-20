#ifndef _SFEM_H
#define _SFEM_H

#include <iostream>
#include <fstream>

class Sfem  
{
protected:
  int P;
  int L;
  int ndim;
  int output_order;
  int* nonzindex; // input non-zero block index
  int* nnzblkindex; // output non-zero block index
//  int nosamp_deletelater; 
private:
  int nsample_output;
public:
  Sfem() { ndim = 0; nsample_output =0; alphamat=0; Gauss = false; psisqr = 0; isreduced = false;};
  ~Sfem() {};
  int getL(){ return L;}
  int getP(){return P;}
  int no_terms(int x, int y);	
  void setOrder(int output_order);
  void computeLP();
  void updatendim() {ndim++ ; }
  int getndim() { return ndim; }
  int getoutput_order() { return output_order; }
  void setnsamp_out(int nso) {nsample_output =  nso; } // need to put a switch to print a warning message if called more than once
  int getnsamp_out() {return nsample_output;}
  void init_XiPsi();
  void genXi(int seed);
  void genXiPsi(int seed);
  void build_psi();
  void copyXi(double *xisource) {xi=xisource;}
  double *psisqr;
  void makealpha();
  void build_psisqr();
  double* getxi() {return xi;}
  double* getpsi() {return psi;}
  double* getPsiSqr() {return psisqr;}
  std::ifstream readfile;
  bool Gauss;              
  bool isreduced;
  int* getnonzindex() {return nonzindex;}
  void computeNnzBlocks(double* bln);
  int* getNnzBlocks() {return nnzblkindex;}
  void setreduced() {isreduced = true;}
protected:
  int **alphamat;
  double *xi;
  double *psi;
  int zufall_(int nn, double *a);
  int zufalli_(int seed_local);
  int normalen_(int nn, double *x);
  int normal00_();
  int nchooser(int n, int r);
};


#endif
