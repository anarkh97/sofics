#ifndef _DECOMP_H_
#define _DECOMP_H_

#include <Utils.d/Connectivity.h>

class Decomposition {
  public:
     char * esname;  // Name of the elementset for which this is a decomposition
     int nsub;       // Number of subdomains
     int * pele;     // pointer into eln
     int * eln;      // List of elements in each subdomain

     enum Alg { Greedy , Experimental } ;

     Decomposition() {}
     Decomposition(char *ename, int ns, Alg);
     Decomposition(int ns, int *p, int *l) { esname = (char *) "noname";
         nsub = ns; pele = p; eln = l; }
     ~Decomposition() { delete [] pele; delete [] eln; }

     void setName(char *name) { esname = name; }

     void greedyBuild(int, Connectivity *etoe, Connectivity *ntoe, 
                      Connectivity*eton);
     void cpuGreedyBuild(int,Connectivity *etoe, Connectivity *ntoe, 
                         Connectivity*eton, long *sizes, int numSub);
     void newGreedyBuild(int nsub, Connectivity *ntoe,
                         Connectivity*eton);

     void outputDump(FILE *,int);
     int num(int iSub) const { return pele[iSub+1]-pele[iSub]; }
     int *operator[](int iSub) { return eln+pele[iSub]; }
};

#endif
