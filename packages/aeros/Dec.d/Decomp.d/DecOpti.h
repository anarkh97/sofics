#ifndef _DECOPTI_H_
#define _DECOPTI_H_

class Connectivity;
class Decomposition;
class Elemset;

class DecOpti {
   Decomposition *dec;
   Connectivity *eToN, *nToE, *eToE;
   int numEle;       // number of elements
   int numSub;       // number of subdomains
   int *subAssign;   // to which sub an element belongs
   int *firstElem;   // first element in a sub
   int *nextElem;    // link to next element in subdomain
   int *previousElem;
   int *subWeight;
   int *subMask;
   int *weight;
   int *mask;
   int *subNode;

   int improve;
 public:
   DecOpti(Connectivity *eton, Connectivity *ntoe, Connectivity *etoe, 
           Decomposition *dec);
   void minInterface(int sub);
   Decomposition *minInterface();
   int findBestSub(int elem);
   void reAssign(int elem, int sub);
   void assign(int elem, int sub);
   Decomposition *faceComponents(Elemset *es);
};

Decomposition * improveDec(Decomposition *dec);

Decomposition *
improveDec(Elemset *es,
           Connectivity *eton, Connectivity *ntoe, Connectivity *etoe,
           Decomposition *dec);

#endif
