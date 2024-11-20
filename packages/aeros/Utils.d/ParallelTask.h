#ifndef PARALLELTASK_H_
#define PARALLELTASK_H_

class Connectivity;
class EqNumberer;
class Elemset;

class ParallelTask {
   int numSub;
 public:

   ParallelTask(int _numSub) { numSub = _numSub; }

   int findProfileSize(Connectivity **nodeToNode,EqNumberer **eqNums,
                       int *sizes);

   long findProfileSize(Connectivity *stoe,Elemset &allElements,
                       long *sizes);

   void findSize(int iSub, Connectivity **nodeToNode,
                EqNumberer **eqNums,int *sizes);

   void findSizes(int iSub, long *sizes, 
                  Connectivity *stoe, Elemset *allElements );

};

#endif

