#ifndef _MATRIXTIMERS_H_
#define _MATRIXTIMERS_H_
#include <cstdio>


class MatrixTimers {
 public:
   double factor;
   double assemble;
   double formTime;
   double readTime;
   double readDecomp;
   double setUpDataTime;
   double makeInterface;
   double distributeBCs;
   double makeConnectivity;
   double renumbering;
   double createDofs;
   double makeSubDomains;
   double makeInternalInfo;
   double constructTime;
   double receiveFluidTime;
   double sendFluidTime;
   double formRhs;
   double updateState;

   long memoryParse;
   long memorySetUp;
   long memoryForm;
   long memorySolve;

   long memoryInterface;
   long memoryDistBC;
   long memoryConnect;
   long memoryCPUMAP;
   long memorySubdomain;
   long memoryInternal;
   long memoryMPCToNode;
   long memoryMPCToSub;
   long memoryDistMPC;
   double distributeMPCs;
   long memoryDecomp;

   long memoryElemToNode;
   long memorySubToElem;
   long memoryNodeToSub;
   long memorySubToNode;
   long memoryNodeToElem;
   long memoryElemToSub;

   MatrixTimers() { distributeMPCs = 0.0; memoryMPCToNode = 0; memoryMPCToSub = 0; 
                    factor = 0.0; assemble = 0.0; readTime = 0.0;
                    readDecomp = 0.0;
                    formTime = 0.0; setUpDataTime = 0.0; 
                    makeInterface = 0.0; distributeBCs = 0.0; 
                    makeConnectivity = 0.0; makeSubDomains = 0.0; 
                    makeInternalInfo = 0.0; createDofs = 0.0;
                    renumbering = 0.0; constructTime = 0.0;
                    receiveFluidTime = 0.0; sendFluidTime = 0.0;
                    formRhs = 0.0; updateState = 0.0;
                    memoryParse = 0; memorySetUp = 0; memoryForm = 0;
                    memorySolve = 0; memoryInterface = 0; memoryDistBC = 0;
                    memoryConnect = 0; memorySubdomain = 0; memoryInternal = 0;
                    memoryDecomp = 0; memoryElemToNode = 0; memorySubToElem = 0;
                    memoryNodeToSub = 0; memorySubToNode = 0; memoryCPUMAP = 0;
                    memoryNodeToElem = 0; memoryElemToSub = 0; 
                  }
};

#endif
