#include <Element.d/Element.h>
#include <Dec.d/Decomp.d/DecOpti.h>
#include <Utils.d/Connectivity.h>
#include <Dec.d/Decomp.d/Decomp.h>
#include <Threads.d/PHelper.h>
#include <cstdio>
#include <iostream>

void
DecOpti::minInterface(int sub)
{
// First simple optimization phase: move elements to minimize interface
  int iElem, iNode;
  int elem, node;
// Loop over all the subdomain elements and build the list of nodes for
// this subdomain and count their weight.
  int nodeIndex = 0;
  // loop over the elements
  elem = firstElem[sub];
  while(elem >= 0) {
     for(iNode = 0; iNode < eToN->num(elem); ++iNode) {
       node = (*eToN)[elem][iNode];
       if(mask[node] != sub) {
          mask[node] = sub;
          weight[node] = 1;
          subNode[nodeIndex++] = node;
       } else
         weight[node] += 1;
     }
   elem = nextElem[elem];
  }

  // Loop over the subdomain nodes and examine the elements connected to
  // a node with weight of 1.
  for(iNode = 0; iNode < nodeIndex; ++iNode) {
    node = subNode[iNode];
    if(weight[node] == 1 && nToE->num(node) != 1) {
      // examine all the elements attached to this node
      for(iElem = 0; iElem < nToE->num(node); ++iElem) {
        elem = (*nToE)[node][iElem];
        if(subAssign[elem] != sub) continue;
        int newSub = findBestSub(elem);
        if(newSub == sub) continue;
        // if we move it, update all weights
        else {
          // fprintf(stderr, "Reassign %d to %d\n",elem+1, newSub+1);
          reAssign(elem, newSub);
          for(int iNode = 0; iNode < eToN->num(elem); ++iNode) {
            node = (*eToN)[elem][iNode];
            weight[node] -= 1;
            // see if we reduced the interface size
            if(weight[node] == 0 && nToE->num(node) > 1)
               improve++;
          }
        }
      }
    }
  }
}
  


int
DecOpti::findBestSub(int elem)
{
 int newSub = subAssign[elem];
 int maxND = 0;
 int iElem, iNode;
 // this can be replaced by a double loop if eToE is unavailable
 for(iElem = 0; iElem < eToE->num(elem); ++iElem) {
   subWeight[subAssign[(*eToE)[elem][iElem]]] = 0;
   subMask[subAssign[(*eToE)[elem][iElem]]] = -1;
 }

 for(iNode = 0; iNode < eToN->num(elem); ++iNode) {
   int node = (*eToN)[elem][iNode];
   for(iElem = 0; iElem < nToE->num(node); ++iElem) {
     int elemI = (*nToE)[node][iElem];
     if(elemI != elem) {
       int elemSub = subAssign[elemI];
       if(subMask[elemSub] != node) {
         subMask[elemSub] = node;
         subWeight[elemSub] += 1;
         if(subWeight[elemSub] > maxND ||
             (subWeight[elemSub] == maxND && elemSub == subAssign[elem])) {
             newSub = elemSub;
             maxND = subWeight[elemSub];
         } 
       }
     }
   }
 }
 return newSub;
}

void
DecOpti::reAssign(int elem, int sub)
{
 if(previousElem[elem] >= 0)
   nextElem[previousElem[elem]] = nextElem[elem];
 else {
   firstElem[subAssign[elem]] = nextElem[elem]; 
 }
 if(nextElem[elem] >= 0)
   previousElem[nextElem[elem]] = previousElem[elem];
 if(firstElem[sub] >= 0)
   previousElem[firstElem[sub]] = elem;
 nextElem[elem] = firstElem[sub];
 previousElem[elem] = -1;
 firstElem[sub] = elem;
 subAssign[elem] = sub;
}


DecOpti::DecOpti(Connectivity *eton, Connectivity *ntoe, Connectivity *etoe, 
                 Decomposition *_dec)
{
 int iNode, iEle, iSub;
 dec = _dec;

 eToN = eton;
 nToE = ntoe;
 eToE = etoe;

 numEle = eToN->csize();
 numSub = dec->nsub;
 int numNodes = ntoe->csize();
 subAssign = new int[numEle];
 firstElem = new int[numSub];
 nextElem = new int[numEle];
 previousElem = new int[numEle];

 subWeight = new int[numSub];
 subMask = new int[numSub];
 subNode = new int[numNodes];

 weight = new int[numNodes];
 mask = new int[numNodes];
 for(iNode = 0; iNode < numNodes; ++iNode)
   mask[iNode] = -1;

 for(iEle = 0; iEle < numEle; ++iEle)
   subAssign[iEle] = -1;

 for(iSub = 0; iSub < numSub; ++iSub) {
   firstElem[iSub] = -1;
   for(iEle = dec->pele[iSub]; iEle < dec->pele[iSub+1]; ++iEle) {
     int elem = dec->eln[iEle];
     if(firstElem[iSub] >= 0)
       previousElem[firstElem[iSub]] = elem;
     nextElem[elem] = firstElem[iSub];
     previousElem[elem] = -1;
     firstElem[iSub] = elem;
     subAssign[elem] = iSub;
   }
 }
 improve = 0;
}

Decomposition *
DecOpti::minInterface()
{

 // KHP: MODIFICATION
 int iSub;
 for(iSub = 0; iSub < numSub; ++iSub)
   minInterface(iSub);
 //execParal(numSub, this, &DecOpti::minInterface);

 // fprintf(stderr," ... Total improvement %d\n",improve);

 Decomposition *newDec = new Decomposition;
 newDec->esname = dec->esname;
 newDec->nsub = numSub;

 return newDec;
}

Decomposition *
improveDec(Elemset *es, Connectivity *eton, Connectivity *ntoe, 
           Connectivity *etoe, Decomposition *dec)
{
  std::cout << " *** WARNING Improving the decomposition in this way is no longer guaranteed (facelist has not been ported to FEM)" << std::endl;
  DecOpti dopt(eton,ntoe,etoe,dec);
  Decomposition *newDec = dopt.minInterface();
  newDec = dopt.faceComponents(es);
  newDec->esname = dec->esname;
  newDec->setName((char*)"OptDec");
  return newDec;
}

#include <Driver.d/PolygonSet.h>

Decomposition *
DecOpti::faceComponents(Elemset *es)
{
  PolygonSet pset;
  pset.isID = 1;
  int *pc = new int[numEle+1];
  int totFaces = 0;
  int iEle;
  for(iEle = 0; iEle < numEle; ++iEle) {
    pc[iEle] = totFaces;
    if((*es)[iEle]) totFaces +=  (*es)[iEle]->facelist(pset);
  }
  pc[numEle] = totFaces;

  int *fconnect = new int[totFaces];
  for(iEle = 0; iEle < numEle; ++iEle) {
     if((*es)[iEle]) (*es)[iEle]->facelist(pset, fconnect+pc[iEle]);
  }

 Connectivity eToF(numEle,pc,fconnect);
 Connectivity *fToE = eToF.alloc_reverse();
 Connectivity *eFe = eToF.transcon(fToE);

 //int exactNumEle = eToN->numNonZeroP();
 // Now loop over the subdomains
 int *xComp = new int[numEle+1];
 int *eleMask = new int[numEle];
 int *eleList = new int[numEle];
 for(iEle = 0; iEle < numEle; ++iEle)
   eleMask[iEle] = 1;
 int numComp = 0;
 xComp[0] = 0;
 iEle = 0;
 int iSub;
 for(iSub = 0; iSub < numSub; ++iSub) {
   int ele = firstElem[iSub];
   while(ele >= 0) {
     if(eleMask[ele]) {
       eleMask[ele] = 0;
       eleList[iEle] = ele;
       int jEle = iEle;
       iEle = iEle+1;
       for(; jEle < iEle; ++jEle) {
          int elemJ = eleList[jEle];
          for(int kEle = 0; kEle < eFe->num(elemJ); ++kEle) {
            int elemK = (*eFe)[elemJ][kEle];
            if(subAssign[elemK] != iSub) continue;
            if(eleMask[elemK]) {
               eleMask[elemK] = 0;
               eleList[iEle++] = elemK;
            }
          }
       }
       xComp[numComp+1] = iEle;
       numComp++;
     }
     ele = nextElem[ele];
   }
 }

 /* Here is an optional phase that reattaches 'small' components to bigger
    subdomains */
 // Find the small components
 int small = (int) 0.2*numEle/numSub;
 int *subNum = new int[numComp];
 int *finalP = new int[numComp+1]; //slightly over dimensioned

 int sNum=0;
 for(iSub = 0; iSub < numComp; ++iSub) {
   if(xComp[iSub+1]-xComp[iSub] < small)
     subNum[iSub] = -1;
   else
    {
     subNum[iSub] = sNum;
     finalP[sNum+1] = xComp[iSub+1]-xComp[iSub];
     sNum++;
    }
   int cNum = subNum[iSub];
   for(iEle = xComp[iSub]; iEle < xComp[iSub+1]; ++iEle)
     subAssign[eleList[iEle]] = cNum;
 }

 int *nsubMask = new int[sNum];
 int *subCount = new int[sNum];
 for(iSub = 0; iSub < sNum; ++iSub)
   nsubMask[iSub] = -1;
 for(iSub = 0; iSub < numComp; ++iSub) {
   int maxCount =0;
   int bestSub = -1;
   if(subNum[iSub] < 0) {
     for(iEle= xComp[iSub]; iEle < xComp[iSub+1]; ++iEle) {
        int eleI = eleList[iEle];
        // Loop over the neighbors and count the number of faces they share
        for(int jEle = 0; jEle < eFe->num(eleI); ++jEle) {
          int eleJ = (*eFe)[eleI][jEle];
          if(subAssign[eleJ] < 0) continue;
          int subN = subAssign[eleJ];
          if(nsubMask[subN] != iSub) {
            nsubMask[subN] = iSub;
            subCount[subN] = 1;
          }
          else
            subCount[subN]++;
          if(subCount[subN] > maxCount) {
            maxCount = subCount[subN];
            bestSub = subN;
          }
        }
       
     } 
     if(bestSub < 0) {
        bestSub = sNum++;
        finalP[sNum+1] = 0;
     }
     subNum[iSub] = bestSub;
     finalP[bestSub+1] += xComp[iSub+1]-xComp[iSub];
   }
 }

 int tot = 0;
 finalP[0] = 0;
 for(iSub = 0; iSub < sNum; ++iSub) {
   int pt = tot;
   tot += finalP[iSub+1];
   finalP[iSub+1] = pt;
 }

 int *finalList = new int[numEle];
 for(iSub = 0; iSub < numComp; ++iSub) {
   int sub = subNum[iSub];
   for(iEle = xComp[iSub]; iEle < xComp[iSub+1]; ++iEle) {
     finalList[finalP[sub+1]] = eleList[iEle];
     finalP[sub+1] ++;
   }
 }

 Decomposition *newDec = new Decomposition;
 newDec->pele = finalP;
 newDec->eln = finalList;
 newDec->nsub = sNum;
 return newDec;
}
