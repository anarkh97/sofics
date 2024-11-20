#include <iostream>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <cmath>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Utils.d/DistHelper.h>

// TG : TEMPORARY FIX !!
#ifdef WINDOWS
#define srand48 srand
#define drand48 rand
#endif

#include <vector>
#include <algorithm>
#include "MultiFront.h"
#include <Dec.d/Decomp.d/Decomp.h>
#include <Utils.d/dofset.h>

extern bool nosa;
extern int verboseFlag;
extern bool allowMechanisms;

#include <Element.d/Element.h>
class DecGluedElement: public Element {
    int *nn;
    int numN;
    Element *e1, *e2;
public:
    int ie1, ie2;
    DecGluedElement(Element *_e1, int _ie1, Element *_e2, int _ie2) {
        e1 = _e1; ie1 = _ie1; e2 = _e2; ie2 = _ie2;
        int* n1 = e1->nodes();
        int* n2 = e2->nodes();
        nn = new int[e1->numNodes()+e2->numNodes()];
        int c=0;
        for(int i=0;i<e1->numNodes();i++) nn[c++] = n1[i];
        for(int i=0;i<e2->numNodes();i++) nn[c++] = n2[i];
        delete[] n1;
        delete[] n2;
        std::sort(nn,nn+e1->numNodes()+e2->numNodes());
        int *p = std::unique(nn,nn+e1->numNodes()+e2->numNodes());
        numN = static_cast<int>(p - nn);
    }
    ~DecGluedElement() override { delete[] nn; }

    int getElementType() const override { return -1; }
    Category getCategory() const override { return Category::Structural; }

    void renum(const int *) override {
        fprintf(stderr,"DecGluedElement::renum not implemented\n");
    }
    void renum(EleRenumMap&) override {
        fprintf(stderr,"DecGluedElement::renum not implemented\n");
    }
    int* dofs(DofSetArray &, int *p) const override {
        fprintf(stderr,"DecGluedElement::dofs not implemented\n");
        return 0;
    }
    void markDofs(DofSetArray &) const override {
        fprintf(stderr,"DecGluedElement::markDofs not implemented\n");
    }
    int numDofs() const override{
        fprintf(stderr,"DecGluedElement::numDofs not implemented\n");
        return 0;
    }
    int numNodes() const override {
        return numN;
    }
    int* nodes(int *p) const override {
        if(p == 0) p = new int[numN];
        for(int i=0;i<numN;i++) p[i] = nn[i];
        return p;
    }
    PrioInfo examine(int sub, MultiFront *mf) override {
        PrioInfo p1 = e1->examine(sub,mf);
        PrioInfo p2 = e2->examine(sub,mf);
        PrioInfo res;
        res.isReady = (p1.isReady || p2.isReady);
        if(!res.isReady) return res;
        if (p1.isReady && p2.isReady) {
            res.priority = (p1.priority < p2.priority)? p1.priority:p2.priority;
        } else if (p1.isReady) {
            res.priority =  p1.priority;
        } else {
            res.priority =  p2.priority;
        }
        return res;
    }
//        bool hasRot() const override { return e1->hasRot() || e2->hasRot(); }

};

OrderList::OrderList(int l)
{
    length = l;
    list = new int[length];
    infList = new int[length];
    begin = end = 0;
    infBegin = infEnd = 0;
}

OrderList::~OrderList()
{
    delete [] list;
    delete [] infList;
}

void
OrderList::add(int n)
{
    list[end] = n;
    end ++;
    if(end >= length)
        end -= length;
}

void
OrderList::addInfinity(int n)
{
    infList[infEnd] = n;
    infEnd ++;
    if(infEnd >= length)
        infEnd -= length;
}

int
OrderList::pop()
{
    if(infBegin != infEnd) {
        int l = infList[infBegin];
        infBegin++;
        if(infBegin >= length)
            infBegin -= length;
        return l;
    }
    if(end == begin) return -1;
    int l = list[begin];
    begin++;
    if(begin >= length)
        begin -= length;
    return l;
}

void
OrderList::clear()
{
    begin = end = infBegin = infEnd = 0;
}

int
OrderList::popInfinity()
{
    if(infBegin != infEnd) {
        int l = infList[infBegin];
        infBegin++;
        if(infBegin >= length)
            infBegin -= length;
        return l;
    }
    return -1;
}

BoundList::BoundList(int mSize)
{
    maxIndex = mSize;
    listIndex = new int[maxIndex];
    list = new int[maxIndex];
    for(int i =0; i < maxIndex; ++i)
        listIndex[i] = -1;
    nObj = 0;
}

BoundList::~BoundList()
{
    delete [] listIndex;
    delete [] list;
}

void
BoundList::add(int x)
{
    listIndex[x] = nObj;
    list[nObj] = x;
    nObj += 1;
}

void
BoundList::remove(int x)
{
    if(listIndex[x] < 0) return;
    listIndex[ list[nObj-1] ] = listIndex[x];
    list[ listIndex[x] ] = list[nObj-1];
    listIndex[x] = -1;
    nObj -= 1;
}

// The simulated anealing rule 
class SARule
{
    double ee[1000];
public:
    SARule();
    bool accept(double v, double T) {
/*
    if(nosa) return true;
    if(v < 0 || T == 0.0) return true;
*/
        if(nosa || (T == 0.0)) return (v < 0);  // PJSA DEBUG
        if(v < 0) return true;  // PJSA DEBUG
        int index = int(200*v/T);
        if(index >= 1000) return false;
        return ee[index] > drand48();
    }
};

SARule::SARule() {
    for(int i = 0; i < 1000; ++i)
        ee[i] = exp(-i/200.0);
}

MultiFront::MultiFront(Elemset *eset, CoordSet *cs, bool have_fsi, bool _fsGlFlag)
    : nw(0), nSubPerNode(0), boundIndex(0), boundNode(0), elemPerSub(0), subWeight(0), arInfo(0), nodeMask(0), tag(0)
#ifdef NOTMPL
    , prio(eset->size())
#endif
{
    fsGlFlag = _fsGlFlag;
    iterationN=0;
    elems = eset;
    fsGlued_eset = 0;
    nds = cs;
    eToN = new Connectivity(eset->asSet());
    nToE = eToN->alloc_reverse();
    numEle = elems->size();
    if(verboseFlag) filePrint(stderr, " ... Mesh Contains %d Elements and %d Nodes ...\n", numEle, nToE->csize());

// RT: Glue fluid and solid wet interface elements into elements of two
    fsGluedCounter = -1;
    bool reallyHaveFSI = false;
    if (have_fsi) {
        for(int iEle = 0; iEle < numEle; ++iEle)
            if ((*elems)[iEle] && (*elems)[iEle]->isFsiElement()) {
                reallyHaveFSI = true;
                break;
            }
    }

//int(fsGlFlag));
    if (reallyHaveFSI) fsGlFlag = false;
    if (fsGlFlag) {
        fsGlued_eset = new Elemset();
        fsGluedCounter = 0;
        int *tags = new int[numEle];
        for(int iEle = 0; iEle < numEle; ++iEle) {
            tags[iEle] = -1;
            if ((*elems)[iEle]) fsGlued_eset->elemadd(iEle,(*elems)[iEle]);
        }

        for(int iEle = 0; iEle < numEle; ++iEle) {
            if ( (*elems)[iEle] && (*elems)[iEle]->nDecFaces()>0 ) {
                int eventTag = getEventTag();
                for(int iNode = 0; iNode < eToN->num(iEle); ++iNode) {
                    int node = (*eToN)[iEle][iNode];
                    for(int iEle2=0;iEle2 < nToE->num(node); ++iEle2) {
                        int thisElem = (*nToE)[node][iEle2];
                        if (tags[thisElem]==eventTag) continue;
                        tags[thisElem] = eventTag;
                        if ((*elems)[thisElem] && (*elems)[thisElem]->nDecFaces()>0) {
                            bool isF = (*elems)[iEle]->isFluidElement();
                            bool isF2 = (*elems)[thisElem]->isFluidElement();
                            if ( (isF2 && !isF) || (isF && !isF2) ) {
                                for(int iFace=0; iFace<(*elems)[iEle]->nDecFaces();iFace++) {
                                    int nn[64];
                                    int nFaceN = (*elems)[iEle]->getDecFace(iFace,nn);
                                    std::sort(nn,nn+nFaceN);
                                    for(int iFace2=0; iFace2<(*elems)[thisElem]->nDecFaces();
                                        iFace2++) {
                                        int nn2[64];
                                        int nFaceN2 = (*elems)[thisElem]->getDecFace(iFace2,nn2);
                                        if (nFaceN!=nFaceN2) break;
                                        std::sort(nn2,nn2+nFaceN2);
                                        bool isEqual = true;
                                        int i=0,j=0;
                                        while (i<nFaceN) {
                                            while(j<nFaceN2) {
                                                if (nn[i]!=nn2[j]) j++;
                                                else break;
                                            }
                                            if (j==nFaceN2) { isEqual = false; break; }
                                            i++;
                                        }
                                        if (isEqual) if ((*fsGlued_eset)[iEle]!=0)  {
                                                DecGluedElement *dec_e =
                                                    new DecGluedElement((*elems)[iEle],iEle,
                                                                        (*elems)[thisElem],thisElem);
                                                fsGlued_eset->remove(iEle);
                                                fsGlued_eset->remove(thisElem);
                                                fsGlued_eset->elemadd(numEle+fsGluedCounter,dec_e);
                                                fsGluedCounter++;
                                                //                   fsAffinity.insert(make_pair(iEle,thisElem));
                                            }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        delete[] tags;
        elems = fsGlued_eset;
        delete eToN;
        eToN = new Connectivity(elems->asSet());
        delete nToE;
        nToE = eToN->alloc_reverse();
        numEle = elems->last();
//   if(verboseFlag) filePrint(stderr, " ... Glued Mesh Contains %d Elements and %d Nodes ...\n", numEle, nToE->csize());

    }


    int nrnodes = 0;
    for(int x = 0; x < nToE->csize(); ++x)
        if(nToE->num(x) != 0)
            nrnodes++;
    int iEle, i;
    assignedSubD = new int[numEle];
    elemHasRot   = new bool[numEle];
    // mark all of the elements as not assigned yet
    for(iEle = 0; iEle < numEle; ++iEle) {
        assignedSubD[iEle] = -1;
        if((*elems)[iEle])
            elemHasRot[iEle] = ((*elems)[iEle]->hasRot() || allowMechanisms);
    }
    int numNode = nToE->csize();
    // Now we create a flag array used for checking during an operation
    // if a node or element has already been visited. Since it can be used
    // for nodes or elements, it must be sized by the max.
    int flagSize = (numNode > numEle) ? numNode: numEle;
    flag = new int[flagSize];
    for(i = 0; i < flagSize; ++i)
        flag[i] = -1;
    elemCG = 0;




// JLchange: determine which elements are connected with a fsi node 
    isFsiConnected   = new bool[numEle];
    for(iEle = 0; iEle < numEle; ++iEle)
        if((*elems)[iEle]) isFsiConnected[iEle] = false;
    if (have_fsi) {
        for(iEle = 0; iEle < numEle; ++iEle)
            if ((*elems)[iEle] && (*elems)[iEle]->isFsiElement()) {
                int fluidNode = (*elems)[iEle]->fsiFluidNode();
                int strutNode = (*elems)[iEle]->fsiStrutNode();
                for(int iele = 0; iele < nToE->num(fluidNode); ++iele) {
                    int thisElem = (*nToE)[fluidNode][iele];
                    if ((*elems)[thisElem]) isFsiConnected[thisElem] = true;
                }
                for(int iele = 0; iele < nToE->num(strutNode); ++iele) {
                    int thisElem = (*nToE)[strutNode][iele];
                    if ((*elems)[thisElem]) isFsiConnected[thisElem] = true;
                }
            }

        // Add one such isFsiConnected element to the pop.InfList, so the decomposition will start
        // from this element.
        for(iEle = 0; iEle < numEle; ++iEle)
            if ((*elems)[iEle] && isFsiConnected[iEle]) {
                prio.addInfinity(iEle);
                break;
            }
    }

}

MultiFront::~MultiFront()
{
    if(nw) delete [] nw;
    if(nSubPerNode) delete [] nSubPerNode;
    if(boundIndex) delete [] boundIndex;
    if(boundNode) delete [] boundNode;
    if(elemPerSub) delete[] elemPerSub;
    if(subWeight) delete [] subWeight;
    if(arInfo) delete [] arInfo;
    if(nodeMask) delete [] nodeMask;
    if(eToN) delete eToN;
    if(nToE) delete nToE;
    if(assignedSubD) delete [] assignedSubD;
    if(elemHasRot) delete [] elemHasRot;
    if(flag) delete [] flag;
    if(isFsiConnected) delete [] isFsiConnected;
    if(elemCG) delete [] elemCG;
}

void
MultiFront::redoConnect()
{
    int i;
    delete nToE;
    delete eToN;
    eToN = new Connectivity(elems->asSet());
    nToE = eToN->alloc_reverse();
    numEle = elems->size();
    if(verboseFlag) filePrint(stderr, " ... Mesh Contains %d Elements and %d Nodes ...\n", numEle, nToE->csize());
    delete [] flag;
    int numNode = nToE->csize();
    // Now we create a flag array used for checking during an operation
    // if a node or element has already been visited. Since it can be used
    // for nodes or elements, it must be sized by the max.
    int flagSize = (numNode > numEle) ? numNode: numEle;
    flag = new int[flagSize];
    for(i = 0; i < flagSize; ++i)
        flag[i] = -1;
}

void
MultiFront::updateElement(int cur_elem, int sub, int elem)
{
    PrioInfo ePrio = (*elems)[elem]->examine(sub, this);
    // don't do anything if this element is not ready
// JLchange 
// if(ePrio.isReady == 0) return;
    if ((ePrio.isReady == 0) && !(isFsiConnected[elem])) return;
    // look if this element is in the priority list
#ifdef NOTMPL
    // does this element have infinite priority?
// JLchange 
    if(!(isFsiConnected[elem]) && (ePrio.priority <= -100)) {
        //prio.addInfinity(elem);
        prio.add(elem);
        return;
    }
    // if this element is not yet in the boundary list, add it.
// JLchange:
// prio.add(elem);

    if (isFsiConnected[elem]) prio.addInfinity(elem);
    else prio.add(elem);
#else
    PrioState newState;
 newState.prio = ePrio.priority;
 map<int,PrioState>::iterator it = lastPrio.find(elem);
 if(it != lastPrio.end()) {
   // create the newState with an apperance time unchanged
   newState.time = it->second.time;
   // remove the old priority of this element
   prio.erase(it->second);
   // check if this element is assigned to someone else
   // if so, totally discard it from the lists
   if(assignedSubD[elem] >= 0) {
      lastPrio.erase(it);
      return;
   }
   // insert the element with the new priority
   prio[newState] = elem;
   it->second = newState;
 } else {
   if(assignedSubD[elem] >= 0) return;
   newState.time = time++;
   prio[newState] = elem;
   lastPrio[elem] = newState;
 }
#endif
}

void
MultiFront::addNodesToSub(int cur_elem, int sub, int numNewNodes, int *newNodes)
{
    int iElem, iNode;

    int eventTag = getEventTag();
    for(iNode = 0; iNode < numNewNodes; ++iNode)
        for(iElem = 0; iElem < nToE->num(newNodes[iNode]); ++iElem) {
            int elem = (*nToE)[newNodes[iNode]][iElem];
            if(flag[ elem ] != eventTag && assignedSubD[elem] < 0)
                updateElement(cur_elem, sub, elem);
            flag[elem] = eventTag;
        }
}

void
MultiFront::removeBoundNode(int node)
{
    if(boundIndex[node] >= 0) {
        int replNode = boundNode[ boundIndex[node] ] = boundNode[ --numBoundNodes ];
        // Note no check on the following line is necessary as the only possible
        // problem is if numBoundNodes goes down to zero and in that case
        // this line has no effect
        boundIndex[ replNode ] = boundIndex[ node ];
        boundIndex[node] = -1;
    }
}

void
MultiFront::addBoundNode(int node)
{
    if(boundIndex[node] < 0) {
        boundIndex[node] = numBoundNodes;
        boundNode[ numBoundNodes++ ] = node;
    }
}

int ttt = 0;

void
MultiFront::addElemToSub(int sub, int elem)
{
    int iNode;
    assignedSubD[elem] = sub;
    count++;
    // loop over the connected nodes and update their weight
    // also make up the list of nodes that appear for the first time to
    // be touched by this subdomain
    int numNewNodes = 0;
    // Note: We could get rid of this dbg_alloca by passing an eventTag to the
    // routine examining the elements connected to a node and doing a second
    // loop
    int *newNodes = (int *) dbg_alloca(sizeof(int) * eToN->num(elem));

    for(iNode = 0; iNode < eToN->num(elem); ++iNode) {
        int node = (*eToN)[elem][iNode];
        if((*elems)[elem]->isStart()) {
            --nodeMask[node];
        }
        if(nodeMask[node] == 0) { // remove the node from the interface list
            removeBoundNode(node);
        }
        if(nodeMask[node] < 0) std::cerr << "Bad bad bad" << node << "\n";
        // increment the weight
        incrWeight(sub, node, elemHasRot[elem] ? 1 : 0);

        if(weight(sub, node) == 1 && nodeMask[node] > 0) { // this is the first time we see this node
            newNodes[numNewNodes++] = node;
            addBoundNode(node);
        }
    }
    if(numNewNodes > 0)
        addNodesToSub(elem, sub, numNewNodes, newNodes);
}

void
MultiFront::initDec(int numSub)
{
    // keep the target number of subdomain for the improvement phase
    desiredNSub = numSub;
    nw = new NodeWeight[ nToE->numConnect() ];
    nSubPerNode = new int[nToE->csize()];
    boundIndex = new int[nToE->csize()];
    boundNode = new int[nToE->csize()];
    for(int iNode = 0; iNode < nToE->csize(); ++iNode) {
        boundIndex[iNode] = -1;
        nSubPerNode[iNode] = 0;
    }
    numBoundNodes = 0;
    time = 0;
    tag  = 0;
}

Decomposition *
MultiFront::decompose(int numSub, bool have_fsi)
{
    double   t1 = getTime();
    //double tDec = t1;
    long     m1 = memoryUsed();
    //double mDec = m1;
    initDec(numSub);
    int numNodes = nToE->csize();
    nodeMask = new int[numNodes];
    int iNode, iElem;

    if(verboseFlag) filePrint(stderr," ... Starting Initial Mesh Decomposition ...\n");
    for(iNode = 0; iNode < numNodes; ++iNode)  {
        nodeMask[iNode] = 0;
        for(iElem = 0; iElem < nToE->num(iNode); ++iElem)
            if((*elems)[(*nToE)[iNode][iElem] ]->isStart())
                nodeMask[iNode]++;
    }
    // First compute the total weight
    double totWeight = 0;
    for(iElem = 0; iElem < numEle; ++iElem)
        if((*elems)[iElem])
            totWeight += (*elems)[iElem]->weight();
    double remWeight = totWeight;
    double avgWeight = totWeight/numSub;

    int nStartElem = 0;
    for(iElem = 0; iElem < numEle; ++ iElem)
        if((*elems)[iElem] && (*elems)[iElem]->isStart())
            nStartElem++;
    int startElem = -1;
    int nRemain = 0, nRotEl = 0;
    for(iElem = 0; iElem < numEle; ++ iElem)
        if((*elems)[iElem]) {
            nRemain++;
            if((*elems)[iElem]->hasRot() || allowMechanisms) nRotEl++;
        }
    int iSub = 0;
    //int targetNumEle = nRemain/(numSub);
    // fprintf(stderr, " ... Mesh Contains %d Elements, of Which %d Have Rotational dofs ...\n",
    //         nRemain, nRotEl);

    while(nRemain > 0) {
// JLchange: to start a new element, select a boundary node and an element connected to that 
// boundary node. To start from one fsi element, we can choose the start element from the 
// current priority list, if it is not empty. Also see the comment later.  
        // Pick starting element
        startElem = -1;
        int minDeg = numEle+1;
        int minNode =-1;

// JLchange: add the following lines, such that if prio.infList is not empty then pop one element out 
// and use it to start forming the new subdomain 
        // loop to find an unassigned element from the prio.InfList
        if (have_fsi)
            while (!(prio.empty())) {
                int elem = prio.pop();
                if (assignedSubD[elem] < 0) {
                    startElem = elem;
                    break;
                }
            }

        for(iNode = 0; iNode < numBoundNodes; ++iNode)
            if(nodeMask[boundNode[iNode]] > 0 && nodeMask[boundNode[iNode]] < minDeg) {
                minDeg = nodeMask[boundNode[iNode]];
                minNode = boundNode[iNode];
            }
        if ((minNode >= 0) && (startElem < 0)) {
            for(iElem = 0; iElem < nToE->num(minNode); ++iElem) {
                startElem = (*nToE)[minNode][iElem];
                if((*elems)[ startElem ]->isStart()&& assignedSubD[startElem] < 0)
                    break;
            }
            if(iElem == nToE->num(minNode)) {
                startElem = -1;
            }
        } else {
            // std::cerr << " ...     No Connection! numBoundNodes = " << numBoundNodes << "\n";
        }
        if(startElem < 0) {
            int minDeg = numEle+1;
            int minNode = -1;
            for(iNode = 0; iNode < numNodes; ++iNode) {
                if(nodeMask[iNode] > 0 && nodeMask[iNode] < minDeg) {

                    minDeg = nodeMask[iNode];
                    minNode = iNode;
                }
            }
            if(minNode < 0) break;
            for(iElem = 0; iElem < nToE->num(minNode); ++iElem) {
                startElem = (*nToE)[minNode][iElem];
                if((*elems)[ startElem ]->isStart()&& assignedSubD[startElem] < 0)
                    break;
            }
            if(iElem == nToE->num(minNode)) {
                break;
            }
        }
        prio.clear();
        count = 0;
        // add it to this subdomain and remove it from the lists
        addElemToSub(iSub,startElem);
        double targWeight = avgWeight;
        double curWeight = 0.0;

        for(; curWeight < targWeight && !prio.empty(); ) {
            // loop over the active subdomains
            if(prio.empty()) break;
            // pick the highest priority element in the list
            int elem;
            bool isEmpty = false;
            // loop to find an unassigned element;
            do {
                elem = prio.pop();
                // check if this subdomain has become inactive
                isEmpty = prio.empty();
            } while(isEmpty == false && assignedSubD[elem] >= 0);
            if(isEmpty && assignedSubD[elem] >= 0)
                break;
            // add it to this subdomain and remove it from the lists
            if(assignedSubD[elem] < 0) {
                double thisElemWeight = (*elems)[elem]->weight();
                if(curWeight+0.5*thisElemWeight <= targWeight) {
                    addElemToSub(iSub,elem);
                    curWeight += thisElemWeight;
                    remWeight -= thisElemWeight;
                } else
                    break;
            }
        }
/* 
JLchange: if add all elements connected with fsi to priority list, the following few line 
have the risk of putting all such element in the first subdomain. subdomain zero. 
Maybe we can put a weight limit here and if it is over the weight limit then start a new 
subdomain, from a fsi node/element in the current priority list. 
*/
// JLchange 
/*
   int infElem;
   while((infElem = prio.popInfinity()) >= 0) {
     if(assignedSubD[infElem] < 0) {
       addElemToSub(iSub, infElem);
       remWeight -= (*elems)[infElem]->weight();
     }
   }
*/

        int infElem;
        while (((infElem = prio.popInfinity()) >= 0) && (curWeight < targWeight)) {
            if(assignedSubD[infElem] < 0) {
                double thisElemWeight = (*elems)[infElem]->weight();
                addElemToSub(iSub, infElem);
                curWeight += thisElemWeight;
                remWeight -= thisElemWeight;
            }
        }

        iSub++;
        nRemain -= count;
    }
    int iEle;
    for(iEle = 0; iEle < numEle; ++iEle) {
        if((*elems)[iEle] && assignedSubD[iEle] < 0)
            assignedSubD[iEle] = iSub++;
    }
    Decomposition *dec = new Decomposition;
    dec->nsub = iSub;
    dec->pele = new int[dec->nsub+1];
    dec->eln = new int[numEle];

    numSub = iSub;
    for(iSub = 0; iSub <= numSub; ++iSub)
        dec->pele[iSub] = 0;
    for(iEle = 0; iEle < numEle; ++iEle)
        if((*elems)[iEle]) {
            if(assignedSubD[iEle] < 0) {
                std::cerr << "Element " << iEle+1 << "Is not assigned. it will not be in"
                                                     "the decomposition\n" ;
                continue;
            }
            dec->pele[assignedSubD[iEle]]++;
        }

    for(iSub = 0; iSub < numSub; ++iSub)
        dec->pele[iSub+1] += dec->pele[iSub];

    for(iEle = 0; iEle < numEle; ++iEle)
        if((*elems)[iEle]) {
            if(assignedSubD[iEle] < 0) continue;
            dec->eln[ --dec->pele[assignedSubD[iEle]] ] = iEle;
        }

    // filePrint(stderr," ... Generated Initial Mesh Decomposition In %14.5f sec and %14.3f Mb ...\n",
    //           (getTime() - t1)/1000.0, (memoryUsed() - m1)/(1024.0*1024.0));

    // Now apply improvements to the decomposition
    if(verboseFlag) filePrint(stderr, " ...     Initial Number of Subdomains = %d ...\n", dec->nsub);
    dec = improveDec(dec, 0.2);
    rebuildInfo(dec);

    //optimization
    if(true) {
        if(verboseFlag) filePrint(stderr," ... Optimizing Initial Mesh Decomposition with %d subs ...\n", dec->nsub);
        t1 = getTime();
        m1 = memoryUsed();
        dec = optimize(dec);

        if(verboseFlag) filePrint(stderr," ...     Updated Number of Subdomains = %d ...\n", dec->nsub);
        // filePrint(stderr," ... Optimized Initial Mesh Decomposition In %14.5f sec and %14.3f Mb\n",
        //           (getTime() - t1)/1000.0, (memoryUsed() - m1)/(1024.0*1024.0));

        if(verboseFlag) filePrint(stderr," ... Checking Components            ...\n");
        t1 = getTime();
        m1 = memoryUsed();
        dec = checkComponents(dec);

        if(verboseFlag) filePrint(stderr, " ...     Updated Number of Subdomains = %d ...\n", dec->nsub);
        // filePrint(stderr," ... Checked Optimized Mesh Decomposition In %14.5f sec and %14.3f Mb\n",
        //           (getTime() - t1)/1000.0,(memoryUsed() - m1)/(1024.0*1024.0));

        if(verboseFlag) filePrint(stderr," ... Re-tying Subdomains            ...\n");
        t1 = getTime();
        m1 = memoryUsed();
        rebuildInfo(dec);

        dec = improveDec(dec,0.1);
        rebuildInfo(dec);
        // filePrint(stderr," ... Re-tied Checked Mesh Decomposition In %14.5f sec and %14.3f Mb\n",
        // (getTime() - t1)/1000.0,(memoryUsed() - m1)/(1024.0*1024.0));
    }

// JLchange: 
    // Update fsi elements such that each fsi element is connected to a structural element in a subdomain
    if (have_fsi) {
        dec = updateFsiConnection(dec);
        dec = eliminateEmpty(dec);
        rebuildInfo(dec);
// dec = removeFsiElements(dec); // PJSA 7-31-06
    }

    // filePrint(stderr," ... Total Time for Mesh Decomposition Is %14.5f sec and %14.3f Mb\n",
    //          (getTime() - tDec)/1000.0, (memoryUsed()-mDec)/(1024.0*1024.0));


    if (fsGluedCounter>0) {

        Decomposition *origDec = dec;
        dec = new Decomposition;
        dec->nsub = origDec->nsub;
        dec->pele = new int[dec->nsub+1];
        dec->eln = new int[numEle-fsGluedCounter];

        for(iSub = 0; iSub <= dec->nsub; ++iSub)
            dec->pele[iSub] = 0;

        for(iSub = 0; iSub < dec->nsub; ++iSub) {
            for(int i = origDec->pele[iSub]; i< origDec->pele[iSub+1]; i++) {
                dec->pele[iSub+1]++;
                if (origDec->eln[i] >= numEle-fsGluedCounter) dec->pele[iSub+1]++;
            }
        }

        for(iSub = 0; iSub < dec->nsub; ++iSub) {
            dec->pele[iSub+1] += dec->pele[iSub];
        }

        for(iSub = 0; iSub < dec->nsub; ++iSub) {
            int c = dec->pele[iSub];
            for(int i = origDec->pele[iSub]; i< origDec->pele[iSub+1]; i++) {
                int iEle = origDec->eln[i];
                if (iEle >= numEle-fsGluedCounter) {
                    DecGluedElement *pe = dynamic_cast<DecGluedElement*>( (*elems)[iEle] );
                    dec->eln[c] = pe->ie1;
                    c++;
                    dec->eln[c] = pe->ie2;
                    c++;
                } else {
                    dec->eln[c] = iEle;
                    c++;
                }
            }
        }

        delete origDec;
    }
    if (fsGluedCounter>-1) {
// RT: Where to delete elems = fsGlued_eset;
        delete fsGlued_eset;
    }

    return dec;
}

int
MultiFront::weight(int sub, int node)
{
    int ip;
    for(ip = nToE->offset(node); ip < nToE->offset(node+1); ++ip) {
        if(nw[ip].subd == sub) return nw[ip].weight;
        if(nw[ip].subd < 0) {
            return 0;
        }
    }
    return 0;
}

int
MultiFront::numSubForNode(int node)
{
    return nSubPerNode[node];
}

int
MultiFront::incrWeight(int sub, int node, int rmode)
{
    int ip;
    int initP = nToE->offset(node);
    int finalP = initP+nSubPerNode[node];
    for(ip = initP; ip < finalP; ++ip) {
        if(nw[ip].subd == sub) {
            nw[ip].weight++;
            nw[ip].rotWeight += rmode;
            return nw[ip].weight;
        }
    }
    if(ip < nToE->offset(node+1)) {
        nw[ip].subd = sub;
        nw[ip].weight = 1;
        nw[ip].rotWeight = rmode;
        nSubPerNode[node]++;
        return 1;
    }
    std::cerr << "Error in node weight, subdomain count exceeds maximum " << nSubPerNode[node] << " vs "<< nToE->offset(node+1)-nToE->offset(node);
    throw "Error in node weight, subdomain count exceeds maximum";
    return 0;
}

int
MultiFront::decWeight(int sub, int node, int rmode)
{
    int ip;
    int initP = nToE->offset(node);
    int finalP = initP+nSubPerNode[node];
    for(ip = initP; ip < finalP; ++ip) {
        if(nw[ip].subd == sub) {
            nw[ip].weight--;
            nw[ip].rotWeight -= rmode;
            if(nw[ip].weight == 0) { // eliminate this subdomain from the list
                nw[ip] = nw[finalP-1];
                nSubPerNode[node]--;
                return 0;
            }
            return nw[ip].weight;
        }
    }
    std::cerr << "Error in dec node weight, subdomain was not found\n";
    throw "Error in dec node weight, subdomain was not found";
}

int
MultiFront::numConnectedSub(int node)
{
    return nSubPerNode[node];
}


int
MultiFront::weight(int node)
{
    return nToE->num(node);
}

Decomposition *
MultiFront::improveDec(Decomposition *origDec, double coef)
{
    computeBalCost(origDec);
    int iSub, iEle, jSub, iNode;
    // total number of elements in the decomposition
    // int totNEle = origDec->pele[origDec->nsub];
    // int targetElePerSub = totNEle/desiredNSub;
    // threshhold for trying to reattach a subdomain
    double threshhold = coef*(origDec->nsub)*avgWeight/desiredNSub;

    bool hasRemap = false;
    int *subRemap = new int[origDec->nsub];
    int *subMask = new int[origDec->nsub];
    for(iSub = 0; iSub < origDec->nsub; ++iSub)
    {
        subMask[iSub] = -1;
        subRemap[iSub] = iSub; // initialize the remap to identity
    }
    int *subList = new int[origDec->nsub];
    int totSubs = 0;

    for(iSub = 0; iSub < origDec->nsub; ++iSub)
        if(subWeight[iSub] < threshhold)
        {
            // We try to find the best subto reattach this subdomain
            // first locate all the subdomains that touch this one
            int nCandidate = 0;
            for(iEle = 0; iEle < origDec->num(iSub); ++iEle)
            {
                int ele = (*origDec)[iSub][iEle];
                for(iNode = 0; iNode < eToN->num(ele); ++iNode)
                {
                    int node = (*eToN)[ele][iNode];
                    int ip;
                    int initP = nToE->offset(node);
                    int finalP = initP+nSubPerNode[node];
                    for(ip = initP; ip < finalP; ++ip) {
                        if(nw[ip].subd != iSub && subMask[nw[ip].subd] != iSub)
                        {
                            subMask[nw[ip].subd] = iSub;
                            subList[nCandidate++] = nw[ip].subd;
                        }
                    }
                }
            }
            int bestCandidate = iSub;
            int bestPrio = 0;
            for(jSub = 0; jSub < nCandidate; ++jSub)
            {
                int curPrio = 0;
                // make sure we remap to a lower number!!
                for(iEle = 0; iEle < origDec->num(iSub); ++iEle)
                {
                    int ele = (*origDec)[iSub][iEle];
                    // check if this element provides a bond to the considered subdomain
                    if((*elems)[ele]->examine(subList[jSub], this).isReady)
                    {
                        curPrio++;
                    }
                }
                if(curPrio > bestPrio)
                {
                    bestPrio = curPrio;
                    bestCandidate = subList[jSub];
                }
            }
            if(bestCandidate != iSub)
            {
                hasRemap = true;
                subRemap[iSub] = subRemap[bestCandidate];
                totSubs--;
            }
        }

    if(hasRemap)
    {
        int *finalRemap = new int[origDec->nsub];
        bool *isInspected = new bool[origDec->nsub];
        int *subSequence = new int[origDec->nsub];
        for(iSub = 0; iSub < origDec->nsub; ++iSub)
            isInspected[iSub] = false;
        totSubs = 0;
        // the follwing loop is correct, assuming subRemap always maps to a
        // lower subdomain number
        for(iSub = 0; iSub < origDec->nsub; ++iSub)
        {
            // check if this subdomain has already been inspected below
            if(isInspected[iSub])
                continue;

            if(subRemap[iSub] == iSub)
            {
                finalRemap[iSub] = totSubs++;
                isInspected[iSub] = true;
            } else {
                int numFollow = 1;
                subSequence[0] = iSub;
                int curSub = subRemap[iSub];
                while(isInspected[curSub] == false) {
                    if(subRemap[curSub] == curSub) {
                        finalRemap[curSub] = totSubs++;
                        isInspected[curSub] = true;
                    } else {
                        subSequence[numFollow++] = curSub;
                        curSub =  subRemap[curSub];
                    }
                }
                for(; numFollow-- > 0; ) {
                    finalRemap[subSequence[numFollow]] =
                        finalRemap[subRemap[subSequence[numFollow]]];
                    isInspected[subSequence[numFollow]] = true;
                }
            }
        }
        int *origRemap = subRemap;
        subRemap = finalRemap;
        // apply the remapping to the node and element data structures
        int *subP = new int[totSubs+1];
        for(iSub = 0; iSub < totSubs+1; ++iSub)
            subP[iSub] = 0;
        for(iSub = 0; iSub < origDec->nsub; ++iSub)
            subP[subRemap[iSub]] += origDec->num(iSub);
        for(iSub = 0; iSub < totSubs; ++iSub)
            subP[iSub+1] += subP[iSub];
        int *elL = new int[subP[totSubs]];
        for(iSub = 0; iSub < origDec->nsub; ++iSub)
        {
            int mps = subRemap[iSub];
            for(iEle = 0; iEle < origDec->num(iSub); ++iEle)
                elL[--subP[mps]] = (*origDec)[iSub][iEle];
        }
        delete [] origRemap;
        delete origDec;
        origDec = new Decomposition(totSubs, subP, elL);
    }
    delete [] subRemap;
    delete [] subMask;
    delete [] subList;
    return origDec;
}

void
MultiFront::rebuildInfo(Decomposition *dec)
{
    int iSub, iEle, iNode, i;
    int np = nToE->numConnect();
    NodeWeight initW;
    for(i = 0; i < np; ++i)
        nw[i] = initW;

    for(i = 0; i < nToE->csize(); ++i)
        nSubPerNode[i] = 0;
    for(iNode = 0; iNode < nToE->csize(); ++iNode)
        boundIndex[iNode] = -1;
    numBoundNodes = 0;
    if(elemPerSub) delete [] elemPerSub;
    elemPerSub = new int[dec->nsub];
    for(iSub = 0; iSub < dec->nsub; ++iSub) {
        elemPerSub[iSub] = dec->num(iSub);
        for(iEle = 0; iEle < dec->num(iSub); ++ iEle) {
            int ele = (*dec)[iSub][iEle];
            assignedSubD[ele] = iSub;
            for(iNode = 0; iNode < eToN->num(ele); ++iNode)
                incrWeight(iSub, (*eToN)[ele][iNode], elemHasRot[ele] ? 1 : 0);
        }
    }

    for(iNode = 0; iNode < nToE->csize(); ++iNode) {
        for(iEle = 0; iEle < nToE->num(iNode); ++iEle)
            if(assignedSubD[ (*nToE)[iNode][iEle] ] !=
               assignedSubD[ (*nToE)[iNode][0] ]) {
                addBoundNode(iNode);
                break;
            }
    }
}

Decomposition *
MultiFront::updateFsiConnection(Decomposition *origDec)
{
    int iEle, iele, iSub;
    int numSub = origDec->nsub;

    for(iEle = 0; iEle < numEle; ++iEle)
        if ((*elems)[iEle] && (*elems)[iEle]->isFsiElement()) {
            // Fid the structral node attached to this Fsi element
            //int fluidNode = (*elems)[iEle]->fsiFluidNode();
            int strutNode = (*elems)[iEle]->fsiStrutNode();
            // find if the structure node is in a structral element inside this assigned subdomain
            int thisSub = assignedSubD[iEle];
            int targetSub;
            bool strutEleFound = false;
            for(iele = 0; iele < nToE->num(strutNode); ++iele) {
                int thisElem = (*nToE)[strutNode][iele];
                if (!((*elems)[thisElem]->isFsiElement())) {
                    targetSub = assignedSubD[thisElem];
                    if (thisSub == targetSub) strutEleFound = true;
                }
            }
            if (strutEleFound == false)
                assignedSubD[iEle] = targetSub;
        }

    for(iSub = 0; iSub <= numSub; ++iSub)
        origDec->pele[iSub] = 0;
    for(iEle = 0; iEle < numEle; ++iEle)
        if((*elems)[iEle]) {
            if(assignedSubD[iEle] < 0) {
                std::cerr << "Element " << iEle+1 << "Is not assigned. it will not be in"
                                                     "the decomposition\n";
                continue;
            }
            origDec->pele[assignedSubD[iEle]]++;
        }

    for(iSub = 0; iSub < numSub; ++iSub)
        origDec->pele[iSub+1] += origDec->pele[iSub];

    for(iEle = 0; iEle < numEle; ++iEle)
        if((*elems)[iEle]) {
            if(assignedSubD[iEle] < 0)
                continue;
            origDec->eln[ --origDec->pele[assignedSubD[iEle]] ] = iEle;
        }

    return origDec;
}

Decomposition *
MultiFront::eliminateEmpty(Decomposition *origDec)
{
    int iEle, iSub;
    int numSub = origDec->nsub;
    int realNumSub = 0;
    for(int iSub = 0; iSub < numSub; ++iSub)
        if (origDec->num(iSub) > 0) realNumSub++;

    if (realNumSub == numSub) return origDec;
    else {
        int *subP = new int[realNumSub+1];
        for(iSub = 0; iSub < realNumSub+1; iSub++)
            subP[iSub] = 0;
        int index = 0;
        for (iSub = 0; iSub < origDec->nsub; ++iSub)
            if (origDec->num(iSub) > 0) subP[++index] = origDec->num(iSub);
        for(iSub = 0; iSub < realNumSub; ++iSub)
            subP[iSub+1] += subP[iSub];
        int *elL = new int[subP[realNumSub]];
        for(iEle = 0; iEle < numEle; iEle++)
            if((*elems)[iEle]) {
                if(assignedSubD[iEle] < 0) continue;
                elL[iEle] = origDec->eln[iEle];
            }

        delete [] origDec->pele;
        delete [] origDec->eln;
        delete origDec;
        origDec = new Decomposition(realNumSub, subP, elL);
        return origDec;
    }
}

Decomposition *
MultiFront::removeFsiElements(Decomposition *dec)
{
    filePrint(stderr, " ... Removing Fsi Elements \n");
    int nsubs = 0;
    int *ptr = new int[dec->nsub+1];
    int *list = new int[numEle];

    ptr[0] = 0;
    int total = 0;
    for(int iSub = 0; iSub < dec->nsub; ++iSub) {
        int count = 0;
        for(int iEle = 0; iEle < dec->num(iSub); ++ iEle) {
            int ele = (*dec)[iSub][iEle];
            if(!(*elems)[ele]->isFsiElement()) {
                count++;
                list[total++] = ele;
            }
        }
        if(count>0) {
            nsubs++;
            ptr[nsubs] = ptr[nsubs-1]+count;
        }
    }
    delete dec;
    dec = new Decomposition(nsubs,ptr,list);
    return dec;
}

double
MultiFront::buildARInfo(Decomposition *dec)
{
    bool newAlloc = (elemCG == 0);
    if(newAlloc)
        elemCG = new double[ numEle ][3];
    int iSub, iEle;
    double cost = 0;
    for(iSub = 0; iSub < dec->nsub; ++iSub) {
        arInfo[iSub].ncgx = arInfo[iSub].ncgy = arInfo[iSub].ncgz = 0.0;
        arInfo[iSub].nxs = arInfo[iSub].nys = arInfo[iSub].nzs = 0.0;
        for(iEle = 0; iEle < dec->num(iSub); ++iEle) {
            int el = (*dec)[iSub][iEle];
            if(newAlloc)
                (*elems)[el]->getCG(*nds, elemCG[el][0], elemCG[el][1], elemCG[el][2]);
            arInfo[iSub].ncgx += elemCG[el][0];
            arInfo[iSub].ncgy += elemCG[el][1];
            arInfo[iSub].ncgz += elemCG[el][2];
            arInfo[iSub].nxs += elemCG[el][0]*elemCG[el][0];
            arInfo[iSub].nys += elemCG[el][1]*elemCG[el][1];
            arInfo[iSub].nzs += elemCG[el][2]*elemCG[el][2];
        }
        arInfo[iSub].nNd = dec->num(iSub);
        cost +=
            (arInfo[iSub].nxs + arInfo[iSub].nys + arInfo[iSub].nzs -
             (arInfo[iSub].ncgx*arInfo[iSub].ncgx +
              arInfo[iSub].ncgy*arInfo[iSub].ncgy +
              arInfo[iSub].ncgz*arInfo[iSub].ncgz)/arInfo[iSub].nNd)/arInfo[iSub].nNd;
    }
    return cost;
}

Decomposition *
MultiFront::checkComponents(Decomposition *dec)
{
    int iSub, iEle, iNode, i;
    int np = nToE->numConnect();
    NodeWeight initW;
    for(i = 0; i < np; ++i)
        nw[i] = initW;

    for(iNode = 0; iNode < nToE->csize(); ++iNode) {
        nSubPerNode[iNode] = 0;
        boundIndex[iNode] = -1;
        nodeMask[iNode] = 0;
        for(iEle = 0; iEle < nToE->num(iNode); ++iEle)
            if((*elems)[(*nToE)[iNode][iEle] ]->isStart())
                nodeMask[iNode]++;
    }
    numBoundNodes = 0;

    int *oldSub = new int[numEle];
    for(iEle = 0; iEle < numEle; ++iEle) {
        oldSub[iEle] = assignedSubD[iEle];
        assignedSubD[iEle] = -1;
    }
    int cSub = 0;
    for(iSub =0; iSub < dec->nsub; ++iSub) {
        for(iEle = 0; iEle < dec->num(iSub); ++iEle)
            if(assignedSubD[(*dec)[iSub][iEle] ] < 0) {
                int startElem = (*dec)[iSub][iEle];
                if((*elems)[ startElem ]->isStart() == false) continue;
                prio.clear();
#ifndef NOTMPL
                lastPrio.clear();
#endif
                addElemToSub(cSub, startElem);
                while(prio.empty() == false) {
#ifdef NOTMPL
                    int elem = prio.pop();
                    if(oldSub[elem] == iSub && assignedSubD[elem] < 0)
                        addElemToSub(cSub, elem);
#else
                    map<PrioState,int>::iterator it = prio.begin();
         int elem = it->second;
         lastPrio.erase(elem);
         prio.erase(it);
         if(oldSub[elem] == iSub && assignedSubD[elem] < 0)
           addElemToSub(cSub, elem);
#endif
                }
                cSub++;
            }
    }
    for(iEle = 0; iEle < numEle; ++iEle) // check for rogue non start elements
        if(oldSub[iEle] >= 0 && assignedSubD[iEle] < 0)
            assignedSubD[iEle]=cSub++;
    if(cSub != dec->nsub) {
        delete [] dec->pele;
        dec->pele = new int[cSub+1];
        dec->nsub = cSub;
        for(iSub = 0; iSub <= cSub; ++iSub)
            dec->pele[iSub] = 0;
        for(iEle = 0; iEle < numEle; ++iEle)
            if((*elems)[iEle]) {
                if(assignedSubD[iEle] < 0) {
                    std::cerr << "Element " << iEle+1 << "Is not assigned."
                                                         "it will not be in the decomposition\n" ;
                    continue;
                }
                dec->pele[assignedSubD[iEle]]++;
                if(assignedSubD[iEle] >= cSub || assignedSubD[iEle] < 0)
                    fprintf(stderr, "Bug %d %d\n", assignedSubD[iEle], cSub);
            }

        for(iSub = 0; iSub < cSub; ++iSub)
            dec->pele[iSub+1] += dec->pele[iSub];

        for(iEle = 0; iEle < numEle; ++iEle)
            if((*elems)[iEle]) {
                if(assignedSubD[iEle] < 0) continue;
                dec->eln[ --dec->pele[assignedSubD[iEle]] ] = iEle;
            }
    }
    delete [] oldSub;
    return dec;
}

double
MultiFront::balFct(double ideal, double invIdeal, int num) {
    double diffRatio = (num-ideal)*invIdeal;
    double ratio = num*invIdeal;
    double cost = 1.0 + 16.0*diffRatio*diffRatio;
#ifdef NONLINEAR
    if(diffRatio > 0)
    cost += 128*diffRatio*diffRatio*diffRatio*diffRatio;
#endif
    if(ratio < 0.5)
        cost += 1/ratio-2.0;
    return cost;
}

// double oc1 = 0.4, oc2 = 1.0, oc3 = 0.2;
double oc1 = 0.4, oc2 = 1.5, oc3 = 0.2;
// Balance, AR, bound
double dc1=0, dc2=0, dc3=0;
double deltaBound, deltaAR, deltaBal;

Decomposition *
MultiFront::optimize(Decomposition *origDec)
{
    int iNode, iEle, iSub;

    int elemFlag = getEventTag();
    BoundList boundList(numEle);
    // number of boundary nodes touching to this element.
    int *nBNode = new int[numEle];
    //int numNodes = nToE->csize();
    int numSub = origDec->nsub;
    double xDiff = 0;
    for(iEle = 0; iEle < numEle; ++iEle)
        nBNode[iEle] = 0;
    double balCost = 0;
    // build the balance cost
    balCost = computeBalCost(origDec);
    // balCost cannot be zero, because it starts with a one, so it can
    // only be zero if there are no subdomains
    double balCoef = oc1*1.0/balCost;

    //fprintf(stderr, "Balance cost: %e, numBoundNodes %d\n", balCost, numBoundNodes);

    // Compute the cost function and normalize it
    // if no boundary nodes are found, we cannot improve the decomposition
    if(numBoundNodes == 0) { delete [] nBNode; return origDec; }
    boundCoef = oc3/numBoundNodes;
    arInfo =  new ARInfo[numSub];
    double arCost = buildARInfo(origDec);
    // a zero arCost can truly only happen if each subdomain is a single element
    if(arCost == 0) arCoef = oc2;
    else arCoef = oc2/arCost;

    // Loop through the nodes to gather the boundary elements
    for(iNode = 0; iNode < numBoundNodes; ++iNode) {
        int node = boundNode[iNode];
        for(iEle = 0; iEle < nToE->num(node); ++iEle) {
            int elem = (*nToE)[node][iEle];
            nBNode[elem] += 1;
            if(flag[elem] != elemFlag) {
                boundList.add(elem);
                flag[elem] = elemFlag;
            }
        }
    }

    // For SA compute a reasonable starting temperature based on average cost
    double T0 = 0.0;
    int n = 0;
    for(iEle = 0; iEle < boundList.size(); ++iEle) {
        int elem = boundList[iEle];
        ARInfo resAR;
        int oldSub = assignedSubD[elem];
        ARInfo curAR = arInfo[oldSub];
        int newSub;
        if(curAR.nNd == 1) continue;
        double dc = findBestDelta(elem, newSub, resAR, curAR);
        if(newSub < 0) continue;
        dc += balCoef*(computeBalDiff(oldSub,-1,elem)
                       + computeBalDiff(newSub,+1,elem));
        n++;
        T0 += dc >= 0 ? dc : -dc;
    }
    if(n > 0) T0 = T0/n;
    else { delete [] nBNode; return origDec; }
    // Iterate through the boundary elements to try to flip them to another
    // subdomain.
    int nIter = 80*boundList.size();
    int idleCount = 0;
    // adjust the starting temperature to go toward zero when the number
    // of boundary elements increases
    double coef = double(nIter)/1e7;
    if(coef < 1) coef = 1.0;
    if(coef > 1) coef = exp(-10*(coef-1)*(coef-1));
    coef = 1.0;
    double T = T0*coef;
    double rT = 1.0-1.0/8;// 1.0-1.0/20;
    int tChange = boundList.size();
    if(tChange > 150000 || nosa ) {
        if(verboseFlag) filePrint(stderr, " ...     Using a Deterministic Scheme ...\n");
        T = 0.0;
    }
    else {
        if(verboseFlag) filePrint(stderr, " ...     Using a Simulated Annealing Scheme ...\n");
        srand48(45302172); // will make sure that XPost's decomposer and fem's dec agree ! Thomas
    }

    SARule sa;
    for(int iter = 0; iter < nIter; ++iter) {
        if((iter+1)%tChange == 0) T *= rT;
        if(idleCount == boundList.size()) break;
        idleCount++;
        // randomly pick an element
        int elem = boundList[iter % boundList.size()];
        ARInfo resAR;
        int oldSub = assignedSubD[elem];
        // avoid obtaining an empty subdomain
        if(elemPerSub[oldSub] == 1) continue;

        ARInfo curAR = arInfo[oldSub];
        int newSub;
        double dc = findBestDelta(elem, newSub, resAR, curAR);
        if(newSub < 0) continue;
        double tmp;
        dc += (tmp = balCoef*(computeBalDiff(oldSub,-1,elem)
                              + computeBalDiff(newSub,+1,elem)));
        deltaBal = tmp;
        if(sa.accept(dc,T)) {  // perform the swap
            double thisWeight = (*elems)[elem]->weight();
            subWeight[oldSub] -= thisWeight;
            subWeight[newSub] += thisWeight;
            dc1 += deltaBal;
            dc2 += deltaAR;
            dc3 += deltaBound;
            idleCount = 0;
            xDiff += dc;
            // remove the element from its currently assigned subdomain
            arInfo[oldSub] = curAR;
            elemPerSub[oldSub] -= 1;
            for(iNode = 0; iNode < eToN->num(elem); ++iNode) {
                int node = (*eToN)[elem][iNode];
                if(decWeight(oldSub, node, elemHasRot[elem] ? 1 : 0) == 0) {
                    // check if this node is disappearing from the boundary
                    if(numConnectedSub(node) == 1) {
                        removeBoundNode(node);
                        for(iEle = 0; iEle < nToE->num(node); ++iEle) {
                            int eleI = (*nToE)[node][iEle];
                            nBNode[eleI]--;
                            if(nBNode[eleI] == 0) {
                                // remove this element from the boundary list
                                boundList.remove(eleI);
                            }
                        }
                    }
                }
            }
            // add the element to the new subdomain
            assignedSubD[elem] = newSub;
            elemPerSub[newSub] += 1;
            arInfo[newSub] = resAR;
            for(iNode = 0; iNode < eToN->num(elem); ++iNode) {
                int node = (*eToN)[elem][iNode];
                if(incrWeight(newSub, node, elemHasRot[elem] ? 1 : 0) == 1) {
                    // check if this node is new on the boundary
                    if(numConnectedSub(node) == 2) {
                        addBoundNode(node);
                        for(iEle = 0; iEle < nToE->num(node); ++iEle) {
                            int eleI = (*nToE)[node][iEle];
                            nBNode[eleI]++;
                            if(nBNode[eleI] == 1) {
                                // add this element to the boundary list
                                boundList.add(eleI);
                            }
                        }
                    }
                }
            }
        } // end of swaping !
    }
    for(iSub = 0; iSub <= numSub; ++iSub)
        origDec->pele[iSub] = 0;
    for(iEle = 0; iEle < numEle; ++iEle)
        if((*elems)[iEle]) {
            if(assignedSubD[iEle] < 0) {
                std::cerr << "Element " << iEle+1 << "Is not assigned. it will not be in"
                                                     "the decomposition\n";
                continue;
            }
            origDec->pele[assignedSubD[iEle]]++;
        }

    for(iSub = 0; iSub < numSub; ++iSub)
        origDec->pele[iSub+1] += origDec->pele[iSub];

    for(iEle = 0; iEle < numEle; ++iEle)
        if((*elems)[iEle]) {
            if(assignedSubD[iEle] < 0)
                continue;
            origDec->eln[ --origDec->pele[assignedSubD[iEle]] ] = iEle;
        }
    arCost = buildARInfo(origDec);

    delete [] nBNode;
    return origDec;
}

double
MultiFront::add(ARInfo &ari, int elem)
{
    double cCost = (ari.nxs + ari.nys + ari.nzs -
                    (ari.ncgx*ari.ncgx +
                     ari.ncgy*ari.ncgy +
                     ari.ncgz*ari.ncgz)/ari.nNd)/ari.nNd;
    ari.nxs += elemCG[elem][0]*elemCG[elem][0];
    ari.nys += elemCG[elem][1]*elemCG[elem][1];
    ari.nzs += elemCG[elem][2]*elemCG[elem][2];
    ari.ncgx += elemCG[elem][0];
    ari.ncgy += elemCG[elem][1];
    ari.ncgz += elemCG[elem][2];
    ari.nNd += 1;
    double nCost = (ari.nxs + ari.nys + ari.nzs -
                    (ari.ncgx*ari.ncgx +
                     ari.ncgy*ari.ncgy +
                     ari.ncgz*ari.ncgz)/ari.nNd)/ari.nNd;
    return nCost-cCost;
}

double
MultiFront::remove(ARInfo &ari, int elem)
{
    double cCost = (ari.nxs + ari.nys + ari.nzs -
                    (ari.ncgx*ari.ncgx +
                     ari.ncgy*ari.ncgy +
                     ari.ncgz*ari.ncgz)/ari.nNd)/ari.nNd;
    ari.nxs -= elemCG[elem][0]*elemCG[elem][0];
    ari.nys -= elemCG[elem][1]*elemCG[elem][1];
    ari.nzs -= elemCG[elem][2]*elemCG[elem][2];
    ari.ncgx -= elemCG[elem][0];
    ari.ncgy -= elemCG[elem][1];
    ari.ncgz -= elemCG[elem][2];
    ari.nNd -= 1;
    double nCost = (ari.nxs + ari.nys + ari.nzs -
                    (ari.ncgx*ari.ncgx +
                     ari.ncgy*ari.ncgy +
                     ari.ncgz*ari.ncgz)/ari.nNd)/ari.nNd;
    return nCost-cCost;
}

double
MultiFront::findBestDelta(int elem, int &bestSub, ARInfo &resAR, ARInfo &curAR)
{
    int iNode, jNode, ip;
    int tag = getEventTag();
    // cost benefit of removing this element from its current subdomain
    double leaveCostDiff = 0;
    int curSub = assignedSubD[elem];
// curAR = arInfo[curSub];
    double minDiff =-10000000;
    bestSub = -1;
    deltaBound = 0.0;
    deltaAR = 0.0;
#ifndef NODEBASEDARCOST
    double tmp;
    leaveCostDiff += (tmp = arCoef*remove(curAR, elem));
    deltaAR += tmp;
#endif
    double addAR, addBound, minAR = 0, minBound = 0;
    for(iNode = 0; iNode < eToN->num(elem); ++iNode) {
        int node = (*eToN)[elem][iNode];
        // compute the cost of removing the element from its current subdomain
        int curWeight = weight(curSub, node);
        if(curWeight == 1) {
            if(nSubPerNode[node] > 1) {
                leaveCostDiff -= boundCoef;
                deltaBound -= boundCoef;
            }
#ifdef NODEBASEDARCOST
            double tmp;
      leaveCostDiff += (tmp = arCoef*remove(curAR, node));
      deltaAR += tmp;
#endif
        }
        // look for neighbors and compute the cost of adding it to them
        int initP = nToE->offset(node);
        int finalP = initP+nSubPerNode[node];
        for(ip = initP; ip < finalP; ++ip) {
            int candidate = nw[ip].subd;
            if(candidate == curSub) continue;
            addAR = addBound = 0;
            if(flag[candidate] != tag) {
                flag[candidate] = tag;
                // check if this subdomain can accept this element
                if((*elems)[elem]->examine(candidate, this).isReady == false) {
                    continue;
                }

                ARInfo candARI = arInfo[candidate];
                double addDiff = 0;
                for(jNode = 0; jNode < eToN->num(elem); ++jNode) {
                    int nodeJ = (*eToN)[elem][jNode];
                    if(weight(candidate, nodeJ) == 0) {
                        if(nToE->num(nodeJ) > 1) { // there are other subdomains touching
                            addDiff += boundCoef;  // this node
                            addBound += boundCoef;
                        }
#ifdef NODEBASEDARCOST
                        double tmp;
              addDiff += (tmp =  arCoef*add(candARI, nodeJ));
              addAR += tmp;
#endif
                    }
                }
#ifndef NODEBASEDARCOST
                addDiff += (tmp = arCoef*add(candARI, elem));
                addAR += tmp;
#endif
                if(bestSub < 0 || addDiff < minDiff) {
                    bestSub = candidate;
                    minDiff = addDiff;
                    resAR = candARI;
                    minAR = addAR; minBound = addBound;
                }
            }
        }
    }

    deltaBound += minBound;
    deltaAR += minAR;
    if(minDiff == -10000000) return leaveCostDiff; else // PJSA DEBUG
        return leaveCostDiff + minDiff;
}

void
MultiFront::memEstimate(Decomposition *dec, int dofsPerNode, double memory[5],
                        FILE *memFile)
{
    memory[1] = memory[2] = memory[3] = memory[4] = 0.0;
    double *memSize = new double[dec->nsub];
    double minSize = -1;
    double maxSize = 0;
    double totSize = 0.0;
    int iNode, jNode, iEle, iSub;
    int numNodes = nToE->csize();
    double cost = 0;
    int dualSize = 0;
    rebuildInfo(dec);
    int repNodes = numNodes;
    // Node coordinates cost:
    cost += 3*sizeof(double)*numNodes;
    for(iNode = 0; iNode < numBoundNodes; ++iNode) {
        int node = boundNode[iNode];
        int ns = numSubForNode(node);
        dualSize += (ns*(ns-1))/2;
        repNodes += ns-1;
    }
    // Add reorthogonalization cost for 100 vectors:
    cost += sizeof(double)*2*dofsPerNode*dualSize*100;
    memory[4] = sizeof(double)*2*dofsPerNode*dualSize*100;
    // GRRBM memory estimate:
    cost += 4*sizeof(double)*dofsPerNode*repNodes*6;
    for(iNode = 0; iNode < numNodes; ++iNode)
        nodeMask[iNode] = -1;
    int matMem = 0;
    // make the sub to sub connectivity:
    int *subMask = new int[dec->nsub];
    int *ccount = new int[dec->nsub+1];
    for(iSub = 0; iSub < dec->nsub; ++iSub) {
        ccount[iSub] = 0;
        subMask[iSub] = -1;
    }
    for(iSub = 0; iSub < dec->nsub; ++iSub)
    {
        for(iEle = 0; iEle < dec->num(iSub); ++iEle)
        {
            int elem = (*dec)[iSub][iEle];
            for(iNode = 0; iNode < eToN->num(elem); ++iNode)
            {
                int node = (*eToN)[elem][iNode];
                for(int ip = nToE->offset(node); ip < nToE->offset(node+1)
                                                 && nw[ip].subd < 0; ++ip)
                {
                    if(subMask[ nw[ip].subd ] != iSub)
                    {
                        subMask[ nw[ip].subd ] = iSub;
                        ccount[iSub] ++;
                    }
                }
            }
        }
    }
    // pass two
    ccount[dec->nsub] = 0;
    for(iSub = 0; iSub < dec->nsub; ++iSub)
    {
        ccount[iSub+1] += ccount[iSub];
        subMask[iSub] = -1;
    }
    int *subPtr = new int[ccount[dec->nsub]];
    for(iSub = 0; iSub < dec->nsub; ++iSub)
    {
        for(iEle = 0; iEle < dec->num(iSub); ++iEle)
        {
            int elem = (*dec)[iSub][iEle];
            for(iNode = 0; iNode < eToN->num(elem); ++iNode)
            {
                int node = (*eToN)[elem][iNode];
                for(int ip = nToE->offset(node); ip < nToE->offset(node+1)
                                                 && nw[ip].subd < 0; ++ip)
                {
                    if(subMask[ nw[ip].subd ] != iSub)
                    {
                        subMask[ nw[ip].subd ] = iSub;
                        subPtr[--ccount[iSub]] = nw[ip].subd;
                    }
                }
            }
        }
    }
    Connectivity *subToSub = new Connectivity(dec->nsub, ccount, subPtr);
    compStruct subRenum = subToSub->renumByComponent(1);
    //int profileSize = subToSub->findProfileSize(subRenum.renum);
    int profileSize = 0;

    memory[3] = dofsPerNode*dofsPerNode*sizeof(double)*profileSize
                +sizeof(int)*(subToSub->numConnect() + dec->nsub);
    cost += memory[3];

    int *nodeMap = new int[numNodes];
    int *locToGl = new int[numNodes];
    int *locDofsPerNode = new int[numNodes];
    for(iSub = 0; iSub < dec->nsub; ++iSub)
    {
        int kbbSize = 0;
        int ndRenum = 0;
        int nNd = 0;
        int *ptr = new int[dec->num(iSub) + 1];
        ptr[0] = 0;
        for(iEle = 0; iEle < dec->num(iSub); ++iEle) {
            int elem = (*dec)[iSub][iEle];
            nNd += eToN->num(elem);
        }
        int *target = new int[nNd];
        int tgIndex = 0;
        for(iEle = 0; iEle < dec->num(iSub); ++iEle)
        {
            int elem = (*dec)[iSub][iEle];
            for(iNode = 0; iNode < eToN->num(elem); ++iNode)
            {
                int node = (*eToN)[elem][iNode];
                if(node < 0) fprintf(stderr, " *** ERROR: node < 0\n");
                if(nodeMask[node] != iSub)
                {
                    nodeMask[node] = iSub;
                    locToGl[ndRenum] = node;
                    locDofsPerNode[ndRenum] = (rotWeight(iSub, node) > 0) ? 6 : 3;
                    nodeMap[node] = ndRenum++;
                }
                target[tgIndex++] = nodeMap[node];
            }
            ptr[iEle+1] = tgIndex;
        }
        Connectivity eToN(dec->num(iSub), ptr, target);
        Connectivity *nToE = eToN.alloc_reverse();
        Connectivity *nToN = nToE->transcon(&eToN);
        // Find Kbb cost:
        for(iNode = 0; iNode <ndRenum; ++iNode)
        {
            int nodeI = locToGl[iNode];
            if(numSubForNode(nodeI) == 1) continue;
            for(jNode = 0; jNode < nToN->num(iNode); ++jNode)
            {
                int nodeJ = locToGl[(*nToN)[iNode][jNode]];
                if(numSubForNode(nodeJ) > 1) kbbSize++;
            }
        }

        memory[1] += (sizeof(double)+sizeof(int))*kbbSize*dofsPerNode*dofsPerNode;
        cost+=(sizeof(double)+sizeof(int))*kbbSize*dofsPerNode*dofsPerNode;

        compStruct renum = nToN->renumByComponent(1);

        DofSetArray dsa(ndRenum, locDofsPerNode, renum.renum);
        //int profileSize = nToN->findProfileSize(renum.renum);

        long profileSize = nToN->findProfileSize(&dsa);

        memSize[iSub] = 8.0*profileSize/(1024*1024);
        if(minSize < 0 || memSize[iSub] < minSize)
            minSize = memSize[iSub];
        if(memSize[iSub] > maxSize)
            maxSize = memSize[iSub];
        totSize += memSize[iSub];

        //matMem += profileSize*sizeof(double)*dofsPerNode*dofsPerNode;
        // matrix integer arrays costs
        cost += sizeof(int)*3*nToE->csize();
        //cost += sizeof(double)*dofsPerNode*dofsPerNode*profileSize;
        cost += sizeof(double)*profileSize;
        // connectivities costs including for master processing
        cost += (3*(eToN.numConnect()+eToN.csize()+nToE->numConnect()+nToE->csize())+
                 nToN->csize()+nToN->numConnect()+3)*sizeof(int);
        // element costs:
        cost += eToN.csize()*(46*sizeof(int)+sizeof(void *));
        // renum array cost and dsa, c_dsa
        cost += 4*sizeof(int)*nToE->csize();
        if(nToE)
            delete nToE;
        if(nToN)
            delete nToN;
        // if(renum.renum)
        //delete [] renum.renum;
    }

    delete [] nodeMap;
    memory[2] = matMem;
    memory[1] += matMem;
    memory[0] = cost/(1024*1024.0);
    memory[1] /= (1024*1024.0);
    memory[2] /= (1024*1024.0);
    memory[3] /= (1024*1024.0);
    memory[4] /= (1024*1024.0);
    // nomre_sous_domaine  min  ave  max  max/ave
    double aveMem = totSize/dec->nsub;
    fprintf(memFile, "# Num_Sub   Min_Mem   Ave_Mem   Max_Mem   LBF_Mem = Max_Mem/Ave_Mem\n");
    fprintf(memFile, "# %d  %f %f %f %f\n", dec->nsub, minSize,
            aveMem, maxSize, maxSize/aveMem);
    for(iSub = 0; iSub < dec->nsub; ++iSub)
        fprintf(memFile, "%f\n", memSize[iSub]);

    delete[] memSize;
}
int
MultiFront::rotWeight(int sub, int node)
{
    int ip;
    int initP = nToE->offset(node);
    int finalP = initP+nSubPerNode[node];
    for(ip = initP; ip < finalP; ++ip) {
        if(nw[ip].subd == sub) return nw[ip].rotWeight;
    }
    return 0;
}

int
MultiFront::bestSubFor(int nNd, int *nd)
{
    int iNode;
    for(iNode = 0; iNode < nNd; ++iNode)
        if(nToE->num(nd[iNode]) > 0)
            return assignedSubD[ (*nToE)[nd[iNode]][0] ];
    return -1;
}

double
MultiFront::computeBalCost(Decomposition *dec)
{
    int iEle,iSub, numSub = dec->nsub;
    double totWeight = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub)
        for(iEle = 0; iEle < dec->num(iSub); ++iEle)
            totWeight += (*elems)[ (*dec)[iSub][iEle] ]->weight();
    avgWeight = totWeight/numSub;
    invAvgWeight = numSub/totWeight;
    if(subWeight) delete [] subWeight;
    subWeight = new double[numSub];
    double cost = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub) {
        subWeight[iSub] = 0;
        for(iEle = 0; iEle < dec->num(iSub); ++iEle)
            subWeight[iSub] += (*elems)[ (*dec)[iSub][iEle] ]->weight();
        cost += balFct(avgWeight, invAvgWeight, int(subWeight[iSub]));

    }
    return cost;
}

double
MultiFront::computeBalDiff(int iSub, int sign, int ele)
{
    double cost = 0.0;
    cost -= balFct(avgWeight, invAvgWeight, int(subWeight[iSub]));
    double newWeight = subWeight[iSub]+ sign*(*elems)[ele]->weight();
    cost += balFct(avgWeight, invAvgWeight, int(newWeight));
    return cost;
}

