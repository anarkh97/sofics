#ifndef _MULTI_FRONT_H_
#define _MULTI_FRONT_H_

#define NOTMPL
#include <map>

#include <cstdio>

#include <Element.d/Element.h>
#include <Utils.d/Connectivity.h>                                               

struct PrioState { 
   int prio;
   int time;
   bool operator < (const PrioState &o) const
     { return(prio < o.prio || (prio == o.prio && time < o.time)); }
};

struct NodeWeight {
    int subd;
    int weight;
    int rotWeight;
    NodeWeight() { subd = -1; weight = 0; rotWeight = 0;}
};

class Decomposition;

class BoundList {
    int maxIndex;
    int nObj;
    int *listIndex;
    int *list;
 public:
    BoundList(int size);
    ~BoundList();
    void add(int n);
    void remove(int n);
    int size() { return nObj; }
    int operator[](int i) { return list[i]; }
    bool isPresent(int i) { return listIndex[i] >= 0; }
};

class OrderList {
    int begin, end, length;
    int infBegin, infEnd;
    int *listIndex;
    int *list;
    int *infList;
  public:
    OrderList(int s);
    ~OrderList();
    int pop();
//    void remove(int n);
    void add(int n);
    void addInfinity(int n);
    int popInfinity();
    void clear();
    bool empty() { return begin == end && infBegin == infEnd; }
};

class MultiFront {

        Elemset *elems;
	CoordSet *nds;
        double (*elemCG)[3];
	Connectivity *nToE;
	Connectivity *eToN;
        bool *elemHasRot;
        bool *isFsiConnected;
        int *assignedSubD;
        int *flag; // for flagging a single visit of nodes or elements
        int *nodeMask, *boundNode, *boundIndex;
	int numBoundNodes;
	int count;
        int tag;
        int numEle;
        int time;
        int desiredNSub;
	int iterationN;
#ifdef NOTMPL
	OrderList prio;
#else
        // priority of elements in subdomains
        map<int,PrioState> lastPrio;
        // now sorted per priority
        map<PrioState, int> prio;
#endif
        // data structure for node weights:
        NodeWeight *nw;
        int *nSubPerNode;
        int *elemPerSub;
        double *subWeight;
        double invAvgWeight;
        double avgWeight;

        double arCoef, boundCoef, elemCoef;

        struct ARInfo {
           double ncgx, ncgy, ncgz;
           double nxs, nys, nzs;
           int nNd;
           ARInfo() { ncgx = ncgy = ncgz = nxs = nys = nzs = 0.0; nNd = 0; }
        };

        ARInfo * arInfo;

//        std::multimap<int,int> fsAffinity;
        int fsGluedCounter;
        Elemset* fsGlued_eset;
        int fsGlFlag;

        double add(ARInfo &, int node);
        double remove(ARInfo &, int node);
        void rebuildInfo(Decomposition *dec);
        void remapAll(int *subRemap);
        int getEventTag() { return tag++; }
        void addNodesToSub(int,int, int, int *);
        void addElemToSub(int, int);
        void updateElement(int,int, int);
        void initDec(int);
	void addBoundNode(int);
	void removeBoundNode(int);
        double findBestDelta(int elem, int &sub, ARInfo &res, ARInfo &cur);
        double buildARInfo(Decomposition *);
        double balFct(double bal, double invbal, int nEle);
        double computeBalCost(Decomposition *);
        Decomposition *optimize(Decomposition *);
        Decomposition *improveDec(Decomposition *, double threshCoef);
	Decomposition *checkComponents(Decomposition *dec);
        Decomposition *removeFsiElements(Decomposition *);
	
// JLchange
        Decomposition *updateFsiConnection(Decomposition *); 
        Decomposition *eliminateEmpty(Decomposition *); 

	int numSubForNode(int node);
  public:
        MultiFront(Elemset *eSet, CoordSet *nds = 0, bool have_fsi = false, bool _fsGlFlag = false);
        ~MultiFront();
        Decomposition * decompose(int nsub, bool have_fsi = false);
        // Weight of a node in a subdomain
        int weight(int sub, int node);
        int rotWeight(int sub, int node);
        int incrWeight(int, int, int rMode=0);
        int decWeight(int, int, int rMode=0);
        int numConnectedSub(int node);
        // Total weight of a node
        int weight(int);
        Connectivity *getNToE() { return nToE; }
	void memEstimate(Decomposition *, int dofsPerNode, 
              double me[5], FILE *);
        int bestSubFor(int nNd, int *nd) ;

        void redoConnect();
        double computeBalDiff(int iSub, int sign, int ele);
};
#endif
