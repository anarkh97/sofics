#include <Hetero.d/FlExchange.h>
#include <Hetero.d/FilteredFile.h>
#include <Utils.d/dofset.h>
#include <Element.d/State.h>

#include <Comm.d/Communicator.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/SolverInfo.h>
#include <Corotational.d/GeomState.h>

#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>

// std
#include <cstdlib>
#include <algorithm>
#include <map>
using std::map;
using std::pair;
typedef pair<int, pair<double, double> >  locoord;
//         elem id      xi1     xi2

#define DETERMINISTIC_FLEXCHANGER

const char *RECEIVE_LIST_KW = "RCVF";
const char *SUBDOMAIN_KW = "SUBD";
const char *SEND_LIST_KW = "SNDF";

extern Communicator *structCom, *fluidCom, *heatStructCom;
extern int verboseFlag;
extern GeoSource *geoSource;
extern SolverInfo &solInfo;
extern Domain *domain;

FlExchanger::FlExchanger(CoordSet &_cs, Elemset &_eset, SurfaceEntity *_surface, DofSetArray *_dsa,
                         OutputInfo *_oinfo, bool _wCracking)
  : cs(_cs), eset(_eset)
{
  dsa     = _dsa;
  oinfo   = _oinfo;
  tmpDisp = 0;
  tmpVel  = 0;
  useFaceElem = true;
  wCracking = _wCracking;
  sentInitialCracking = false;

  if(wCracking) {
    // The deep copy of the surface entity is required because face elements may be removed from
    // the surface due to element deletion. However, this class currently requires access to the
    // original surface topology.
    surface = new SurfaceEntity(*_surface);
    int nElems = (surface->GetIsShellFace() && domain->tdenforceFlag()) ? surface->nFaceElements()/2 : surface->nFaceElements();
    faceElemToNode = new Connectivity(SetAccess<FaceElemSet>{*surface->GetPtrFaceElemSet(), nElems});
    nodeToFaceElem = faceElemToNode->alloc_reverse();
  }
  else {
    surface = _surface;
    faceElemToNode = 0;
    nodeToFaceElem = 0;
  }

  idSendTo = 0;
  nbSendTo = 0;
  consOrigin = 0;
  localF = 0;
  buffer = 0;
  sndTable = 0;
}

FlExchanger::FlExchanger(CoordSet &_cs, Elemset &_eset, DofSetArray *_dsa,
                         OutputInfo *_oinfo) : cs(_cs), eset(_eset), surface(0)
{
  dsa      = _dsa;
  oinfo    = _oinfo;
  tmpDisp  = 0;
  tmpVel   = 0;
  useFaceElem = false;
  wCracking = false;
  sentInitialCracking = false;
  faceElemToNode = 0;
  nodeToFaceElem = 0;
  idSendTo = 0;
  nbSendTo = 0;
  consOrigin = 0;
  localF = 0;
  buffer = 0;
  sndTable = 0;
}

FlExchanger::~FlExchanger()
{
  if(tmpDisp) {
    delete tmpDisp;
  }

  if(tmpVel) {
    delete tmpVel;
  }

  if(wCracking && surface) {
    delete surface;
  }

  if(faceElemToNode) {
    delete faceElemToNode;
  }

  if(nodeToFaceElem) {
    delete nodeToFaceElem;
  }

  if(idSendTo) {
    delete [] idSendTo;
  }

  if(nbSendTo) {
    delete [] nbSendTo;
  }

  if(consOrigin) {
    delete [] consOrigin;
  }

  if(localF) {
    delete [] localF;
  }

  if(buffer) {
    delete [] buffer;
  }

  if(sndTable) {
    delete [] sndTable;
  }
}

double
FlExchanger::getFluidLoad(Vector &force, int tIndex, double time,
                          double alphaf, int &iscollocated, GeomState *geomState)
{
  aforce[0] = aforce[1] = aforce[2] = 0.0;
  FaceElemSet *feset;

  if(useFaceElem) {
    feset = &(surface->GetFaceElemSet());
  }

  Element     *thisElement;
  FaceElement *thisFaceElem;

#ifdef DETERMINISTIC_FLEXCHANGER
  // Message-passing programming models are by default nondeterministic: the arrival order of messages
  // sent from two processes, A and B, to a third process, C, is not defined. It is the programmer's
  // responsibility to ensure that a computation is deterministic when this is required. 
  std::vector<double> *buffers = new std::vector<double>[nbrReceivingFromMe];
  std::map<int,int> cpupos;
  for(int i = 0; i < nbrReceivingFromMe; i++) {
    int tag = FLTOSTMT + ((rcvParity > 0) ? 1 : 0);
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
    cpupos[rInfo.cpu] = i;
    buffers[i].insert(buffers[i].begin(), buffer, buffer + nbSendTo[consOrigin[rInfo.cpu]]*3);
  }

  for(std::map<int,int>::iterator it = cpupos.begin(); it != cpupos.end(); ++it) {
    double *buffer = buffers[it->second].data();
    int fromNd = it->first;
    int origin = consOrigin[fromNd];
#else
  for(int i = 0; i < nbrReceivingFromMe; i++) {
    int tag = FLTOSTMT + ((rcvParity > 0) ? 1 : 0);
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
    int fromNd = rInfo.cpu;
    int origin = consOrigin[fromNd];
#endif

    for(int j = 0; j < nbSendTo[origin]; ++j) {
      int nDof;

      if(!useFaceElem) {
        thisElement = eset[sndTable[origin][j].elemNum];
        thisElement->getFlLoad(cs, sndTable[origin][j], buffer + 3 * j, localF, geomState);
        nDof = thisElement->numDofs();
        transformVector(localF, thisElement);
      }
      else {
        thisFaceElem = (*feset)[sndTable[origin][j].elemNum];
        thisFaceElem->getFlLoad(sndTable[origin][j], buffer + 3 * j, localF);
        nDof = thisFaceElem->numDofs();
        transformVector(localF, thisFaceElem);
      }

      int *dof = sndTable[origin][j].dofs;

      for(int iDof = 0; iDof < nDof; ++iDof)
        if(dof[iDof] >= 0) {
          force[dof[iDof]] += localF[iDof];
        }

      aforce[0] += buffer[3 * j];
      aforce[1] += buffer[3 * j + 1];
      aforce[2] += buffer[3 * j + 2];
    }
  }

  flipRcvParity();

  if(oinfo) {
    if(tIndex % oinfo->interval == 0 && oinfo->filptr != NULL) {
      fprintf(oinfo->filptr, "%e   ", time);
      fprintf(oinfo->filptr, "%e %e %e\n", aforce[0], aforce[1], aforce[2]);
      fflush(oinfo->filptr);
    }
  }

#ifdef DETERMINISTIC_FLEXCHANGER
  delete [] buffers;
#endif

  // ML & KP For 'corrected' aeroelastic force
  iscollocated = (isCollocated) ? 1 : 0;
  return (time + alphaf * dt);
}

void
FlExchanger::sendDisplacements(State &state, int tag, GeomState *geomState)
{
  static int it = 0;

  if(tmpDisp == 0) {
    tmpDisp = new Vector(state.getDisp());
  }

  if(tmpVel == 0) {
    tmpVel = new Vector(state.getVeloc());
  }

  *tmpDisp = state.getDisp();
  tmpDisp->linAdd(dt * alpha[0], state.getVeloc(), dt * alpha[1], state.getPrevVeloc());

  *tmpVel = state.getVeloc();
  if (it > 0) {
    (*tmpVel) *= (1.0+alphasv);
    tmpVel->linAdd(-alphasv, state.getPrevVeloc());
  }

  State newState(state, *tmpDisp);
  if(verboseFlag) fprintf(stderr, "Disp Norm %e Veloc Norm %e\n", tmpDisp->norm(), state.getVeloc().norm());

  FaceElemSet *feset;
  int         *fnId;

  if(useFaceElem) {
    feset = &(surface->GetFaceElemSet());
    fnId  = surface->GetPtrGlNodeIds();
  }

  Element     *thisElement;
  FaceElement *thisFaceElem;
  double xxx = 0, yyy = 0;
  int pos = 0;

  for(int i = 0; i < nbrReceivingFromMe; i++) {
    int fluidNode = idSendTo[i];
    int origPos = pos;

    // compute displacements for each match point
    for(int j = 0; j < nbSendTo[i]; ++j) {
      if(!useFaceElem) {
        thisElement = eset[sndTable[i][j].elemNum];
        thisElement->computeDisp(cs, newState, sndTable[i][j], buffer + pos, geomState);
      }
      else {
        thisFaceElem = (*feset)[sndTable[i][j].elemNum];
        thisFaceElem->computeDisp(cs, newState, sndTable[i][j], buffer + pos, geomState, fnId);
      }

      // PJSA 4/14/2014: Now that the predictor is applied on the structure side for A6 and
      // A7, we need to compensate for the legacy predictor which is applied on the fluid side
      // (see MatchNodeSet::getDisplacement in MatchNodeCore.C).
      if(algnum == 6) {
        for(int k = 0; k < 3; ++k) {
          buffer[pos + k] -= 0.5 * dt * buffer[pos + k + 3];
        }
      }
      else if(algnum == 7) {
        for(int k = 0; k < 3; ++k) {
          buffer[pos + k] -= dt * buffer[pos + k + 3];
        }
      }

      pos += 6;
    }

    if(tag < 0) {
      tag = STTOFLMT + ((sndParity > 0) ? 1 : 0);
    }

    for(int ip = origPos; ip < pos; ip += 6) {
      xxx += buffer[ip + 0] * buffer[ip + 0] + buffer[ip + 1] * buffer[ip + 1] +
             buffer[ip + 2] * buffer[ip + 2];
      yyy += buffer[ip + 3] * buffer[ip + 3] + buffer[ip + 4] * buffer[ip + 4] +
             buffer[ip + 5] * buffer[ip + 5];
    }

    fluidCom->sendTo(fluidNode, tag, buffer + origPos, pos - origPos);
  }

  fluidCom->waitForAllReq();
  if (verboseFlag) fprintf(stderr, "Sending %e and %e to %d CPUs\n", xxx, yyy, nbrReceivingFromMe);
  flipSndParity();
}

void
FlExchanger::sendModeShapes(int numModes, int numN, double(**v)[6],
                            State &st, double ampFactor, int numIDis6, BCond* iDis6)
{
  int iMode;

  for(iMode = 0; iMode < numModes; ++iMode) {
    // Create the state
    Vector &d = st.getDisp();
    d.zero();
    int i, iNode, dof;
    for(i = 0; i < numIDis6; ++i) {
      dof = dsa->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
      if(dof >= 0)
        d[dof] += iDis6[i].val;
    }  

    for(iNode = 0; iNode < numN; ++iNode) {
      dof = dsa->locate(iNode, DofSet::Xdisp);

      if(dof >= 0) {
        d[dof] += v[iMode][iNode][0] * ampFactor;
      }

      dof = dsa->locate(iNode, DofSet::Ydisp);

      if(dof >= 0) {
        d[dof] += v[iMode][iNode][1] * ampFactor;
      }

      dof = dsa->locate(iNode, DofSet::Zdisp);

      if(dof >= 0) {
        d[dof] += v[iMode][iNode][2] * ampFactor;
      }

      dof = dsa->locate(iNode, DofSet::Xrot);

      if(dof >= 0) {
        d[dof] += v[iMode][iNode][3] * ampFactor;
      }

      dof = dsa->locate(iNode, DofSet::Yrot);

      if(dof >= 0) {
        d[dof] += v[iMode][iNode][4] * ampFactor;
      }

      dof = dsa->locate(iNode, DofSet::Zrot);

      if(dof >= 0) {
        d[dof] += v[iMode][iNode][5] * ampFactor;
      }
    }

    sendDisplacements(st, 1201 + iMode);
  }
}

void
FlExchanger::waitOnSend()
{
  fluidCom->waitForAllReq();
}

void
FlExchanger::sendParam(int _algnum, double step, double totaltime,
                       int rstinc, int _isCollocated, double _a[2], double _b)
{
  int TNd  = 0;
  int thisNode = structCom->myID();
  double buffer[5];
  buffer[0] = (double)((_algnum == 20 && solInfo.dyna3d_compat) ? 22 : _algnum);
  buffer[1] = (_algnum == 5) ? step / 2 : step;
  buffer[2] = totaltime;
  buffer[3] = (double) rstinc;
  buffer[4] = (double) 2; // Used to be AeroScheme now always conservative
  int msglen = 5;
  int tag = 3000;

  if(thisNode == 0) {
    fluidCom->sendTo(TNd, tag, buffer, msglen);
    fluidCom->waitForAllReq();
  }

  algnum = _algnum;
  isCollocated = _isCollocated;
  dt = step;
  /* PJSA 4/14/2014: see comments in Parser.d/p.y under AeroInfo
    if (algnum > 4)
      alpha[0] = alpha[1] = 0;
    else  {
      alpha[0] = _a[0];
      alpha[1] = _a[1];
    } */
  alpha[0] = _a[0];
  alpha[1] = _a[1];

  alphasv = _b;
}

void
FlExchanger::sendSubcyclingInfo(int sub)
{
  double buffer = (double)sub;
  fluidCom->sendTo(0, 777/*tag*/, &buffer, 1);
  fluidCom->waitForAllReq();
}

void
FlExchanger::sendTempParam(int _algnum, double step, double totaltime,
                           int rstinc, double alphat[2])
{
  int TNd  = 0;
  int thisNode;
  thisNode = structCom->myID();
  double buffer[5];
  buffer[0] = (double) _algnum;
  buffer[1] = step;
  buffer[2] = totaltime;
  buffer[3] = (double) rstinc;
  buffer[4] = (double) 2; // Used to be AeroScheme now always conservative
  int msglen = 5;
  int tag = 3000;

  if(thisNode == 0) {
    fluidCom->sendTo(TNd, tag, buffer, msglen);
    fluidCom->waitForAllReq();
  }

  algnum = _algnum;
  dtemp = step;
  alph[0] = alphat[0];
  alph[1] = alphat[1];
}

void
FlExchanger::sendModeFreq(double *modfrq, int nummod)
{
  int TNd  = 0;
  int ThisNode = structCom->myID();
  int Type = 1200;
  double *sBuffer = new double[1 + nummod];
  sBuffer[0]   = (double) nummod;
  int mod;

  for(mod = 0; mod < nummod; mod++) {
    sBuffer[1 + mod] = modfrq[mod];
  }

  int Pos  = 1 + nummod;

  if(ThisNode == 0) {
    fluidCom->sendTo(TNd, Type, sBuffer, Pos);
    fluidCom->waitForAllReq();
  }

  delete [] sBuffer;
}

//KW: send the embedded wet surface to fluid
void
FlExchanger::sendEmbeddedWetSurface()
{
  if(structCom->myID() != 0) {
    return;  /* do nothing */
  }

  if(!useFaceElem) {
    fprintf(stderr, "ERROR: Embedded Wet Surface undefined! Aborting...\n");
    exit(-1);
  }

  // info about surface
  FaceElemSet &feset = surface->GetFaceElemSet();
  CoordSet   &fnodes = surface->GetNodeSet();
  int          *fnId = surface->GetPtrGlNodeIds();
  map<int, int> *g2l = surface->GetPtrGlToLlNodeMap();
  int         nNodes = fnodes.size();
  int         nElems = (surface->GetIsShellFace() && domain->tdenforceFlag()) ? surface->nFaceElements()/2 : surface->nFaceElements();
  bool    renumbered = surface->IsRenumbered();
  int eType;

  if(wCracking) {
    // for element deletion, always use quads (convert triangles to degenerate quads)
    eType = 4;
  }
  else {
    // otherwise, convert triangles to degenerate quads only if there is at least one quad face.
    eType = 3;

    for(int i = 0; i < nElems; ++i) {
      if(feset[i]->GetFaceElemType() == FaceElement::QUADFACEL4) {
        eType = 4;
        break;
      }
    }
  }

  // data preparation
  int    buf[4] = {eType, (int)wCracking, nNodes, nElems};
  double nodes[3 * nNodes];
  int    nne = (eType == 3) ? 3 : 4; // number of nodes per element
  int    *elems = new int[nne * nElems];

  for(int i = 0; i < nNodes; i++) {
    nodes[3 * i]   = cs[fnId[i]]->x;
    nodes[3 * i + 1] = cs[fnId[i]]->y;
    nodes[3 * i + 2] = cs[fnId[i]]->z;
  }

  if(renumbered)
    for(int i = 0; i < nElems; i++) {
      FaceElement *ele = feset[i];

      for(int j = 0; j < ele->nNodes(); ++j) {
        elems[nne * i + j] = ele->GetNode(j);
      }

      if(ele->GetFaceElemType() == FaceElement::TRIFACEL3 && eType == 4) { // degenerate quad
        elems[nne * i + 3] = ele->GetNode(2);
      }
    }
  else
    for(int i = 0; i < nElems; i++) {
      FaceElement *ele = feset[i];

      for(int j = 0; j < ele->nNodes(); ++j) {
        elems[nne * i + j] = (*g2l)[ele->GetNode(j)];
      }

      if(ele->GetFaceElemType() == FaceElement::TRIFACEL3 && eType == 4) { // degenerate quad
        elems[nne * i + 3] = (*g2l)[ele->GetNode(2)];
      }
    }

  // send the package sizes
  fluidCom->sendTo(0, 555/*tag*/, buf, 4);
  fluidCom->waitForAllReq();

  if(solInfo.dyna3d_compat) {
    fluidCom->sendTo(0, 999/*tag*/, buf + 2, 2);
    fluidCom->waitForAllReq();
  }

  // send the node set
  fluidCom->sendTo(0, 666/*tag*/, nodes, nNodes * 3);
  fluidCom->waitForAllReq();
  // send the element set
  fluidCom->sendTo(0, 888/*tag*/, elems, nElems * nne);
  fluidCom->waitForAllReq();
  delete [] elems;
}

void
FlExchanger::read(int myNode, const char *inputFileName)
{
  // This function will either (a) read the matcher file or (b) query the embedded surface directly,
  // and then do some pre-processing to initialize various data structures.

  FilteredFile *ffile;
  char *cLine;
  FaceElemSet *feset;
  CoordSet *fnodes;
  int *fnId;
  map<int, locoord> exy;
  int numSnd, numRcv;

  if(!useFaceElem) {

    ffile = new FilteredFile(inputFileName);

    // Locate where the data pertaining to this node starts
    while(1) {
      if(ffile->findToken(SUBDOMAIN_KW) < 0) {
        fprintf(stderr, "Token %s not found. Aborting\n", SUBDOMAIN_KW);
        exit(1);
      }

      cLine = ffile->getLineAfterToken();
      fflush(stdout);
      int nodeNum;
      sscanf(cLine, "%d", &nodeNum);

      if(nodeNum - 1 == myNode) {
        break;
      }

      ffile->findToken("END");
    }

    // Let's look for the receive list
    if(ffile->findToken(RECEIVE_LIST_KW) < 0) {
      fprintf(stderr, "Token %s not found. Aborting\n", RECEIVE_LIST_KW);
      exit(1);
    }

    cLine = ffile->getLineAfterToken();
    sscanf(cLine, "%d", &numSnd);
  }
  else {
     // info about surface
     feset = surface->GetPtrFaceElemSet();
     fnodes = surface->GetPtrNodeSet();
     fnId = surface->GetPtrGlNodeIds();
     exy = feset->computeNodeLocalCoords(fnId, fnodes->size());

     numSnd = 0;
  }

  int *rcvcomid  = new int[numSnd];
  int *rcvcomlen = new int[numSnd];
  int **rcvEleList = new int *[numSnd];
  int **rcvGaussList = new int *[numSnd];
  int actualSenders = 0; // actual number of fluid nodes sending to me
  int maxSender = 0; // number associated with the highest actual sender
  int maxPRec = 0;   // number associated with the highest fluid process that receives from us

  for(int sender = 0; sender < numSnd; ++sender) {
    if(!useFaceElem) {
      cLine = ffile->getLine();
      sscanf(cLine, "%d%d", rcvcomid + sender, rcvcomlen + sender);
    }

    if(rcvcomlen[sender] > 0) {
      rcvEleList[sender]   = new int[rcvcomlen[sender]];
      rcvGaussList[sender] = new int[rcvcomlen[sender]];
      ++actualSenders;

      if(rcvcomid[sender] > maxSender) {
        maxSender = rcvcomid[sender];
      }
    }
    else {
      rcvEleList[sender]   = 0;
      rcvGaussList[sender] = 0;
    }

    for(int j = 0; j < rcvcomlen[sender]; ++j) {
      if(!useFaceElem) {
        cLine = ffile->getLine();
        sscanf(cLine, "%d%d", rcvEleList[sender] + j, rcvGaussList[sender] + j);
      }
    }
  }

  if(!useFaceElem) {
    // Let's look for the send list
    if(ffile->findToken(SEND_LIST_KW) < 0) {
      fprintf(stderr, "Token %s not found. Aborting\n", SEND_LIST_KW);
      exit(1);
    }

    cLine = ffile->getLineAfterToken();
    sscanf(cLine, "%d", &numRcv);
  }
  else {
    numRcv = 1;
  }

  int *sndcomid  = new int[numRcv];
  int *sndcomlen = new int[numRcv];
  int actualReceivers = 0;
  InterpPoint **interpPoints  = new InterpPoint *[numRcv];

  for(int receiver = 0; receiver < numRcv; ++ receiver) {
    if(!useFaceElem) {
      cLine = ffile->getLine();
      sscanf(cLine, "%d%d", sndcomid + receiver, sndcomlen + receiver);
    }
    else {
      sndcomid[receiver] = 1;
      sndcomlen[receiver] = fnodes->size();
    }

    if(sndcomid[receiver] > maxPRec) {
      maxPRec = sndcomid[receiver];
    }

    if(sndcomlen[receiver] > 0) {
      ++actualReceivers;
      interpPoints[receiver] = new InterpPoint[sndcomlen[receiver]];
    }
    else {
      interpPoints[receiver] = 0;
    }

    for(int j = 0; j < sndcomlen[receiver]; ++j) {
      if(!useFaceElem) {
        cLine = ffile->getLine();
        sscanf(cLine, "%d%lf%lf%lf%lf%lf", &interpPoints[receiver][j].elemNum,
               interpPoints[receiver][j].xy, interpPoints[receiver][j].xy + 1,
               interpPoints[receiver][j].gap, interpPoints[receiver][j].gap + 1,
               interpPoints[receiver][j].gap + 2);
        // map global elem number to packed elemset numbering (note: numbering by convention starts at 0)
        interpPoints[receiver][j].elemNum = geoSource->glToPackElem(interpPoints[receiver][j].elemNum-1);
      }
      else {
        map<int, locoord>::iterator it = exy.find(j);
        interpPoints[receiver][j].elemNum = (it->second).first;
        interpPoints[receiver][j].xy[0] = (it->second).second.first;
        interpPoints[receiver][j].xy[1] = (it->second).second.second;
        interpPoints[receiver][j].gap[0] = 0.0;
        interpPoints[receiver][j].gap[1] = 0.0;
        interpPoints[receiver][j].gap[2] = 0.0;
      }
    }
  }

  // Now create the compacted tables
  //  *** First the receive table is compacted to the actual receivers
  //  *** and we create the element wet mask
  nbrSendingToMe = actualSenders;
  bufferLen = 0;

  for(int sender = 0; sender < numSnd; ++sender) {
    if(rcvcomlen[sender] > 0) {

      if(bufferLen < 2 + NBPRESSDATAMAX * rcvcomlen[sender]) {
        bufferLen = 2 + NBPRESSDATAMAX * rcvcomlen[sender];
      }
    }
  }

  // *** Make the send tables
  nbrReceivingFromMe = actualReceivers;
  sndTable = new InterpPoint*[nbrReceivingFromMe];
  idSendTo = new int[nbrReceivingFromMe];
  nbSendTo = new int[nbrReceivingFromMe];
  int rId = 0;
  int mysize = 0;
  consOrigin = new int[maxPRec + 1];

  for(int receiver = 0; receiver < numRcv; ++receiver) {
    mysize += sndcomlen[receiver];

    if(sndcomlen[receiver] > 0) {
      if(bufferLen < 2 + 6 * sndcomlen[receiver]) {
        bufferLen = 2 + 6 * sndcomlen[receiver];
      }

      sndTable[rId] = interpPoints[receiver];
      nbSendTo[rId] = sndcomlen[receiver];
      idSendTo[rId] = sndcomid[receiver] - 1;
      consOrigin[sndcomid[receiver] - 1] = rId;
      rId++;
    }
  }

  buffer = new double[mysize * 6];

  // Last step: Create the dof nDof arrays and allocate localF
  int totalNDof = 0;
  int maxNDof = 0;
  int nDof;

  for(int i = 0; i < nbrReceivingFromMe; i++) {
    for(int j = 0; j < nbSendTo[i]; ++j) {
      if(!useFaceElem) {
        Element *thisElement = eset[sndTable[i][j].elemNum];
        nDof = thisElement->numDofs();
      }
      else {
        FaceElement *thisElement = (*feset)[sndTable[i][j].elemNum];
        nDof = thisElement->numDofs();
      }
      totalNDof += nDof;

      if(nDof > maxNDof) {
        maxNDof = nDof;
      }
    }
  }

  localF = new double[maxNDof];
  int *array = new int[totalNDof];

  for(int i = 0; i < nbrReceivingFromMe; i++) {
    for(int j = 0; j < nbSendTo[i]; ++j) {
      if(!useFaceElem) {
        Element *thisElement = eset[sndTable[i][j].elemNum];
        nDof = thisElement->numDofs();
        sndTable[i][j].dofs = array;
        thisElement->dofs(*dsa, array);
      }
      else {
        FaceElement *thisElement = (*feset)[sndTable[i][j].elemNum];
        nDof = thisElement->numDofs();
        sndTable[i][j].dofs = array;
        thisElement->dofs(*dsa, array, fnId);
      }
      array += nDof;
    }
  }

  if(!useFaceElem) {
    delete ffile;
    for(int sender = 0; sender < numSnd; ++sender) {
      if(rcvcomlen[sender] > 0) {
        delete [] rcvEleList[sender];
        delete [] rcvGaussList[sender];
      }
    }
  }

  delete [] rcvcomid;
  delete [] rcvcomlen;
  delete [] rcvEleList;
  delete [] rcvGaussList;
  delete [] sndcomid;
  delete [] sndcomlen;
  delete [] interpPoints;
}

// Temperature routines

void
FlExchanger::sendTemperature(State &state)
{
  /*  nbrReceivingFromMe = number of fluid subdomains to which a structural
                           subdomain has to send values to.
      nbSendTo[i] = number of fluid nodes of fluid subdomain i communicating
                    to a structural subdomain
      buffer[0]   = fluid subdomain number
      buffer[1]   = number of information contained in one node  */
  if(tmpDisp == 0) {
    tmpDisp = new Vector(state.getDisp());
  }

  *tmpDisp = state.getDisp();
  tmpDisp->linAdd(dtemp * alph[0], state.getVeloc(), dtemp * alph[1], state.getPrevVeloc());
  State newState(state, *tmpDisp);
  int i, j;
  int pos = 0;

  for(i = 0; i < nbrReceivingFromMe; i++) {
    int origPos = pos;

    for(j = 0; j < nbSendTo[i]; ++j) {
      Element *thisElement = eset[sndTable[i][j].elemNum];
      thisElement->computeTemp(cs, newState, sndTable[i][j].xy, buffer + pos);
      pos += 1;
    }

    int tag = STTOFLHEAT;
    int fluidNode  = idSendTo[i];
    fluidCom->sendTo(fluidNode, tag, buffer + origPos, pos - origPos);
  }

  waitOnSend();
}

void
FlExchanger::sendStrucTemp(Vector &tempsent)
{
  // Sends temperature to mechanical structure
  // tempsent are NODAL temperatures
  int tag = STTOSTMT;
  int zero = 0;
  int i;
  int pos = tempsent.size();

  for(i = 0; i < pos; i++) {
    buff[i] = tempsent[i];
  }

  heatStructCom->sendTo(zero, tag, buff, pos);
  heatStructCom->waitForAllReq();
}

void
FlExchanger::getStrucTemp(double *temprcvd)
{
  // Receives temperature from Thermal Elements
  // temprcvd are NODAL temperatures
  int tag = STTOSTMT;
  RecInfo rInfo = heatStructCom->recFrom(tag, buff, buffLen);
  int fromNd = rInfo.cpu;
  int rsize = rInfo.len;

  for(int i = 0; i < rsize; i++) {
    temprcvd[i] = buff[i] ;
  }
}

void
FlExchanger::sendHeatSource(Vector &heatsent)
{
  // Sends heat source from mechanical structure to thermal
  // heatsent are NODAL heat sources
  int tag = STTOSTHEAT;
  int zero = 0;
  int i;
  int pos = heatsent.size();

  for(i = 0; i < pos; i++) {
    buff[i] = heatsent[i];
  }

  heatStructCom->sendTo(zero, tag, buff, pos);
}

void
FlExchanger::getHeatSource(double *heatrcvd)
{
  // Receives heat sources from Structural Elements
  // heatrcvd are NODAL heat sources
  int tag = STTOSTHEAT;
  int rsize;
  int fromNd;
  RecInfo rInfo = heatStructCom->recFrom(tag, buff, buffLen);
  fromNd = rInfo.cpu;
  rsize = rInfo.len;
  int i;

  for(i = 0; i < rsize; i++) {
    heatrcvd[i] = buff[i] ;
  }
}

double
FlExchanger::getFluidFlux(Vector &flux, double time, double &bflux)
{
  // sndTable's extensions dofs, xy, elemNum are declared in FlExchange.h
  // in structure InterpPoint.
  // nbrReceivingFromMe =  number of fluid subdomains
  bflux = 0.;

  for(int i = 0; i < nbrReceivingFromMe; i++) {
    int tag = FLTOSTHEAT;
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
    int fromNd = rInfo.cpu;
    int rsize = rInfo.len;
    int origin = consOrigin[fromNd];

    // Loop Over wet points of each fluid subdomain

    for(int j = 0; j < nbSendTo[origin]; ++j) {
      Element *thisElement = eset[sndTable[origin][j].elemNum];
      thisElement->getFlFlux(sndTable[origin][j].xy, buffer + j, localF);
      int nDof = thisElement->numDofs();
      int *dof = sndTable[origin][j].dofs;

      // Summing the fluxes in each structural wet node

      for(int iDof = 0; iDof < nDof; ++iDof)
        if(dof[iDof] >= 0) {
          flux[dof[iDof]] += localF[iDof];
          bflux += localF[iDof];
        }
    }
  }

  return bflux;
}

void
FlExchanger::thermoread(int &bLen)
{
  // Initialize the buffer
  buffLen = bLen;
  buff = new double[buffLen];
}

void
FlExchanger::printreceiving()
{
  fprintf(stderr, " ::::::: nbrReceivingFromMe = %d\n", nbrReceivingFromMe);
}

void
FlExchanger::sendRelativeResidual(double relres)
{
  int returnFlag = 0;
  int FldNd = 0;
  int tag;
  int thisNode = structCom->myID();
  double buffer[1];
  buffer[0] = relres;
  int msglen = 1;

  if(thisNode == 0) {
    tag = STRELRESFL;
    fluidCom->sendTo(FldNd, tag, buffer, msglen);
  }
}

void
FlExchanger::sendNumParam(int numParam, int actvar, double steadyTol)
{
  int returnFlag = 0;
  int FldNd = 0;
  int tag;
  int thisNode = structCom->myID();
  double buffer[3];
  buffer[0] = (double) numParam;
  buffer[1] = (double) actvar;
  buffer[2] = steadyTol;
  int msglen = 3;

  if(thisNode == 0) {
    tag = STNUMPAFL;
    fluidCom->sendTo(FldNd, tag, buffer, msglen);
  }
}

void
FlExchanger::getNumParam(int &numParam)
{
  int tag;
  int thisNode = structCom->myID();
  double buffer[1];
  int msglen = 1;

  if(thisNode == 0) {
    tag =  FLNUMPAST;
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, msglen);
    numParam = (int) buffer[0];
  }
}

int
FlExchanger::cmdCom(int commandFlag)
{
  int returnFlag = 0;
  int FldNd  = 0;
  int tag;
  int thisNode = structCom->myID();
  double buffer[1];
  buffer[0] = (double) commandFlag;
  int msglen = 1;

  if(thisNode == 0) {
    tag = STCMDMSG;
    fluidCom->sendTo(FldNd, tag, buffer, msglen);
    tag =  FLCMDMSG;
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, msglen);
    returnFlag = (int) buffer[0];
  }

  return returnFlag;
}

int
FlExchanger::cmdComHeat(int commandFlag)
{
  int returnFlag = 0;
  int FldNd  = 0;
  int thisNode;
  int tag;
  thisNode = structCom->myID();
  double buffer[1];
  buffer[0] = (double) commandFlag;
  int msglen = 1;

  if(thisNode == 0) {
    tag = STCMDMSG;
    fluidCom->sendTo(FldNd, tag, buffer, msglen);
    fluidCom->waitForAllReq();
    tag =  FLCMDMSG;
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, msglen);
    returnFlag = (int) buffer[0];
  }

  return returnFlag;
}

// This routine negotiate with the fluid codes which match points go where
void
FlExchanger::negotiate()
{
  int thisNode = structCom->myID();
  int totSize = 0;
  int numFl = 0;
  numFl = fluidCom->remoteSize();
  int iFluid;
  int *flSize = new int[numFl];
  int totmatch = 0;

  for(iFluid = 0; iFluid < numFl; ++iFluid) {
    int tag = FL_NEGOT;
    int nFlMatched;
    int bufferLen = 1;
    int fromNd;
    RecInfo rInfo = fluidCom->recFrom(tag, &nFlMatched, bufferLen);
    fromNd = rInfo.cpu;
    flSize[fromNd] = nFlMatched;
    totmatch += nFlMatched;
  }

  if(totmatch == 0) {
    fprintf(stderr, " *** WARNING: by-passing negotiate step\n");
    fflush(stderr);
    return;
  }

  if(consOrigin) {
    delete [] consOrigin;
  }

  if(nbSendTo) {
    delete [] nbSendTo;
  }

  if(idSendTo) {
    delete [] idSendTo;
  }

  consOrigin = new int[numFl];
  nbSendTo = new int[numFl]; // should be actual number of senders
  idSendTo = new int[numFl]; // should be actual number of senders
  int nSender = 0;

  for(iFluid = 0; iFluid < numFl; ++iFluid) {
    if(flSize[iFluid] > 0) {
      idSendTo[nSender] = iFluid;
      nbSendTo[nSender] = flSize[iFluid];
      consOrigin[iFluid] = nSender;
      totSize += flSize[iFluid];
      nSender++;
    }
  }

  nbrSendingToMe = nbrReceivingFromMe = nSender;
  int *index = new int[totSize];
  InterpPoint **allPt = sndTable;
  sndTable = new InterpPoint *[nSender];

  for(iFluid = 0; iFluid < nSender; ++iFluid) {
    int tag = FL_NEGOT + 1;
    int bufferLen = totSize;
    int rsize, fromNd;
    RecInfo rInfo = fluidCom->recFrom(tag, index, bufferLen);
    fromNd = rInfo.cpu;
    rsize = rInfo.len;
    int sender = consOrigin[fromNd];
    sndTable[sender] = new InterpPoint[nbSendTo[sender]];

    for(int ipt = 0; ipt < nbSendTo[sender]; ++ipt) {
      sndTable[sender][ipt] = allPt[0][index[ipt]];
      index[ipt] = ipt;
    }

    tag = FL_NEGOT;
    int num = nbSendTo[sender];
    int one = 1;
    fluidCom->sendTo(fromNd, tag, &num, one);
    tag = FL_NEGOT + 1;
    fluidCom->sendTo(fromNd, tag, index, nbSendTo[sender]);
    // To make sure we can reuse the buffers
    fluidCom->waitForAllReq();
  }

  delete [] flSize;
  delete [] index;
  delete [] allPt;
  delete [] buffer;
  buffer = new double[6 * totSize];
}

void
FlExchanger::sendNoStructure()
{
  if(!sentInitialCracking) { // this means there is no initial cracking.
    sentInitialCracking = true;
  }

  int buf[4] = {0, 0, 0, 0};
  fluidCom->sendTo(0, 22/*tag*/, buf, 4);
  fluidCom->waitForAllReq();
}

void
FlExchanger::sendNewStructure(std::set<int> &newDeletedElements)
{
  std::set<int> newDeletedFaceElements;
  int fnodes[12];

  for(std::set<int>::iterator it = newDeletedElements.begin(); it != newDeletedElements.end(); ++it) {
    Element *ele = eset[geoSource->glToPackElem(*it)];
    int *enodes = ele->nodes();
    std::map<int, int> *GlToLlNodeMap = surface->GetPtrGlToLlNodeMap();

    for(int iNode = 0; iNode < ele->numNodes(); ++iNode) {
      std::map<int, int>::iterator it2 = GlToLlNodeMap->find(enodes[iNode]);

      if(it2 == GlToLlNodeMap->end()) {
        continue;
      }

      int *GlNodeIds = surface->GetPtrGlNodeIds();

      for(int j = 0; j < nodeToFaceElem->num(it2->second); j++) { // loop over the face elements connected to the iNode-th node
        int k = (*nodeToFaceElem)[it2->second][j];

        if(newDeletedFaceElements.find(k) != newDeletedFaceElements.end()) {
          continue;
        }

        FaceElement *faceEl = surface->GetFaceElemSet()[k];

        if(faceEl && (faceEl->nNodes() <= ele->numNodes())) {
          faceEl->GetNodes(fnodes, GlNodeIds);
#if (((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && HAS_CXX11_ALL_OF && HAS_CXX11_LAMBDA)
          if(std::all_of(fnodes, fnodes + faceEl->nNodes(),
          [&](int i) {
          return (std::find(enodes, enodes + ele->numNodes(), i) != enodes + ele->numNodes());
          })) {
            //std::cerr << "face element " << k+1 << " on embedded surface is associated with deleted element " << *it+1 << std::endl;
            newDeletedFaceElements.insert(k);
            break;
          }
#else
          std::cerr << " *** ERROR: C++11 support required in FlExchanger::sendNewStructure().\n";
          exit(-1);
#endif
        }
      }
    }

    delete [] enodes;
  }

  if(newDeletedFaceElements.empty()) {
    // in this case, none of the deleted elements are on the wet interface
    sendNoStructure();
  }
  else {
    int numConnUpdated = newDeletedFaceElements.size();
    int *newConn = new int[5 * numConnUpdated];
    int numLvlUpdated2 = newDeletedFaceElements.size();
    int *lvlsetElemNum = new int[numLvlUpdated2];
    double *lvlsets = new double[4 * numLvlUpdated2];
    int numNewNodes = 0;
    int *phantom_nodes_indices = NULL;
    double *phantom_nodes_xyz0 = NULL;
    int *new2old = NULL;
    int i;
    std::set<int>::iterator it;

    for(it = newDeletedFaceElements.begin(), i = 0; it != newDeletedFaceElements.end(); ++it, ++i) {
      newConn[5 * i + 0] = lvlsetElemNum[i] = *it;

      for(int j = 0; j < faceElemToNode->num(*it); ++j) {
        newConn[5 * i + 1 + j] = (*faceElemToNode)[*it][j];
        lvlsets[4 * i + j] = -1;
      }

      if(faceElemToNode->num(*it) == 3) { // convert 3-node triangle to degenerate quad
        newConn[5 * i + 1 + 3] = (*faceElemToNode)[*it][2];
        lvlsets[4 * i + 3] = -1;
      }
    }

    // send cracking "stats".
    int buf[4] = {1, numConnUpdated, numLvlUpdated2, numNewNodes};
    fluidCom->sendTo(0, 22/*tag*/, buf, 4);
    fluidCom->waitForAllReq();

    // for initial cracking, send phantom nodes.
    if(!sentInitialCracking) {
      fluidCom->sendTo(0, 55/*tag*/, phantom_nodes_xyz0, numNewNodes * 3);
      fluidCom->waitForAllReq();
      sentInitialCracking = true;
    }

    // send the integer package (topology + phi indices + new2old)
    int pack1[numConnUpdated * 5 + numLvlUpdated2 + 2 * numNewNodes];

    for(int i = 0; i < numConnUpdated; i++) {
      pack1[5 * i] = newConn[5 * i];

      for(int j = 1; j < 5; j++) {
        pack1[5 * i + j] = newConn[5 * i + j];
      }
    }

    for(int i = 0; i < numLvlUpdated2; i++) {
      pack1[5 * numConnUpdated + i] = lvlsetElemNum[i];
    }

    for(int i = 0; i < numNewNodes; i++) {
      pack1[5 * numConnUpdated + numLvlUpdated2 + 2 * i]   = phantom_nodes_indices[i];
      pack1[5 * numConnUpdated + numLvlUpdated2 + 2 * i + 1] = new2old[i];
    }

    fluidCom->sendTo(0, 33/*tag*/, pack1, numConnUpdated * 5 + numLvlUpdated2 + 2 * numNewNodes);
    fluidCom->waitForAllReq();
    // send phi
    fluidCom->sendTo(0, 44/*tag*/, lvlsets, numLvlUpdated2 * 4);
    fluidCom->waitForAllReq();
    delete [] newConn;
    delete [] lvlsetElemNum;
    delete [] lvlsets;
  }
}

void
FlExchanger::transformVector(double *localF, Element *ele)
{
  if(!solInfo.basicDofCoords) {
    int *nn = ele->nodes();

    for(int k = 0; k < ele->numNodes(); ++k)
      if(NFrameData *cd = cs.dofFrame(nn[k])) {
        if(ele->hasRot()) {
          cd->transformVector6(localF + 6 * k);
        }
        else {
          cd->transformVector3(localF + 3 * k);
        }
      }

    delete [] nn;
  }
}

void
FlExchanger::transformVector(double *localF, FaceElement *ele)
{
  if(!solInfo.basicDofCoords) {
    for(int k = 0; k < ele->nNodes(); ++k) {
      int glNode = (surface->IsRenumbered()) ? surface->GetPtrLlToGlNodeMap()[ele->GetNode(k)] : ele->GetNode(k);

      if(NFrameData *cd = cs.dofFrame(glNode)) {
        cd->transformVector3(localF + 3 * k);
      }
    }
  }
}
