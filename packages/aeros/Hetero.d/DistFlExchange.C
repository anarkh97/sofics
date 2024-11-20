#include <Hetero.d/DistFlExchange.h>
#include <Hetero.d/FilteredFile.h>
#include <Utils.d/dofset.h>
#include <Element.d/State.h>

#include <Comm.d/Communicator.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/SolverInfo.h>
#include <Corotational.d/DistrGeomState.h>

#include <Driver.d/SubDomain.h>
#include <Utils.d/DistHelper.h>
#include <Feti.d/DistrVector.h>
#include <Driver.d/SysState.h>

#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>

// std
#include <map>
#include <set>
using std::map;

#define DETERMINISTIC_DISTFLEXCHANGER

extern Communicator *structCom, *fluidCom, *heatStructCom;
extern int verboseFlag;
extern GeoSource *geoSource;
extern SolverInfo &solInfo;
extern Domain *domain;

DistFlExchanger::DistFlExchanger(CoordSet **_cs, Elemset **_eset, SurfaceEntity *_surface, CoordSet *_globalCoords,
                                 Connectivity *_nodeToElem, Connectivity *_elemToSub, SubDomain **_sd,
                                 DofSetArray **_cdsa, DofSetArray **_dsa, OutputInfo *_oinfo, bool _wCracking)
  : cs(_cs), eset(_eset), globalCoords(_globalCoords), nodeToElem(_nodeToElem),
    elemToSub(_elemToSub), sd(_sd), cdsa(_cdsa), dsa(_dsa), oinfo(_oinfo),
    tmpDisp(0), dsp(0), vel(0), acc(0), pVel(0)
{
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

  buff = 0;
  idSendTo = 0;
  nbSendTo = 0;
  consOrigin = 0;
  localF = 0;
  buffer = 0;
  sndTable = 0;
}

DistFlExchanger::DistFlExchanger(CoordSet **_cs, Elemset **_eset,
                                 DofSetArray **_cdsa, DofSetArray **_dsa,
                                 OutputInfo *_oinfo)
  : cs(_cs), eset(_eset), cdsa(_cdsa), dsa(_dsa), oinfo(_oinfo),
    tmpDisp(0), dsp(0), vel(0), acc(0), pVel(0)
{
  useFaceElem = false;
  wCracking = false;
  sentInitialCracking = false;
  surface = 0;
  faceElemToNode = 0;
  nodeToFaceElem = 0;

  buff = 0;
  idSendTo = 0;
  nbSendTo = 0;
  consOrigin = 0;
  localF = 0;
  buffer = 0;
  sndTable = 0;
}

DistFlExchanger::~DistFlExchanger()
{
  if(useFaceElem) {
    for(int i = 0; i < geoSource->getCpuToSub()->num(structCom->myID()); ++i) {
      delete [] fnId2[i];
    }

    delete [] fnId2;
  }
  delete [] cs;
  delete [] eset;
  delete [] cdsa;
  delete [] dsa;

  if(tmpDisp) {
    delete tmpDisp;
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
    for(int i=0; i<nSender; ++i) if(sndTable[i]) { delete [] sndTable[i][0].dofs; delete [] sndTable[i]; }
    delete [] sndTable;
  }

  if(buff) {
    delete [] buff;
  }

  if(dsp) {
    delete [] dsp;
  }

  if(vel) {
    delete [] vel;
  }

  if(acc) {
    delete [] acc;
  }

  if(pVel) {
    delete [] pVel;
  }
}

double
DistFlExchanger::getFluidLoad(DistrVector &force, int tIndex, double time,
                              double alphaf, int &iscollocated, DistrGeomState *distrGeomState)
{
  aforce[0] = aforce[1] = aforce[2] = 0.0;
  FaceElemSet *feset;

  if(useFaceElem) {
    feset = &(surface->GetFaceElemSet());
  }

  Element     *thisElement;
  FaceElement *thisFaceElem;

#ifdef DETERMINISTIC_DISTFLEXCHANGER
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
      int locSub = sndTable[origin][j].subNumber;
      GeomState *geomState = (distrGeomState) ? (*distrGeomState)[locSub] : 0;
      int nDof;

      if(!useFaceElem) {
        thisElement = (*eset[locSub])[sndTable[origin][j].elemNum];
        thisElement->getFlLoad(*(cs[locSub]), sndTable[origin][j], buffer + 3 * j, localF, geomState);
        nDof = thisElement->numDofs();
        transformVector(localF, thisElement, locSub);
      }
      else {
        thisFaceElem = (*feset)[sndTable[origin][j].elemNum];
        thisFaceElem->getFlLoad(sndTable[origin][j], buffer + 3 * j, localF);
        nDof = thisFaceElem->numDofs();
        transformVector(localF, thisFaceElem, locSub);
      }

      int *dof = sndTable[origin][j].dofs;

      for(int iDof = 0; iDof < nDof; ++iDof)
        if(dof[iDof] >= 0) {
          force.subData(locSub)[dof[iDof]] += localF[iDof];
        }

      aforce[0] += buffer[3 * j];
      aforce[1] += buffer[3 * j + 1];
      aforce[2] += buffer[3 * j + 2];
    }
  }

  flipRcvParity();

  if(oinfo) {
    if(tIndex % oinfo->interval == 0 && oinfo->filptr != NULL) {
      structCom->reduce(3, aforce);
      filePrint(oinfo->filptr, "%e   ", time);
      filePrint(oinfo->filptr, "%e %e %e\n", aforce[0], aforce[1], aforce[2]);
      fflush(oinfo->filptr);
    }
  }

#ifdef DETERMINISTIC_DISTFLEXCHANGER
  delete [] buffers;
#endif

  // ML & KP For 'corrected' aeroelastic force
  iscollocated = (isCollocated) ? 1 : 0;
  return (time + alphaf * dt);
}

void
DistFlExchanger::sendDisplacements(SysState<DistrVector> &state,
                                   double **usrDefDisps, double **usrDefVels, int tag,
                                   DistrGeomState *distrGeomState)
{
  auto cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

  if(tmpDisp == 0) {
    tmpDisp = new DistrVector(state.getDisp());
  }

  *tmpDisp = state.getDisp();
  tmpDisp->linAdd(dt * alpha[0], state.getVeloc(), dt * alpha[1], state.getPrevVeloc());
  //if(verboseFlag)
  //  filePrint(stderr, "Disp Norm %e Veloc Norm %e\n", tmpDisp->norm(), state.getVeloc().norm());
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

  // get the velocities, accelerations, and previous velocities
  if(dsp == 0)  {
    dsp = new Vector[numSub];
    vel = new Vector[numSub];
    acc = new Vector[numSub];
    pVel = new Vector[numSub];
  }

  for(int iSub = 0; iSub < numSub; iSub++) {
    dsp[iSub].setData(tmpDisp->subData(iSub), tmpDisp->subLen(iSub), false);
    vel[iSub].setData(state.getVeloc().subData(iSub),
                      state.getVeloc().subLen(iSub), false);
    acc[iSub].setData(state.getAccel().subData(iSub),
                      state.getAccel().subLen(iSub), false);
    pVel[iSub].setData(state.getPrevVeloc().subData(iSub),
                       state.getPrevVeloc().subLen(iSub), false);
  }

  for(int i = 0; i < nbrReceivingFromMe; i++) {
    int fluidNode = idSendTo[i];
    int origin = consOrigin[fluidNode];
    int origPos = pos;

    // compute displacements for each match point
    for(int j = 0; j < nbSendTo[origin]; ++j) {
      int locSub = sndTable[origin][j].subNumber;
      State localState(cdsa[locSub], dsa[locSub], usrDefDisps[locSub],
                       usrDefVels[locSub], dsp[locSub],
                       vel[locSub], acc[locSub], pVel[locSub], cs[locSub]);
      GeomState *geomState = (distrGeomState) ? (*distrGeomState)[locSub] : 0;

      if(!useFaceElem) {
        thisElement = (*eset[locSub])[sndTable[origin][j].elemNum];
        thisElement->computeDisp(*cs[locSub], localState, sndTable[origin][j], buffer + pos, geomState);
      }
      else {
        thisFaceElem = (*feset)[sndTable[origin][j].elemNum];
        thisFaceElem->computeDisp(*cs[locSub], localState, sndTable[origin][j], buffer + pos, geomState, fnId2[locSub]); // TODO fnId2 must map from embedded surface to locSub node numbers
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
  //if (verboseFlag)
  //  fprintf(stderr, "Sending %e and %e to %d CPUs\n", xxx, yyy, nbrReceivingFromMe);
  flipSndParity();
}

void
DistFlExchanger::sendModeShapes(int numModes, int numN, double(**v)[6],
                                SysState<DistrVector> &st, double ampFactor)
{
  std::cerr << "ERROR: DistFlExchanger::sendModeShapes is not implemented\n";
  exit(-1);
}

void DistFlExchanger::sendParam(int _algnum, double step, double totaltime,
                                int rstinc, int _isCollocated, double _a[2])
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
}

void
DistFlExchanger::sendSubcyclingInfo(int sub)
{
  double buffer = (double)sub;
  fluidCom->sendTo(0, 777/*tag*/, &buffer, 1);
  fluidCom->waitForAllReq();
}

void
DistFlExchanger::sendTempParam(int _algnum, double step, double totaltime,
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
DistFlExchanger::sendModeFreq(double *modfrq, int nummod)
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
DistFlExchanger::sendEmbeddedWetSurface()
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
    nodes[3 * i]   = (*globalCoords)[fnId[i]]->x; // TODO do it without globalCoords
    nodes[3 * i + 1] = (*globalCoords)[fnId[i]]->y;
    nodes[3 * i + 2] = (*globalCoords)[fnId[i]]->z;
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

MatchMap *
DistFlExchanger::getMatchData()
{
  auto cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);
  MatchMap *globMatches = new MatchMap();
  // create global sub to local sub map
  int totSub = cpuToSub->numConnect();
  int *glToLocSubMap = new int[totSub];

  for(int iSub = 0; iSub < totSub; iSub++) {
    glToLocSubMap[iSub] = -1;
  }

  for(int iSub = 0; iSub < numSub; iSub++) {
    glToLocSubMap[(*cpuToSub)[myCPU][iSub] ] = iSub;
  }

  if(!useFaceElem) {
    for(int iSub = 0; iSub < numSub; iSub++) {
      int glSub = (*cpuToSub)[myCPU][iSub];
      int locSub = glToLocSubMap[glSub];

      if(locSub < 0) {
        fprintf(stderr, "*** ERROR: Bad Global to Local Sub Map\n");
      }

      MatchData *matchData = geoSource->getMatchData(locSub);
      int *numMatch = geoSource->getNumMatchData();

      for(int iMatch = 0; iMatch < numMatch[locSub]; iMatch++) {
        int gPos = matchData[iMatch].glPos;
        InterpPoint iPoint;
        iPoint.subNumber = iSub;
        iPoint.elemNum = matchData[iMatch].elemNum;
        iPoint.xy[0] = matchData[iMatch].xi;
        iPoint.xy[1] = matchData[iMatch].eta;
        iPoint.gap[0] = 0.0;
        iPoint.gap[1] = 0.0;
        iPoint.gap[2] = 0.0;
        (*globMatches)[gPos] = iPoint;
      }
    }
  }
  else {
    // info about surface
    FaceElemSet &feset = surface->GetFaceElemSet();
    CoordSet   &fnodes = surface->GetNodeSet();
    int          *fnId = surface->GetPtrGlNodeIds();
    // assign face elements to a subdomain
    int totele = elemToSub->csize();
    int *eleTouch = new int[totele];
    int *eleCount = new int[totele];

    for(int j = 0; j < totele; j++) {
      eleTouch[j] = -1;
    }

    int numFaceEle = (surface->GetIsShellFace() && domain->tdenforceFlag()) ? feset.last()/2 : feset.last();
    int *faceElemToSub = new int[numFaceEle];

    for(int j = 0; j < numFaceEle; ++j) {
      int glEle = feset[j]->findEle(nodeToElem, eleTouch, eleCount, j, fnId);
      int glSub = (*elemToSub)[glEle][0];
      faceElemToSub[j] = glToLocSubMap[glSub];
    }

    delete [] eleTouch;
    delete [] eleCount;
    // find node to elem connectivity and local coords.
    map<int, locoord> exy = feset.computeNodeLocalCoords(fnId, fnodes.size());
    map<int, locoord>::iterator it;

    for(int j = 0; j < fnodes.size(); ++j) {
      it = exy.find(j);

      int faceElemNum = (it->second).first;
      int subNumber = faceElemToSub[faceElemNum];

      if(subNumber > -1) {
        InterpPoint iPoint;
        iPoint.subNumber = subNumber;
        iPoint.elemNum = faceElemNum;
        iPoint.xy[0] = (it->second).second.first;
        iPoint.xy[1] = (it->second).second.second;
        iPoint.gap[0] = 0.0;
        iPoint.gap[1] = 0.0;
        iPoint.gap[2] = 0.0;
        (*globMatches)[j] = iPoint;
      }
    }

    delete [] faceElemToSub;
    fnId2 = new int *[numSub];

    for(int i = 0; i < numSub; ++i) {
      fnId2[i] = new int[fnodes.size()];

      for(int j = 0; j < fnodes.size(); ++j) {
        fnId2[i][j] = sd[i]->globalToLocal(fnId[j]);
      }
    }
  }

  delete [] glToLocSubMap;
  return globMatches;
}

void
DistFlExchanger::sendTemperature(SysState<DistrVector> &dState)
{
  auto cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

  if(tmpDisp == 0) {
    tmpDisp = new DistrVector(dState.getDisp());
  }

  *tmpDisp = dState.getDisp();
  tmpDisp->linAdd(dtemp * alpha[0], dState.getVeloc(), dtemp * alpha[1], dState.getPrevVeloc());
  double xxx = 0, yyy = 0;
  int pos = 0;

  // get the temperature, velocities and previous velocities
  if(dsp == 0)  {
    dsp = new Vector [numSub];
    vel = new Vector [numSub];
    pVel = new Vector [numSub];
  }

  for(int iSub = 0; iSub < numSub; iSub++)  {
    dsp[iSub].setData(tmpDisp->subData(iSub), tmpDisp->subLen(iSub));
    vel[iSub].setData(dState.getVeloc().subData(iSub),
                      dState.getVeloc().subLen(iSub));
    pVel[iSub].setData(dState.getPrevVeloc().subData(iSub),
                       dState.getPrevVeloc().subLen(iSub));
  }

  for(int iNeigh = 0; iNeigh < nbrReceivingFromMe; iNeigh++) {
    // get fluid mpi to send disps to
    int fluidNode = idSendTo[iNeigh];
    int origin = consOrigin[fluidNode];
    int origPos = pos;

    // compute displacements for each match point
    for(int iData = 0; iData < nbSendTo[origin]; ++iData) {
      int locSub = sndTable[origin][iData].subNumber;
      Element *thisElement = (*eset[locSub])[sndTable[origin][iData].elemNum];
      State localState(cdsa[locSub], dsa[locSub], (double *) 0 /*usrDefDisps[locSub]*/,
                       dsp[locSub], vel[locSub], pVel[locSub]);
      thisElement->computeTemp(*cs[locSub], localState,
                               sndTable[origin][iData].xy, buffer + pos);
      pos += 1;
    }

    int tag = STTOFLHEAT;
    fluidCom->sendTo(fluidNode, tag, buffer + origPos, pos - origPos);
  }

  fluidCom->waitForAllReq();
}

void
DistFlExchanger::sendStrucTemp(DistrVector &tempsent)
{
  // for now we assume that the structure and thermal models have the same mesh and the same decomposition
  auto cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);
  // Sends temperature to mechanical structure
  // for now we assume that the structure and thermal models have the same mesh and the same decomposition
  int thisNode;
  thisNode = structCom->myID();
  int tag = STTOSTMT;
  int zero = 0;
  int i;
  int pos = tempsent.size();

  for(i = 0; i < pos; i++) {
    buff[i] = tempsent[i];
  }

  heatStructCom->sendTo(thisNode, tag, buff, pos);
}

void
DistFlExchanger::getStrucTemp(double *temprcvd)
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

double
DistFlExchanger::getFluidFlux(DistrVector &flux, double time, double &bflux)
{
  for(int i = 0; i < nbrReceivingFromMe; i++) {
    int tag = FLTOSTHEAT;
    RecInfo rInfo = fluidCom->recFrom(tag, buffer, bufferLen);
    int fromNd = rInfo.cpu;
    int origin = consOrigin[fromNd];

    // Loop Over wet points of each fluid subdomain

    for(int j = 0; j < nbSendTo[origin]; ++j) {
      int locSub = sndTable[origin][j].subNumber;
      Element *thisElement = (*eset[locSub])[sndTable[origin][j].elemNum];
      thisElement->getFlFlux(sndTable[origin][j].xy, buffer + j, localF);
      int nDof = thisElement->numDofs();
      int *dof = sndTable[origin][j].dofs;

      // Summing the fluxes in each structural wet node

      for(int iDof = 0; iDof < nDof; ++iDof)
        if(dof[iDof] >= 0) {
          flux.subData(locSub)[dof[iDof]] += localF[iDof];
          bflux += localF[iDof];
        }
    }
  }

  return bflux;
}

void
DistFlExchanger::thermoread(int bLen)
{
  // Initialize the buffer
  buffLen = bLen;
  buff = new double[buffLen];
}

int
DistFlExchanger::cmdCom(int commandFlag)
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
  structCom->broadcast(1, &returnFlag);

  return returnFlag;
}

// This routine negotiate with the fluid codes which match points go where
void DistFlExchanger::negotiate()
{
  FaceElemSet *feset;

  if(useFaceElem) {
    feset = &(surface->GetFaceElemSet());
    sendEmbeddedWetSurface();
  }

  int numFl = fluidCom->remoteSize();  // number of fluid mpi processes
  int *flSize = new int[numFl];  // number of matches per fluid mpi
  int iFluid;

  for(iFluid = 0; iFluid < numFl; ++iFluid) {
    int tag = FL_NEGOT;
    int nFlMatched;
    int buffLen = 1;
    RecInfo rInfo = fluidCom->recFrom(tag, &nFlMatched, buffLen);
    int fromNd = rInfo.cpu;
    flSize[fromNd] = nFlMatched;
  }

  // allocate with number of fluid cpu's doing a negotiate
  consOrigin = new int[numFl];
  nbSendTo = new int[numFl];
  nSender = 0;  // number of fluid mpi's sending data
  int maxSize = 0;

  for(iFluid = 0; iFluid < numFl; ++iFluid) {
    if(flSize[iFluid] > 0) {
      nbSendTo[nSender] = flSize[iFluid];
      consOrigin[iFluid] = nSender;

      if(flSize[iFluid] > maxSize) {
        maxSize = flSize[iFluid];
      }

      nSender++;
    }
  }

  MatchMap *globMatches = getMatchData();
  // # of fluid mpi's receiving data from this mpi
  nbrReceivingFromMe = 0;
  int *index = new int[maxSize];
  sndTable = new InterpPoint *[nSender];
  // receive matches from each sender
  int bufferLen = maxSize;
  int totSize = 0;  // total number of matches received by structure
  int maxNDof = 0;  // to set size of local force vector

  for(iFluid = 0; iFluid < nSender; ++iFluid) {
    int tag = FL_NEGOT + 1;
    RecInfo rInfo = fluidCom->recFrom(tag, index, bufferLen);
    int fromNd = rInfo.cpu;
    int sender = consOrigin[fromNd]; // gives fluid mpi #
    // determine how many matches are in this mpi process
    int numHits = 0;
    int iMatch;

    for(iMatch = 0; iMatch < nbSendTo[sender]; iMatch++)
      if(globMatches->find(index[iMatch]) != globMatches->end()) {
        numHits++;
      }

    totSize += numHits;
    int *matchList;  // list of matches
    int totalNDof = 0;

    if(numHits > 0) {
      matchList = new int[numHits];
      sndTable[sender] = new InterpPoint[numHits];
      // create list of matches in this mpi
      int iHit = 0;

      for(iMatch = 0; iMatch < nbSendTo[sender]; iMatch++)
        if(globMatches->find(index[iMatch]) != globMatches->end()) {
          matchList[iHit++] = iMatch;
        }

      // reset nbSendTo to actual number of matches in this mpi
      nbSendTo[sender] = numHits;

      // assign match data to the sndTable
      for(int ipt = 0; ipt < nbSendTo[sender]; ++ipt) {
        int arrayPos = index[matchList[ipt]];
        sndTable[sender][ipt].subNumber = (*globMatches)[arrayPos].subNumber;
        sndTable[sender][ipt].elemNum = (*globMatches)[arrayPos].elemNum;
        sndTable[sender][ipt].xy[0] = (*globMatches)[arrayPos].xy[0];
        sndTable[sender][ipt].xy[1] = (*globMatches)[arrayPos].xy[1];
        sndTable[sender][ipt].gap[0] = 0;
        sndTable[sender][ipt].gap[1] = 0;
        sndTable[sender][ipt].gap[2] = 0;
        // assign dofs to element
        int locSub = sndTable[sender][ipt].subNumber;
        int locElem = sndTable[sender][ipt].elemNum;
        int nDof;

        if(!useFaceElem) {
          Element *thisElement = (*eset[locSub])[locElem];

          if(thisElement == NULL) {
            std::cerr << " ERROR: DistFlExchanger::negotiate() cpu " << structCom->myID() << ", locSub = "
                      << locSub << ", locElem = " << locElem << std::endl;
          }

          nDof = thisElement->numDofs();
        }
        else {
          FaceElement *thisFaceElement = (*feset)[locElem];
          nDof = thisFaceElement->numDofs();
        }

        totalNDof += nDof;

        if(nDof > maxNDof) {
          maxNDof = nDof;
        }
      }

      nbrReceivingFromMe++;
    }
    else  {
      sndTable[sender] = 0;
      matchList = 0;
      nbSendTo[sender] = 0;
    }

    // allocate the local dofset array
    int *array = new int[totalNDof];

    // populate the dofs in the send table
    for(int iData = 0; iData < numHits; iData++)  {
      int locSub = sndTable[sender][iData].subNumber;
      int locElem = sndTable[sender][iData].elemNum;
      sndTable[sender][iData].dofs = array;

      if(!useFaceElem) {
        Element *thisElement = (*eset[locSub])[locElem];
        thisElement->dofs(*(cdsa[locSub]), array);
        array += thisElement->numDofs();
      }
      else {
        FaceElement *thisFaceElement = (*feset)[locElem];
        thisFaceElement->dofs(*(cdsa[locSub]), array, fnId2[locSub]); // TODO fnId2 must map from embedded surface to locSub node numbers
        array += thisFaceElement->numDofs();
      }
    }

    tag = FL_NEGOT;
    int one = 1;
    fluidCom->sendTo(fromNd, tag, &numHits, one);

    if(numHits > 0)  {
      tag = FL_NEGOT + 1;
      fluidCom->sendTo(fromNd, tag, matchList, numHits);
    }

    // To make sure we can reuse the buffers
    fluidCom->waitForAllReq();
    if(matchList) delete [] matchList;
  }

  /*for(int i=0; i<structCom->numCPUs(); ++i) {
    if(structCom->myID() == i) {
      std::cerr << "here in DistFlExchanger::negotiate(), i = " << i << ", numFl = " << numFl << ", nSender = " << nSender
                << ", nbrReceivingFromMe = " << nbrReceivingFromMe << std::endl;
    }
    structCom->sync();
  }*/

  // allocate the size of local force vector
  localF = new double[maxNDof];
  // create list of fluid mpi's to send to
  idSendTo = new int[nbrReceivingFromMe];
  int count = 0;

  for(iFluid = 0; iFluid < numFl; iFluid++)
    if(flSize[iFluid] > 0 && sndTable[consOrigin[iFluid]] != 0) {
      idSendTo[count++] = iFluid;
    }

  delete [] flSize;
  delete [] index;
  delete globMatches;
  buffer = new double[6 * totSize];
  setBufferLength(6 * totSize);
}

void
DistFlExchanger::sendNoStructure()
{
  if(structCom->myID() != 0) {
    return;  /* do nothing */
  }

  if(!sentInitialCracking) { //this means there is no initial cracking.
    sentInitialCracking = true;
  }

  int buf[4] = {0, 0, 0, 0};
  fluidCom->sendTo(0, 22/*tag*/, buf, 4);
  fluidCom->waitForAllReq();
}

void
DistFlExchanger::sendNewStructure()
{
  auto cpuToSub = geoSource->getCpuToSub();
  int myCPU = structCom->myID();
  int numSub = cpuToSub->num(myCPU);

  std::set<int> localDeletedElements, newDeletedFaceElements;
  int fnodes[12];

  for(int iSub = 0; iSub < numSub; iSub++) {
    std::set<int> &newDeletedElements = sd[iSub]->getNewDeletedElements();
    for(std::set<int>::iterator it = newDeletedElements.begin(); it != newDeletedElements.end(); ++it) {
      Element *ele = (*eset[iSub])[sd[iSub]->globalToLocalElem(geoSource->glToPackElem(*it))];
      int *enodes = ele->nodes();
      for(int iNode=0; iNode<ele->numNodes(); ++iNode) enodes[iNode] = sd[iSub]->localToGlobal(enodes[iNode]);
      std::map<int, int> *GlToLlNodeMap = surface->GetPtrGlToLlNodeMap();

      for(int iNode = 0; iNode < ele->numNodes(); ++iNode) {
        std::map<int, int>::iterator it2 = GlToLlNodeMap->find(enodes[iNode]);

        if(it2 == GlToLlNodeMap->end()) {
          continue;
        }

        int *GlNodeIds = surface->GetPtrGlNodeIds();

        for(int j = 0; j < nodeToFaceElem->num(it2->second); j++) { // loop over the face elements connected to the iNode-th node
          int k = (*nodeToFaceElem)[it2->second][j];

          if(localDeletedElements.find(k) != localDeletedElements.end()) {
            continue;
          }

          FaceElement *faceEl = surface->GetFaceElemSet()[k];

          if(faceEl && (faceEl->nNodes() <= ele->numNodes())) {
            faceEl->GetNodes(fnodes, GlNodeIds);
#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && HAS_CXX11_ALL_OF && HAS_CXX11_LAMBDA
            if(std::all_of(fnodes, fnodes + faceEl->nNodes(),
            [&](int i) {
            return (std::find(enodes, enodes + ele->numNodes(), i) != enodes + ele->numNodes());
            })) {
              //std::cerr << "face element " << k+1 << " on embedded surface is associated with deleted element " << *it+1 << std::endl;
              localDeletedElements.insert(k);
              break;
            }
#else
            std::cerr << " *** ERROR: C++11 support required in DistFlExchanger::sendNewStructure().\n";
            exit(-1);
#endif
          }
        }
      }

      delete [] enodes;
    }
  }

  // now reduce over all of the mpi processes
  int localCount = localDeletedElements.size();
  int *recvbuf = new int[structCom->numCPUs()];
  structCom->gather(&localCount, 1, recvbuf, 1);
  int globalCount;
  if(structCom->myID() == 0) {
    globalCount = 0;
    for(int i=0; i<structCom->numCPUs(); ++i) globalCount += recvbuf[i];
  }
  structCom->broadcast(1, &globalCount);
  if(globalCount > 0) {
    int *sendbuf2 = new int[localCount];
    if(localCount > 0) {
      int i=0;
      for(std::set<int>::iterator it = localDeletedElements.begin(); it != localDeletedElements.end(); ++it, ++i) {
        sendbuf2[i] = *it;
      }
    }
    int *recvbuf2, *displs;
    if(structCom->myID() == 0) {
      recvbuf2 = new int[globalCount];
      displs = new int[structCom->numCPUs()];
      displs[0] = 0;
      for(int i=1; i<structCom->numCPUs(); ++i) {
        displs[i] = displs[i-1] + recvbuf[i-1];
      }
    }
    structCom->gatherv(sendbuf2, localCount, recvbuf2, recvbuf, displs);
    delete [] sendbuf2;
    if(structCom->myID() == 0) {
      for(int i=0; i<globalCount; ++i) {
        newDeletedFaceElements.insert(recvbuf2[i]);
      }
      delete [] recvbuf2;
      delete [] displs;
    }
  }
  delete [] recvbuf;

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
DistFlExchanger::transformVector(double *localF, Element *ele, int locSub)
{
  if(!solInfo.basicDofCoords) {
    int *nn = ele->nodes();

    for(int k = 0; k < ele->numNodes(); ++k)
      if(NFrameData *cd = cs[locSub]->dofFrame(nn[k])) {
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
DistFlExchanger::transformVector(double *localF, FaceElement *ele, int locSub)
{
  if(!solInfo.basicDofCoords) {
    for(int k = 0; k < ele->nNodes(); ++k) {
      int glNode = (surface->IsRenumbered()) ? surface->GetPtrLlToGlNodeMap()[ele->GetNode(k)] : ele->GetNode(k);

      if(NFrameData *cd = cs[locSub]->dofFrame(sd[locSub]->globalToLocal(glNode))) {
        cd->transformVector3(localF + 3 * k);
      }
    }
  }
}

