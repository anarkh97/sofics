#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>
#include <Driver.d/Domain.h>
#include <Driver.d/Mpc.h>
#include <Element.d/MpcElement.d/MpcElement.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Helm.d/ARubberF.h>
#include <Element.d/ElemAccess.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/Connectivity.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Memory.h>
#include <Utils.d/ModeData.h>
#include <Dec.d/Decomp.d/Decomp.h>
#include <Driver.d/MultiFront.h>
#include <Driver.d/Header.h>
#include <Corotational.d/DistrGeomState.h>
#ifndef WINDOWS
#include <dlfcn.h>
#endif
#include <map>
#include <utility>
#include <list>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <iterator>
#include <string>
#include <Sfem.d/Sfem.h>
#ifdef USE_EIGEN3
#include <Math.d/rref.h>
#include <Eigen/Core>
#endif

EFrameData null_eframe;
NFrameData null_nframe;
StructProp null_sprops;
Attrib null_attrib;
BCond null_bcond;
CoefData null_coef;
OutputInfo emptyInfo;
extern ModeData modeData;
extern bool callDec;
extern bool exitAfterDec;
extern bool estFlag;
extern bool weightOutFlag;
extern bool trivialFlag;
extern bool randomShuffle;
extern bool allowMechanisms;
extern bool useScotch;
extern int verboseFlag;
#ifndef SALINAS
extern Sfem *sfem;
#endif
//----------------------------------------------------------------------

template< class RandomIt >
void random_shuffle( RandomIt first, RandomIt last )
{
	typename std::iterator_traits<RandomIt>::difference_type i, n;
	n = last - first;
	for (i = n-1; i > 0; --i) {
		using std::swap;
		swap(first[i], first[std::rand() % (i+1)]);
	}
}

GeoSource::GeoSource(int iniSize) : oinfo(emptyInfo, iniSize), nodes(iniSize*16), elemSet(iniSize*16),
   layInfo(0, iniSize), coefData(0, iniSize), layMat(0, iniSize), efd(null_eframe, iniSize), csfd(null_eframe, iniSize), cframes(0, iniSize), nfd(null_nframe, iniSize),
   elementLumpingWeights_(1)
{
  decJustCalled=false;
  exitAfterDec=false;
  nGlobNodes = 0;
  numOutInfo = 0;
  outLimit = -1;
  numNodalOutput = 0;
  curCluster = -1;

  nElem = 0;
  nElemFluid = 0;
  nMpcElem = 0;
  phantomFlag = 0;  // 0 => no phantom elems
  numClusters = 0;

  elemTypeNumNodesMap = 0;
  na = 0;
  namax = 0;
  numEframes = 0;
  numNframes = 0;
  numCframes = 0;
  numCSframes = 0;
  numLayInfo = 0;
  numCoefData = 0;
  numLayMat = 0;
  prsflg = 0;

  constpflg = 1;
  constqflg = 1;
  fixedEndM = 1;
  maxGlobNode = 0;
  maxClusNode = 0;
  numProps = 0;

  // set file names to NULL
  conName = NULL;
  geoName = NULL;
  decName = NULL;
  mapName = (char*) "CPUMAP";  // default cpu map
  matchName = NULL;

  // initialize bc's
  numTextDirichlet = 0;
  numTextNeuman = 0;
  numDirichlet = 0;
  numDirichletFluid = 0;
  numNeuman = 0;
  numNeumanModal = 0;
  numIDis = 0;
  numIDisModal = 0;
  numIDis6 = 0;
  numIVel = 0;
  numIVelModal = 0;
  numDampedModes = 0;
  numComplexDirichlet = 0;
  numComplexNeuman = 0;
  numSurfaceDirichlet = 0;
  numSurfaceNeuman = 0;
  numSurfacePressure = 0;
  numSurfaceConstraint = 0;
  numConstraintElementsIeq = 0;
  lmpcflag = true;

  // PITA
  // Initial seed conditions
  PitaIDis6 = PitaIVel6 = 0;
  numPitaIDis6   = 0;
  numTSPitaIDis6 = 0;
  numPitaIVel6   = 0;
  numTSPitaIVel6 = 0;

  subToCPU = 0;
  cpuToCPU = 0;
  subToNode = 0;
  subToElem = 0;
  unsortedSubToElem = 0;
  subToSub = 0;

  // init match data
  numMatchData = 0;
  matchData = 0;
  numGapVecs = 0;
  gapVec = 0;

  cinfo = new ControlInfo;
  claw = 0;

  isShift = false;
  shiftV = 0.0;

  dbc = 0;
  dbcFluid = 0;
  nbc = 0;
  nbcModal = 0;
  textDBC = dbc = textNBC = nbc = iDis = iDisModal = iDis6 = iVel = iVelModal = modalDamping = 0;
  //cvbc = 0; rdbc = 0; iTemp = 0;
  cdbc = cnbc = 0;
  surface_dbc = surface_nbc = 0;
  surface_pres = 0;
  surface_cfe = 0;

  maxattrib = -1;
  optDec = 0;
  optDecCopy = 0;

  numInternalNodes = 0;
  allNumClusElems = 0;
  binaryInput = false;
  binaryInputControlLeft = false;
  binaryOutput = false;
/*
#ifdef DISTRIBUTED
  binaryOutput = true;
#else
  binaryOutput = false;
#endif
*/
  clusToSub = 0;

  mratio = 1.0; // consistent mass matrix

  elemSet.setMyData(true);

  localIndex_ = 0;
}

//----------------------------------------------------------------------
void GeoSource::cleanUp()
{
  elemSet.deleteElems();
  nElem = 0;
  nElemFluid = 0;
}
//----------------------------------------------------------------------

GeoSource::~GeoSource()
{
  /* not fully implemented */
  if(cinfo) delete cinfo;
  if(unsortedSubToElem) delete unsortedSubToElem;
  if(subToCPU) delete [] subToCPU;
  // claw is deleted by domain
  for(int iInfo = 0; iInfo < numOutInfo; ++iInfo) {
	if(oinfo[iInfo].filptr) fclose(oinfo[iInfo].filptr);
	if(oinfo[iInfo].binfilptr) delete oinfo[iInfo].binfilptr;
  }
  if(surface_dbc) delete [] surface_dbc;
  if(surface_nbc) delete [] surface_nbc;
  if(surface_pres) delete [] surface_pres;
  if(surface_cfe) delete [] surface_cfe;
  if(optDec) delete optDec;
  if(optDecCopy) delete optDecCopy;
}

//----------------------------------------------------------------------

int GeoSource::addNode(int nd, double xyz[3], int cp, int cd)
{
  nodes.nodeadd(nd, xyz, cp, cd);
  nGlobNodes++;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addElem(int en, int type, int nn, int *nodeNumbers)
{
#ifndef SALINAS
  elemSet.elemadd(en, type, nn, nodeNumbers);
#else
  std::cerr << "*** ERROR: GeoSource::addElem(...) not included in Salinas library \n";
#endif


  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addMat(int nmat, const StructProp &p)
{
  if (numProps < nmat+1) // attempt to get numProps right -- Julien & Thomas
	numProps = nmat+1;

  sProps[nmat] = p;

  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addLay(int l, LayInfo *linfo)
{
  if (numLayInfo <= l)
	numLayInfo = l+1;

  layInfo[l] = linfo;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addCoefInfo(int cin, CoefData &cdata)  {

  if (numCoefData <= cin)
	numCoefData = cin+1;

  CoefData *cd  = new CoefData(cdata);
  coefData[cin] = cd;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addLayMat(int m, double *d)
{
  if (numLayMat <= m)
	numLayMat = m+1;

  LayMat *lm = new LayMat(m, d);
  layMat[m] = lm;
  return 0;
}

//----------------------------------------------------------------------

int GeoSource::setAttrib(int n, int a, int ca, int cfrm, double ctheta)
{
  if(a>maxattrib) maxattrib = a;
  attrib[na].nele = n;
  attrib[na].attr = a;
  attrib[na].cmp_attr = ca;
  attrib[na].cmp_frm  = cfrm;
  attrib[na].cmp_theta = ctheta;
  na++;
  atoe[a].elems.push_back(n);
  return 0;
}

//----------------------------------------------------------------------

void GeoSource::setMortarAttrib(int n, int a)
{
  if(a>maxattrib) maxattrib = a;
  mortar_attrib[n] = a;
}

//----------------------------------------------------------------------

int GeoSource::setFrame(int el, double *data)
{
  efd[numEframes].elnum = el;
  int i,j;
  for(i = 0; i < 3; ++i)
   for(j = 0; j < 3; ++j)
	 efd[numEframes].frame[i][j] = data[i*3+j];

  numEframes++;

  return 0;
}

//----------------------------------------------------------------------

int GeoSource::setNodalFrame(int id, double *origin, double *data, int type)
{
  int i,j;
  for(i = 0; i < 3; ++i) {
   nfd[id].origin[i] = origin[i];
   for(j = 0; j < 3; ++j)
	 nfd[id].frame[i][j] = data[i*3+j];
  }
  nfd[id].type = (NFrameData::FrameType) type;

  return 0;
}

//----------------------------------------------------------------------

int GeoSource::setCSFrame(int el, double *data)
{
  csfd[numCSframes].elnum = el;
  int i,j;
  for(i = 0; i < 3; ++i)
   for(j = 0; j < 3; ++j)
	 csfd[numCSframes].frame[i][j] = data[i*3+j];

  numCSframes++;

  return 0;
}

//----------------------------------------------------------------------

int GeoSource::addCFrame(int fn, double *f)  {

  if (fn >= numCframes)
	numCframes = fn+1;

  cframes[fn] = new double[9];

  for (int i = 0; i < 9; ++i)
	cframes[fn][i] = f[i];

  return 0;
}

//----------------------------------------------------------------------
bool GeoSource::checkLMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc)
{
  if(verboseFlag && numLMPC && domain->solInfo().dbccheck) 
	filePrint(stderr," ... Checking for MPCs involving constrained DOFs ...\n");
  int count = 0; // number of inequality constraints to be enforced with lagrange multipliers
  for(int i=0; i < numLMPC; ++i) {
	if(lmpc[i]->type == 1 && ((lmpc[i]->lagrangeMult == 1)
	   || (lmpc[i]->lagrangeMult == -1 && domain->solInfo().lagrangeMult))) count++;
	if(!domain->solInfo().dbccheck) continue;
	for(int j=0; j < lmpc[i]->nterms; ++j) {
	  int mpc_node = lmpc[i]->terms[j].nnum;
	  int mpc_dof = lmpc[i]->terms[j].dofnum;
	  for(int k=0; k<numDirichlet; ++k) {
		if(dbc[k].type == BCond::Usdd) continue; // the following treatment is not appropriate for USDD's
		int dbc_node = dbc[k].nnum;
		int dbc_dof = dbc[k].dofnum;
		if((dbc_node == mpc_node) && (dbc_dof == mpc_dof)) {
		  if(!lmpc[i]->isComplex) {
			lmpc[i]->rhs.r_value -= lmpc[i]->terms[j].coef.r_value * dbc[k].val;
			lmpc[i]->terms[j].coef.r_value = 0;
		  }
		  else {
			lmpc[i]->rhs.c_value -= lmpc[i]->terms[j].coef.c_value * dbc[k].val;
			lmpc[i]->terms[j].coef.c_value = 0;
		  }
		  lmpc[i]->original_rhs = lmpc[i]->rhs;
		}
	  }
	  for(int k=0; k<numComplexDirichlet; ++k) {
		int cdbc_node = cdbc[k].nnum;
		int cdbc_dof = cdbc[k].dofnum;
		if((cdbc_node == mpc_node) && (cdbc_dof == mpc_dof)) {
		  if(!lmpc[i]->isComplex) {
			lmpc[i]->rhs.r_value -= lmpc[i]->terms[j].coef.r_value * cdbc[k].reval;
			lmpc[i]->terms[j].coef.r_value = 0;
		  }
		  else {
			lmpc[i]->rhs.c_value -= lmpc[i]->terms[j].coef.c_value * DComplex(cdbc[k].reval, cdbc[k].imval);
			lmpc[i]->terms[j].coef.c_value = 0;
		  }
		  lmpc[i]->original_rhs = lmpc[i]->rhs;
		}
	  }
	 // note: could also eliminate the term from the mpc to simplify further instead of setting coef to zero
	}
  }
  return (count > 0);
}

void GeoSource::transformLMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc)
{
#ifdef USE_EIGEN3
  // transform LMPCs from basic to DOF_FRA coordinates, if necessary
  for(int i = 0; i < numLMPC; ++i) {
	if(lmpc[i]->getSource() == mpc::Lmpc || lmpc[i]->getSource() == mpc::NodalContact) continue;

	// first, fill the lists of term indices corresponding to nodes with DOF_FRA
	std::map<int, std::list<int> > lists;
	for(int j=0; j<lmpc[i]->nterms; ++j) {
	  int nnum = lmpc[i]->terms[j].nnum;
	  if(nodes[nnum]->cd != 0) lists[nnum].push_back(j);
	}

	if(lists.size() == 0) continue;

	std::vector<LMPCTerm> newterms;
	for(std::map<int, std::list<int> >::iterator it = lists.begin(); it != lists.end(); ++it) {
	  int nnum = it->first;

	  // get the original coefficients
	  Eigen::Vector3d tcoefs, rcoefs; tcoefs.setZero(); rcoefs.setZero();
	  for(std::list<int>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
		int dofnum = lmpc[i]->terms[*it2].dofnum;
		switch(dofnum) {
		  case 0: case 1: case 2:
			tcoefs[dofnum] = lmpc[i]->terms[*it2].coef.r_value;
			break;
		  case 3: case 4: case 5:
			rcoefs[dofnum-3] = lmpc[i]->terms[*it2].coef.r_value;
			break;
		}
	  }

	  // apply the transformation
	  int cd = nodes[nnum]->cd;
	  Eigen::Matrix3d T;
	  T << nfd[cd].frame[0][0], nfd[cd].frame[0][1], nfd[cd].frame[0][2],
		   nfd[cd].frame[1][0], nfd[cd].frame[1][1], nfd[cd].frame[1][2],
		   nfd[cd].frame[2][0], nfd[cd].frame[2][1], nfd[cd].frame[2][2];

	  if(!(tcoefs.array() == 0).all()) { // translation dofs
		tcoefs = (T*tcoefs).eval();
		for(int j=0; j<3; ++j) newterms.push_back(LMPCTerm(nnum,j,tcoefs[j]));
	  }
	  if(!(rcoefs.array() == 0).all()) { // rotation dofs
		rcoefs = (T*rcoefs).eval();
		for(int j=0; j<3; ++j) newterms.push_back(LMPCTerm(nnum,j+3,rcoefs[j]));
	  }
	}

	// remove the old terms and insert new terms lmpc[i]
	std::vector<LMPCTerm>::iterator it = lmpc[i]->terms.begin();
	while(it != lmpc[i]->terms.end()) {
	  if(nodes[it->nnum]->cd != 0)
		it = lmpc[i]->terms.erase(it);
	  else
		++it;
	}
	lmpc[i]->terms.insert(lmpc[i]->terms.end(), newterms.begin(), newterms.end());
	lmpc[i]->removeNullTerms();
  }
#endif
}

void GeoSource::addMpcElements(int numLMPC, ResizeArray<LMPCons *> &lmpc)
{
  int nEle = elemSet.last();
  if(numLMPC) {
	for(int i = 0; i < numLMPC; ++i) {
	  elemSet.mpcelemadd(nEle, lmpc[i], (domain->solInfo().isNonLin() || domain->solInfo().gepsFlg == 1));
	  // if constraint options have been set and they are different from the defaults
	  // then create a StructProp and set attribute for this mpc element
	  if((lmpc[i]->lagrangeMult != -1) &&
		 (lmpc[i]->lagrangeMult != domain->solInfo().lagrangeMult || lmpc[i]->penalty != domain->solInfo().penalty)) {
		int a = maxattrib + 1;
		StructProp p;
		p.lagrangeMult = lmpc[i]->lagrangeMult;
		p.initialPenalty = p.penalty = lmpc[i]->penalty;
		p.type = StructProp::Constraint;
		addMat(a, p);
		setAttrib(nEle, a);
	  }
	  nEle++;
	  nMpcElem++;
	}
	//std::cerr << " ... Converted " << numLMPC << " LMPCs to constraint elements ...\n";
  }
  if(!(domain->solInfo().rbmflg == 1 && domain->solInfo().grbm_use_lmpc)) { // still needed for GRBM
	for(int i=0; i<numLMPC; ++i) if(lmpc[i]) delete lmpc[i];
	lmpc.deleteArray();
	lmpc.restartArray();
	domain->setNumLMPC(0);
  }
}

void GeoSource::addMpcElementsIeq(int numLMPC, ResizeArray<LMPCons *> &lmpc)
{
  if(numLMPC) {
	for(int i = 0; i < numLMPC; ++i) {
	  if(lmpc[i]->type == 1) {
		packedEsetConstraintElementIeq->mpcelemadd(numConstraintElementsIeq, lmpc[i], (domain->solInfo().isNonLin() || domain->solInfo().gepsFlg == 1));
		StructProp *p = new StructProp;
		p->lagrangeMult = (lmpc[i]->lagrangeMult == -1) ? domain->solInfo().lagrangeMult : lmpc[i]->lagrangeMult;
		p->initialPenalty = p->penalty = (lmpc[i]->lagrangeMult == -1) ? domain->solInfo().penalty : lmpc[i]->penalty;
		p->type = StructProp::Constraint;
		(*packedEsetConstraintElementIeq)[numConstraintElementsIeq]->setProp(p,true);
		numConstraintElementsIeq++;
	  }
	}
  }
}

void GeoSource::UpdateContactSurfaceElements(DistrGeomState *geomState, std::map<std::pair<int,int>,double> &mu)
{
  SolverInfo &sinfo = domain->solInfo();
  ResizeArray<LMPCons *> &lmpc = *domain->getLMPC();
  int numLMPC = domain->getNumLMPC();

  // copy the Lagrange multipliers from geomState
  mu.clear();
  for(int i = 0; i < numLMPC; ++i) {
	if(lmpc[i]->getSource() == mpc::ContactSurfaces) {
	  if(sProps[mortar_attrib[lmpc[i]->id.first]].lagrangeMult)
		mu.insert(std::pair<std::pair<int,int>,double>(lmpc[i]->id, 0.0));
	}
  }
  geomState->getMultipliers(mu);

  int count = 0;
  int nEle = elemSet.last();
  int count1 = 0;
  int nNode = nodes.size();
  for(int i = 0; i < numLMPC; ++i) {
	if(lmpc[i]->getSource() == mpc::ContactSurfaces) {
	  if(count < contactSurfElems.size()) { // replace
		//cerr << "replacing element " << contactSurfElems[count] << " with lmpc " << i << std::endl;
		elemSet.deleteElem(contactSurfElems[count]);
		elemSet.mpcelemadd(contactSurfElems[count], lmpc[i]); // replace
		elemSet[contactSurfElems[count]]->setProp(&sProps[mortar_attrib[lmpc[i]->id.first]]);
		if(elemSet[contactSurfElems[count]]->numInternalNodes() == 1) { // i.e. lagrange multiplier
		  int in[1] = { nNode++ };
		  elemSet[contactSurfElems[count]]->setInternalNodes(in);
		}
		count1++;
	  }
	  else { // new
		//cerr << "adding lmpc " << i << " to elemset at index " << nEle << std::endl;
		elemSet.mpcelemadd(nEle, lmpc[i]); // new
		elemSet[nEle]->setProp(&sProps[mortar_attrib[lmpc[i]->id.first]]);
		if(elementLumpingWeights_.size() > 0) {
		  for(int locMesh = 0; locMesh < elementLumpingWeights_.size(); locMesh++){
			//fprintf(stderr,"adding %d to weight maps\n", nEle);
			elementLumpingWeights_[locMesh].insert(std::make_pair(nEle,1.0));
		  }
		}
		if(elemSet[nEle]->numInternalNodes() == 1) {
		  int in[1] = { nNode++ };
		  elemSet[nEle]->setInternalNodes(in);
		}
		contactSurfElems.push_back(nEle);
		nEle++;
	  }
	  count++;
	}
  }
  int count2 = 0;
  while(count < contactSurfElems.size()) {
	//cerr << "deleting elemset " << contactSurfElems.back() << std::endl;
	if(elementLumpingWeights_.size() > 0) {
	  for(int locMesh = 0; locMesh < elementLumpingWeights_.size(); locMesh++){
		//fprintf(stderr,"removing %d from weight maps\n", contactSurfElems.back());
		ElementWeightMap::iterator wIt = elementLumpingWeights_[locMesh].find(contactSurfElems.back());
		elementLumpingWeights_[locMesh].erase(wIt);
	  }
	}
	elemSet.deleteElem(contactSurfElems.back());
	contactSurfElems.pop_back();
	count2++;
  }
  //cerr << "replaced " << count1 << " and added " << count-count1 << " new elements while removing " << count2 << std::endl;
  nElem = elemSet.last();

  if(optDecCopy == 0) {
	optDecCopy = optDec;
  }
  else {
	delete optDec;
  }
  modifyDecomposition(nElem-contactSurfElems.size());
}

void
GeoSource::initializeParameters()
{
  SPropContainer::iterator it = sProps.begin();
  while(it != sProps.end()) {
	StructProp* p = &(it->second);
	p->penalty = p->initialPenalty;
	it++;
  }
}

void
GeoSource::updateParameters()
{
  SPropContainer::iterator it = sProps.begin();
  while(it != sProps.end()) {
	StructProp* p = &(it->second);
	p->penalty *= domain->solInfo().penalty_beta;
	it++;
  }
}

// Order the terms in MPCs so that the first term (slave) can be directly written in terms of the others (master)
void GeoSource::makeDirectMPCs(int &numLMPC, ResizeArray<LMPCons *> &lmpc)
{
  if(numLMPC) {
	if(verboseFlag) std::cerr << " ... Using direct elimination method for " << numLMPC << " constraints ...\n";

	using std::map;
	using std::pair;
	using std::vector;
	// Create a unique integer ID for every DOF involved in an MPC and count
	// in how many MPC each DOF appears.
	int nID = 0;
	map<pair<int,int>, int> dofID;

	std::set<pair<int,int> > dispBC;
	for(int i = 0; i < numDirichlet; ++i)
	  dispBC.insert(pair<int, int>(dbc[i].nnum, dbc[i].dofnum));

	for(int i = 0; i < numLMPC; ++i) {
	  // First flush the MPC from any zero terms
	  lmpc[i]->removeNullTerms();
	  // Second flush the MPC from any terms which have a boundary condition
	  vector<LMPCTerm>::iterator it = lmpc[i]->terms.begin();
	  while(it != lmpc[i]->terms.end()) {
		pair<int,int> p(it->nnum, it->dofnum);
		if(dispBC.find(pair<int,int>(it->nnum, it->dofnum)) != dispBC.end())
		  it = lmpc[i]->terms.erase(it);
		else
		  ++it;
	  }
	  // now assign a unique number to each mpc dof
	  lmpc[i]->nterms = lmpc[i]->terms.size();
	  for(int j = 0; j < lmpc[i]->nterms; ++j) {
		pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
		map<pair<int,int>, int>::iterator it = dofID.find(p);
		if(it == dofID.end())
		  dofID[p] = nID++;
	  }
	}

	// Obtain the MPC to DOF connectivity
	SetAccess<LMPCons> lmpcAccess(numLMPC, lmpc.data(), dofID);
	Connectivity lmpcToDof(lmpcAccess);
	Connectivity *dofToLMPC = lmpcToDof.alloc_reverse();
	//std::cerr << "Number of DOFs in MPCS: " << dofToLMPC->csize() << std::endl;
	Connectivity *lmpcToLmpc = lmpcToDof.transcon(dofToLMPC);
	compStruct renumb = lmpcToLmpc->renumByComponent(1);
	delete lmpcToLmpc;
	//std::cerr << "Number of components = " << renumb.numComp << std::endl;

	// Determine for each MPC which DOF will be slave
	vector<int> mpcSlaveDOF(numLMPC);
	for(int i=0; i<numLMPC;++i) mpcSlaveDOF[i] = -1;
	vector<int> dofSlaveOf(dofToLMPC->csize(), -1);
	for(int i=0; i<dofToLMPC->csize(); ++i) dofSlaveOf[i] = -1;
	int nMPCtoView = numLMPC;
	for(int i = 0; i < dofToLMPC->csize(); ++i)
	  if(dofToLMPC->num(i) == 1) {
		int j = (*dofToLMPC)[i][0];
		if(mpcSlaveDOF[j] == -1) {
		  mpcSlaveDOF[j] = i;
		  nMPCtoView--;
		  // swap the slave term to index 0
		  for(int k = 0; k < lmpc[j]->nterms; ++k) {
			pair<int,int> p(lmpc[j]->terms[k].nnum, lmpc[j]->terms[k].dofnum);
			if(dofID[p] == i) {
			  if(k != 0) {
				LMPCTerm t = lmpc[j]->terms[k];
				lmpc[j]->terms[k] = lmpc[j]->terms[0];
				lmpc[j]->terms[0] = t;
			  }
			}
		  }
		}
	  }
	if(nMPCtoView > 0) {
	  //std::cerr << "Could not find a slave for each MPC. Number of MPCs remaining: " << nMPCtoView << std::endl;
	  int *order = new int[numLMPC];
	  for(int i = 0; i < numLMPC; ++i)
		order[i] = -1;
	  for(int i = 0; i < numLMPC; ++i)
		if(renumb.renum[i] >= 0)
		  order[renumb.renum[i]] = i;

	  for(int i = 0; i < renumb.numComp; ++i) {
		int nMpcToViewI = renumb.xcomp[i+1] - renumb.xcomp[i];
		for(int j = renumb.xcomp[i]; j < renumb.xcomp[i+1]; ++j) {
		  int k = order[j];
		  if(mpcSlaveDOF[k] > -1) nMpcToViewI--;
		}
		if(nMpcToViewI > 0) {
		  //std::cerr << " Group " << i << " has " << nMpcToViewI << " MPCs remaining \n";
		  ResizeArray<LMPCons* > lmpcI(0,renumb.xcomp[i+1]-renumb.xcomp[i]);
		  int l = 0;
		  for(int j = renumb.xcomp[i]; j < renumb.xcomp[i+1]; ++j) {
			int k = order[j];
			lmpcI[l++] = lmpc[k];
		  }
		  reduceMPCs(l, lmpcI);
		}
	  }
	  delete [] order;
	}
	//else std::cerr << "Found a Slave for each MPC!"<< std::endl;

	// compress the mpcs by removing the redundant ones
	int j = 0;
	for(int i = 0; i < numLMPC; ++i) {
	  if(lmpc[i]->nterms != 0) {
		if(i != j) lmpc[j] = lmpc[i];
		j++;
	  }
	  else delete lmpc[i];
	}
	if(verboseFlag && j < numLMPC) std::cerr << " ... Found " << numLMPC-j << " redundant constraints ...\n";
	numLMPC = j;
	delete dofToLMPC;
	renumb.clearMemory();
	if(domain->probType() == SolverInfo::Modal) {
	  for(int i = 0; i < numLMPC; ++i) lmpc[i]->rhs.r_value = 0;
	}
  }
}

int
GeoSource::reduceMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc)
{
#ifdef USE_EIGEN3
  using std::map;
  using std::pair;
  // create a unique integer ID for every DOF involved in the MPCs
  int nID = 0;
  map<pair<int,int>, int> dofID;
  for(int i = 0; i < numLMPC; ++i) {
	for(int j = 0; j < 1; ++j) { // number any slaves first
	  pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
	  map<pair<int,int>, int>::iterator it = dofID.find(p);
	  if(it == dofID.end())
		dofID[p] = nID++;
	}
	for(int j = 1; j < lmpc[i]->nterms; ++j) {
	  pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
	  map<pair<int,int>, int>::iterator it = dofID.find(p);
	  if(it == dofID.end())
		dofID[p] = nID++;
	}

  }

  // Obtain the MPC to DOF connectivity
  SetAccess<LMPCons> lmpcAccess(numLMPC, lmpc.data(), dofID);
  Connectivity lmpcToDof(lmpcAccess);
  Connectivity *dofToLMPC = lmpcToDof.alloc_reverse();

  std::vector<int> *term2col = new std::vector<int>[numLMPC];
  std::vector<pair<int,int> > col2pair(dofToLMPC->csize());
  for(int i = 0; i < numLMPC; ++i) {
	for(int j = 0; j < lmpc[i]->nterms; ++j) {
	  pair<int, int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
	  map<pair<int,int>, int>::iterator it = dofID.find(p);
	  term2col[i].push_back(it->second); col2pair[it->second] = p;
	}
  }
  // copy lmpc coefficients into a dense matrix
  int optc = 1; // set to 1 to deal with non-zero rhs
  int n = numLMPC, m = dofToLMPC->csize();
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> c(n, m+optc);
  c.setZero();
  for(int i = 0; i < n; ++i) {
	for(int j = 0; j < lmpc[i]->nterms; ++j) {
	  c(i,term2col[i][j]) = lmpc[i]->terms[j].coef.r_value;
	}
	if(optc) c(i,m) = lmpc[i]->rhs.r_value;
  }
  delete [] term2col;

  /*double t = -getTime();
  std::cerr << " ... Converting " << numLMPC << " LMPCs to Reduced Row Echelon Form";*/
  int *colmap = new int[c.cols()];
  for(int i = 0; i < c.cols(); ++i) colmap[i] = i;
  int rank = rowEchelon<double, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >(c, true, NULL, colmap, optc, domain->solInfo().mpcDirectTol, domain->solInfo().usePrescribedThreshold);
  //cerr << "took " << (t += getTime())/1000. << " seconds ...\n";
  //if(rank != numLMPC)
  //  std::cerr << "rowEchelon detected " << numLMPC-rank << " redundant constraints\n";

  // copy the coefficients of the rref matrix into the lmpc data structure 
  double tol1 = domain->solInfo().coefFilterTol*std::numeric_limits<double>::epsilon(),
		 tol2 = domain->solInfo().rhsZeroTol   *std::numeric_limits<double>::epsilon();
  for(int i = 0; i < n; ++i) {
	lmpc[i]->terms.clear();
	lmpc[i]->nterms = 0;
	lmpc[i]->rhs.r_value = 0;
	if(i >= rank) {
	  if(std::fabs(c(i,m)) > domain->solInfo().inconsistentTol)
		std::cerr << "warning: inconsistent constraint detected (" << c(i,m) << ")\n";
	  continue;
	}
	for(int j = i; j < m; ++j) {
	  if(j > i && j < rank) continue; // for reduced row echelon form these terms are zero by definition
	  if(std::fabs(c(i,j)) > tol1) {  // filter out the very small coefficients
									  // (note: use tol1 == 0) to keep them all
		LMPCTerm t(col2pair[colmap[j]].first, col2pair[colmap[j]].second, c(i,j));
		lmpc[i]->terms.push_back(t);
		lmpc[i]->nterms++;
	  }
	}
	if(optc) {
	  if(colmap[m] != m) std::cerr << "error: mpc rhs was pivoted\n"; // this should not happen
	  if(std::fabs(c(i,m)) > tol2) // set the rhs to exactly zero if it is already very close
		lmpc[i]->rhs.r_value = c(i,m);
	}
  }

  delete dofToLMPC;
  delete [] colmap;
  return rank;
#else
  std::cerr << " *** ERROR: GeoSource::reduceMPCs requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
  return 0;
#endif
}

//----------------------------------------------------------------------

void GeoSource::addFsiElements(int numFSI, ResizeArray<LMPCons *> &fsi)
{
#ifndef SALINAS
  if(numFSI) filePrint(stderr," ... Adding FSI Elements \n");
  int nEle = elemSet.last();
  int i, j;
// JLchange:
//  for(i = 0; i < numFSI; ++i) {
//    elemSet.fsielemadd(nEle, fsi[i]);
//    nEle++;
//  }
// Instead of adding one fsi to ONE element, add each term (between one fluid node and one
// structure node) to a single element.
  for(i = 0; i < numFSI; ++i) {
	int fluidNode = fsi[i]->lmpcnum;
	for(j=0; j< (fsi[i])->nterms; j++) {
	  LMPCons *thisFsi = new LMPCons(fluidNode, 0.0);
	  LMPCTerm thisLmpcTerm((fsi[i])->terms[j], 1.0);
	  thisFsi->addterm(&(thisLmpcTerm));
	  elemSet.fsielemadd(nEle, thisFsi);
	  nEle++;
	}
  }

//  fsi.deleteArray(); domain->setNumFSI(0); // DEBUG fsi_element
#else
  std::cerr << "*** ERROR: GeoSource::addFsiElements(...) not included in Salinas library \n";
#endif
}

//----------------------------------------------------------------------
void GeoSource::duplicateFilesForPita(int localNumSlices, const int* sliceRankSet)
{
  int initialNumFiles = numOutInfo;
  std::vector<char*> initialFileNameSet(initialNumFiles);
  OutputInfo oI;

  if (localNumSlices > 0)
  {
	for (int j = 0; j < initialNumFiles; ++j)
	{
	  oinfo[j].timeSliceRank = sliceRankSet[0];
	  initialFileNameSet[j] = oinfo[j].filename;
	}
	timeSliceOutputFiles.insert(std::make_pair(sliceRankSet[0], std::pair<int, int>(0, initialNumFiles)));
  }

  // If necessary, create a whole new set of output files
  // for each time-slice beyond the first on the local CPU
  // and initialize the timeSliceRank tags
  for (int i = 1; i < localNumSlices; ++i)
  {
	int firstRequest = numOutInfo;
	for (int j = 0; j < initialNumFiles; ++j)
	{
	  oI.copyParam(oinfo[j]);
	  oI.timeSliceRank = sliceRankSet[i];
	  addOutput(oI);
	}
	int lastRequest = numOutInfo;
	timeSliceOutputFiles.insert(std::make_pair(sliceRankSet[i], std::make_pair(firstRequest, lastRequest)));
  }

  // Append the time-slice rank to the output file names
  for (int i = 0; i < numOutInfo; ++i)
  {
	std::stringstream s;
	s << oinfo[i].filename << '.' << oinfo[i].timeSliceRank;
	const std::string & newFileNameString = s.str();

	int newFileNameLength = newFileNameString.size();
	char* newFileName = (char*) malloc(sizeof(char) * (newFileNameLength + 1));
	newFileNameString.copy(newFileName, newFileNameLength);
	newFileName[newFileNameLength] = '\0';

	oinfo[i].filename = newFileName;
  }

  for (int i = 0; i < initialNumFiles; ++i)
  {
	if (initialFileNameSet[i] != 0)
	{
	  free(initialFileNameSet[i]);
	  initialFileNameSet[i] = 0;
	}
  }
}

void GeoSource::transformCoords()
{
#ifdef USE_EIGEN3
  // transform node coordinates from POS_FRM to basic coordinates
  // see: http://en.wikipedia.org/wiki/List_of_canonical_coordinate_transformations#3-Dimensional
  using std::cos;
  using std::sin;
  SolverInfo &sinfo = domain->solInfo();
  int lastNode = numNodes = nodes.size();
  for(int i=0; i<lastNode; ++i) {  // loop over all nodes
	if(nodes[i] == NULL) continue; // if node pointer is empty, skip
	int cp = nodes[i]->cp;         // get coordinate frame type
	if(cp == 0) {                  // if eulerian scale and continue
	  nodes[i]->x *= sinfo.xScaleFactor;
	  nodes[i]->y *= sinfo.yScaleFactor;
	  nodes[i]->z *= sinfo.zScaleFactor;
	  continue;
	}

	// otherwise apply appropriate transformation
	Eigen::Vector3d v;
	switch(nfd[cp].type) {
	  default:
	  case NFrameData::Rectangular: {
		v << nodes[i]->x, nodes[i]->y, nodes[i]->z;
	  } break;
	  case NFrameData::Cylindrical: {
		double rho   = nodes[i]->x;
		double theta = nodes[i]->y;
		v << rho*cos(theta), rho*sin(theta), nodes[i]->z;
	  } break;
	  case NFrameData::Spherical: {
		double rho   = nodes[i]->x;
		double theta = nodes[i]->y;
		double phi   = nodes[i]->z;
		v << rho*sin(theta)*cos(phi), rho*sin(theta)*sin(phi), rho*cos(theta);
	  } break;
	}

	Eigen::Matrix3d T;
	T << nfd[cp].frame[0][0], nfd[cp].frame[0][1], nfd[cp].frame[0][2],
		 nfd[cp].frame[1][0], nfd[cp].frame[1][1], nfd[cp].frame[1][2],
		 nfd[cp].frame[2][0], nfd[cp].frame[2][1], nfd[cp].frame[2][2];

	v = (T.transpose()*v).eval();

	nodes[i]->x = sinfo.xScaleFactor*(v[0] + nfd[cp].origin[0]);
	nodes[i]->y = sinfo.yScaleFactor*(v[1] + nfd[cp].origin[1]);
	nodes[i]->z = sinfo.zScaleFactor*(v[2] + nfd[cp].origin[2]);
  }
#endif
}

void GeoSource::setNewCoords(std::string nodeFile)
{
  // this function is for setting new xyz locations for a mesh 
  // during an HROM training with training vectors from multiple
  // simulations in which the mesh deformation is a parameter

  // open file stream to file containing nodes numbers and locations
  std::ifstream newNodes; 
  newNodes.open(nodeFile.c_str());
  int lastNode = numNodes = nodes.size();
  
  // now loop over nodes
  for(int i=0; i<lastNode; ++i) {
	 if(nodes[i] == NULL) continue; // if node pointer is empty, skip
	 int nnum;
	 double x,y,z;
	 newNodes >> nnum; newNodes >> x; newNodes >> y; newNodes >> z;
	 nnum--;
	 assert(nnum == i);
	 nodes[i]->x = x;
	 nodes[i]->y = y;
	 nodes[i]->z = z;
  } 
}

void GeoSource::setUpData(int topFlag)
{
  using std::map;
  using std::list;
  int lastNode = numNodes = nodes.size();
  const int nMaxEle = elemSet.last();

  // Check that the selected load cases have been defined
  domain->checkCases();

  // Set up element pressure load
  BlastLoading::BlastData *conwep = (domain->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  for(ElemPressureContainer::iterator i = eleprs.begin(); i != eleprs.end(); ++i) {
	PressureBCond &pbc = *i;
	int elemNum = i->elnum;
	if(elemSet[elemNum]) {
	  pbc.mftt = domain->getMFTT(pbc.loadsetid);
	  pbc.conwep = conwep;
	  pbc.loadfactor = domain->getLoadFactor(pbc.loadsetid);
	  if(pbc.face == -1)
		elemSet[elemNum]->setPressure(&pbc);
	  else {
		int *nodes = new int[elemSet[elemNum]->numNodes()];
		int nNodes = elemSet[elemNum]->getFace(pbc.face, nodes);
		if(nNodes > 0) {
		  int type;
		  switch(nNodes) {
			case 3 : type = 15; break;
			case 4 : type = 16; break;
			case 6 : type = 17; break;
			case 8 : type = 18; break;
			case 9 : type = 19; break;
			case 10 : type = 21; break;
			case 12 : type = 20; break;
		  }
		  domain->addNeumElem(-1, type, pbc.val, nNodes, nodes, &pbc);
		  domain->neum[domain->numNeum-1]->setAdjElementIndex(elemNum);
		  delete [] nodes;
		}
		else
		  filePrint(stderr, " *** WARNING: Pressure was found for unsupported element %d\n", elemNum+1);
	  }
	}
	else
	  filePrint(stderr, " *** WARNING: Pressure was found for non-existent element %d\n", elemNum+1);
  }

  // Set up element preload
  for(ElemPreloadContainer::iterator i = eleprl.begin(); i != eleprl.end(); ++i) {
	int elemNum = i->first;
	if(elemSet[elemNum])
	  elemSet[elemNum]->setPreLoad(i->second);
	else
	  filePrint(stderr, " *** WARNING: Preload was found for non-existent element %d\n", elemNum+1);
  }

  // Set up element frames
  for (int iFrame = 0; iFrame < numEframes; iFrame++)  {
	Element *ele = elemSet[efd[iFrame].elnum];
	if(ele == 0) {
	  filePrint(stderr, " *** WARNING: Frame was found for non-existent element %d \n", efd[iFrame].elnum+1);
	}
	else {
	  ele->setFrame(&(efd[iFrame].frame));
	}
  }
  for (int i = 0; i < nMaxEle; ++i)
	if(elemSet[i]) elemSet[i]->buildFrame(nodes);

  // Set up element attributes
  SolverInfo &sinfo = domain->solInfo();
/*if((na == 0) && (sinfo.probType != SolverInfo::Top) && (sinfo.probType != SolverInfo::Decomp)) {
	filePrint(stderr," **************************************\n");
	filePrint(stderr," *** ERROR: ATTRIBUTES not defined  ***\n");
	filePrint(stderr," **************************************\n");
	exit(-1);
  }*/

  // check for elements with no attribute, and add dummy properties in certain cases
  bool *hasAttr = new bool[nMaxEle];
  for(int i = 0; i < nMaxEle; ++i) hasAttr[i] = false;
  for(std::map<int,Attrib>::iterator it = attrib.begin(); it != attrib.end(); ++it) {
	Attrib &attrib_i = it->second;
	if(attrib_i.nele < nMaxEle) {
	  hasAttr[attrib_i.nele] = true;
	  elemSet[attrib_i.nele]->setElementAttribute(attrib_i.attr);
	}
  }
  int dattr, dattr2;
  bool hasAddedDummy = false, hasAddedDummy2 = false, hasMultiplier = false;
  for(int i = 0; i < nMaxEle; ++i) {
	if(elemSet[i] && !hasAttr[i]) {
	  if(!elemSet[i]->isRigidElement() && !hasAddedDummy) {
		dattr = ++maxattrib;
		addMat(dattr, StructProp());
		hasAddedDummy = true;
	  }
	  else if(elemSet[i]->isRigidElement() && !hasAddedDummy2) {
		// for rigid elements the default density needs to be changed to zero
		dattr2 = ++maxattrib;
		StructProp p;
		p.rho = 0.0;
		addMat(dattr2, p);
		hasAddedDummy2 = true;
	  }
	  if(sinfo.probType == SolverInfo::Top || sinfo.probType == SolverInfo::Decomp || elemSet[i]->isConstraintElement()
		 || elemSet[i]->isSloshingElement()) {
		if(!elemSet[i]->isRigidElement()) setAttrib(i, dattr);
		else setAttrib(i, dattr2);
	  }
	  else filePrint(stderr, " *** WARNING: Element %d has no attribute defined\n", i+1);
	}
  }
  delete [] hasAttr;

  // add properties for mortar conditions
  bool printOne = true, printTwo = true;
  bool checkConstraintMethod = (sinfo.probType != SolverInfo::Top && sinfo.probType != SolverInfo::Decomp &&
								sinfo.probType != SolverInfo::PodRomOffline && !domain->solInfo().ROMPostProcess);
  for(int i=0; i<domain->GetnMortarConds(); ++i) {
	MortarHandler* CurrentMortarCond = domain->GetMortarCond(i);
	ConstraintOptions *copt = CurrentMortarCond->GetConstraintOptions();
	if(copt) {
	  int a = maxattrib + 1;
	  StructProp p;
	  p.lagrangeMult = copt->lagrangeMult;
	  p.initialPenalty = p.penalty = copt->penalty;
	  p.type = StructProp::Constraint;
	  addMat(a, p);
	  setMortarAttrib(domain->GetMortarCond(i)->GetId(), a);
	}
	else {
	  if(!hasAddedDummy) {
		dattr = maxattrib + 1;
		addMat(dattr, StructProp());
		hasAddedDummy = true;
	  }
	  setMortarAttrib(domain->GetMortarCond(i)->GetId(), dattr);
	}
	// check constraint method for tied/contact surfaces in cases where mpc elements are not initially generated
	// i.e. nonlinear contactsurfaces with "tdenforce off", and linear/nonlinear tied/contactsurfaces with "tdenforce on".
	if(checkConstraintMethod) {
	  bool lagrangeMult = (copt) ? copt->lagrangeMult : sinfo.lagrangeMult;
	  double penalty = (copt) ? copt->penalty : sinfo.penalty;
	  if(printOne && !domain->tdenforceFlag() && CurrentMortarCond->GetInteractionType() == MortarHandler::CTC &&
		 sinfo.isNonLin() && sinfo.newmarkBeta == 0 && (sinfo.mpcDirect || lagrangeMult == true || penalty == 0)) {
		filePrint(stderr, "\x1B[31m *** WARNING: Constraint method is not \n"
						  "     supported for explicit dynamics.  \x1B[0m\n");
		printOne = false;
	  }
	  else if(printOne && domain->tdenforceFlag() && (sinfo.mpcDirect || (lagrangeMult == false && penalty == 0) ||
		 (lagrangeMult == true && penalty != 0))) {
		filePrint(stderr, "\x1B[31m *** WARNING: Constraint method is not \n"
						  "     supported for explicit dynamics.  \x1B[0m\n");
		printOne = false;
	  }
	}
  }

  num_arubber = 0;
  // assign default properties
  SPropContainer::iterator it = sProps.begin();
  while(it != sProps.end()) {
	StructProp* p = &(it->second);
// RT : 02/01/2013 changing syntax to AMAT
//    if(p->soundSpeed == 1.0)
//      p->soundSpeed = omega()/complex<double>(p->kappaHelm, p->kappaHelmImag);
	complex<double> ka = (shiftVal() >= 0) ? omega()/p->soundSpeed : sqrt(complex<double>(shiftVal()))/p->soundSpeed;
	p->kappaHelm = real(ka);
	p->kappaHelmImag = imag(ka);
	domain->updateSDETAF(p,omega()+1e-9);
	domain->updateRUBDAFT(p,omega()+1e-9);
	if(p->E0!=0.0 || p->mu0!=0.0) num_arubber++;
	if(p->type != StructProp::Constraint) {
	  p->lagrangeMult = (sinfo.mpcDirect) ? false : sinfo.lagrangeMult;
	  p->initialPenalty = p->penalty = (sinfo.mpcDirect) ? 0.0 : sinfo.penalty;
	}
	p->constraint_hess = sinfo.constraint_hess;
	p->constraint_hess_eps = sinfo.constraint_hess_eps;
	it++;
  }
  if (sinfo.doFreqSweep) {
	domain->solInfo().curSweepParam = 0;
	domain->setFrequencySet(0);
  }

  // Set up material properties
  // & compute average E and nu (for coupled_dph)
  // & set the mpc data
  int structure_element_count = 0;
  int fluid_element_count = 0;
  global_average_E = 0.0;
  global_average_nu = 0.0;
  global_average_rhof = 0.0;
  int solverClass = sinfo.classifySolver();
  for(std::map<int,Attrib>::iterator it = attrib.begin(); it != attrib.end(); ++it) {
	Attrib &attrib_i = it->second;
	Element *ele = elemSet[ attrib_i.nele ];

	// Check if element exists
	if(ele == 0) {
	  filePrint(stderr, " *** WARNING: Attribute was found for non existent element %d\n", attrib_i.nele+1);
	  continue;
	}
	if(attrib_i.attr < -1) { // phantom elements
	  phantomFlag = 1;
      if(ele->isConstraintElement()) {
        filePrint(stderr, " *** ERROR: Phantom attribute was found for rigid element %d\n", attrib_i.nele+1);
        exit(-1);
      }
	  ele->setProp(0);
	}
	else {
	  SPropContainer::iterator it = sProps.find(attrib_i.attr);
	  if(it == sProps.end()) {
        if(topFlag != 2 && topFlag != 7 && !(callDec && exitAfterDec)) 
		  filePrint(stderr, " *** WARNING: The material for element %d does not exist\n", attrib_i.nele+1);
	  }
	  else {
		StructProp *prop = &(it->second);
		ele->setProp(prop);

		// compute global average structural and fluid properties
		if(!ele->isConstraintElement()) {
		  if(! dynamic_cast<HelmElement *>(elemSet[attrib_i.nele])) { // not a fluid element
			global_average_E += prop->E;
			global_average_nu += prop->nu;
			structure_element_count++;
		  } else {
			global_average_rhof += prop->rho;
			fluid_element_count++;
		  }
		}
		else if (prop->relop != 0) numConstraintElementsIeq++;

		// check constraint method
		if(checkConstraintMethod && ele->isMpcElement()) {
		  if(printOne && sinfo.newmarkBeta == 0 && (prop->lagrangeMult == true || prop->penalty == 0)) {
			// for explicit dynamics, the penalty constraint method with non-zero penalty parameter must be used
			// with the exception of tied/contact surfaces using "tdenforce on". In the latter case no elements
			// are created so this error message will not be encountered.
			filePrint(stderr, "\x1B[31m *** WARNING: Constraint method is not \n"
							  "     supported for explicit dynamics.  \x1B[0m\n");
			printOne = false;
		  }
		  else if(printTwo && sinfo.newmarkBeta != 0 && sinfo.mpcDirect == 0 && prop->lagrangeMult == true && prop->penalty == 0 && solverClass == 0) {
			// for statics, impe, eigen and implicit dynamics, certain solvers cannot be used with the multipliers
			// constraint method. This check could/should be more specific by considering inequality constraints.
			filePrint(stderr, "\x1B[31m *** WARNING: Equation solver unsuited \n"
							  "     for multipliers constraint method.\x1B[0m\n");
			printTwo = false;
		  }
		}
	  }
	}

	if(attrib_i.cmp_attr >= 0) {
	  if(coefData[attrib_i.cmp_attr] != 0) {
		int type = (coefData[attrib_i.cmp_attr]->coefFlag) ? 5 : 1;
		if(attrib_i.cmp_frm > -1) { // cframe
		  ele->setCompositeData(type, 0, 0, coefData[attrib_i.cmp_attr]->values(),
								cframes[attrib_i.cmp_frm]);
		}
		else { // ctheta
		  ele->setCompositeData2(type, 0, 0, coefData[attrib_i.cmp_attr]->values(),
								 nodes, attrib_i.cmp_theta);
		}
	  }
	  else {
		LayInfo *li = layInfo[attrib_i.cmp_attr];
		if(li == 0) {
		  filePrint(stderr, " *** WARNING: Attribute found that refers to nonexistant composite data: %d\n", attrib_i.cmp_attr+1);
		  continue;
		}
		// Set up layer material properties if necessary
		for(int k=0; k<li->nLayers(); ++k) {
		  int mid = li->getLayerMaterialId(k);
		  if(mid > -1) {
			LayMat *lmk = layMat[mid];
			li->setLayerMaterialProperties(k,lmk->data);
		  }
		}
		// type is 3 for LAYC 2 for LAYN
		int type = 3 - li->getType();
		if(attrib_i.cmp_frm > -1) { // cframe
		  ele->setCompositeData(type, li->nLayers(), li->values(), 0,
								cframes[attrib_i.cmp_frm]);
		}
		else { // ctheta
		  ele->setCompositeData2(type, li->nLayers(), li->values(), 0,
								 nodes, attrib_i.cmp_theta);
		}
	  }
	}
	else if(attrib_i.cmp_frm > -1) { // cframe
	  ele->setCompositeData(1, 0, 0, 0, cframes[attrib_i.cmp_frm]);
	}
    else if(attrib_i.cmp_frm == -2) { // ctheta
      ele->setCompositeData2(1, 0, 0, 0, nodes, attrib_i.cmp_theta);
    }
  }
  if(structure_element_count > 0) {
	global_average_E /= double(structure_element_count);
	global_average_nu /= double(structure_element_count);
  }
  if(fluid_element_count > 0) {
	global_average_rhof /= double(fluid_element_count);
  }

  // Set up beam element offsets
  for(std::vector<OffsetData>::iterator offIt = offsets.begin(); offIt != offsets.end(); ++offIt) {
	for(int i = offIt->first; i <= offIt->last; ++i) {
	  if(elemSet[i] == 0) {
	filePrint(stderr, " *** WARNING: Setting up offset on non-existent element %d ***\n", i+1);
	  } else
	elemSet[i]->setOffset(offIt->o);
	}
  }

  // Set up material data in elements
  std::pair<int, ResizeArray<MFTTData*>* > *ysst = domain->getYSST();
  for(int i = 0; i < ysst->first; ++i) {
	for(map<int, NLMaterial *>::iterator matIter = materials.begin(); matIter != materials.end(); ++matIter)
	  matIter->second->setSDProps((*ysst->second)[i]);
  }
  delete ysst;
  std::pair<int, ResizeArray<MFTTData*>* > *yssrt = domain->getYSSRT();
  for(int i = 0; i < yssrt->first; ++i) {
	for(map<int, NLMaterial *>::iterator matIter = materials.begin(); matIter != materials.end(); ++matIter)
	  matIter->second->setSRDProps((*yssrt->second)[i]);
  }
  delete yssrt;
  std::pair<int, ResizeArray<MFTTData*>* > *ss1dt = domain->getSS1DT();
  for(int i = 0; i < ss1dt->first; ++i) {
	for(map<int, NLMaterial *>::iterator matIter = materials.begin(); matIter != materials.end(); ++matIter)
	  matIter->second->setS1DProps((*ss1dt->second)[i]);
  }
  delete ss1dt;
  std::pair<int, ResizeArray<SS2DTData*>* > *ss2dt = domain->getSS2DT();
  for(int i = 0; i < ss2dt->first; ++i) {
	for(map<int, NLMaterial *>::iterator matIter = materials.begin(); matIter != materials.end(); ++matIter)
	  matIter->second->setS2DProps((*ss2dt->second)[i]);
  }
  delete ss2dt;
  std::pair<int, ResizeArray<MFTTData*>* > *ymst = domain->getYMST();
  for(int i = 0; i < ymst->first; ++i) {
    for(map<int, NLMaterial *>::iterator matIter = materials.begin(); matIter != materials.end(); ++matIter)
      matIter->second->setEDProps((*ymst->second)[i]);
  }
  delete ymst;
  map<int,int>::iterator uIter = matUsage.begin();
  while(uIter != matUsage.end()) {
	int elemNum = uIter->first;
	int matNum = uIter->second;
	Element *ele = elemSet[ elemNum ];

	// Check if element exists
	if (ele == 0) {
	  filePrint(stderr, " *** WARNING: Material was found for non existent element %d \n", elemNum+1);
	  uIter++; // go to the next mapping
	  continue;
	}
	map<int, NLMaterial *>::iterator matIter = materials.find(matNum);
	if(matIter == materials.end()) {
	  filePrint(stderr, " *** WARNING: Non Existent material (%d) was assigned to element %d \n", matNum+1, elemNum+1);
	} else
	  ele->setMaterial(matIter->second);
	uIter++;
  }
/*
  SPropContainer::iterator it = sProps.begin();
  while(it != sProps.end()) {
	StructProp* p = &(it->second);
	if(p->soundSpeed == 1.0)
	  p->soundSpeed = omega()/complex<double>(p->kappaHelm, p->kappaHelmImag);
	it++;
  }
*/

  // Set up the internal nodes
  for (int iElem = 0; iElem < nMaxEle; ++iElem) {
	Element *ele = elemSet[iElem];
	if(ele == 0) continue;
	int nIN = ele->numInternalNodes();
	if(nIN > 0) {
	  int *nn = new int[nIN];
	  for(int i = 0; i < nIN; ++i)
		nn[i] = lastNode++;
	  ele->setInternalNodes(nn);
	  delete [] nn;
	}
	localToGlobalElementsNumber.push_back(iElem);
  }
  numInternalNodes = lastNode-numNodes;

  // preprocess the surface node groups
  ResizeArray<SurfaceEntity*> *SurfEntities = domain->viewSurfEntities();
  for(map<int, list<int> >::iterator i = surfaceGroup.begin(); i != surfaceGroup.end(); ++i) {
	for(list<int>::iterator j = i->second.begin(); j != i->second.end(); ++j) {
	  for(int k = 0; k < domain->getNumSurfs(); k++) {
		if((*SurfEntities)[k]->ID()-1 == *j) {
		  for(int l = 0; l < (*SurfEntities)[k]->GetnNodes(); ++l)
			nodeGroup[i->first].insert((*SurfEntities)[k]->GetPtrGlNodeIds()[l]);
		}
	  }
	}
  }

  for (int iOut = 0; iOut < numOutInfo; iOut++) {

	switch (oinfo[iOut].type) {
	  case OutputInfo::Statevector :
		if(verboseFlag) filePrint(stderr, " ... Saving state snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().statevectPodRom = true;
		domain->solInfo().skipState = oinfo[iOut].interval;
		domain->solInfo().statePodRomFile.push_back(std::string(oinfo[iOut].filename));
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::Residual :
		if(verboseFlag) filePrint(stderr, " ... Saving residual snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().residvectPodRom = true;
		domain->solInfo().skipResidual = oinfo[iOut].interval;
		domain->solInfo().residualPodRomFile = oinfo[iOut].filename;
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::Jacobian :
		if(verboseFlag) filePrint(stderr, " ... Saving jacobian snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().jacobvectPodRom = true;
		domain->solInfo().skipJacobian = oinfo[iOut].interval;
		domain->solInfo().jacobianPodRomFile = oinfo[iOut].filename;
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::RobData :
		if(verboseFlag) filePrint(stderr, " ... Reduced Order Basis Construction: saving to %s ...\n",
								  oinfo[iOut].filename);
		domain->solInfo().SVDoutput = oinfo[iOut].filename;
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::SampleMesh :
		if(verboseFlag) filePrint(stderr, " ... Computing Hyper-Reduction Coefficients: saving to %s ...\n",
								  oinfo[iOut].filename);
		domain->solInfo().reducedMeshFile = oinfo[iOut].filename;
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::Accelvector :
		if(verboseFlag) filePrint(stderr, " ... Saving acceleration snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().accelvectPodRom = true;
		domain->solInfo().skipAccel = oinfo[iOut].interval;
		domain->solInfo().accelPodRomFile.push_back(std::string(oinfo[iOut].filename));
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::Velocvector :
		if(verboseFlag) filePrint(stderr, " ... Saving velocity snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().velocvectPodRom = true;
		domain->solInfo().skipVeloc = oinfo[iOut].interval;
		domain->solInfo().velocPodRomFile.push_back(std::string(oinfo[iOut].filename));
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::InternalStateVar :
		if(verboseFlag) filePrint(stderr, " ... Saving internal state variables snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().isvPodRom = true;
		domain->solInfo().skipInternalStateVar = oinfo[iOut].interval;
		domain->solInfo().isvPodRomFile = oinfo[iOut].filename;
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::DualStateVar :
		if(verboseFlag) filePrint(stderr, " ... Saving dual state variables snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().dsvPodRom = true;
		domain->solInfo().skipDualStateVar = oinfo[iOut].interval;
		domain->solInfo().dsvPodRomFile.push_back(oinfo[iOut].filename);
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::MuStateVar :
		if(verboseFlag) filePrint(stderr, " ... Saving Lagrange Multiplier snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().muvPodRom = true;
		domain->solInfo().skipMuStateVar = oinfo[iOut].interval;
		domain->solInfo().muvPodRomFile.push_back(oinfo[iOut].filename);
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::Forcevector :
		if(verboseFlag) filePrint(stderr, " ... Saving force snapshots every %d time steps to %s ...\n",
								  oinfo[iOut].interval, oinfo[iOut].filename);
		domain->solInfo().activatePodRom = true;
		domain->solInfo().snapshotsPodRom = true;
		domain->solInfo().forcevectPodRom = true;
		domain->solInfo().skipForce = oinfo[iOut].interval;
		domain->solInfo().forcePodRomFile = oinfo[iOut].filename;
		oinfo[iOut].PodRomfile = true;
		break;
	  case OutputInfo::Constraintvector :
		if(verboseFlag) filePrint(stderr," ... Saving constraint snapshots to %s ...\n", oinfo[iOut].filename);
		domain->solInfo().ConstraintBasisPod = true;
		domain->solInfo().constraintSnapshotFile = oinfo[iOut].filename;
	  break;
	  case OutputInfo::Constraintviolation :
		if(verboseFlag) filePrint(stderr," ... Saving constraint violations to %s ...\n", oinfo[iOut].filename);
		domain->solInfo().constraintViolationFile = oinfo[iOut].filename;
	  break;
	  default :
		break;
	}

	if (oinfo[iOut].groupNumber > 0)  {
	  if (topFlag < 0 && nodeGroup.find(oinfo[iOut].groupNumber) == nodeGroup.end())
		filePrint(stderr, " *** WARNING: Requested group output id not found: %d\n", oinfo[iOut].groupNumber);
	}
  }
}

//----------------------------------------------------------------

int GeoSource::getNodes(CoordSet &nds)  {

  nds = nodes;

  // the routines calling this expect to know the total number of nodes
  // including internally created ones.
  return numNodes+numInternalNodes;
}

CoordSet& GeoSource::GetNodes() { return nodes; }

//----------------------------------------------------------------

/* If nElems and elemList exist, then we are getting Elems for the
   binary case.  In this case, we will not need to create a glToPck
   since this is only used when the decomposition from text input
   Also, the packed ElemSet will have sorted the phantom and real
   elements, so that the phantom elements will be at the end of the list */

int GeoSource::getElems(Elemset &packedEset, int nElems, int *elemList)
{
  SolverInfo &sinfo = domain->solInfo();
  bool flagCEIeq = (!sinfo.readInDualROB.empty() || sinfo.ConstraintBasisPod);

  if(sinfo.HEV) { packedEsetFluid = new Elemset(); nElemFluid = 0; }
  if(flagCEIeq) { packedEsetConstraintElementIeq = new Elemset(); numConstraintElementsIeq = 0; lmpcflag = true; }

  int iEle, numele;

  if(nElems) numele = nElems;
  else numele = elemSet.last();

  int packFlag = 0;
  if(!nElems) packFlag = 1;

  nElem = 0; // counting the phantom elements as well
  int nPhantoms = 0;

  // add real elements to list
  glToPckElems.clear();
  for(iEle = 0; iEle < numele; ++iEle) {
	Element *ele = elemSet[iEle];
	if(ele) {
	  if(!ele->isPhantomElement()) {
		if((!sinfo.HEV || !ele->isHEVFluidElement()) && (!flagCEIeq || !ele->isConstraintElementIeq())) {
		  packedEset.elemadd(nElem, ele);
		  if(packFlag)
			glToPckElems[iEle] = nElem;
		  nElem++;
		}
		else if(sinfo.HEV && ele->isHEVFluidElement()) {
		  packedEsetFluid->elemadd(nElemFluid, ele);
		  nElemFluid++;
		}
		else if(flagCEIeq && ele->isConstraintElementIeq()) {
		  MpcElement *mpcele = dynamic_cast<MpcElement*>(ele);
		  if(mpcele->functionType() != MpcElement::LINEAR) lmpcflag = false;
		  packedEsetConstraintElementIeq->elemadd(numConstraintElementsIeq, ele);
		  numConstraintElementsIeq++;
		}
	  }
	  else
		nPhantoms++;
	}
  }

  if(flagCEIeq) {
	// convert inequality LMPCs to mpc elements
	addMpcElementsIeq(domain->getNumLMPC(), *domain->getLMPC());
  }

  // add sommer elements
  if(!domain->getSowering()) {
	if (sinfo.ATDARBFlag != -2.0) {
	  for (int i = 0; i < domain->numSommer; i++) {
		packedEset.elemadd(nElem,domain->sommer[i]);
		if(packFlag)
		  glToPckElems[numele+i] = nElem;
		nElem++;
	  }
	}
  }

  int nRealElems = nElem;

  // set number of real elements
  packedEset.setEmax(nElem);

  if(sinfo.HEV) {
	packedEsetFluid->setEmax(nElemFluid);
  }

  // add phantom elements to list
  iEle = 0;
  int nPhants = 0;
  while (nPhants < nPhantoms) {
	Element *ele = elemSet[iEle];
	if(ele) {
	  if(ele->isPhantomElement()) {
		packedEset.elemadd(nElem, ele);
		if(packFlag)
		  glToPckElems[iEle] = nElem;
		nElem++;
		nPhants++;
	  }
	}
	iEle++;
  }

  return nRealElems;
}

void GeoSource::setElemTypeMap()
{
  // build elem type to num nodes map
  elemTypeNumNodesMap = new int[89];
  for (int i = 0; i < 89; i++)
	elemTypeNumNodesMap[i] = -1;

  elemTypeNumNodesMap[1] = 2;
  elemTypeNumNodesMap[2] = 4;
  elemTypeNumNodesMap[3] = 4;
  elemTypeNumNodesMap[4] = 4;
  elemTypeNumNodesMap[6] = 2;
  elemTypeNumNodesMap[7] = 3;
  elemTypeNumNodesMap[8] = 3;
  elemTypeNumNodesMap[9] = 2;
  elemTypeNumNodesMap[10] = 4;
  elemTypeNumNodesMap[11] = 1;
  elemTypeNumNodesMap[17] = 8;
  elemTypeNumNodesMap[18] = 4;
  elemTypeNumNodesMap[19] = 3;
  elemTypeNumNodesMap[20] = 3;
  elemTypeNumNodesMap[21] = 2;
  elemTypeNumNodesMap[22] = 2;
  elemTypeNumNodesMap[23] = 4;
  elemTypeNumNodesMap[24] = 5;
  elemTypeNumNodesMap[25] = 10;
  elemTypeNumNodesMap[30] = 4;
  elemTypeNumNodesMap[31] = 4;
  elemTypeNumNodesMap[35] = 3;
  elemTypeNumNodesMap[36] = 3;
  elemTypeNumNodesMap[40] = 4;
  elemTypeNumNodesMap[41] = 4;
  elemTypeNumNodesMap[52] = 8;
  elemTypeNumNodesMap[60] = 4;
  elemTypeNumNodesMap[61] = 3;
  elemTypeNumNodesMap[88] = 4; //HB

  //return elemTypeNumNodesMap;
}

void GeoSource::setElementPressure(PressureBCond& pbc)
{
 prsflg = 1;
 eleprs.push_back(pbc);
}

void GeoSource::setElementPreLoad(int elemNum, double _preload)
{
 std::vector<double> preload; 
 preload.push_back(_preload);
 eleprl.push_back(std::pair<int,std::vector<double> >(elemNum,preload));
}

void GeoSource::setElementPreLoad(int elemNum, double _preload[3])
{
 std::vector<double> preload;
 for(int i=0; i<3; ++i) preload.push_back(_preload[i]);
 eleprl.push_back(std::pair<int,std::vector<double> >(elemNum,preload));
}

void GeoSource::setConsistentPFlag(int _constpflg)
{
 constpflg = _constpflg;
}

void GeoSource::setConsistentQFlag(int _constqflg, int _fixedEndM)
{
 constqflg = _constqflg;
 fixedEndM = _fixedEndM;
}

int GeoSource::getTextDirichletBC(BCond *&bc)
{
  bc = textDBC;
  return numTextDirichlet;
}

int GeoSource::getTextNeumanBC(BCond *&bc)
{
  bc = textNBC;
  return numTextNeuman;
}

int GeoSource::getDirichletBC(BCond *&bc)
{
  bc = dbc;
  return numDirichlet;
}

int GeoSource::getDirichletBCFluid(BCond *&bc)
{
  if (dbcFluid) {
	bc = dbcFluid;
	return numDirichletFluid;
  }
  else {
   bc = 0;
   return 0;
  }
}

int GeoSource::getNeumanBC(BCond *&bc)
{
  bc = nbc;
  return numNeuman;
}

int GeoSource::getNeumanBCModal(BCond *&bc)
{
  bc = nbcModal;
  return numNeumanModal;
}

int GeoSource::getIDis(BCond *&bc)
{
  bc = iDis;
  return numIDis;
}

int GeoSource::getIDisModal(BCond *&bc)
{
  bc = iDisModal;
  return numIDisModal;
}

int GeoSource::getIDis6(BCond *&bc)
{
  bc = iDis6;
  return numIDis6;
}

int GeoSource::getIVel(BCond *&bc)
{
  bc = iVel;
  return numIVel;
}

int GeoSource::getIVelModal(BCond *&bc)
{
  bc = iVelModal;
  return numIVelModal;
}

int GeoSource::getModalDamping(BCond *&damping)
{
  damping = modalDamping;
  return numDampedModes;
}

int GeoSource::getSurfaceDirichletBC(BCond *&bc)
{
  bc = surface_dbc;
  return numSurfaceDirichlet;
}

int GeoSource::getSurfaceNeumanBC(BCond *&bc)
{
  bc = surface_nbc;
  return numSurfaceNeuman;
}

int GeoSource::getSurfacePressure(PressureBCond *&bc)
{
  bc = surface_pres;
  return numSurfacePressure;
}

int GeoSource::getSurfaceConstraint(BCond *&bc)
{
  bc = surface_cfe;
  return numSurfaceConstraint;
}

void GeoSource::computeGlobalNumElements()
{
  // PJSA: note timing file will print nElem
  // need do a global max across all clusters
#ifdef DISTRIBUTED
  for(int i=0; i<numClusters; ++i)
	allNumClusElems[i] = structCom->globalMax(allNumClusElems[i]);  // find a better way!!
#endif
  nElemAllClusters = 0;
  for(int i=0; i<numClusters; ++i) nElemAllClusters += allNumClusElems[i];
  delete [] allNumClusElems; allNumClusElems = 0;
}

void GeoSource::setTextBC()
{
  // keep track of bc's defined in text file for addition to binary data
  numTextDirichlet = numDirichlet;
  if (numTextDirichlet)
	textDBC = dbc;
  numTextNeuman = numNeuman;
  if (numTextNeuman)
	textNBC = nbc;
}

void GeoSource::applyAuxData(int *cl2LocElem, int *cl2LocNode,
							 int minElemNum, int maxElemNum)
{
  // renumber elements in this subdomain
  for (int iElem = 0; iElem < nElem; iElem++)
	elemSet[iElem]->renum(cl2LocNode);

  // Set up element attributes

  // set material properties for element
  std::map<int, StructProp *> usedMat;
  for (int iAttr = minElemNum; iAttr < maxElemNum; iAttr++)  {

	int locElemNum = cl2LocElem[ attrib[iAttr].nele ];

	if (locElemNum >= 0)  {
	  if (attrib[iAttr].attr >= 0)  {
		StructProp *prop;
		std::map<int, StructProp *>::iterator it = usedMat.find(attrib[iAttr].attr);
		if(it == usedMat.end())  {
		  // prop = new StructProp(sProps[attrib[iAttr].attr]);
		  prop = &sProps[attrib[iAttr].attr]; // PJSA: for helmsweep, helmCoef is updated in sProps[0] so don't copy!!
											  // check with Thuan why he decided to make a copy of the material here
		  usedMat[attrib[iAttr].attr] = prop;
		}
		else
		  prop = it->second;

		elemSet[locElemNum]->setProp(prop);
	  }
	  else  {
		phantomFlag = 1;
		elemSet[locElemNum]->setProp(0);
	  }
	}

	if(attrib[iAttr].cmp_attr >= 0) {
	  if(coefData[attrib[iAttr].cmp_attr] != 0) {
		if(attrib[iAttr].cmp_frm > -1) { // user input cframe
		  elemSet[locElemNum]->setCompositeData(1, 0, 0, coefData[attrib[iAttr].cmp_attr]->values(),
												cframes[attrib[iAttr].cmp_frm]);
		}
		else { // user input ctheta
		  elemSet[locElemNum]->setCompositeData2(1, 0, 0, coefData[attrib[iAttr].cmp_attr]->values(),
												 nodes, attrib[iAttr].cmp_theta);
		}
	  }
	  else  {
		if (layInfo[attrib[iAttr].cmp_attr] == 0) {
		  fprintf(stderr," *** WARNING: Attribute found that refers to"
						 " nonexistant composite data: %d"
						 " for elemNum %d \n",
						   attrib[iAttr].cmp_attr+1, iAttr);
		  continue;
		}
	  }

	  // type is 3 for LAYC 2 for LAYN
	  int type = 3-layInfo[attrib[iAttr].cmp_attr]->getType();
	  if(attrib[iAttr].cmp_frm > -1) { // user input cframe
		elemSet[locElemNum]->setCompositeData(type, layInfo[attrib[iAttr].cmp_attr]->nLayers(), layInfo[attrib[iAttr].cmp_attr]->values(), 0,
											  cframes[attrib[iAttr].cmp_frm]);
	  }
	  else { // user input ctheta
		elemSet[locElemNum]->setCompositeData2(type, layInfo[attrib[iAttr].cmp_attr]->nLayers(), layInfo[attrib[iAttr].cmp_attr]->values(), 0,
											   nodes, attrib[iAttr].cmp_theta);
	  }
	}
  }

  // set up element frames
  for(int iFrame = 0; iFrame < numEframes; iFrame++)  {

	// get local element number
	if(efd[iFrame].elnum >= maxElemNum)
	  continue;

	int locElemNum =  cl2LocElem[ efd[iFrame].elnum ];

	// check if eframe is in this subdomain
	if(locElemNum >= 0)
	  elemSet[locElemNum]->setFrame(&(efd[iFrame].frame));

  }
}

int GeoSource::getSubCtrl(BCond *glCtrlData, int numGlData,
			   BCond *&locCtrlData, int glSub, int *&locToGlUserData)
{
  int i, iSub;

  // Form NodeToSub connectivity
  Connectivity *nodeToSub = subToNode->alloc_reverse();

  // count number of data in this subdomain
  int numLocCtrlData = 0;
  for (i = 0; i < numGlData; i++)  {
	for (iSub = 0; iSub < nodeToSub->num(glCtrlData[i].nnum); iSub++)
	  if ( (*nodeToSub)[glCtrlData[i].nnum][iSub] == glSub )  {
		numLocCtrlData++;
		break;
	  }
  }

  if (numLocCtrlData)  {

	// allocate arrays to number of data in this subdomain
	locCtrlData = new BCond[numLocCtrlData];
	locToGlUserData = new int[numLocCtrlData];

	// populate arrays
	int count = 0;
	for (i = 0; i < numGlData; i++)  {
	  for (iSub = 0; iSub < nodeToSub->num(glCtrlData[i].nnum); iSub++)
		if ( (*nodeToSub)[glCtrlData[i].nnum][iSub] == glSub )  {

		  locCtrlData[count] = glCtrlData[i];
		  locToGlUserData[count] = i;

		  // increment count
		  count++;
		}
	}
  }

  return numLocCtrlData;
}

//---------------------------------------------------------------------

int GeoSource::getSubCtrl(BCond *glCtrlData, int numGlData,
						  BCond *&locCtrlData, int *gl2ClMap, int *cl2LocMap,
						  int *&locToGlUserData)
{
  int i;

  // count number of data in this subdomain
  int numLocCtrlData = 0;
  for (i = 0; i < numGlData; i++)  {
	if (glCtrlData[i].nnum <= maxGlobNode)  {
	  int clusterNum = gl2ClMap[glCtrlData[i].nnum];
	  if (clusterNum < maxClusNode && clusterNum >= 0)
		if (cl2LocMap[clusterNum] >= 0)
		  numLocCtrlData++;
	}
  }

  if (numLocCtrlData)  {

	// allocate arrays to number of data in this subdomain
	locCtrlData = new BCond[numLocCtrlData];
	locToGlUserData = new int[numLocCtrlData];

	// populate arrays
	int count = 0;
	for (i = 0; i < numGlData; i++)  {

	  if (glCtrlData[i].nnum <= maxGlobNode)  {
		int clusterNodeNum = gl2ClMap[glCtrlData[i].nnum];
		if (clusterNodeNum < maxClusNode && clusterNodeNum >= 0)  {
		  if (cl2LocMap[clusterNodeNum] >= 0)  {

			locCtrlData[count] = glCtrlData[i];
			locToGlUserData[count] = i;

			// renumber to local node number
			locCtrlData[count].nnum = cl2LocMap[clusterNodeNum];

			// increment count
			count++;
		  }
		}
	  }
	}
  }
  return numLocCtrlData;
}

//--------------------------------------------------------------

void GeoSource::augmentBC(int numTextBC, BCond *textBC, BCond *&binBC,
						  int &numBinBC)  {

  BCond *tmpBC = 0;
  if (numBinBC)
	tmpBC = binBC;

  binBC = new BCond[numTextBC + numBinBC];

  int iBC;
  for (iBC = 0; iBC < numBinBC; iBC++)
	binBC[iBC] = tmpBC[iBC];

  for (iBC = 0; iBC < numTextBC; iBC++)
	binBC[numBinBC+iBC] = textBC[iBC];

  numBinBC += numTextBC;
}


//-------------------------------------------------------------------

int GeoSource::getBC(BCond *bc, int numBC, int *cl2LocNodeMap,
					 BCond *&subBC, int *gl2ClNodeMap){

  int numLocBC = 0;
  int bcNodeNum;

  // check if mapping exists for bc nodes
  int iBC;
  for (iBC = 0; iBC < numBC; iBC++)  {

	if (gl2ClNodeMap)  {
	  if (bc[iBC].nnum <= maxGlobNode)
		bcNodeNum = gl2ClNodeMap[bc[iBC].nnum];
	  else
		bcNodeNum = -1;
	}
	else
	  bcNodeNum = bc[iBC].nnum;

	if (bcNodeNum >= 0  && bcNodeNum < maxClusNode)
	  if (cl2LocNodeMap[ bcNodeNum ] >= 0)
		numLocBC++;
  }

  // set BC
  subBC = new BCond[numLocBC];

  int count = 0;
  for (iBC = 0; iBC < numBC; iBC++)  {

	if (gl2ClNodeMap)  {
	  if (bc[iBC].nnum <= maxGlobNode)
		bcNodeNum = gl2ClNodeMap[bc[iBC].nnum];
	  else
		bcNodeNum = -1;
	}
	else
	  bcNodeNum = bc[iBC].nnum;

	if (bcNodeNum >= 0 && bcNodeNum < maxClusNode)  {

	  if (cl2LocNodeMap[bcNodeNum] >= 0)  {

		  subBC[count] = bc[iBC];

		  // renumber to local node number
		  subBC[count].nnum = cl2LocNodeMap[bcNodeNum];

		  // increment count
		  count++;
	  }
	}
  }

  return numLocBC;
}

//--------------------------------------------------------------

void GeoSource::cleanAuxData()  {

  //layInfo.deleteArray();
  layInfo.restartArray();
  numLayInfo = 0;

  //coefData.deleteArray();
  coefData.restartArray();
  numCoefData = 0;

  numProps = 0;
  //sProps.deleteArray();
  //sProps.restartArray();

  na = 0;
  namax = 0;
  //attrib.deleteArray();
  //attrib.restartArray();

  numEframes = 0;

  cframes.restartArray();
  numCframes = 0;

  if (dbc)  {
	delete [] dbc;
	dbc = 0;
  }

  if (nbc)  {
	delete [] nbc;
	nbc = 0;
  }

  if (iDis)  {
	delete [] iDis;
	iDis = 0;
  }

  if (iVel) {
	delete [] iVel;
	iVel = 0;
  }

  if (iVelModal) {
	delete [] iVelModal;
	iVelModal = 0;
  }

}

//--------------------------------------------------------------

int GeoSource::readRanges(BinFileHandler &file, int &numRanges,
						  int (*&ranges)[2])
{
  // read number of ranges
  file.read(&numRanges, 1);

  // read number of singles
  int numSingles;
  file.read(&numSingles, 1);

  // allocate memory for ranges
  ranges = new int[numRanges+numSingles][2];

  // read in ranges
  file.read(reinterpret_cast<int *>(ranges), 2*numRanges);

  // read the singles for this subdomain
  int *singles = new int[numSingles];
  file.read(singles, numSingles);

  // compute number of values in subdomain
  int numValues = 0;
  for (int iRange = 0; iRange < numRanges; iRange++)
	numValues += ranges[iRange][1] - ranges[iRange][0] + 1;

  numValues += numSingles;

  // convert singles to range format
  for (int iSingle = 0; iSingle < numSingles; iSingle++)  {

	ranges[numRanges+iSingle][0] = singles[iSingle];
	ranges[numRanges+iSingle][1] = singles[iSingle];
  }
  numRanges += numSingles;

  delete [] singles;
  return numValues;
}

//------------------------------------------------------------

void GeoSource::outputNodeScalars(int fileNum, double *data,
				  int outputSize, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time != -1.0) {
	if(outputSize == 1)
	  fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
	else
	  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
  }

  if (oinfo[fileNum].groupNumber > 0)  {

	int group = oinfo[fileNum].groupNumber;
	std::set<int>::iterator it = nodeGroup[group].begin();

	while (it != nodeGroup[group].end() )  {
	  int inode = (domain->outFlag == 1) ? domain->nodeTable[*it]-1 : *it;
	  filePrint(oinfo[fileNum].filptr, " %d % *.*E % *.*E % *.*E % *.*E\n", *it+1, w,p,
				nodes[*it]->x, w,p, nodes[*it]->y, w,p, nodes[*it]->z, w,p, data[inode]);
	  it++;
	}

  } else {
	if(outputSize == 1) {
		fprintf(oinfo[fileNum].filptr," % *.*E\n", w, p, data[0]);
	} else {
	  for (int i = 0; i < outputSize; i++) {
		filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i]);
	  }
	}
  }
  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputNodeScalars(int fileNum, DComplex *data, int outputSize, double time) {

  int i;
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;


  switch(oinfo[fileNum].complexouttype) {
	default:
	case OutputInfo::realimag :
	  // print real part, or in the case group option both the real and imaginary parts
	  // RT: 12/10/2016 - Why only real part??
	  // RT: 12/10/2016 - Why do fprintf and filePrint alternate???
	  if(time != -1.0) {
		if(outputSize == 1)
		  fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  if (oinfo[fileNum].groupNumber > 0)  {

		int group = oinfo[fileNum].groupNumber;
		std::set<int>::iterator it = nodeGroup[group].begin();

		while (it != nodeGroup[group].end() )  {
		  int inode = (domain->outFlag == 1) ? domain->nodeTable[*it]-1 : *it;
		  filePrint(oinfo[fileNum].filptr, " %d % *.*E % *.*E % *.*E % *.*E % *.*E\n", *it+1, w, p,
					nodes[*it]->x, w,p, nodes[*it]->y, w,p, nodes[*it]->z,
					w,p, data[inode].real(), w, p, data[inode].imag());
		  it++;
		}
	  }
	  else  {
		for(i = 0; i < outputSize; i++) {
		  if(outputSize == 1) {
			fprintf(oinfo[fileNum].filptr," % *.*E", w, p, data[i].real());
		  }
		  else
			filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i].real());
		}
		// print imaginary part
		if(time != -1.0) {
		  if(outputSize != 1)
			 filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
		}
		for (i = 0; i < outputSize; i++) {
		// RT: 12/10/2016 - changed filePrint to fprintf
		  if (outputSize == 1)
			fprintf(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i].imag());
		  else
			filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, data[i].imag());
		}
	  }
	  break;
	case OutputInfo::modulusphase :
	  // print modulus
	  if(time != -1.0) {
		if(outputSize == 1)
		  fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  if (oinfo[fileNum].groupNumber > 0)  {

		int group = oinfo[fileNum].groupNumber;
		std::set<int>::iterator it = nodeGroup[group].begin();

		while (it != nodeGroup[group].end() )  {
		  int inode = (domain->outFlag == 1) ? domain->nodeTable[*it]-1 : *it;
		  filePrint(oinfo[fileNum].filptr, " %d % *.*E % *.*E % *.*E % *.*E % *.*E\n", *it+1, w, p,
					nodes[*it]->x, w,p, nodes[*it]->y, w,p, nodes[*it]->z,
					w, p, std::abs(data[inode]), w, p, arg(data[inode]));
		  it++;
		}
	  }
	  else  {
		for(i = 0; i < outputSize; i++) {
		  if (outputSize == 1)
			fprintf(oinfo[fileNum].filptr," % *.*E", w, p, std::abs(data[i]));
		  else
			filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, std::abs(data[i]));
		}
		// print phase part
		if (time != -1.0) {
		  if(outputSize != 1) filePrint(oinfo[fileNum].filptr,"  % *.*E  \n", w, p, time);
		}
		for (i = 0; i < outputSize; i++)
		  filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, arg(data[i]));
	  }
	  break;
	case OutputInfo::animate :
	  if(outputSize != 1) {
		double phi = 0;
		double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
		for(i=0; i<oinfo[fileNum].ncomplexout; ++i) {
		  filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,phi);
		  for(int j = 0; j < outputSize; j++) {
			double proj = std::abs(data[j])*cos(arg(data[j])-phi);
			filePrint(oinfo[fileNum].filptr," % *.*E\n", w, p, proj);
		  }
		  phi += incr;
		}
	  }
	  else std::cerr << " *** WARNING: animate not supported for single-node or nodal group output \n";
	  break;
  }

  fflush(oinfo[fileNum].filptr);
}


//------------------------------------------------------------

void GeoSource::outputElemVectors(int fileNum, double *data,
								  int outputSize, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time != -1.0) {
	if(outputSize == 1)
	  fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
	else
	  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
  }

  for(int i = 0; i < outputSize; i++) {
	if(outputSize == 1)
	  fprintf(oinfo[fileNum].filptr," % *.*E % *.*E\n",
			  w, p, data[2*i], w, p, data[2*i+1]);
	else
	  filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
				w, p, data[2*i], w, p, data[2*i+1]);
  }

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputElemVectors(int fileNum, DComplex *data,
								  int outputSize, double time)
{
  int i;
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  switch(oinfo[fileNum].complexouttype) {
	default:
	case OutputInfo::realimag :
	  // output real part or both real & imag parts in the case of single node output
	  if(time != -1.0) {
		if(outputSize == 1)
		  fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  for(i = 0; i < outputSize; i++) {
		if(outputSize == 1)
		  fprintf(oinfo[fileNum].filptr," % *.*E % *.*E  % *.*E % *.*E\n",
				  w, p, data[2*i].real(), w, p, data[2*i+1].real(),
				  w, p, data[2*i].imag(), w, p, data[2*i+1].imag());
		else
		  filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
					w, p, data[2*i].real(), w, p, data[2*i+1].real());
	  }
	  // output imaginary part
	  if(outputSize != 1) {
		if(time != -1.0) {
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
		}
		for(i = 0; i < outputSize; i++)
		  filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
					w, p, data[2*i].imag(), w, p, data[2*i+1].imag());
	  }
	  break;
   case OutputInfo::modulusphase :
	  // output modulus part or both modulus & phase in the case of single node output
	  if(time != -1.0) {
		if(outputSize == 1)
		  fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  for(i = 0; i < outputSize; i++) {
		if(outputSize == 1)
		  fprintf(oinfo[fileNum].filptr," % *.*E % *.*E  % *.*E % *.*E\n",
				  w, p, std::abs(data[2*i]), w, p, std::abs(data[2*i+1]),
				  w, p, arg(data[2*i]), w, p, arg(data[2*i+1]));
		else
		  filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
					w, p, std::abs(data[2*i]), w, p, std::abs(data[2*i+1]));
	  }
	  // output phase part
	  if(outputSize != 1) {
		if(time != -1.0) {
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
		}
		for(i = 0; i < outputSize; i++)
		  filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
					w, p, arg(data[2*i]), w, p, arg(data[2*i+1]));
	  }
	  break;
	case OutputInfo::animate :
	  if(outputSize != 1) {
		double phi = 0;
		double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
		for(i=0; i<oinfo[fileNum].ncomplexout; ++i) {
		  filePrint(oinfo[fileNum].filptr,"  % *.*E  \n",w,p,phi);
		  for(int j = 0; j < outputSize; j++) {
			double proj[2];
			for(int k=0; k<2; ++k)
			  proj[k] = std::abs(data[2*j+k])*cos(arg(data[2*j+k])-phi);
			filePrint(oinfo[fileNum].filptr," % *.*E % *.*E\n",
					  w, p, proj[0], w, p, proj[1]);
		  }
		  phi += incr;
		}
	  }
	  else std::cerr << " *** WARNING: animate not supported for single-node output \n";
	  break;
  }

  fflush(oinfo[fileNum].filptr);
}

//-----------------------------------------------------------------

// now same function exists not reading from file
void GeoSource::getTextDecomp(bool sowering)
{
	MatrixTimers &mt = domain->getTimers();
	startTimerMemory(mt.readDecomp, mt.memoryDecomp);

	// Get Decomposition File pointer
	FILE *f = cinfo->decPtr;

	if(f == 0)
		f = fopen("DECOMPOSITION","r");

	if(f == 0) {
		filePrint(stderr," **************************************************\n");
		filePrint(stderr," *** ERROR: DECOMPOSITION file does not exist   ***\n");
		filePrint(stderr," ***        Please provide a DECOMPOSITION file ***\n");
		filePrint(stderr," **************************************************\n");
		exit(-1);
	}

	// get decomposition
	if(subToElem) {
		fprintf(stderr," ... Already read Decomposition\n");
		return;
	}
	mt.memorySubToElem -= memoryUsed();
	// Allocate memory for subdomain to element connectivity
	int *connect = new int[nElem];

	// Get the number of subdomains in Decomposition file
	int numSub;
	int error = fscanf(f,"%d",&numSub);

	// Decomposition file error checking
	if(error == 0) {
		char s1[14],s2[40],s3[4],s4[40];
		int error = fscanf(f,"%s%s%s%s",s1,s2,s3,s4);

		// Get the number of subdomains in Decomposition file
		error = fscanf(f,"%d",&numSub);
	}

	int *cx = new int[numSub+1];

	int curEle = 0;
	int isub;
	for (isub = 0; isub < numSub; ++isub) {
		int nele;
		int n = fscanf(f,"%d",&nele);
		cx[isub] = curEle;
		if(curEle + nele > numElem()) {
			fprintf(stderr," *** ERROR: This decomposition contains more elements "
			               "than the original mesh:\n");
			fprintf(stderr," *** %d vs %d\n", curEle + nele, nElem);
			exit(1);
		}

		int iele;
		for (iele = 0; iele < nele; ++iele) {
			int n = fscanf(f,"%d",connect+curEle);
			connect[curEle] -= 01;
			curEle++;
		}
	}

	cx[numSub] = curEle;
	// PHIL FIX THIS
	/*int iele = 0;
	int cEle = 0;
	for(isub = 0; isub < numSub; ++isub) {
	  for(; iele < cx[isub+1]; iele++)
		if(glToPckElems.find(connect[iele]) != glToPckElems.end())
		  connect[cEle++] =  connect[iele];
	  cx[isub+1] = cEle;
	  }*/

	if(domain->solInfo().isNonLin() && domain->GetnContactSurfacePairs() && !domain->tdenforceFlag()) {
		optDec = new Decomposition();
		optDec->nsub = numSub;
		optDec->pele = new int[numSub+1];
		optDec->eln  = new int[curEle];
		for(int i=0; i<=numSub; i++) optDec->pele[i] = cx[i];
		for(int i=0; i<curEle; i++) optDec->eln[i] = connect[i];
	}

	subToElem = new Connectivity(numSub,cx,connect);

	subToElem->renumberTargets(glToPckElems);

	mt.memorySubToElem += memoryUsed();

	stopTimerMemory(mt.readDecomp, mt.memoryDecomp);
#ifdef DISTRIBUTED
	if(!binaryInput) {
		if(conName) {
			BinFileHandler connectivityFile(conName, "rb");
			clusToSub = std::make_unique<Connectivity>(connectivityFile, true);
			Connectivity *subToClusConnect = clusToSub->alloc_reverse();
			subToClus.resize(numSub);
			for(int i = 0; i < numSub; i++)
				subToClus[i] = (*subToClusConnect)[i][0];
			delete subToClusConnect;
			numClusters = clusToSub->csize();
		}
		else {
			subToClus.resize(numSub);
			for(int i = 0; i < numSub; i++)
				subToClus[i] = 0;
			numClusters = 1;
			int *ptr = new int[2];
			int *target = new int[numSub];
			for(int i = 0; i < numSub; i++)
				target[i] = i;
			ptr[0] = 0;
			ptr[1] = numSub;
			clusToSub = std::make_unique<Connectivity>(1,ptr,target);
			numClusNodes = nGlobNodes;
			numClusElems = nElem; //HB: not sure this is be always correct (i.e. phantoms els ...)
		}
	}
#endif
}

void GeoSource::setNumNodalOutput()
{
  // for single node output
  numNodalOutput = 0;
  for(int iInfo = 0; iInfo < numOutInfo; iInfo++)
	if(oinfo[iInfo].nodeNumber != -1)
	  numNodalOutput++;
  if(numNodalOutput)  {
	outputNodes.resize(numNodalOutput);
	outNodeIndex.resize(numNodalOutput);
	int count = 0;
	for(int iInfo = 0; iInfo < numOutInfo; iInfo++)
	  if(oinfo[iInfo].nodeNumber != -1)  {
		outputNodes[count] = oinfo[iInfo].nodeNumber;
		outNodeIndex[count] = iInfo;
		count++;
	  }
  }
}

//----------------------------------------------------------------------

int GeoSource::setDirichlet(int _numDirichlet, BCond *_dbc)
{
  if(domain->solInfo().dmpc) {
	domain->addDirichletLMPCs(_numDirichlet, _dbc);
	return 0;
  }

  if(dbc) {

	// Allocate memory for correct number of dbc
	BCond *nd = new BCond[numDirichlet+_numDirichlet];

	// copy old dbc
	int i;
	for(i = 0; i < numDirichlet; ++i)
	   nd[i] = dbc[i];

	// copy new dbc
	for(i = 0; i<_numDirichlet; ++i)
	  nd[i+numDirichlet] = _dbc[i];

	// set correct number of dbc
	numDirichlet += _numDirichlet;

	// delete old array of dbc
	delete [] dbc;

	// set new pointer to correct number of dbc
	dbc = nd;

  }

  else {
	numDirichlet = _numDirichlet;
	dbc          = _dbc;
  }

  return 0;
}


//-------------------------------------------------------------------
void GeoSource::convertHEVDirToHelmDir()
{
  if(dbcFluid) {

	// Allocate memory for correct number of dbc
	BCond *nd = new BCond[numDirichlet+numDirichletFluid];

	// copy old dbc
	int i;
	for(i = 0; i < numDirichlet; ++i)
	   nd[i] = dbc[i];

	// copy new dbc and modify the dofs
	for(i = 0; i<numDirichletFluid; ++i) {
	  nd[i+numDirichlet] = dbcFluid[i];
	  nd[i+numDirichlet].dofnum = 8-1;
	}

	// set correct number of dbc
	numDirichlet += numDirichletFluid;

	// delete old array of dbc
	if (dbc) delete [] dbc;
	delete[] dbcFluid;

	// set new pointer to correct number of dbc
	dbc = nd;
	numDirichletFluid = 0;
  }
}

//-------------------------------------------------------------------
int GeoSource::setDirichletFluid(int _numDirichletFluid, BCond *_dbcFluid)
{
  if(dbcFluid) {

	// Allocate memory for correct number of dbcFluid
	BCond *ndFluid = new BCond[numDirichletFluid+_numDirichletFluid];

	// copy old dbcFluid
	int i;
	for(i = 0; i < numDirichletFluid; ++i)
	   ndFluid[i] = dbcFluid[i];

	// copy new dbcFluid
	for(i = 0; i<_numDirichletFluid; ++i)
	  ndFluid[i+numDirichletFluid] = _dbcFluid[i];

	// set correct number of dbcFluid
	numDirichletFluid += _numDirichletFluid;

	// delete old array of dbcFluid
	delete [] dbcFluid;

	// set new pointer to correct number of dbcFluid
	dbcFluid = ndFluid;

  }

  else {
	numDirichletFluid = _numDirichletFluid;
	dbcFluid          = _dbcFluid;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setNeuman(int _numNeuman, BCond *_nbc)
{
  if (nbc) {

	// Allocate memory for correct number of nbc
	BCond *nd = new BCond[numNeuman+_numNeuman];

	// copy old nbc
	int i;
	for (i = 0; i < numNeuman; ++i)
	  nd[i] = nbc[i];

	// copy new nbc
	for (i = 0; i < _numNeuman; ++i)
	  nd[i+numNeuman] = _nbc[i];

	// set correct number of nbc
	numNeuman += _numNeuman;

	// delete old array of nbc
	delete [] nbc;

	// set new pointer to correct number of nbc
	nbc = nd;

  }
  else  {
	numNeuman = _numNeuman;
	nbc       = _nbc;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setNeumanModal(int _numNeumanModal, BCond *_nbcModal)
{
  if (nbcModal) {
	// Allocate memory for correct number of nbcModal
	BCond *nd = new BCond[numNeumanModal+_numNeumanModal];

	// copy old nbcModal
	for (int i = 0; i < numNeumanModal; ++i)
	  nd[i] = nbcModal[i];

	// copy new nbcModal
	for (int i = 0; i < _numNeumanModal; ++i)
	  nd[i+numNeumanModal] = _nbcModal[i];

	// set correct number of nbcModal
	numNeumanModal += _numNeumanModal;

	// delete old array of nbcModal
	delete [] nbcModal;

	// set new pointer to correct number of nbcModal
	nbcModal = nd;
  }
  else {
	numNeumanModal = _numNeumanModal;
	nbcModal    = _nbcModal;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIDis6(int _numIDis6, BCond *_iDis6)
{
/*
  if (iDis6) {

	// Allocate memory for correct number of iDis6
	BCond *nd = new BCond[numIDis6+_numIDis6];

	// copy old iDis6
	int i;
	for (i = 0; i < numIDis6; ++i)
	  nd[i] = iDis6[i];

	// copy new iDis6
	for ( i = 0; i < _numIDis6; ++i)
	  nd[i+numIDis6] = _iDis6[i];

	// set correct number of iDis6
	numIDis6 += _numIDis6;

	// delete old array of iDis6
	delete [] iDis6;

	// set new pointer to correct number of iDis6
	iDis6 = nd;

  }
  else  {
*/

	numIDis6 = _numIDis6;
	iDis6    = _iDis6;
//  }
  return 0;
}

//-------------------------------------------------------------------
int GeoSource::setPitaIDis6(int n, BCond *i, int numTSPitaIDis6_)
{
  // numTSPitaIDis6 corresponds to the number of time-slices and also
  // to the number of displacement vectors read from the input file.
  // numPitaIDis6 corresponds to the number of boundary conditions for
  // each vector, that is 6 * #dof.
  // These values will be used to check that the all necessary seed
  // values have been specified
  numPitaIDis6   = n / numTSPitaIDis6_;
  numTSPitaIDis6 = numTSPitaIDis6_;
  PitaIDis6 = i;
  return 0;
}

//-------------------------------------------------------------------
int GeoSource::setPitaIVel6(int n, BCond *i, int numTSPitaIVel6_)
{
  // Same principle as setPitaIDis6()
  numPitaIVel6   = n / numTSPitaIVel6_;
  numTSPitaIVel6 = numTSPitaIVel6_;
  PitaIVel6 = i;
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIDis(int _numIDis, BCond *_iDis)
{
  if (iDis) {
	// Allocate memory for correct number of iDis
	BCond *nd = new BCond[numIDis+_numIDis];

	// copy old iDis
	for (int i = 0; i < numIDis; ++i)
	  nd[i] = iDis[i];

	// copy new iDis
	for (int i = 0; i < _numIDis; ++i)
	  nd[i+numIDis] = _iDis[i];

	// set correct number of iDis
	numIDis += _numIDis;

	// delete old array of iDis
	delete [] iDis;

	// set new pointer to correct number of iDis
	iDis = nd;
  }
  else {
	numIDis = _numIDis;
	iDis    = _iDis;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIDisModal(int _numIDisModal, BCond *_iDisModal)
{
  if (iDisModal) {
	// Allocate memory for correct number of iDisModal
	BCond *nd = new BCond[numIDisModal+_numIDisModal];

	// copy old iDisModal
	for (int i = 0; i < numIDisModal; ++i)
	  nd[i] = iDisModal[i];

	// copy new iDisModal
	for (int i = 0; i < _numIDisModal; ++i)
	  nd[i+numIDisModal] = _iDisModal[i];

	// set correct number of iDisModal
	numIDisModal += _numIDisModal;

	// delete old array of iDisModal
	delete [] iDisModal;

	// set new pointer to correct number of iDisModal
	iDisModal = nd;
  }
  else {
	numIDisModal = _numIDisModal;
	iDisModal    = _iDisModal;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIVel(int _numIVel, BCond *_iVel)
{
  if (iVel) {

	// Allocate memory for correct number of iVel
	BCond *nd = new BCond[numIVel+_numIVel];

	// copy old iVel
	int i;
	for (i = 0; i < numIVel; ++i)
	  nd[i] = iVel[i];

	// copy new iVel
	for ( i = 0; i < _numIVel; ++i)
	  nd[i+numIVel] = _iVel[i];

	// set correct number of iVel
	numIVel += _numIVel;

	// delete old array of iVel
	delete [] iVel;

	// set new pointer to correct number of iVel
	iVel = nd;

  }
  else  {
	numIVel = _numIVel;
	iVel    = _iVel;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::setIVelModal(int _numIVelModal, BCond *_iVelModal)
{
  if (iVelModal) {
	// Allocate memory for correct number of iVelModal
	BCond *nd = new BCond[numIVelModal+_numIVelModal];

	// copy old iVelModal
	int i;
	for (i = 0; i < numIVelModal; ++i)
	  nd[i] = iVelModal[i];

	// copy new iVelModal
	for ( i = 0; i < _numIVelModal; ++i)
	  nd[i+numIVelModal] = _iVelModal[i];

	// set correct number of iVelModal
	numIVelModal += _numIVelModal;

	// delete old array of iVelModal
	delete [] iVelModal;

	// set new pointer to correct number of iVelModal
	iVelModal = nd;

  }
  else  {
	numIVelModal = _numIVelModal;
	iVelModal    = _iVelModal;
  }
  return 0;
}

//-------------------------------------------------------------------
/*
int GeoSource::setITemp(int _numITemp, BCond *_iTemp)
{
  if (iTemp) {

	// Allocate memory for correct number of iTemp
	BCond *nd = new BCond[numITemp+_numITemp];

	// copy old iTemp
	int i;
	for (i = 0; i < numITemp; ++i)
	  nd[i] = iTemp[i];

	// copy new iTemp
	for ( i = 0; i < _numITemp; ++i)
	  nd[i+numITemp] = _iTemp[i];

	// set correct number of iTemp
	numITemp += _numITemp;

	// delete old array of iTemp
	delete [] iTemp;

	// set new pointer to correct number of iTemp
	iTemp = nd;

  }
  else  {
	numITemp = _numITemp;
	iTemp = _iTemp;
  }
  return 0;
}
*/
//----------------------------------------------------------------

int GeoSource::setModalDamping(int _numDampedModes, BCond *_modalDamping)
{
  if(modalDamping) {

	// Allocate memory for correct number of damped modes
	BCond *nd = new BCond[numDampedModes + _numDampedModes];

	// copy old modalDamping
	int i;
	for(i = 0; i < numDampedModes; ++i)
	  nd[i] = modalDamping[i];

	for(i = 0; i < _numDampedModes; ++i)
	  nd[i+numDampedModes] = _modalDamping[i];

	numDampedModes += _numDampedModes;

	delete[] modalDamping;
	modalDamping = nd;
  }
  else {

	numDampedModes = _numDampedModes;
	modalDamping = _modalDamping;
  }
  return 0;
}

//-------------------------------------------------------------------

int GeoSource::addSurfaceDirichlet(int _numSurfaceDirichlet, BCond *_surface_dbc)
{
  if(surface_dbc) {

	// Allocate memory for correct number of dbc
	BCond *nd = new BCond[numSurfaceDirichlet+_numSurfaceDirichlet];

	// copy old dbc
	int i;
	for(i = 0; i < numSurfaceDirichlet; ++i)
	   nd[i] = surface_dbc[i];

	// copy new dbc
	for(i = 0; i<_numSurfaceDirichlet; ++i)
	  nd[i+numSurfaceDirichlet] = _surface_dbc[i];

	// set correct number of dbc
	numSurfaceDirichlet += _numSurfaceDirichlet;

	// delete old array of dbc
	delete [] surface_dbc;

	// set new pointer to correct number of dbc
	surface_dbc = nd;

  }

  else {
	numSurfaceDirichlet = _numSurfaceDirichlet;
	surface_dbc          = _surface_dbc;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::addSurfaceNeuman(int _numSurfaceNeuman, BCond *_surface_nbc)
{
  if(surface_nbc) {

	// Allocate memory for correct number of nbc
	BCond *nd = new BCond[numSurfaceNeuman+_numSurfaceNeuman];

	// copy old nbc
	int i;
	for(i = 0; i < numSurfaceNeuman; ++i)
	   nd[i] = surface_nbc[i];

	// copy new nbc
	for(i = 0; i<_numSurfaceNeuman; ++i)
	  nd[i+numSurfaceNeuman] = _surface_nbc[i];

	// set correct number of nbc
	numSurfaceNeuman += _numSurfaceNeuman;

	// delete old array of nbc
	delete [] surface_nbc;

	// set new pointer to correct number of nbc
	surface_nbc = nd;

  }

  else {
	numSurfaceNeuman = _numSurfaceNeuman;
	surface_nbc          = _surface_nbc;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::addSurfacePressure(int _numSurfacePressure, PressureBCond *_surface_pres)
{
  if(surface_pres) {

	// Allocate memory for correct number of pres
	PressureBCond *nd = new PressureBCond[numSurfacePressure+_numSurfacePressure];

	// copy old pres
	int i;
	for(i = 0; i < numSurfacePressure; ++i)
	   nd[i] = surface_pres[i];

	// copy new pres
	for(i = 0; i<_numSurfacePressure; ++i)
	  nd[i+numSurfacePressure] = _surface_pres[i];

	// set correct number of pres
	numSurfacePressure += _numSurfacePressure;

	// delete old array of pres
	delete [] surface_pres;

	// set new pointer to correct number of pres
	surface_pres = nd;

  }

  else {
	numSurfacePressure = _numSurfacePressure;
	surface_pres       = _surface_pres;
  }

  return 0;
}

//-------------------------------------------------------------------

int GeoSource::addSurfaceConstraint(int _numSurfaceConstraint, BCond *_surface_cfe)
{
  if(surface_cfe) {

	// Allocate memory for correct number of cfe
	BCond *nd = new BCond[numSurfaceConstraint+_numSurfaceConstraint];

	// copy old cfe
	int i;
	for(i = 0; i < numSurfaceConstraint; ++i)
	   nd[i] = surface_cfe[i];

	// copy new cfe
	for(i = 0; i<_numSurfaceConstraint; ++i)
	  nd[i+numSurfaceConstraint] = _surface_cfe[i];

	// set correct number of cfe
	numSurfaceConstraint += _numSurfaceConstraint;

	// delete old array of cfe
	delete [] surface_cfe;

	// set new pointer to correct number of cfe
	surface_cfe = nd;

  }

  else {
	numSurfaceConstraint = _numSurfaceConstraint;
	surface_cfe          = _surface_cfe;
  }

  return 0;
}

//-------------------------------------------------------------------

void GeoSource::readMatchInfo(BinFileHandler &matchFile,
		int (*matchRanges)[2], int numMatchRanges, int subNum,
		int *clusToLocElem, int myID)
{
  // get TOC
  BinFileHandler::OffType tocLoc;
  matchFile.read(&tocLoc, 1);
  matchFile.seek(tocLoc);
  int numClusMatches;
  matchFile.read(&numClusMatches, 1);
  BinFileHandler::OffType dataStart;
  matchFile.read(&dataStart, 1);

  // get matches according to ranges
  for (int iR = 0; iR < numMatchRanges; iR++)  {

	// get number of data in range
	int nMatch = matchRanges[iR][1] - matchRanges[iR][0] + 1;

	// seek to correct file position
	matchFile.seek(dataStart + matchRanges[iR][0] * ( 2*sizeof(int) + 5*sizeof(double) ));

	for (int iMatch = 0; iMatch < nMatch; iMatch++)  {

	  int tmpData[2];  // cluster elem #, global match file position
	  matchFile.read(tmpData, 2);

	  //matchData[subNum][iMatch].elemNum = clusToLocElem[tmpData[0]];
	  //this modification is necessary because of the sort of subToElem in DecDomain
	  int glSubNum = (*cpuToSub)[myID][subNum];
	  int glElemNum = (*unsortedSubToElem)[glSubNum][clusToLocElem[tmpData[0]]];
	  matchData[subNum][iMatch].elemNum = subToElem->cOffset(glSubNum,glElemNum);
	  matchData[subNum][iMatch].glPos = tmpData[1];

	  // read natual coordinates of match
	  double coords[2];
	  matchFile.read(coords, 2);

	  matchData[subNum][iMatch].xi = coords[0];
	  matchData[subNum][iMatch].eta = coords[1];

	  // read gap vector
	  matchFile.read(gapVec[subNum][iMatch], 3);
	}
  }
}

//-------------------------------------------------------------------

int GeoSource::getCPUMap(FILE *f, Connectivity *subToSub, int glNumSub, int numCPU)
{
	int numSub = glNumSub;
	if(subToSub)
		numSub = subToSub->csize();
	int totSub = numSub;

	if(!cpuToSub) {
		if(f == 0) { // Trivial map

#ifdef USE_SCOTCH
			if(subToSub) {
	if(verboseFlag) filePrint(stderr, " ... Making CPU Map using SCOTCH, numCPU = %d ...\n", numCPU);
	Connectivity *graph = subToSub->modifyAlt(); // scotch doesn't allow loops
	cpuToSub = std::shared_ptr<Connectivity>(graph->SCOTCH_graphPart(numCPU));
	delete graph;
	  } else
#endif
			{
				if(verboseFlag) filePrint(stderr, " ... Making Trivial CPU Map, numCPU = %d ... \n", numCPU);
				int *cx  = new int[numCPU+1];
				int *connect = new int[totSub];
				if(subToCPU) delete [] subToCPU;
				subToCPU = new int[totSub];

				int curSub = 0;
				int icpu;
				int nSubPerCPU = totSub/numCPU;
				int remain = totSub%numCPU;
				for (icpu = 0; icpu < numCPU; ++icpu) {
					cx[icpu] = curSub;
					int iSub;
					int subForCPU = (icpu < remain) ? nSubPerCPU+1 : nSubPerCPU;
					for(iSub = 0; iSub < subForCPU; ++iSub) {
						connect[curSub] = curSub;
						subToCPU[curSub] = icpu;
						curSub++;
					}
				}
				cx[numCPU] = curSub;
				cpuToSub = std::make_shared<Connectivity>(numCPU, cx, connect);
			}
		} else {
			cpuToSub = std::make_shared<Connectivity>(f,totSub);
			numCPU = cpuToSub->csize();
			filePrint(stderr, " ... Reading CPU Map from file %s, numCPU = %d ... \n", mapName, numCPU);
			if(numCPU != structCom->numCPUs()) {
				fprintf(stderr, " *** ERROR: CPUMAP file is for %d MPI processes\n", numCPU);
				exit(-1);
			}
		}
	}

	Connectivity subDomainToCPU = cpuToSub->reverse();
	if(cpuToCPU) delete cpuToCPU;
	cpuToCPU = cpuToSub->transcon(&subDomainToCPU);

	if(domain->solInfo().aeroFlag >= 0 || domain->solInfo().aeroheatFlag >= 0) {
		int numLocSub = 0;
#ifdef USE_MPI
		int myID = structCom->myID();
#else
		int myID = 0;
#endif
		numLocSub = cpuToSub->num(myID);

		// allocate for match data arrays
		typedef double (*gVec)[3];
		numMatchData = new int[numLocSub];
		matchData = new MatchData *[numLocSub];
		numGapVecs = new int[numLocSub];
		gapVec = new gVec[numLocSub];

		// PJSA 02-04-2010
		if(matchName != NULL) {
			BinFileHandler connectivityFile(conName, "rb");
			Connectivity clusToSub(connectivityFile, true);
			Connectivity subToClus = clusToSub.reverse();

			// build global to cluster subdomain map
			int *gl2ClSubMap = new int[totSub];
			for (int iClus = 0; iClus < clusToSub.csize(); iClus++)  {
				int clusNum = 0;
				for (int iSub = 0; iSub < clusToSub.num(iClus); iSub++)
					gl2ClSubMap[ clusToSub[iClus][iSub] ] = clusNum++;
			}

			for(int locSub = 0; locSub < numLocSub; ++locSub) {

				int glSub = (*cpuToSub)[myID][locSub];
				int clusNum = subToClus[glSub][0];
				char fullDecName[128];
				const char *suffix = computeClusterSuffix(clusNum + 1, clusToSub.csize());
				sprintf(fullDecName, "%s%s", decName, suffix);
				BinFileHandler decFile(fullDecName, "rb");

				// read some stuff from decFile which can now be discarded
				int numClusSub;
				decFile.read(&numClusSub, 1);
				BinFileHandler::OffType curLoc = decFile.tell();
				int clusSub = gl2ClSubMap[glSub];
				decFile.seek(curLoc + sizeof(BinFileHandler::OffType) * clusSub);
				BinFileHandler::OffType infoLoc;
				decFile.read(&infoLoc, 1);
				decFile.seek(infoLoc);
				int (*nodeRanges)[2] = 0;
				int numNodeRanges;
				int numLocNodes = readRanges(decFile, numNodeRanges, nodeRanges);
				int (*elemRanges)[2] = 0;
				int numElemRanges;
				int numLocElems = readRanges(decFile, numElemRanges, elemRanges);

				int minElemNum = elemRanges[0][0];
				int maxElemNum = 0;
				for(int iR = 0; iR < numElemRanges; iR++)  {
					if(elemRanges[iR][0] < minElemNum)
						minElemNum = elemRanges[iR][0];
					if(elemRanges[iR][1] > maxElemNum)
						maxElemNum = elemRanges[iR][1];
				}
				maxElemNum++;  // for easy allocation
				int iElem;
				int *cl2LocElem = new int[maxElemNum];
				for(iElem = 0; iElem < maxElemNum; iElem++)
					cl2LocElem[iElem] = -1;
				iElem = 0;
				for(int iR = 0; iR < numElemRanges; ++iR)
					for(int cElem = elemRanges[iR][0]; cElem <= elemRanges[iR][1]; ++cElem) {
						cl2LocElem[cElem] = iElem++;
					}
				if(nodeRanges) delete [] nodeRanges;
				if(elemRanges) delete [] elemRanges;

				int nConnects;
				decFile.read(&nConnects, 1);
				if(nConnects > 0) {
					int *connectedDomain = new int[nConnects];
					decFile.read(connectedDomain, nConnects);
					int size, numtarget;
					decFile.read(&size, 1);
					decFile.read(&numtarget, 1);
					int *pointer = new int[size+1];
					decFile.read(pointer, size+1);
					if(numtarget > 0) {
						int *target = new int[numtarget];
						decFile.read(target, numtarget);
						delete [] target;
					}
					delete [] pointer;
					delete [] connectedDomain;
				}

				// now we are at the right place in decFile to read matcher stuff
				int numMatchRanges;
				int (*matchRanges)[2] = 0;
				numMatchData[locSub] = readRanges(decFile, numMatchRanges, matchRanges);
				matchData[locSub] = new MatchData[numMatchData[locSub]];
				gapVec[locSub] = new double[numMatchData[locSub]][3];
				if(numMatchRanges) {
					char fullMatchName[128];
					sprintf(fullMatchName, "%s%s", matchName, suffix);
					BinFileHandler matchFile(fullMatchName, "rb");
					readMatchInfo(matchFile, matchRanges, numMatchRanges, locSub, cl2LocElem, myID);
				}
				if(matchRanges) delete [] matchRanges;
				delete [] cl2LocElem;
				delete [] suffix;
			}
			delete [] gl2ClSubMap;
			delete unsortedSubToElem;
			unsortedSubToElem = 0;
		}
/*
	else {
	  fprintf(stderr,"*** ERROR: Binary Match File not specified\n");
	  exit (-1);
	}
*/
	}

	return numCPU;
}

//-----------------------------------------------------------------------

void GeoSource::deleteMatchArrays(int numLocSub) {

  // de-allocation for match data arrays
  if(numMatchData) {
	delete [] numMatchData;
	numMatchData = 0;
  }
  if(matchData) {
	if(matchName != NULL) for(int i=0; i<numLocSub; ++i) delete [] matchData[i];
	delete [] matchData;
	matchData = 0;
  }
  if(numGapVecs) {
	delete [] numGapVecs;
	numGapVecs = 0;
  }
  if(gapVec) {
	if(matchName != NULL) for(int i=0; i<numLocSub; ++i) delete [] gapVec[i];
	delete [] gapVec;
	gapVec = 0;
  }
}

//-----------------------------------------------------------------------
void GeoSource::addOutput(OutputInfo &outputInfo) {

  oinfo[numOutInfo++] = outputInfo;
  if(outputInfo.type == OutputInfo::Farfield ||
	 outputInfo.type == OutputInfo::Kirchhoff)
	domain->solInfo().farfield = true;
}

//-----------------------------------------------------------------------
std::pair<int, int>
GeoSource::getTimeSliceOutputFileIndices(int timeSliceRank) {
  std::map<int, std::pair<int, int> >::const_iterator it = timeSliceOutputFiles.find(timeSliceRank);
  return it != timeSliceOutputFiles.end() ? it->second : std::pair<int, int>(0, 0);
}

//--------------------------------------------------------------------
void GeoSource::openOutputFilesForPita(int sliceRank)
{
  std::pair<int, int> indices = getTimeSliceOutputFileIndices(sliceRank);
  for(int iInfo = indices.first; iInfo < indices.second; ++iInfo) {
   if(!oinfo[iInfo].PodRomfile) {
   if (oinfo[iInfo].interval != 0) {
	  char *fileName = oinfo[iInfo].filename;
	  if (strlen(cinfo->outputExt) != 0) {
		int len1 = strlen(fileName);
		int len2 = strlen(cinfo->outputExt);
		char *nfn = new char[len1+len2+1];
		strcpy(nfn, fileName);
		strcat(nfn,cinfo->outputExt);
		fileName = nfn;
	  }

	  if((oinfo[iInfo].filptr= fopen(fileName,"w")) == (FILE *) 0 ) {
		fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n",oinfo[iInfo].filename );
		exit(0);
	  }
	  outputHeader(iInfo);
	  fflush(oinfo[iInfo].filptr);
	}
  }
 }
}


//--------------------------------------------------------------------
void GeoSource::openOutputFiles(int *outNodes, int *outIndex, int numOuts)
{
  int iInfo;

  if(numOuts == 0) { // open all output files and write their corresponding TOPDOM/DEC header
	for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
	 if(!oinfo[iInfo].PodRomfile) {
	  if(oinfo[iInfo].interval != 0) {
		char *fileName = oinfo[iInfo].filename;
		if (strlen(cinfo->outputExt) != 0) {
		  int len1 = strlen(fileName);
		  int len2 = strlen(cinfo->outputExt);
		  char *nfn = new char[len1+len2+1];
		  strcpy(nfn, fileName);
		  strcat(nfn,cinfo->outputExt);
		  fileName = nfn;
		}

		if((oinfo[iInfo].filptr= fopen(fileName,"w")) == (FILE *) 0 ) {
		  fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n",oinfo[iInfo].filename );
		  exit(0);
		}
		outputHeader(iInfo);
		fflush(oinfo[iInfo].filptr);
	  }
	 }
	}
  }
  else { // open selected output files
	for(int iOut = 0; iOut < numOuts; iOut++)  {
	  iInfo = outIndex[iOut];
	  if(!oinfo[iInfo].PodRomfile) {
	   if(oinfo[iInfo].interval != 0) {
		char *fileName = oinfo[iInfo].filename;
		if(strlen(cinfo->outputExt) != 0) {
		  int len1 = strlen(fileName);
		  int len2 = strlen(cinfo->outputExt);
		  char *nfn = new char[len1+len2+1];
		  strcpy(nfn, fileName);
		  strcat(nfn,cinfo->outputExt);
		  fileName = nfn;
		}

		if((oinfo[iInfo].filptr= fopen(fileName,"w")) == (FILE *) 0 ) {
		  fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n",oinfo[iInfo].filename );
		  exit(0);
		}
		outputHeader(iInfo);
		fflush(oinfo[iInfo].filptr);
	  }
	}
   }
  }
}

//--------------------------------------------------------------------
void GeoSource::openSensorOutputFiles()
{
  int iInfo;

  // open all single node output files and write their corresponding TOPDOM/DEC header
  for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
	if(!(oinfo[iInfo].nodeNumber == -1)) {
	  if(oinfo[iInfo].interval != 0) {
		char *fileName = oinfo[iInfo].filename;
		if (strlen(cinfo->outputExt) != 0) {
		  int len1 = strlen(fileName);
		  int len2 = strlen(cinfo->outputExt);
		  char *nfn = new char[len1+len2+1];
		  strcpy(nfn, fileName);
		  strcat(nfn,cinfo->outputExt);
		  fileName = nfn;
		}

		if((oinfo[iInfo].filptr = fopen(fileName,"w")) == (FILE *) 0) {
		  fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n", oinfo[iInfo].filename);
		  exit(0);
		}
		outputHeader(iInfo);
		fflush(oinfo[iInfo].filptr);
	  }
	}
  }
}

//--------------------------------------------------------------------
void GeoSource::closeOutputFiles()
{
  for(int i = 0; i < numOutInfo; ++i) {
   if(!oinfo[i].PodRomfile) {
	this->closeOutputFileImpl(i);
   }
  }
}

//--------------------------------------------------------------------
void GeoSource::closeOutputFilesForPita(int sliceRank)
{
  std::pair<int, int> indices = getTimeSliceOutputFileIndices(sliceRank);
  for (int i = indices.first; i < indices.second; ++i) {
   if(!oinfo[i].PodRomfile) {
	this->closeOutputFileImpl(i);
   }
  }
}

//--------------------------------------------------------------------
void GeoSource::closeOutputFileImpl(int fileIndex)
{
  if((oinfo[fileIndex].interval != 0) && oinfo[fileIndex].filptr) {
	fclose(oinfo[fileIndex].filptr);
	oinfo[fileIndex].filptr = 0;
  }
}

//--------------------------------------------------------------------

void GeoSource::outputHeader(int fileNumber)
{
  // only one node is requested for output,
  const int& averageFlg = oinfo[fileNumber].averageFlg;
  std::string s = (averageFlg == -1 || averageFlg == 0) ? "element" : "node";
  if(oinfo[fileNumber].nodeNumber != -1) {
	if(domain->isComplex()) {
	  switch(oinfo[fileNumber].complexouttype) {
		case(OutputInfo::realimag) :
          fprintf(oinfo[fileNumber].filptr, "# %s %d (real and imaginary parts)\n", s.c_str(), oinfo[fileNumber].nodeNumber+1);
		  break;
		case(OutputInfo::modulusphase) :
          fprintf(oinfo[fileNumber].filptr, "# %s %d (complex modulus and phase)\n", s.c_str(), oinfo[fileNumber].nodeNumber+1);
		  break;
	  }
	}
	else
      fprintf(oinfo[fileNumber].filptr, "# %s %d\n", s.c_str(), oinfo[fileNumber].nodeNumber+1);
	return;
  }

  // get data description
  char headDescrip[200];
  getHeaderDescription(headDescrip, fileNumber);
  fprintf(oinfo[fileNumber].filptr, "%s", headDescrip);
}

//---------------------------------------------------------------

void GeoSource::setControl(char *_checkfile, char *_nodeSetName,
						   char *_elemSetName, char *_bcondSetName)
{
/*
  cinfo->checkfile   = _checkfile;
  cinfo->nodeSetName = _nodeSetName;
  cinfo->elemSetName = _elemSetName;
  cinfo->bcondSetName = _bcondSetName;
*/
  // PJSA: limit these 4 strings to 31 characters + terminating null-character for xpost compatability
  cinfo->checkfile = new char [32];
  for(int i=0; i<32; ++i) {
	if(i==31) cinfo->checkfile[i] = '\0';
	else cinfo->checkfile[i] = _checkfile[i];
	if(_checkfile[i] == '\0') break;
  }
  cinfo->nodeSetName = new char [32];
  for(int i=0; i<32; ++i) {
	if(i==31) cinfo->nodeSetName[i] = '\0';
	else cinfo->nodeSetName[i] = _nodeSetName[i];
	if(_nodeSetName[i] == '\0') break;
  }
  cinfo->elemSetName = new char [32];
  for(int i=0; i<32; ++i) {
	if(i==31) cinfo->elemSetName[i] = '\0';
	else cinfo->elemSetName[i] = _elemSetName[i];
	 if(_elemSetName[i] == '\0') break;
  }
  if(_bcondSetName != 0) {
	cinfo->bcondSetName = new char [32];
	for(int i=0; i<32; ++i) {
	  if(i==31) cinfo->bcondSetName[i] = '\0';
	  else cinfo->bcondSetName[i] = _bcondSetName[i];
	   if(_bcondSetName[i] == '\0') break;
	}
  }
}

//---------------------------------------------------------------------

void GeoSource::createSingleCpuToSub(int numSub)
{
  MatrixTimers &mt = domain->getTimers();
  mt.memoryCPUMAP -= memoryUsed();
  int *ptr = new int[2];
  int *trg = new int[numSub];
  ptr[0] = 0;
  ptr[1] = numSub;
  for (int iSub = 0; iSub < numSub; ++iSub)
	trg[iSub] = iSub;
  cpuToSub = std::make_shared<Connectivity>(1, ptr, trg);
  mt.memoryCPUMAP += memoryUsed();
}

//---------------------------------------------------------------------

void GeoSource::setControlFile(char *_filename)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->fileName = _filename;
}

//--------------------------------------------------------------------

void GeoSource::setControlRoutine(char *_routinename)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->routineName = _routinename;
}

//--------------------------------------------------------------------

int GeoSource::setSensorLocations(int _numSensor, BCond *_sensor)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->numSensor = _numSensor;
  claw->sensor    = _sensor;
  return 0;
}

//--------------------------------------------------------------------

int GeoSource::setActuatorLocations(int _numActuator, BCond *_actuator)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->numActuator = _numActuator;
  claw->actuator    = _actuator;
  return 0;
}

//--------------------------------------------------------------------

int GeoSource::setUsddLocation(int _numSensor, BCond *_sensor)
{
  if (claw == 0) claw = new ControlLawInfo;

  if (claw->userDisp) {
	// Allocate memory for correct number of usdd
	BCond *nd = new BCond[claw->numUserDisp + _numSensor];

	// copy old usdd
	int i;
	for(i = 0; i < claw->numUserDisp; ++i)
	   nd[i] = claw->userDisp[i];

	// copy new usdd
	for(i = 0; i<_numSensor; ++i) {
	  nd[i+claw->numUserDisp].nnum = _sensor[i].nnum;
	  nd[i+claw->numUserDisp].dofnum = _sensor[i].dofnum;
	  nd[i+claw->numUserDisp].val = _sensor[i].val;
	}

	// set correct number of usdd
	claw->numUserDisp += _numSensor;

	// delete old array of usdd
	delete [] claw->userDisp;

	// set new pointer to correct number of usdd
	claw->userDisp = nd;
  }
  else {
	claw->numUserDisp = _numSensor;

	// need to copy the pointer data
	claw->userDisp    = new BCond[claw->numUserDisp];
	for (int k = 0; k < claw->numUserDisp; k++)  {
	  claw->userDisp[k].nnum = _sensor[k].nnum;
	  claw->userDisp[k].dofnum = _sensor[k].dofnum;
	  claw->userDisp[k].val = _sensor[k].val;
	}
  }

  return 0;
}

//--------------------------------------------------------------------

int GeoSource::setUsdfLocation(int _numActuator, BCond *_actuator)
{
  if(claw == 0) claw = new ControlLawInfo;
  claw->numUserForce = _numActuator;
  claw->userForce    = _actuator;
  return 0;
}

//--------------------------------------------------------------------

void GeoSource::outputEnergy(int fileNum, double time, double W) {

  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;
  
  fprintf(oinfo[fileNum].filptr," %e % *.*E\n", time, w, p, W);

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputEnergyPerAttribute(int fileNum, double time, double* W, int attrib)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;
  
  fprintf(oinfo[fileNum].filptr," %e ", time);
  for(int i=0; i<attrib; i++) fprintf(oinfo[fileNum].filptr, "% *.*E", w, p, W[i]);
  fprintf(oinfo[fileNum].filptr, "\n");
  
  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputEnergy(int fileNum, double time, DComplex W) {

  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  // print real part
  fprintf(oinfo[fileNum].filptr," %e % *.*E\n", time, w, p, W.real());

  // print imaginary part
  fprintf(oinfo[fileNum].filptr," %e % *.*E\n", time, w, p, W.imag());

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputEnergies(int fileNum, double time, double Wext, double Waero,
				   double Wela, double Wkin, double Wdmp, double error) {

  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  fprintf(oinfo[fileNum].filptr," %e % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
		  time, w, p, Wext, w, p, Waero, w, p, Wela, w, p, Wkin, w, p, Wdmp, w, p, error);

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputEnergies(int fileNum, double time, DComplex Wext, DComplex Waero,
							   DComplex Wela, DComplex Wkin, DComplex Wdmp, DComplex error) {

  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  // print real part
  fprintf(oinfo[fileNum].filptr," %e % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
		  time, w, p, Wext.real(), w, p, Waero.real(), w, p, Wela.real(),
				w, p, Wkin.real(), w, p, Wdmp.real(), w, p, error.real());

  // print imaginary part
  fprintf(oinfo[fileNum].filptr," %e % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
		  time, w, p, Wext.imag(), w, p, Waero.imag(), w, p, Wela.imag(),
				w, p, Wkin.imag(), w, p, Wdmp.imag(), w, p, error.imag());

  fflush(oinfo[fileNum].filptr);
}

//--------------------------------------------------------------------

void GeoSource::getHeaderDescription(char *headDescrip, int fileNumber)
{
  int dataType = 0;  // 1 for nodal, 2 for elemental

  // solver information structure
  SolverInfo &sinfo = domain->solInfo();

  char prbType[20];
  if(sinfo.probType == SolverInfo::Static) strcpy(prbType,"Static");
  else if(sinfo.probType == SolverInfo::Dynamic) strcpy(prbType,"Dynam");
  else if(sinfo.probType == SolverInfo::NonLinStatic ||
		  sinfo.probType == SolverInfo::MatNonLinStatic) strcpy(prbType,"NLStatic");
  else if(sinfo.probType == SolverInfo::NonLinDynam ||
		  sinfo.probType == SolverInfo::MatNonLinDynam ||
		  sinfo.probType == SolverInfo::PodRomOffline) strcpy(prbType,"NLDynamic");
  else if(sinfo.probType == SolverInfo::ArcLength) strcpy(prbType,"Arclength");
  else if(sinfo.probType == SolverInfo::TempDynamic) strcpy(prbType,"Temp");
  else if(sinfo.probType == SolverInfo::AxiHelm) strcpy(prbType,"AxiHelm");
  else if(sinfo.buckling) strcpy(prbType,"Buckling");
  else if(sinfo.probType == SolverInfo::Modal) strcpy(prbType,"Modal");
  if(isShifted() && sinfo.probType != SolverInfo::Modal) {
	if(sinfo.doFreqSweep) strcpy(prbType,"FrequencySweep");
	else strcpy(prbType,"FrequencyResponse");
  }
  if(sinfo.isAcousticHelm()) {
	if(sinfo.doFreqSweep) strcpy(prbType,"FAcousticSweep");
	else strcpy(prbType,"FAcoustic");
  }

  if (oinfo[fileNumber].type == OutputInfo::YModulus)
	strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::MDensity)
	strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::Thicknes)
	strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::ShapeAtt)
	strcpy(prbType,"Attributes");

  if (oinfo[fileNumber].type == OutputInfo::ShapeStc)
	strcpy(prbType,"Static");

  int type = oinfo[fileNumber].type;

  // No header for CompositeData
  if (type == OutputInfo::Composit)
	return;

  if (type == OutputInfo::InXForce || type == OutputInfo::InYForce ||
	  type == OutputInfo::InZForce || type == OutputInfo::AXMoment ||
	  type == OutputInfo::AYMoment || type == OutputInfo::AZMoment )  {
	sprintf(headDescrip, header[type], prbType, cinfo->elemSetName, nElem);
	dataType = 2;
	return;
  }

  int avgnum = oinfo[fileNumber].averageFlg;
  if (avgnum > 0)  {
	if (oinfo[fileNumber].groupNumber > 0)  {
	  char ng[80];
	  sprintf(ng, "NODALGROUP %d: NODENUMBER X0  Y0  Z0  RESULT-1 RESULT-2 ... RESULT-N", oinfo[fileNumber].groupNumber);

	  sprintf(headDescrip, header[type], prbType, ng, nodeGroup[oinfo[fileNumber].groupNumber].size());
	}
	else {
	  int numNodesOut = (domain->outFlag) ? nodes.nnz() : numNodes;
	  sprintf(headDescrip, header[type], prbType, cinfo->nodeSetName, numNodesOut);
	}
	dataType = 1;
  }
  else {
	sprintf(headDescrip, ele_header[type], prbType, cinfo->elemSetName);
	dataType = 2;
  }

  return;
}

void GeoSource::computeAndCacheHeaderLength(int fileNum)
{
  char headDescrip[200];
  getHeaderDescriptionAndLength(headDescrip, fileNum);
}

int GeoSource::getHeaderDescriptionAndLength(char *headDescrip, int fileNumber) {
	getHeaderDescription(headDescrip, fileNumber);
	const int len = std::strlen(headDescrip);

	// Cache header length
	headLen.resize(numOutInfo);
	headLen[fileNumber] = len;

	return len;
}

ControlInterface* GeoSource::getUserSuppliedFunction()
{
  ControlInterface *userSupFunc = 0;

#if !defined(WINDOWS) && !defined(SALINAS)
  if(claw) {
	void *handle;
	dlerror(); // forget about the last error
	handle = dlopen(claw->fileName, RTLD_NOW);
	const char *errorMsg;
	if ((errorMsg = dlerror() ) != 0) {
	  fprintf(stderr," *** ERROR: in dynamic loading of %s: %s\n",
			  claw->fileName,errorMsg);
	  exit(-1);
	}

  ControlInterface ** fcp =
	  (ControlInterface **) dlsym(handle, claw->routineName);

  if (fcp == 0) {
	fprintf(stderr," *** ERROR: in dynamic loading of %s: "
					 "control function not found\n",
					 claw->routineName);
	exit(-1);
  }

  userSupFunc = *fcp;
  }
#else
	std::cerr << " (W) : Warning : not supported under windows !" << std::endl;
#endif
  return userSupFunc;
}

//--------------------------------------------------------------------

void GeoSource::outputElemStress(int fileNum, double *stressData,
								 int numOutElems, const std::vector<size_t> &offsets, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  if(time != -1.0) {
	if(numOutElems == 1)
	  fprintf(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
	else
	  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
  }

  for(int i = 0; i < numOutElems; i++)  {
	int numNodes = offsets[i+1] - offsets[i];
	for(int iNode = 0; iNode < numNodes; iNode++)
	  filePrint(oinfo[fileNum].filptr," % *.*E", w, p, stressData[offsets[i]+iNode]);
	filePrint(oinfo[fileNum].filptr,"\n");
  }

  fflush(oinfo[fileNum].filptr);
}

void GeoSource::outputElemStress(int fileNum, DComplex *stressData,
								 int numOutElems, const std::vector<size_t> &offsets, double time)
{
  int w = oinfo[fileNum].width;
  int p = oinfo[fileNum].precision;

  switch(oinfo[fileNum].complexouttype) {
	default:
	case OutputInfo::realimag :
	  // print real part
	  if(time != -1.0) {
		if(numOutElems == 1)
		  filePrint(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  for(int i = 0; i < numOutElems; i++)  {
		int numNodes = offsets[i+1] - offsets[i];
		for(int iNode = 0; iNode < numNodes; iNode++)
		  filePrint(oinfo[fileNum].filptr," % *.*E", w, p, stressData[offsets[i]+iNode].real());
		filePrint(oinfo[fileNum].filptr,"\n");
	  }
	  // print imaginary part
	  if(time != -1.0) {
		if(numOutElems == 1)
		  filePrint(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  for(int i = 0; i < numOutElems; i++)  {
		int numNodes = offsets[i+1] - offsets[i];
		for(int iNode = 0; iNode < numNodes; iNode++)
		  filePrint(oinfo[fileNum].filptr," % *.*E", w, p, stressData[offsets[i]+iNode].imag());
		filePrint(oinfo[fileNum].filptr,"\n");
	  }
	  break;
	case OutputInfo::modulusphase :
	  // print modulus
	  if(time != -1.0) {
		if(numOutElems == 1)
		  filePrint(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  for(int i = 0; i < numOutElems; i++)  {
		int numNodes = offsets[i+1] - offsets[i];
		for(int iNode = 0; iNode < numNodes; iNode++)
		  filePrint(oinfo[fileNum].filptr," % *.*E", w, p, std::abs(stressData[offsets[i]+iNode]));
		filePrint(oinfo[fileNum].filptr,"\n");
	  }
	  // print phase part
	  if(time != -1.0) {
		if(numOutElems == 1)
		  filePrint(oinfo[fileNum].filptr,"  % *.*E  ", w, p, time);
		else
		  filePrint(oinfo[fileNum].filptr,"  % *.*E\n", w, p, time);
	  }
	  for(int i = 0; i < numOutElems; i++)  {
		int numNodes = offsets[i+1] - offsets[i];
		for(int iNode = 0; iNode < numNodes; iNode++)
		  filePrint(oinfo[fileNum].filptr," % *.*E", w, p, arg(stressData[offsets[i]+iNode]));
		filePrint(oinfo[fileNum].filptr,"\n");
	  }
	  break;
	case OutputInfo::animate :
	  double phi = 0;
	  double incr = 2.0*PI/double(oinfo[fileNum].ncomplexout);
	  for(int i=0; i<oinfo[fileNum].ncomplexout; ++i) {
		filePrint(oinfo[fileNum].filptr,"  % *.*E  ",w,p,phi);
		if(numOutElems != 1) filePrint(oinfo[fileNum].filptr,"\n");
		for(int j = 0; j < numOutElems; j++) {
		  int numNodes = offsets[j+1] - offsets[j];
		  for(int iNode = 0; iNode < numNodes; iNode++) {
			double proj = std::abs(stressData[offsets[j]+iNode])*cos(arg(stressData[offsets[j]+iNode])-phi);
			filePrint(oinfo[fileNum].filptr," % *.*E", w, p, proj);
		  }
		  filePrint(oinfo[fileNum].filptr,"\n");
		}
		phi += incr;
	  }
	  break;
  }

  fflush(oinfo[fileNum].filptr);
}

//----------------------------------------------------------------------

int GeoSource::glToPackElem(int i) const
{
  const std::map<int, int>::const_iterator it = glToPckElems.find(i);
  return (it != glToPckElems.end()) ? it->second : -1;
}

//-----------------------------------------------------------------------

void
GeoSource::loadMaterial(const char *matName, const char *fileName)
{
#if !defined(WINDOWS) && !defined(SALINAS)
  void *handle;
  handle = dlopen(fileName, RTLD_NOW);
  const char *errorMsg;
  if ((errorMsg = dlerror() ) != 0) {
	fprintf(stderr," *** ERROR: in dynamic loading of Material file %s: %s\n",
		fileName,errorMsg);
	exit(-1);
  }
  MatLoader fct = (MatLoader) dlsym(handle, "materialLoader");
  if(fct == 0) {
	fprintf(stderr," *** ERROR: in dynamic loading of %s: "
	"materialLoader not found\n", fileName);
	exit(-1);
  }
  fprintf(stderr, "Loaded material is at %p\n", fct);
  userDefMat[matName] = fct;

#else
		std::cerr << " (W) : Warning : not supported under windows/salinas !" << std::endl;
#endif
}

void
GeoSource::addMaterial(int i, const char *matName, DoubleList &args)
{
	auto it = userDefMat.find(matName);
	if(it == userDefMat.end()) {
		fprintf(stderr, "User defined material %s is not found\n", matName);
		exit(-1);
	}
	MatLoader ml = userDefMat[matName];
	materials[i] = (*ml)(args.nval, args.v);
}

void
GeoSource::simpleDecomposition(int numSubdomains, bool estFlag, bool weightOutFlag, bool trivialFlag, bool fsglFlag)
{
 decJustCalled=true;
 double    t1 = getTime();

 int maxEle = elemSet.last();

 if(trivialFlag) {

   int numEle = 0;
   for(int i=0; i<maxEle; ++i)
	 if(elemSet[i]) numEle++;

   optDec = new Decomposition;
   optDec->nsub = numSubdomains;
   optDec->pele = new int[optDec->nsub+1];
   optDec->eln = new int[numEle];

   int div = numEle / optDec->nsub;
   int rem = numEle % optDec->nsub;

   optDec->pele[0] = 0;
   for(int i=0; i<optDec->nsub; i++)
	 optDec->pele[i+1] = optDec->pele[i] + ((i < rem) ? div+1 : div);

   int k = 0;
   for(int i=0; i<maxEle; i++)
	 if(elemSet[i]) optDec->eln[k++] = i;

   if(randomShuffle) { 
	 filePrint(stderr, " ... Trivial Decomposition With Random Shuffle (RAND_MAX = %d) ... \n", RAND_MAX);
	 std::srand(1);
	 random_shuffle(optDec->eln, optDec->eln+k);
   }

   if(verboseFlag)
	 filePrint(stderr, " ... %d Elements Have Been Arranged in %d Subdomains ...\n",
			   maxEle, optDec->nsub);
 }
 else if(useScotch) {
#ifdef USE_SCOTCH
   if(verboseFlag) filePrint(stderr, " ... Decomposing mesh using SCOTCH  ...\n");

   Elemset packedEset;
   std::map<int, int> pckToGlElems;
   int numele = elemSet.last(), nElem = 0;
   for(int iEle = 0; iEle < numele; ++iEle) {
	 Element *ele = elemSet[iEle];
	 if(ele) {
	   packedEset.elemadd(nElem, ele);
	   pckToGlElems[nElem] = iEle;
	   nElem++;
	 }
   }

   Connectivity elemToNode(packedEset.asSet());
   Connectivity nodeToElem = elemToNode.reverse();
   Connectivity elemToElem = elemToNode.transcon(nodeToElem);
   Connectivity *graph = elemToElem.modifyAlt(); // scotch doesn't allow loops
   Connectivity *subToElem = graph->SCOTCH_graphPart(numSubdomains);
   subToElem->renumberTargets(pckToGlElems);
   int *ptr = new int[subToElem->csize()+1];
   for(int i = 0; i < subToElem->csize()+1; ++i) ptr[i] = subToElem->getPointer()[i];
   int *tgt = new int[subToElem->getNumTarget()];
   for(int i = 0; i < subToElem->getNumTarget(); ++i) tgt[i] = subToElem->getTarget()[i];
   optDec = new Decomposition(numSubdomains, ptr, tgt);
   delete subToElem; 
   delete graph; 

   if(verboseFlag)
	 filePrint(stderr, " ... %d Elements Have Been Arranged in %d Subdomains ...\n",
			   maxEle, optDec->nsub);
#else
   std::cerr << " *** WARNING: USE_SCOTCH is not defined in GeoSource::simpleDecomposition\n";
   exit(-1);
#endif
 }
 else {

 elemSet.setWeights();

 Elemset baseSet(maxEle);
 int nSpring = 0, maxSprNodes = 0, nMass = 0;
 int iEle, iSub;
 for(iEle = 0; iEle < maxEle; ++iEle)
   if(elemSet[iEle]) {
	 if(allowMechanisms || (elemSet[iEle]->isSpring() == false && elemSet[iEle]->isMass() == false))
	   baseSet.elemadd(iEle, elemSet[iEle]);
	 else {
	   if(elemSet[iEle]->isSpring()) {
	 nSpring++;
	 if(elemSet[iEle]->numNodes() > maxSprNodes)
	   maxSprNodes = elemSet[iEle]->numNodes();
	   }
	   else nMass++;
	 }
   }
 if(nSpring > 0) filePrint(stderr, " ... This mesh has %d springs ...\n", nSpring);

 MultiFront mf(&baseSet, &nodes, bool(domain->getNumFSI()), fsglFlag);

 // Decompose and optimize the structure into subdomains
 if ( domain->solInfo().isCoupled && (!isFeti(domain->solInfo().solvercntl->type)||
	  (!domain->solInfo().isMatching && domain->solInfo().solvercntl->fetiInfo.fsi_corner != 0) ) )
   optDec = mf.decompose(numSubdomains, bool(domain->getNumFSI()));
 else
   optDec = mf.decompose(numSubdomains);

 if(verboseFlag) filePrint(stderr, " ... %d Elements Have Been Arranged in %d Subdomains and %d Springs ...\n",
		   optDec->pele[optDec->nsub], optDec->nsub, nSpring);

 nSpring += nMass;
 if(nSpring || nMass) {
   int (*springAssign)[2] = new int[nSpring][2];
   int *subIncr = new int[optDec->nsub];
   for(iSub = 0; iSub < optDec->nsub; ++iSub)
	 subIncr[iSub] = 0;
   int *lNd = (int *) dbg_alloca(sizeof(int)*maxSprNodes);
   nSpring = 0;
   for(iEle = 0; iEle < maxEle; ++iEle)
	 if(elemSet[iEle] && (!allowMechanisms && (elemSet[iEle]->isSpring() || elemSet[iEle]->isMass()))) {
	   elemSet[iEle]->nodes(lNd);
	   springAssign[nSpring][0] = iEle;
	   springAssign[nSpring][1] = mf.bestSubFor(elemSet[iEle]->numNodes(), lNd);
	   if(springAssign[nSpring][1] < 0)
	 filePrint(stderr, " *** WARNING: Spring element %d is not connected to any non-spring element\n",
				   springAssign[nSpring][0]+1);
	   else
	 subIncr[springAssign[nSpring][1]]++;
	   nSpring++;
	 }
   int *nDecPtr = new int[optDec->nsub+1];
   int *nDecTarget = new int[optDec->pele[optDec->nsub] + nSpring];
   int ptr = 0;
   for(iSub = 0; iSub < optDec->nsub; ++iSub) {
	 ptr += subIncr[iSub] + optDec->num(iSub);
	 nDecPtr[iSub] = ptr;
   }
   nDecPtr[optDec->nsub] = ptr;
   for(iSub = 0; iSub < optDec->nsub; ++iSub)
	 for(iEle = 0; iEle < optDec->num(iSub); ++iEle)
	   nDecTarget[--nDecPtr[iSub]] = (*optDec)[iSub][iEle];
   for(iEle = 0; iEle < nSpring; ++iEle)
	 nDecTarget[--nDecPtr[springAssign[iEle][1]]] = springAssign[iEle][0];
   delete [] optDec->pele;
   delete [] optDec->eln;
   optDec->pele = nDecPtr;
   optDec->eln  = nDecTarget;

   // fprintf(stderr, "Going to make a check on springs, max = %d\n", maxSprNodes);
   int lastNode = nodes.size();
   bool *ndIsUsed = new bool[lastNode];
   int iNode;
   for(iNode = 0; iNode < lastNode; ++iNode)
	 ndIsUsed[iNode] = false;
   for(iSub = 0; iSub < optDec->nsub; ++iSub) {
	 for(iEle = 0; iEle < optDec->num(iSub); ++iEle) {
	   int nds[128];
	   int enm = (*optDec)[iSub][iEle];
	   if(allowMechanisms || (elemSet[enm]->isSpring() == false && elemSet[enm]->isMass() == false)) {
	 int nn = elemSet[enm]->numNodes();
	 elemSet[enm]->nodes(nds);
	 for(int i = 0; i < nn; ++i)
	   ndIsUsed[nds[i]] = true;
	   }
	 }
	 for(iEle = 0; iEle < optDec->num(iSub); ++iEle) {
	   int nds[128];
	   int enm = (*optDec)[iSub][iEle];
	   if(!allowMechanisms && (elemSet[enm]->isSpring() || elemSet[enm]->isMass())) {
	 elemSet[enm]->nodes(nds);
	 if(ndIsUsed[nds[0]] == false && ndIsUsed[nds[1]] == false)
	   filePrint(stderr, " *** WARNING: Found a badly assigned spring\n");
	   }
	 }
	 for(iEle = 0; iEle < optDec->num(iSub); ++iEle) {
	   int nds[128];
	   int enm = (*optDec)[iSub][iEle];
	   if(allowMechanisms || (elemSet[enm]->isSpring() == false && elemSet[enm]->isMass() == false)) {
	 int nn = elemSet[enm]->numNodes();
	 elemSet[enm]->nodes(nds);
	 for(int i = 0; i < nn; ++i)
	   ndIsUsed[nds[i]] = false;
	   }
	 }
   }
   // filePrint(stderr," ... Assigned Springs to Subdomains In %14.5f sec ...\n", (getTime() - t1)/1000.0);
 }

 if(weightOutFlag){
   FILE *weightFile;
   weightFile = domain->openFile(cinfo->checkfile, ".load");
   filePrint(stderr," ... Evaluating Subdomain Weight Distribution ...\n");
   t1 = getTime();
   double *subW = new double[optDec->nsub];
   double minW, maxW = 0;
   double totW = 0;
   int iij;
   for(iij = 0; iij < optDec->nsub; ++iij) {
	 double trueW = 0;
	 for(iEle = 0; iEle < optDec->num(iij); ++iEle)
	   trueW += elemSet[(*optDec)[iij][iEle]]->trueWeight();
	 subW[iij] = trueW;
	 if(iij == 0 || trueW < minW)
	   minW = trueW;
	 if(trueW > maxW)
	   maxW = trueW;
	 totW += trueW;
   }
   fprintf(weightFile, "# %d %e %e %e %f\n", optDec->nsub, minW, totW/optDec->nsub,
   maxW, maxW*optDec->nsub/totW);
   for(iij = 0; iij < optDec->nsub; ++iij)
	 fprintf(weightFile, "%e\n", subW[iij]);
   fclose(weightFile);
 }

 if(estFlag) {
   double mem[5];
   // Open and output memory file
   filePrint(stderr," ... Estimating Subdomain Memory Requirements ...\n");
   t1 = getTime();
   FILE *memFilePtr = domain->openFile(cinfo->checkfile, ".mem");
   mf.memEstimate(optDec, 3, mem, memFilePtr);
   fclose(memFilePtr);
 }
 }

 // Open optimized decomposition file
 if(verboseFlag) filePrint(stderr," ... Saving Decomposition File      ...\n");
 t1 = getTime();
 FILE *optFilePtr = domain->openFile(cinfo->checkfile, ".optDec");

 // Output optimized decomposition to file
 optDec->setName(cinfo->elemSetName);
#ifndef SALINAS
 optDec->outputDump(optFilePtr, 0);
#endif
 fclose(optFilePtr);
}

void
GeoSource::modifyDecomposition(int maxEleCopy)
{
 int maxEle = elemSet.last();

 optDec = new Decomposition;
 optDec->nsub = optDecCopy->nsub;
 optDec->pele = new int[optDec->nsub+1];
 optDec->eln = new int[maxEle];

 int div = (maxEle-maxEleCopy) / optDec->nsub;
 int rem = (maxEle-maxEleCopy) % optDec->nsub;

 optDec->pele[0] = 0;
 for(int i=0; i<optDec->nsub; i++)
   optDec->pele[i+1] = optDec->pele[i] + (optDecCopy->pele[i+1]-optDecCopy->pele[i]) + ((i < rem) ? div+1 : div);

 int k=0, l=maxEleCopy;
 for(int i=0; i<optDec->nsub; i++) {
   for(int j=optDecCopy->pele[i]; j<optDecCopy->pele[i+1]; ++j) {
	 optDec->eln[k++] = optDecCopy->eln[j];
   }
   for(int j=optDec->pele[i]+(optDecCopy->pele[i+1]-optDecCopy->pele[i]); j<optDec->pele[i+1]; ++j) {
	 optDec->eln[k++] = l++;
   }
 }

 if(verboseFlag)
   filePrint(stderr, " ... %d Elements Have Been Arranged in %d Subdomains ...\n",
			 maxEle, optDec->nsub);

 if(subToElem) { delete subToElem; subToElem = 0; }
 clusToSub.reset();
 if(unsortedSubToElem) { delete unsortedSubToElem; unsortedSubToElem = 0; }
}

void
GeoSource::setExitAfterDec(bool exit)
{
  exitAfterDec=exit;
  if(exit) domain->solInfo().setProbType(SolverInfo::Decomp);
}

#include <Driver.d/Sower.h>
#include <vector>

Connectivity
bindSommerToSolid(const Elemset &elemSet, gsl::span<SommerElement*> sommerElems, bool wetFlag = false) {
	std::vector<std::pair<int,int>> sToE;
	Connectivity nToE = Connectivity{ elemSet.asSet() }.reverse();
	Connectivity sToN = Connectivity::fromElements( sommerElems.size(),
		[&sommerElems] (auto i){ return sommerElems[i] ? sommerElems[i]->getNodeSpan() : gsl::span<const int>{}; });
	Connectivity sToAnyE = sToN.transcon(nToE);
	std::vector<int> anyENodes;
	auto nBaseSommer = wetFlag ? sommerElems.size()/2 : sommerElems.size();
	for (size_t i = 0; i < nBaseSommer; ++i) {
		auto sommerNodes = sommerElems[i]->getNodeSpan();
		bool found = false;
		bool foundWet = false;
		for (auto anyI : sToAnyE[i]) {
			auto ae = elemSet[anyI];
			anyENodes.resize(ae->numNodes());
			ae->nodes(anyENodes.data());
			auto inElement = [ &anyENodes ](auto nd) {
				return std::find(anyENodes.begin(), anyENodes.end(), nd) != anyENodes.end();
			};
			if ( std::all_of(sommerNodes.begin(), sommerNodes.end(), inElement) ) {
				if(wetFlag){
					bool isFluid = ae->isFluidElement();
					auto target = isFluid ? i : i + nBaseSommer;
					if (ae->isFluidElement() && !found) {
						sToE.push_back({target, anyI});
						sommerElems[target]->sFlag = !isFluid;
						sommerElems[target]->soundSpeed = ae->getProperty()->soundSpeed;
					}
					(isFluid ? found : foundWet) = true;
					if(found && foundWet)
						break;
				} else {
					sToE.push_back({i, anyI});
					found = true;
					break;
				}
			}
		}
		if (!found || (wetFlag && !foundWet)) {
			fprintf(stderr,
			        "Connectivity.C: No adjacent element to the Sommerfeld element.\n");
			exit(-1);
		}
	}
	return Connectivity::fromLinkRange(sToE);
}

/**
 * writeDistributedInputFiles will do the first step of the distributed input for fem :
 * it will assign subdomains to nCluster clusters and create a data file for each cluster
 * containing only the data pertinent to that cluster. (all of this is done by creating and initializing a Sower object properly )
 * @param nCluster number of cluster fem will run on
 * @see class Sower Driver.d/Sower.h
*/
void GeoSource::writeDistributedInputFiles(int nCluster, Domain *domain, int nCpu)
{
  if(subToElem == 0)
	getTextDecomp(true); //the option true or false is for sowering or not
						 // if sowering there is no reordering

  FILE *f = fopen(mapName,"r");
  if(!f) {
    getCPUMap(f, (Connectivity *)NULL, subToElem->csize(), nCpu);
    f = fopen(mapName, "w");
    cpuToSub->write(f);
    fclose(f);
  }
  else {
    filePrint(stderr, " ... Reading CPU Map from file %s ... \n", mapName);
    cpuToSub = std::make_shared<Connectivity>(f, subToElem->csize());
    if(cpuToSub->csize() != nCpu) {
      filePrint(stderr, "ERROR: specified CPU Map file %s has incorrect number of CPUs\n", mapName);
      exit(-1);
    }
  }

  // TOPOLOGY (also element data PRESSURE and PRELOAD)
  Sower sower(subToElem, domain->getElementSet(), nCluster, domain->viewSurfEntities(), cpuToSub.get()); //HB
  /* the sower object creates the cluster to element connectivity and from it will
	 create cluster to data connectivities for each data added subsequently.
	 Thus It will finally know which cluster needs what */

  // begin adding data to sower object

  // NODES
  typedef ImplicitConnectivity<Elemset*, EsetGeomAccessor> implicitElem;
  implicitElem *eToN = new implicitElem(&(domain->getElementSet()));
  sower.addParentToChildData<NodesIO, CoordSet*, implicitElem*>(NODES_TYPE, ELEMENTS_TYPE, 0, &nodes, eToN);
  //Connectivity *eToN = new Connectivity(&(domain->getElementSet()));
  //sower.addParentToChildData<NodesIO, CoordSet*, Connectivity*>(NODES_TYPE, ELEMENTS_TYPE, 0, &nodes, eToN);
  delete eToN;

  // ATTRIBUTES
  std::pair<int, std::map<int, Attrib>* > attrPair = std::make_pair(na, &attrib);
  ImplicitConnectivity<std::pair<int,std::map<int, Attrib>* >*, ElemAttrAccessor>
	*EtoAtt = new ImplicitConnectivity<std::pair<int, std::map<int, Attrib>* >*, ElemAttrAccessor>(&attrPair);
  sower.addParentToChildData<AttribIO>(ATTRIBUTES_TYPE, ELEMENTS_TYPE, 0, &attrPair, EtoAtt);

  // NEUMAN BOUNDARY CONDITIONS
  typedef std::pair<int, BCond*> BCPair;
  typedef ImplicitConnectivity<BCPair*, BCDataAccessor> implicitBC;
  BCPair forcesPair = std::make_pair(numNeuman, nbc);
  implicitBC *forcesToNodes = new implicitBC(&forcesPair);
  sower.addChildToParentData<BCDataIO<FORCES_TYPE> >(FORCES_TYPE, NODES_TYPE, 0, &forcesPair, forcesToNodes);
  delete forcesToNodes;

  // DIRICHLET BOUNDARY CONDITIONS
  BCPair dispPair = std::make_pair(numDirichlet, dbc);
  implicitBC *dispToNodes = new implicitBC(&dispPair);
  sower.addChildToParentData<BCDataIO<DISPLACEMENTS_TYPE> >(DISPLACEMENTS_TYPE, NODES_TYPE, 0, &dispPair, dispToNodes);
  delete dispToNodes;

  // MATERIALS
  typedef std::pair<int, SPropContainer* > MATPair;
  typedef ImplicitConnectivity<std::pair<int,std::map<int,Attrib>* >*, MatAttrAccessor> implicitMat;
  MATPair matPair = std::make_pair(numProps, &sProps);
  implicitMat* EleToMat = new implicitMat(&attrPair);
  sower.addParentToChildData<MatIO>(MATERIALS_TYPE, ELEMENTS_TYPE, 0, &matPair, EleToMat);
  delete EleToMat;

  delete EtoAtt; // not used any further

  // TETT
  typedef ImplicitConnectivity<std::pair<int,SPropContainer* >*, CurveMatAccessor> implicitCurve;
  implicitCurve* MatToCurve = new implicitCurve(&matPair);
  sower.addParentToChildData<MFTTDataIO<TETT_TYPE> >(TETT_TYPE, MATERIALS_TYPE, 0, domain->getCTETT(), MatToCurve);
  delete MatToCurve;

  // YMTT
  typedef ImplicitConnectivity<std::pair<int,SPropContainer* >*, CurveYoungMatAccessor> implicitYoungCurve;
  implicitYoungCurve* MatYoungToCurve = new implicitYoungCurve(&matPair);
  sower.addParentToChildData<MFTTDataIO<YMTT_TYPE> >(YMTT_TYPE, MATERIALS_TYPE, 0, domain->getYMTT(), MatYoungToCurve);
  delete MatYoungToCurve;

  // LMPC
  typedef ImplicitConnectivity<std::pair<int,ResizeArray<LMPCons *>* >*, LMPCAccessor> implicitLMPC;
  std::pair<int,ResizeArray<LMPCons *>* > LMPCPair = std::make_pair(domain->getNumLMPC(),domain->getLMPC());
  implicitLMPC* ImplicitLMPC = new implicitLMPC(&LMPCPair);
  sower.addChildToParentData<LMPCIO>(LMPC_TYPE, NODES_TYPE, 0, &LMPCPair, ImplicitLMPC);
  delete ImplicitLMPC;

  // COMPOSITE
  typedef std::pair<int,std::map<int,Attrib>* > AttPair;
  typedef ImplicitConnectivity<AttPair*, CmpAttrAccessor> implicitCMP;
  implicitCMP* EleToCmp =  new implicitCMP(&attrPair);
  typedef std::pair<int, ResizeArray<LayInfo *>* > layInfoPair;
  layInfoPair  lip = std::make_pair(numLayInfo, &layInfo);
  sower.addParentToChildData<CompositeLIO>(COMPOSITEL_TYPE, ELEMENTS_TYPE, 0, &lip, EleToCmp);
  typedef std::pair<int, ResizeArray<CoefData *>* > coefDataPair;
  coefDataPair cdp = std::make_pair(numCoefData, &coefData);
  sower.addParentToChildData<CompositeCIO>(COMPOSITEC_TYPE, ELEMENTS_TYPE, 0, &cdp, EleToCmp);

  // CFRAMES
  typedef ImplicitConnectivity<AttPair*, CmpFrAttrAccessor> implicitCFM;
  implicitCFM* EleToCFM =  new implicitCFM(&attrPair);
  typedef std::pair<int, ResizeArray<double *>* > cFramesPair;
  cFramesPair cfp = std::make_pair(numCframes, &cframes);
  sower.addParentToChildData<CFramesIO>(CFRAMES_TYPE, ELEMENTS_TYPE, 0, &cfp, EleToCFM);

  // BOFFSET
  typedef ImplicitConnectivity<std::vector<OffsetData>*, BoffsetAccessor> implicitBoffset;
  implicitBoffset* offToElem = new implicitBoffset(&offsets);
  sower.addChildToParentData<BoffsetIO>(BOFFSET_TYPE, ELEMENTS_TYPE, 0, &offsets, offToElem);
  delete offToElem;

  // EFRAMES
  typedef std::pair<int,ResizeArray<EFrameData>* > EfPair;
  typedef ImplicitConnectivity<EfPair*, EFrameDataAccessor > implicitEFrame;
  EfPair efPair = std::make_pair(numEframes,&efd);
  implicitEFrame* efToElem = new implicitEFrame(&efPair);
  sower.addChildToParentData<EFrameIO>(EFRAME_TYPE, ELEMENTS_TYPE, 0, &efPair, efToElem);

  // DIMASS
  std::vector<DMassData*> dmv; // what can I do with a chained list !
  DMassData* fdm = domain->getFirstDMassData();
  if(fdm != 0) // some DMASS
	{
	  do
	{
	  dmv.push_back(fdm);
	  fdm = fdm->next;
	}
	  while(fdm->next!= 0);
	  typedef ImplicitConnectivity<std::vector<DMassData* >*, DimassAccessor > implicitDMass;
	  implicitDMass * massToNode = new implicitDMass(&dmv);
	  sower.addChildToParentData<DMassIO>(DIMASS_TYPE, NODES_TYPE, 0, &dmv, massToNode);
	}

  // INITIAL DISPLACEMENTS
  BCPair idispPair = std::make_pair(numIDis, iDis);
  implicitBC *idispToNodes = new implicitBC(&idispPair);
  sower.addChildToParentData<BCDataIO<IDISP_TYPE> >(IDISP_TYPE, NODES_TYPE, 0, &idispPair, idispToNodes);
  delete idispToNodes;

  // INITIAL DISPLACEMENTS (6 column)
  BCPair idisp6Pair = std::make_pair(numIDis6, iDis6);
  implicitBC *idisp6ToNodes = new implicitBC(&idisp6Pair);
  sower.addChildToParentData<BCDataIO<IDISP6_TYPE> >(IDISP6_TYPE, NODES_TYPE, 0, &idisp6Pair, idisp6ToNodes);
  delete idisp6ToNodes;

  // INITIAL VELOCITIES
  BCPair ivelPair = std::make_pair(numIVel, iVel);
  implicitBC *ivelToNodes = new implicitBC(&ivelPair);
  sower.addChildToParentData<BCDataIO<IVEL_TYPE> >(IVEL_TYPE, NODES_TYPE, 0, &ivelPair, ivelToNodes);
  delete ivelToNodes;
/*
  // INITIAL TEMPERATURES
  BCPair itempPair = std::make_pair(numITemp, iTemp);
  implicitBC *itempToNodes = new implicitBC(&itempPair);
  sower.addChildToParentData<BCDataIO<ITEMP_TYPE> >(ITEMP_TYPE, NODES_TYPE, 0, &itempPair, itempToNodes);
  delete itempToNodes;
*/
  // COMPLEX DIRICHLET
  typedef std::pair<int, ComplexBCond*> ComplexBCPair;
  typedef ImplicitConnectivity<ComplexBCPair*, ComplexBCDataAccessor> implicitComplexBC;
  ComplexBCPair cdirPair = std::make_pair(numComplexDirichlet, cdbc);
  implicitComplexBC cdirToNodes(&cdirPair);
  sower.addChildToParentData<ComplexBCDataIO<HDIR_TYPE> >(HDIR_TYPE, NODES_TYPE, 0, &cdirPair, &cdirToNodes);

  // COMPLEX NEUMANN
  ComplexBCPair cneuPair = std::make_pair(numComplexNeuman, cnbc);
  implicitComplexBC cneuToNodes(&cneuPair);
  sower.addChildToParentData<ComplexBCDataIO<HNEU_TYPE> >(HNEU_TYPE, NODES_TYPE, 0, &cneuPair, &cneuToNodes);

  // ATDDNB & HDNB
  typedef std::pair<int, SommerElement **> SommerPair;
  SommerPair dnbPair = std::make_pair(domain->numNeum,(domain->neum).yield());
  Connectivity dnbToE = bindSommerToSolid( domain->getElementSet(), {dnbPair.second, dnbPair.first});
  sower.addChildToParentData<SommerDataIO<DNB_TYPE> >(DNB_TYPE, ELEMENTS_TYPE, 0, &dnbPair, &dnbToE);

  // ATDROB & HSCB
  SommerPair scatPair = std::make_pair(domain->numScatter,(domain->scatter).yield());
  Connectivity scatToE = bindSommerToSolid(domain->getElementSet(), {scatPair.second, scatPair.first});
  sower.addChildToParentData<SommerDataIO<SCAT_TYPE> >(SCAT_TYPE, ELEMENTS_TYPE, 0, &scatPair, &scatToE);

  // ATDARB & HARB
  SommerPair arbPair = std::make_pair(domain->numSommer,(domain->sommer).yield());
  Connectivity arbToE = bindSommerToSolid(domain->getElementSet(), {arbPair.second, arbPair.first});
  sower.addChildToParentData<SommerDataIO<ARB_TYPE> >(ARB_TYPE, ELEMENTS_TYPE, 0, &arbPair, &arbToE);

  // HWIBO
  // RT: Duplicate the wet elements so they can be sent to both sides
  ResizeArray<SommerElement *> wet2(0);
  for(int i=0;i<domain->numWet;i++) wet2[i] = domain->wet[i];
  for(int i=0;i<domain->numWet;i++) wet2[i+domain->numWet] = domain->wet[i]->clone();
 
  SommerPair wetPair = std::make_pair(2*domain->numWet,wet2.yield());
  Connectivity wetToE = bindSommerToSolid(domain->getElementSet(), {wetPair.second,2*domain->numWet},true);
  sower.addChildToParentData<SommerDataIO<WET_TYPE> >(WET_TYPE, ELEMENTS_TYPE, 0, &wetPair, &wetToE);

  if (claw) {
	// distribute control law data -> have to do that for all 4 categories...
	for (int i = 0 ; i < claw->numSensor ; i++) {
	  claw->sensor[i].val=i;
	}
	BCPair sensorPair = std::make_pair(claw->numSensor, claw->sensor);
	implicitBC sensorToNodes(&sensorPair);
	sower.addChildToParentData<BCDataIO<SENSOR_TYPE> >(SENSOR_TYPE, NODES_TYPE, 0, &sensorPair, &sensorToNodes);

	for (int i = 0 ; i < claw->numActuator ; i++) {
	  claw->actuator[i].val=i;
	}
	BCPair actuatorPair = std::make_pair(claw->numActuator, claw->actuator);
	implicitBC actuatorToNodes(&actuatorPair);
	sower.addChildToParentData<BCDataIO<ACTUATOR_TYPE> >(ACTUATOR_TYPE, NODES_TYPE, 0, &actuatorPair, &actuatorToNodes);

	for (int i = 0 ; i < claw->numUserDisp ; i++) {
	  claw->userDisp[i].val=i;
	}
	BCPair usddPair = std::make_pair(claw->numUserDisp, claw->userDisp);
	implicitBC usddToNodes(&usddPair);
	sower.addChildToParentData<BCDataIO<USDD_TYPE> >(USDD_TYPE, NODES_TYPE, 0, &usddPair, &usddToNodes);

	for (int i = 0 ; i < claw->numUserForce ; i++) {
	  claw->userForce[i].val=i;
	}
	BCPair usdfPair = std::make_pair(claw->numUserForce, claw->userForce);
	implicitBC usdfToNodes(&usdfPair);
	sower.addChildToParentData<BCDataIO<USDF_TYPE> >(USDF_TYPE, NODES_TYPE, 0, &usdfPair, &usdfToNodes);
  }

  // done adding data to sower object
#ifdef SOWER_DEBUG
  std::cerr << "  Debug requested for sower\n" << std::endl;
  sower.printDebug();
#endif
  sower.write();
}

std::unique_ptr<Connectivity>
GeoSource::getDecomposition()
{
	if(binaryInput)
		getBinaryDecomp();
	else if(!optDec)
		getTextDecomp();
	else {
		int i;
		int numSub = optDec->nsub;
		// Allocate memory for offsets
		int *cx = new int[numSub+1];
		// Allocate memory for subdomain to element connectivity
		nElem=0;
		for(i=0; i<numSub; i++) nElem += optDec->num(i);
		int *connect = new int[nElem];
		for(i=0; i<=numSub; i++) cx[i] = optDec->pele[i];
		for(i=0; i<nElem; i++) connect[i] = optDec->eln[i];
		subToElem = new Connectivity(numSub,cx,connect);
		subToElem->renumberTargets(glToPckElems);  // required if gaps in element numbering

#ifdef DISTRIBUTED
		if(conName) {
			BinFileHandler connectivityFile(conName, "rb");
			clusToSub = std::make_unique<Connectivity>(connectivityFile, true);
			Connectivity subToClusConnect = clusToSub->reverse();
			subToClus.resize(numSub);
			for(int i = 0; i < numSub; i++)
				subToClus[i] = subToClusConnect[i][0];
			numClusters = clusToSub->csize();
		}
		else {
			subToClus.resize(numSub);
			for(int i = 0; i < numSub; i++)
				subToClus[i] = 0;
			numClusters = 1;
			std::vector<int> target(numSub); std::iota(target.begin(), target.end(), 0);
			std::vector<size_t> ptr(2);
			ptr[0] = 0;
			ptr[1] = numSub;
			clusToSub = std::make_unique<Connectivity>(1,std::move(ptr),std::move(target));
			numClusNodes = nGlobNodes;
			numClusElems = nElem; //HB: not sure this is be always correct (i.e. phantoms els ...)
		}
#endif
	}
	if(matchName != nullptr && !unsortedSubToElem) unsortedSubToElem = new Connectivity(*subToElem);
	return std::make_unique<Connectivity>( *subToElem );
}

void GeoSource::getBinaryDecomp()
{
	if(!subToElem) {
		int myID = structCom->myID();
		std::ostringstream oss;

		FILE *f = fopen(mapName,"r"); // JAT 080515
		if(f == 0) {
			filePrint(stderr, "*** ERROR: Cannot open CPU Map file %s\n", mapName);
			exit(-1);
		}
		cpuToSub = std::make_shared<Connectivity>(f,0);
		fclose(f);
		int numCPU = cpuToSub->csize();
		if(numCPU != structCom->numCPUs()) {
			fprintf(stderr, " *** ERROR: CPU Map file %s is for %d CPUs\n", mapName, numCPU);
			exit(-1);
		}

		int firstSubInCPU = (*cpuToSub)[myID][0]; // JAT 080315
		for(int i = 0; i < cpuToSub->num(myID); ++i) {
			if(subToClus[(*cpuToSub)[myID][i]] != subToClus[firstSubInCPU]) {
				fprintf(stderr, " *** ERROR: Subdomains mapped to CPU %d are not in the same cluster\n", myID);
				exit(-1);
			}
		}
		oss << decomposition_ << subToClus[firstSubInCPU]+1;
		BinFileHandler fp(oss.str().c_str(), "rb");

		ConnectivityT<size_t,int> csubToSub2(fp);
		auto *csubToNode = new ConnectivityT<size_t,gl_node_idx>(fp);
		ConnectivityT<size_t,gl_node_idx> cnodeToNode(fp);
		auto *cnodeToSub = new ConnectivityT<size_t,int>(fp);

		std::map<gl_sub_idx, lc_sub_idx> glToLocSub2;
		for(int j=0; j<csubToSub2.csize(); ++j)
			glToLocSub2.insert(std::pair<int,int>(csubToSub2[j][0], j));
		std::map<gl_node_idx, lc_node_idx> glToLocNode;
		for(int j=0; j<cnodeToNode.csize(); ++j)
			glToLocNode.insert(std::pair<gl_node_idx,int>(cnodeToNode[j][0], j));
		numClusNodes = cnodeToNode.csize();

		auto *subnode = new SparsePairType1(glToLocSub2,csubToNode);
		subToNode_sparse = new SparseConnectivityType1(subnode);

		auto *nodesub = new SparsePairType2(glToLocNode,cnodeToSub);
		nodeToSub_sparse = new SparseConnectivityType2(nodesub);

#ifdef SOWER_DEBUG
		for(int i=0; i<structCom->numCPUs(); ++i) {
	  if(i == structCom->myID()) {
		//std::cerr << "subToElem_sparse = \n"; subToElem_sparse->print();
		//std::cerr << "elemToSub_sparse = \n"; elemToSub_sparse->print();
		std::cerr << "subToNode_sparse = \n"; subToNode_sparse->print();
		std::cerr << "nodeToSub_sparse = \n"; nodeToSub_sparse->print();
	  }
	}
#endif
	}
}

void GeoSource::readGlobalBinaryData()
{
	if(!subToSub || subToClus.size() == 0) {
		BinFileHandler fp2(connectivity_.c_str(), "rb");
		clusToSub =std::make_unique<Connectivity>(fp2);
		numClusters = clusToSub->csize();
		subToClus.resize(clusToSub->getNumTarget());
		for(int k = 0; k < clusToSub->csize(); k++)
			for(int i = 0; i < clusToSub->num(k); i++)
				subToClus[(*clusToSub)[k][i]] = k;
		int numSub;
#ifdef SUBTOSUBINFILE
		subToSub = new Connectivity(fp2);
	numSub = subToSub->csize();
#else
		fp2.read(&numSub, 1);
#endif
#ifdef SOWER_DEBUG
		std::cerr << "*** subToClus, from connectivity binary file: \n";
		for(int i=0; i<clusToSub->getNumTarget(); ++i) std::cerr << subToClus[i] << " "; std::cerr << std::endl;
		if(subToSub) { std::cerr << "*** subToSub, from connectivity binary file: \n"; subToSub->print(); }
#endif

		// build global to cluster subdomain map
		gl2ClSubMap = new int[numSub];
		for(int iClus = 0; iClus < numClusters; iClus++)  {
			int clusNum = 0;
			for(int iSub = 0; iSub < clusToSub->num(iClus); iSub++)
				gl2ClSubMap[ (*clusToSub)[iClus][iSub] ] = clusNum++;
		}

		fp2.read(&nGlobNodes, 1);
		domain->setNumNodes(nGlobNodes);
		int nGlobElems;
		fp2.read(&nGlobElems, 1);
		domain->setNumElements(nGlobElems);
#ifdef SOWER_DEBUG
		std::cerr << "*** global number of nodes = " << nGlobNodes << std::endl;
		std::cerr << "*** global number of elements = " << nGlobElems << std::endl;
#endif
		if(numClusters == 1) {
			numClusNodes = nGlobNodes;
			numClusElems = nGlobElems;
		}
		/*else {
		  filePrint(stdout," *** ERROR: only one cluster is supported \n");
		  exit(-1);
		}*/

#ifdef SOWER_SURFS
		int nGlobSurfs;
		fp2.read(&nGlobSurfs,1);
		domain->setNumSurfs(nGlobSurfs);
#ifdef SOWER_DEBUG
		std::cerr << "*** global number of surfaces = " << nGlobSurfs << std::endl;
#endif
/*
	if(nGlobSurfs > 0) {
	  Sower sower;
	  BinFileHandler *f = sower.openBinaryFile(0);
	  for(int isurf=0; isurf<domain->getNumSurfs(); isurf++) {
		int* dummy= 0;
		SurfaceEntity* surf = sower.read<SurfaceIO>(*f, isurf, dummy, true);
		if(dummy) { delete [] dummy; dummy= 0; } // not used
		domain->AddSurfaceEntity(surf,isurf);
	  }
	}
*/
#endif
	}
}

void GeoSource::computeClusterInfo(int glSub, Connectivity *_subToNode)
{
  int clusNum = subToClus[glSub];
  Connectivity *clusToNode = (_subToNode) ? clusToSub->transcon(_subToNode)
										  : clusToSub->transcon(subToNode);
  numClusNodes = clusToNode->num(clusNum);
  delete clusToNode;
  Connectivity *clusToElem = clusToSub->transcon(subToElem);
  numClusElems = clusToElem->num(clusNum);
  delete clusToElem;
}

void GeoSource::makeEframe(int ele, int refnode, double *d)
{
  d[0] = d[1] = d[2] = 0.0;
  int beam_nodes[2];
  elemSet[ele]->nodes(beam_nodes);
  Node *node1 = nodes[beam_nodes[0]];
  Node *node2 = nodes[refnode-1];
  d[3] = node2->x - node1->x;
  d[4] = node2->y - node1->y;
  d[5] = node2->z - node1->z;
  d[6] = d[7] = d[8] = 0.0;
}


void GeoSource::setAttributeGroup(int a, int g)
{
  group[g].attributes.push_back(a);
}

//-------------------------------------------------------

void GeoSource::setNodeGroup(int nn, int id)  {

  nodeGroup[id].insert(nn);
}

//-------------------------------------------------------

void GeoSource::setSurfaceGroup(int sn, int id)  {

  surfaceGroup[id].push_back(sn);
}

bool GeoSource::isInNodeGroup(int id) {
        bool flag = (nodeGroup.count(id) > 0) ? true : false;
        return flag;
}

bool GeoSource::isInAttrGroup(int id) {
        bool flag = (group.count(id) > 0) ? true : false;
        return flag;
}

void GeoSource::printGroups()
{
  for(std::map<int, Group >::iterator it = group.begin(); it != group.end(); ++it) {
	std::cerr << "group_id = " << it->first+1 << std::endl;
	std::cerr << "  number of attributes in this group =  " << it->second.attributes.size() << std::endl;
	std::cerr << "  number of random properties in this group =  " << it->second.randomProperties.size() << std::endl;
	std::cerr << "  attributes: ";
	for(int i = 0; i < int(it->second.attributes.size()); ++i) std::cerr << it->second.attributes[i]+1 << " ";
	std::cerr << "\n  random properties: ";
	for(int i = 0; i < int(it->second.randomProperties.size()); ++i) {
	  std::cerr << "rprop = " << it->second.randomProperties[i].rprop
		   << ", mean = " << it->second.randomProperties[i].mean
		   << ", std_dev = " << it->second.randomProperties[i].std_dev << std::endl;
	 }
	std::cerr << std::endl;
  }
}

void GeoSource::setGroupRandomProperty(int g, Rprop prop_type, double mean, double std_dev)
{
#ifndef SALINAS
  RandomProperty rp(prop_type, mean, std_dev);
  group[g].randomProperties.push_back(rp);
  sfem->updatendim();
#endif
}

bool GeoSource::elemOutput()
{
  for(int iInfo = 0; iInfo < geoSource->getNumOutInfo(); iInfo++) {
	if(oinfo[iInfo].averageFlg == 0) return true;
	if(oinfo[iInfo].type == OutputInfo::InXForce || oinfo[iInfo].type == OutputInfo::InYForce || oinfo[iInfo].type == OutputInfo::InZForce ||
	   oinfo[iInfo].type == OutputInfo::AXMoment || oinfo[iInfo].type == OutputInfo::AYMoment || oinfo[iInfo].type == OutputInfo::AZMoment)
	  return true;
  }
  return false;
}

bool GeoSource::energiesOutput()
{
  for(int iInfo = 0; iInfo < geoSource->getNumOutInfo(); iInfo++) {
	if(oinfo[iInfo].type == OutputInfo::Energies)
	  return true;
  }
  return false;
}

bool GeoSource::romExtForceOutput()
{
  for(int iInfo = 0; iInfo < geoSource->getNumOutInfo(); iInfo++) {
	if(oinfo[iInfo].type == OutputInfo::RomExtForce || oinfo[iInfo].type == OutputInfo::RomExtForce6)
	  return true;
  }
  return false;
}

bool GeoSource::noOutput(int x, int ndflag)
{
  bool noOut = true;
  for(int i = 0; i < numOutInfo; i++) {
	if(oinfo[i].ndtype != ndflag) continue;
	if(ndflag != 0 && oinfo[i].type != OutputInfo::Disp6DOF && oinfo[i].type != OutputInfo::Displacement) continue;
	if(oinfo[i].interval != 0 && x % oinfo[i].interval == 0) {
	  noOut = false;
	  break;
	}
  }
  return noOut;
}

void GeoSource::setLocalIndex(int j) {
  localIndex_ = j;
}

void GeoSource::setElementLumpingWeight(int iele, double value) {
  if(elementLumpingWeights_.size() < localIndex_+1) elementLumpingWeights_.resize(localIndex_+1);
  elementLumpingWeights_[localIndex_][iele] = value;
}

void GeoSource::pushBackStiffVec(double Kelem) {
  ReducedStiffVec.push_back(Kelem);
}

void GeoSource::pushBackUDEIMVec(double Uelem) {
  UDEIMBasisVec.push_back(Uelem);
}

void GeoSource::pushBackROMLMPCVec(double value) {
  ROMLMPCVec.push_back(value);
}

void GeoSource::setSampleNodesAndSlots(int node, int dof){
  nodeDofSlotPairVec_.push_back(std::make_pair(node,dof));
}

void GeoSource::setSampleElemsAndDOFs(int node, int dof){
  elemDofPairVec_.push_back(std::make_pair(node,dof));
}


void GeoSource::getARubberLambdaMu(double omega, complex<double> *lambda,
								   complex<double> *mu) {

  int i_arubber = 0;
  SPropContainer::iterator it = sProps.begin();
  while(it != sProps.end()) {
	StructProp* p = &(it->second);
	if(p->E0!=0.0 || p->mu0!=0.0) {
	  ARubberF ar(0,omega,
			  p->E0,p->dE,p->mu0,p->dmu,
			  p->eta_E,p->deta_E,p->eta_mu,p->deta_mu);

	  lambda[i_arubber] = ar.d_lambda(0);
	  mu[i_arubber] = ar.d_mu(0);
	  i_arubber++;
	}
	it++;
  }
}


#ifdef SOWER_SURFS 
void GeoSource::readDistributedSurfs(int subNum)
{
  Sower sower;
  auto f = sower.openBinaryFile(subNum);
  //filePrint(stderr," ... Read surfaces data\n");
  for(int isurf=0; isurf<domain->getNumSurfs(); isurf++) {
	//filePrint(stderr," ... Read data of surface %d\n",isurf);
	int* dummy= 0;
	SurfaceEntity* surf = sower.read<SurfaceIO>(*f, isurf, dummy, true);
	if(dummy) { delete [] dummy; dummy= 0; } // not used
	domain->AddSurfaceEntity(surf,isurf);
  }
}
#endif

