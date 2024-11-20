// ------------------------------------------------------------
// HB - 06/30/03 
// ------------------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// STL
#include <algorithm>
#include <vector>
#include <map>

// FEM headers
#include <Utils.d/DistHelper.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/DistrGeomState.h>
#include <Driver.d/SubDomain.h>

#ifdef SOWER_SURFS
#include <Utils.d/BinFileHandler.h>
#endif

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTORS
// -----------------------------------------------------------------------------------------------------
SurfaceEntity::SurfaceEntity():ElemSet()
{
   Initialize();
}

SurfaceEntity::SurfaceEntity(int _Id):ElemSet()
{
   Initialize();
   Id = _Id;
}

SurfaceEntity::SurfaceEntity(const SurfaceEntity &other)
 : Id(other.Id),
   ElemSet(other.ElemSet.size()), 
   nNodes(other.nNodes),
   nVertices(other.nVertices),
   LocalNumbering(other.LocalNumbering),
   nTri3(other.nTri3),
   nTri6(other.nTri6),
   nQuad4(other.nQuad4),
   nQuad8(other.nQuad8),
   ReverseNormals(other.ReverseNormals),
   IsShellFace(other.IsShellFace),
   ShellThickness(other.ShellThickness),
   PreserveOrdering(other.PreserveOrdering)
{
  for(int i=0; i<other.ElemSet.size(); ++i) {
    if(FaceElement *faceEl = other.ElemSet[i]) ElemSet.elemadd(i, faceEl->clone());
    else ElemSet.elemadd(i, (FaceElement*)0);
  }

  NodeCoordMap = (other.NodeCoordMap) ? new std::map<int,Node>(*other.NodeCoordMap) : 0;

  if(other.gNodeIds) {
    gNodeIds = new int[nNodes];
    for(int i=0; i<nNodes; ++i) gNodeIds[i] = other.gNodeIds[i];
  }
  else gNodeIds = 0;

  LlToGlNodeMap = (other.LlToGlNodeMap) ? gNodeIds : 0;
  GlToLlNodeMap = (other.GlToLlNodeMap) ? new std::map<int,int>(*other.GlToLlNodeMap) : 0;

  if(other.NodeSet) {
    int size = other.NodeSet->size();
    NodeSet = new CoordSet(size);
    for(int i=0; i<size; ++i) {
      if(Node *node = (*other.NodeSet)[i]) NodeSet->nodeadd(i, *node);
    }
  }
  else NodeSet = 0;

  if(other.gVertexIds) {
    gVertexIds = new int[nVertices];
    for(int i=0; i<nVertices; ++i) gVertexIds[i] = other.gVertexIds[i];
  }
  else gVertexIds = 0;

  LlToGlVertexMap = (other.LlToGlVertexMap) ? gVertexIds : 0;

  GlToLlVertexMap = (other.GlToLlVertexMap) ? new std::map<int,int>(*other.GlToLlVertexMap) : 0;

  if(other.LlVertexToLlNodeMap) {
    LlVertexToLlNodeMap = new int[nVertices];
    for(int i=0; i<nVertices; ++i) LlVertexToLlNodeMap[i] = other.LlVertexToLlNodeMap[i];
  }
  else LlVertexToLlNodeMap = 0;

  ACMEBlocksMap = (other.ACMEBlocksMap) ? new Connectivity(*other.ACMEBlocksMap) : 0;

#ifdef HB_NODALNORMAL
  if(other.NdNormals) {
    NdNormals = new double[nNodes][3];
    for(int i=0; i<nNodes; ++i) {
      NdNormals[i][0] = other.NdNormals[i][0];
      NdNormals[i][1] = other.NdNormals[i][1];
      NdNormals[i][2] = other.NdNormals[i][2];
    }
  }
  else NdNormals = 0;
#endif
}

// -----------------------------------------------------------------------------------------------------
//                                           DESTRUCTORS 
// -----------------------------------------------------------------------------------------------------
SurfaceEntity::~SurfaceEntity()
{
  Id = -1;
  if(gVertexIds) { delete [] gVertexIds; gVertexIds = 0; nVertices = 0; }
  if(gNodeIds)   { delete [] gNodeIds  ; gNodeIds   = 0; nNodes    = 0; }
  LlToGlVertexMap = 0;
  LlToGlNodeMap   = 0;
  if(GlToLlVertexMap) { delete GlToLlVertexMap; GlToLlVertexMap = 0; }
  if(GlToLlNodeMap  ) { delete GlToLlNodeMap;   GlToLlNodeMap   = 0; }
  if(LlVertexToLlNodeMap) { delete [] LlVertexToLlNodeMap; LlVertexToLlNodeMap = 0; } 
  if(NodeSet) { delete NodeSet; NodeSet = 0; }
  if(ACMEBlocksMap) { delete ACMEBlocksMap; ACMEBlocksMap = 0; }
  nTri3 = nTri6 = nQuad4 = nQuad8 = 0;
  //ElemSet.~FaceElemSet();
}

// -----------------------------------------------------------------------------------------------------
//                                   INITILIZATION & CLEAN/CLEAR METHODS
// -----------------------------------------------------------------------------------------------------
void
SurfaceEntity::Initialize()
{
  Id = -1;
  NodeCoordMap = 0; 
  LocalNumbering = false;
 
  nNodes         = 0;
  gNodeIds       = 0;
  LlToGlNodeMap  = 0;
  GlToLlNodeMap  = 0;

  NodeSet        = 0;
  
  nVertices          = 0;
  gVertexIds         = 0;
  LlToGlVertexMap    = 0;
  GlToLlVertexMap    = 0;
  LlVertexToLlNodeMap= 0;

  nTri3 = nTri6 = 0;
  nQuad4= nQuad8= 0;
  ACMEBlocksMap = 0; 
#ifdef HB_NODALNORMAL
  NdNormals = 0;
#endif
  ReverseNormals = false;
  IsShellFace = false;
  ShellThickness = 0.0;
  PreserveOrdering = false;
}

// -----------------------------------------------------------------------------------------------------
//                                          SETUP & UPDATE METHODS
// -----------------------------------------------------------------------------------------------------
void
SurfaceEntity::SetUpData(CoordSet* cs)
{
  SetUpNodeData();
  SetUpVertexData();
  MakeLlVertexToLlNodeMap();
  MakeACMEBlocksMap();
  if(cs) ExtractNodeSet(*cs);
}

void
SurfaceEntity::SetUpVertexData()
{
  MakeGlVertexIds();
  MakeVertexMaps();
}

void
SurfaceEntity::SetUpNodeData()
{
  MakeGlNodeIds();
  MakeNodeMaps();
}

void
SurfaceEntity::UpdateNodeData(GeomState *geomState)
{
  for(int inode=0; inode<nNodes; inode++) {
    NodeSet->getNode(inode).x = (*geomState)[ gNodeIds[inode] ].x;
    NodeSet->getNode(inode).y = (*geomState)[ gNodeIds[inode] ].y;
    NodeSet->getNode(inode).z = (*geomState)[ gNodeIds[inode] ].z;
  }
}

void
SurfaceEntity::UpdateNodeData(DistrGeomState *geomState, SubDomain **sd)
{
  for(int inode=0; inode<nNodes; inode++) {
    int glNode = gNodeIds[inode];
    for(int j = 0; j < geomState->getNumSub(); ++j) {
      int locNode = sd[j]->globalToLocal(glNode);
      if(locNode > -1) {
        NodeSet->getNode(inode).x = (*(*geomState)[j])[locNode].x;
        NodeSet->getNode(inode).y = (*(*geomState)[j])[locNode].y;
        NodeSet->getNode(inode).z = (*(*geomState)[j])[locNode].z;
      }
    }
  }
}

void
SurfaceEntity::MakeGlVertexIds()
{
  // -> UNDERLAYING LINEAR FACE ELEMENT NODES
  if(gVertexIds) { delete [] gVertexIds; gVertexIds = 0; }

  // Get all the Vertices
  int maxElem = ElemSet.size();
  int nVerts = 0;
  for(int iel=0; iel<maxElem; iel++)
    if(ElemSet[iel]) nVerts += ElemSet[iel]->nVertices();

  int* Vertices = new int[nVerts]; nVerts = 0;
  for(int iel=0; iel<maxElem; iel++) {
    if(ElemSet[iel]) {
      ElemSet[iel]->GetVertices(&Vertices[nVerts]);
      nVerts += ElemSet[iel]->nVertices();
    }
  }

  // Eliminate duplicate vertices
  if(PreserveOrdering) {
    std::set<int> s;
    int nVertices=0;
    for(int i=0; i<nVerts; ++i)
      if(s.find(Vertices[i]) == s.end()) {
        s.insert(Vertices[i]);
        Vertices[nVertices++] = Vertices[i];
      }
  }
  else {
    std::sort(Vertices,Vertices+nVerts);
    nVertices = int(std::unique(Vertices,Vertices+nVerts)-Vertices);
  }

  // Copy into the gVertexIds array
  gVertexIds = new int[nVertices];
  std::copy(Vertices,Vertices+nVertices,gVertexIds);

  if(Vertices) delete [] Vertices;
}

void
SurfaceEntity::MakeGlNodeIds()
{
  if(gNodeIds) { delete [] gNodeIds; gNodeIds = 0; }

  // Get all the Nodes
  int maxElem = ElemSet.size();
  int nNds = 0;
  for(int iel=0; iel<maxElem; iel++) {
    if(ElemSet[iel]) nNds += ElemSet[iel]->nNodes();
  }

  int* Nds = new int[nNds]; nNds = 0;
  for(int iel=0; iel<maxElem; iel++) {
    if(ElemSet[iel]) {
      ElemSet[iel]->GetNodes(&Nds[nNds]);
      nNds += ElemSet[iel]->nNodes();
    }
  }

  // Eliminate duplicate nodes
  if(PreserveOrdering) {
    std::set<int> s;
    nNodes = 0;
    for(int i=0; i<nNds; ++i)
      if(s.find(Nds[i]) == s.end()) {
        s.insert(Nds[i]);
        Nds[nNodes++] = Nds[i];
      }
  }
  else {
    std::sort(Nds,Nds+nNds);
    nNodes = int(std::unique(Nds,Nds+nNds)-Nds);
  }

  // Copy into the gNodeIds array
  gNodeIds = new int[nNodes];
  std::copy(Nds,Nds+nNodes,gNodeIds);

  if(Nds) delete [] Nds;
}

void
SurfaceEntity::MakeNodeMaps()
{
  if(!gNodeIds) MakeGlNodeIds(); 

  LlToGlNodeMap = gNodeIds;

  if(GlToLlNodeMap) { delete GlToLlNodeMap; GlToLlNodeMap = 0; }
  GlToLlNodeMap = new std::map<int,int>(); 
  for(int inode=0; inode<nNodes; inode++)
    (*GlToLlNodeMap)[gNodeIds[inode]] = inode;   
}

void
SurfaceEntity::MakeVertexMaps()
{
  if(!gVertexIds) MakeGlVertexIds(); 

  LlToGlVertexMap = gVertexIds;
  
  if(GlToLlVertexMap) { delete GlToLlVertexMap; GlToLlVertexMap= 0; }
  GlToLlVertexMap = new std::map<int,int>(); 
  for(int ivertex=0; ivertex<nVertices; ivertex++) 
    (*GlToLlVertexMap)[gVertexIds[ivertex]] = ivertex;
}

void
SurfaceEntity::MakeLlVertexToLlNodeMap()
{
  if(!gNodeIds)      MakeGlNodeIds();
  if(!gVertexIds)    MakeGlVertexIds();
  if(!GlToLlNodeMap) MakeNodeMaps();

  if(LlVertexToLlNodeMap) { delete [] LlVertexToLlNodeMap; LlVertexToLlNodeMap = 0; }
  LlVertexToLlNodeMap = new int[nVertices];
  for(int ivertex=0; ivertex<nVertices; ivertex++)
    LlVertexToLlNodeMap[ivertex] = (*GlToLlNodeMap)[LlToGlVertexMap[ivertex]];
}

void
SurfaceEntity::MakeNodeSet(CoordSet& cs)
{
  if(!gNodeIds) MakeGlNodeIds();
  if(NodeSet) { delete NodeSet; NodeSet = 0; }
  NodeSet = new CoordSet(nNodes);

  for(int inode=0; inode<nNodes; inode++) {
    Node &node = cs.getNode(gNodeIds[inode]);
    NodeSet->nodeadd(inode, node); 
  }
}

void
SurfaceEntity::MakeNodeSet(std::map<int,Node>& NodeCoordMap)
{
  if(!gNodeIds) MakeGlNodeIds();
  if(NodeSet) { delete NodeSet; NodeSet = 0; }
  NodeSet = new CoordSet(nNodes);
                                                                                                     
  for(int inode=0; inode<nNodes; inode++){
    std::map<int,Node>::iterator Imap = NodeCoordMap.find(gNodeIds[inode]);
    if(Imap!=NodeCoordMap.end()){
      Node &node = (*Imap).second;
      NodeSet->nodeadd(inode, node);
    } else {
      std::cerr << " ### PB in SurfaceEntity::MakeNodeSet(map<int,Node>&) ###" << std::endl;
      std::cerr << " ### NODE " <<gNodeIds[inode]<<" NOT FOUND IN INPUT NodeCoordMap ###" << std::endl;
    }
  }
}

void
SurfaceEntity::ExtractNodeSet(CoordSet& cs)
{
  MakeNodeSet(cs);
}

void
SurfaceEntity::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  ElemSet.Renumber(OldToNewNodeIds);
}

void
SurfaceEntity::Renumber()
{
  if(!GlToLlNodeMap) MakeNodeMaps(); 
  Renumber((*GlToLlNodeMap));
  LocalNumbering = true;
}

void
SurfaceEntity::MakeACMEBlocksMap()
{
  if(ACMEBlocksMap) delete ACMEBlocksMap;
  
  // Step 1. Count number of Tri3, Tri6, Quad4 & Quad8 face els. in the ElemSet
  //         & store their position (index) in the ElemSet array
  nTri3 = nTri6 = nQuad4 = nQuad8 = 0;
  int nElems = GetnFaceElems();
  std::vector<int> IndexTri3;  IndexTri3.reserve(512); // to avoid too many
  std::vector<int> IndexTri6;  IndexTri6.reserve(512); // automatic memory
  std::vector<int> IndexQuad4; IndexQuad4.reserve(512);// re-allocation
  std::vector<int> IndexQuad8; IndexQuad8.reserve(512);
  for(int iel=0; iel<nElems; iel++){
   int etype = ElemSet[iel]->GetFaceElemType();
   switch(etype)
   {
     case FaceElement::QUADFACEL4:
       IndexQuad4.push_back(iel);
       nQuad4++;
       break;
     case FaceElement::QUADFACEQ8:
     case FaceElement::QUADFACEQ9:
     case FaceElement::QUADFACEC12:
       IndexQuad8.push_back(iel);
       nQuad8++;
       break;
     case FaceElement::TRIFACEL3:
       IndexTri3.push_back(iel);
       nTri3++;
       break;
     case FaceElement::TRIFACEQ6:
     case FaceElement::TRIFACEC10:
       IndexTri6.push_back(iel);
       nTri6++;
       break;
     default:
       std::cerr << "Face element Type " << etype << " is NOT supported/implemented." << std::endl;
       exit(-1);
       return;
   }
  }

  // Step 2. Create & fill the ACMEBlocksMap connectivity:
  //         -> ACMEBlocksMap[0][i] gives the index of the Quad4 els (if any) in the ElemSet array
  //         -> ACMEBlocksMap[1][i] gives the index of the Quad8 els (if any) in the ElemSet array
  //         -> ACMEBlocksMap[2][i] gives the index of the Tri3 els (if any) in the ElemSet array
  //         -> ACMEBlocksMap[3][i] gives the index of the Tri6 els (if any) in the ElemSet array
  int nACMEBlocks = 4; // 1st block: Quad4 
                       // 2nd block: Quad8
                       // 3th block: Tri3 
                       // 4th block: Tri6 
  int nTri = nTri3+nTri6;  // Tri3 + Tri6 face els.
  int nQuad= nQuad4+nQuad8;// Quad4 + Quad8 face els.

  int* ptr = new int[nACMEBlocks+1];
  int* target = new int[nTri+nQuad];
  ptr[0] = 0;
  int count = 0;
  // fill Quad4 indexes
  for(int i=0; i<nQuad4; i++)
    target[count++] = IndexQuad4[i];
  ptr[1] = ptr[0]+nQuad4;
  // fill Quad8 indexes
  for(int i=0; i<nQuad8; i++)
    target[count++] = IndexQuad8[i];
  ptr[2] = ptr[1]+nQuad8;
  // fill Tri3 indexes
  for(int i=0; i<nTri3; i++)
    target[count++] = IndexTri3[i];
  ptr[3] = ptr[2]+nTri3;
  // fill Tri6 indexes
  for(int i=0; i<nTri6; i++)
    target[count++] = IndexTri6[i];
  ptr[4] = ptr[3]+nTri6;
  
  ACMEBlocksMap = new Connectivity(nACMEBlocks,ptr,target);
}

/*
int
SurfaceEntity::GetnFaceElems(int etype)
{
   switch(etype)
   {
     case FaceElement::QUADFACEL4:
       return nQuad4;
       break;
     case FaceElement::QUADFACEQ8:
       return nQuad8;
       break;
     case FaceElement::TRIFACEL3:
       return nTri3;
       break;
     case FaceElement::TRIFACEQ6:
       return nTri6;
       break;
     default:
       std::cerr << "Face element Type " << etype << " is NOT supported/implemented." << std::endl;
       exit(-1);
       return;
   }
}
*/

int 
SurfaceEntity::FillACMEFaceBlocks(int* face_connectivity, std::map<int,int>& OldToNewNodeIds, bool vertexOnly)
{
  if(!ACMEBlocksMap) MakeACMEBlocksMap();
  int offset = 0;
  if(vertexOnly) {
    for(int iBlk=0; iBlk<ACMEBlocksMap->csize(); iBlk++){
      for(int i=0, iend=ACMEBlocksMap->num(iBlk); i<iend; i++){
        int iel = (*ACMEBlocksMap)[iBlk][i];
        ElemSet[iel]->GetVertices(&face_connectivity[offset], OldToNewNodeIds);
        offset += ElemSet[iel]->nVertices();
      }
    }
  } else {
    for(int iBlk=0; iBlk<ACMEBlocksMap->csize(); iBlk++){
      for(int i=0, iend=ACMEBlocksMap->num(iBlk); i<iend; i++){
        int iel = (*ACMEBlocksMap)[iBlk][i];
        ElemSet[iel]->GetNodes(&face_connectivity[offset], OldToNewNodeIds);
        offset += ElemSet[iel]->nNodes();
      }
    }
  }
  /*for(int iBlk=0; iBlk<ACMEBlocksMap->csize(); iBlk++){
    for(int i=0; i<ACMEBlocksMap->num(iBlk); i++){
      int iel = (*ACMEBlocksMap)[iBlk][i];
      if(vertexOnly) {
        ElemSet[iel]->GetVertices(&face_connectivity[offset], OldToNewNodeIds);
        offset += ElemSet[iel]->nVertices();
      } else {
        ElemSet[iel]->GetNodes(&face_connectivity[offset], OldToNewNodeIds);
        offset += ElemSet[iel]->nNodes();
      }
    }
  }*/
  return(offset);
}

void
SurfaceEntity::Reset(CoordSet* cs)
{
  // this function should be called after elements are removed from ElemSet, e.g. due to element deletion
  ElemSet.repack();
  int *GlNodeIds = GetPtrGlNodeIds();
  std::map<int,int> LlToGlNodeMap;
  for(int i=0; i<GetnNodes(); ++i) LlToGlNodeMap[i] = GlNodeIds[i];
  Renumber(LlToGlNodeMap);
  SetUpData(cs);
  Renumber();
}

// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS
// -----------------------------------------------------------------------------------------------------
void 
SurfaceEntity::SetId(int _Id) { Id = _Id; }

void
SurfaceEntity::SetReverseNormals(bool _ReverseNormals) { ReverseNormals = _ReverseNormals; }

void
SurfaceEntity::SetIsShellFace(bool _IsShellFace) { IsShellFace = _IsShellFace; }

void
SurfaceEntity::SetShellThickness(double _ShellThickness) { ShellThickness = _ShellThickness; }

void
SurfaceEntity::SetPreserveOrdering(bool _PreserveOrdering) { PreserveOrdering = _PreserveOrdering; }

void
SurfaceEntity::AddFaceElement(int num, int etype, int nnodes, int* nodes)
{
  ElemSet.elemadd(ElemSet.last(), etype, nnodes, nodes);
}

void
SurfaceEntity::AddFaceElement(FaceElement* FaceElem)
{
  ElemSet.elemadd(ElemSet.last(), FaceElem);
}

void
SurfaceEntity::RemoveFaceElement(int num)
{
  ElemSet.remove(num);
}

void 
SurfaceEntity::SetPtrNodeSet(CoordSet* ndSet) 
{ 
  NodeSet = ndSet; 
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int 
SurfaceEntity::GetId() { return(Id); }

int 
SurfaceEntity::ID() { return(Id); }

bool
SurfaceEntity::GetReverseNormals() { return ReverseNormals; }

bool
SurfaceEntity::GetIsShellFace() { return IsShellFace; }

double
SurfaceEntity::GetShellThickness() { return ShellThickness; }

FaceElemSet*
SurfaceEntity::GetPtrFaceElemSet() { return(&ElemSet); }

FaceElemSet&
SurfaceEntity::GetFaceElemSet() { return(ElemSet); }

CoordSet*
SurfaceEntity::GetPtrNodeSet() { return(NodeSet); }

CoordSet&
SurfaceEntity::GetNodeSet() { return(*NodeSet); }

int
SurfaceEntity::GetnFaceElems() { return(ElemSet.nElems()); }

int
SurfaceEntity::nFaceElements() { return(ElemSet.nElems()); }

int
SurfaceEntity::GetnNodes() { return(nNodes); }

int 
SurfaceEntity::GetnVertices() { return(nVertices); }

bool
SurfaceEntity::IsRenumbered() { return(LocalNumbering); }

int
SurfaceEntity::GetGlNodeId(int inode) { return(gNodeIds[inode]); }

Node&
SurfaceEntity::GetNode(int inode)
{
  return(NodeSet->getNode(inode));
}

int
SurfaceEntity::GetGlVertexId(int ivertex) { return(gVertexIds[ivertex]); }

Node&
SurfaceEntity::GetVertex(int ivertex)
{
  return(NodeSet->getNode(GetLlVertexInLlNode(ivertex)));
}

int
SurfaceEntity::GetLlVertexInLlNode(int ivertex)
{
  return(LlVertexToLlNodeMap[ivertex]);
}

int*
SurfaceEntity::GetPtrGlNodeIds() 
{ 
  if(!gNodeIds) MakeGlNodeIds();
  return(gNodeIds); 
}

int*
SurfaceEntity::GetPtrGlVertexIds() 
{ 
  if(!gVertexIds) MakeGlVertexIds();
  return(gVertexIds); 
}

int*
SurfaceEntity::GetPtrLlToGlNodeMap()
{
  if(!LlToGlNodeMap) MakeNodeMaps();
  return(LlToGlNodeMap);
}

int*
SurfaceEntity::GetPtrLlToGlVertexMap()
{
  if(!LlToGlVertexMap) MakeVertexMaps();
  return(LlToGlVertexMap);
}

std::map<int, int>*
SurfaceEntity::GetPtrGlToLlNodeMap()
{
  if(!GlToLlNodeMap) MakeNodeMaps();
  return(GlToLlNodeMap);
}

std::map<int, int>*
SurfaceEntity::GetPtrGlToLlVertexMap()
{
  if(!GlToLlVertexMap) MakeVertexMaps();
  return(GlToLlVertexMap);
}

Connectivity*
SurfaceEntity::GetPtrACMEBlocksMap()
{
  if(!ACMEBlocksMap) MakeACMEBlocksMap();
  return(ACMEBlocksMap);
}

// -----------------------------------------------------------------------------------------------------
//                                            PRINT METHODS
// -----------------------------------------------------------------------------------------------------
void
SurfaceEntity::PrintFaceElemSet()
{
   ElemSet.print();
}

void 
SurfaceEntity::Print()
{
   filePrint(stderr," ------------------------------- \n");
   filePrint(stderr," * surface ID: %d\n", Id);
   filePrint(stderr," ------------------------------- \n");
   if(LocalNumbering)
     filePrint(stderr," * surface face elem. set in local node numbering: \n");
   else 
     filePrint(stderr," * surface face elem. set: \n");
   PrintFaceElemSet();

   if(nNodes) filePrint(stderr," * surface contains %d nodes\n",nNodes);
   if(gNodeIds){
     if(NodeSet)
       for(int inode=0; inode<nNodes; inode++) {
         filePrint(stderr,"  global Id %6d --- local Id %6d: %e %e %e\n",
                   gNodeIds[inode]+1,inode+1,NodeSet->getNode(inode).x,
                   NodeSet->getNode(inode).y,NodeSet->getNode(inode).z);
       }
     else
       for(int inode=0; inode<nNodes; inode++)
         filePrint(stderr,"  global Id %6d --- local Id %6d\n",gNodeIds[inode]+1,inode+1);
   }
   if(nVertices) filePrint(stderr," * surface contains %d vertices\n",nVertices);
   if(gVertexIds){
     for(int ivertex=0; ivertex<nVertices; ivertex++){
       filePrint(stderr,"  global Id %6d --- local Id %6d\n",gVertexIds[ivertex]+1,ivertex+1);
     }
   }
   if(LlVertexToLlNodeMap){
     filePrint(stderr," * Local vertex to local node map:\n");
     for(int ivertex=0; ivertex<nVertices; ivertex++){
       filePrint(stderr,"  local vertex %6d --- local node %6d\n",ivertex+1,LlVertexToLlNodeMap[ivertex]+1);
     }
   }
   if(ACMEBlocksMap){
     filePrint(stderr," * surface face els set contains %d ACME blocks:\n",ACMEBlocksMap->csize());
     for(int iblk=0; iblk<ACMEBlocksMap->csize(); iblk++){
       filePrint(stderr,"  ACME block %d contains %6d face els.\n",iblk+1,ACMEBlocksMap->num(iblk));
     }
     ACMEBlocksMap->print();
   }
#ifdef HB_NODALNORMAL
   if(NdNormals) {
     filePrint(stderr," * Nodal normals:\n");
     if(gNodeIds) {
       for(int i=0; i<nNodes; i++)
         filePrint(stderr,"  global Id %6d --- local Id %6d, normal = %e %e %e\n",
                   gNodeIds[i]+1,i+1,NdNormals[i][0],NdNormals[i][1],NdNormals[i][2]);
     } else {
       for(int i=0; i<nNodes; i++)
         filePrint(stderr,"  local Id %6d, normal = %e %e %e\n",
                   i+1,NdNormals[i][0],NdNormals[i][1],NdNormals[i][2]);
     }
   }
#endif
   filePrint(stderr," ------------------------------- \n");
}

void
SurfaceEntity::PrintFaceNormal(CoordSet& cs)
{
  double Normal[3];
  double m[2];
  int nElems = GetnFaceElems();
  for(int iel=0; iel<nElems; iel++){
    switch(ElemSet[iel]->GetFaceElemType()) // get center coordinates in paramatric domain
    {
      case FaceElement::QUADFACEL4:
      case FaceElement::QUADFACEQ8:
      case FaceElement::QUADFACEQ9:
      case FaceElement::QUADFACEC12:
        m[0] = m[1] = 0.0;
        break;
      case FaceElement::TRIFACEL3:
      case FaceElement::TRIFACEQ6:
      case FaceElement::TRIFACEC10:
        m[0] = m[1] = 1./3.;
        break;
      default:
        std::cerr << "Face element Type " << ElemSet[iel]->GetFaceElemType() << " is NOT supported/implemented." << std::endl;
        exit(-1);
        return;
    }
    ElemSet[iel]->GetIsoParamMappingNormalAndJacobian(Normal, m, cs);
    filePrint(stderr," - face element %6d, normal = %e, %e, %e\n",iel+1,Normal[0],Normal[1],Normal[2]);
  }
}

#ifdef SOWER_SURFS
void 
SurfaceEntity::WriteSower(BinFileHandler& file)
{
  // Write (global) surface Id
  file.write(&Id, 1);

  // Write faces els (elId in Surface, el. type, el. nodes)
  int nElems = GetnFaceElems();
  file.write(&nElems,1);
  ElemSet.WriteSower(file,LlToGlNodeMap); // using global node Ids
 
  // Write nodes'coordinates
  file.write(&nNodes,1);
  for(int i=0; i<nNodes; i++) {
    file.write(&(gNodeIds[i]),1); // global node Id
    double xyz[3] = {NodeSet->getNode(i).x,NodeSet->getNode(i).y,NodeSet->getNode(i).z};
    file.write(xyz,3); 
  }
}
#endif

// Comute nodal normal using the "Medial-Quadric" approach 
// See "Parallel Feature-Preserving Mesh Smoothing" 
// by X. Jiao & P. J. Alexander in Computational Science & Engineering 2005
// The current implemented below uses a (local) area-based weighting 
// EXPERIMENTAL ...
#ifdef HB_NODALNORMAL
double vect3DNormalize(double* v) 
{
  double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(norm!=0.0) { 
    double inorm = 1.0/norm;
    v[0] *= inorm; v[1] *= inorm; v[2] *= inorm; 
  }
  return(norm);
}
#include <Utils.d/linkfc.h>
extern "C" {
  void _FORTRAN(dspev)(const char& JOBZ, const char& UPLO, const int& N,
                       double* AP, double* W, double* Z, const int& LDZ, double* WORK, int& INFO);
}

double*
SurfaceEntity::ComputeNodalNormals(CoordSet& cs)
{
  //double (*NdNormals)[3] = new double[nNodes][3];
  if(NdNormals) delete [] NdNormals; 
  NdNormals = new double[nNodes][3];
  int* NdCard = new int[nNodes]; // nodal "cardinality" (multiplicity of each node)
  double (*matA)[6] = new double[nNodes][6]; // nodal 3x3 matrix A in packed format (upper part, column-wise)
  double (*rhsb)[3] = new double[nNodes][3]; // 3x1 rhs b at each node
  // Step 0. Initialize 
  for(int i=0; i<nNodes; i++) {
    NdCard[i] = 0;
    NdNormals[i][0] = NdNormals[i][1] = NdNormals[i][2] = 0.0;
    rhsb[i][0] = rhsb[i][1] = rhsb[i][2] = 0.0;
    matA[i][0] = matA[i][1] = matA[i][2] = matA[i][3] = matA[i][4] = matA[i][5] = 0.0;
  }
  // Step 1. Compute 3x3 matrix A & 3x1 rhs b at each node 
  double Nrml[3];
  for(int iel=0, nElems=GetnFaceElems(); iel<nElems; iel++){
    double (*refCoords)[2] = reinterpret_cast<double (*)[2]>(ElemSet[iel]->ViewRefCoords());
    for(int i=0, nElNds=ElemSet[iel]->nNodes(); i<nElNds; i++) {
      int nd = ElemSet[iel]->GetNode(i);
      NdCard[nd]++;
      double jac = ElemSet[iel]->GetIsoParamMappingNormalAndJacobian(Nrml, refCoords[i], cs);
      //printf("surface %d, el %6d, node %2d, jac = %e, normal = %e %e %e\n",Id,iel+1,i+1,jac,Nrml[0],Nrml[1],Nrml[2]);
      matA[nd][0] += jac*Nrml[0]*Nrml[0];
      matA[nd][1] += jac*Nrml[0]*Nrml[1];
      matA[nd][2] += jac*Nrml[1]*Nrml[1];
      matA[nd][3] += jac*Nrml[0]*Nrml[2];
      matA[nd][4] += jac*Nrml[1]*Nrml[2];
      matA[nd][5] += jac*Nrml[2]*Nrml[2];
      rhsb[nd][0] -= jac*Nrml[0];
      rhsb[nd][1] -= jac*Nrml[1];
      rhsb[nd][2] -= jac*Nrml[2];
    }
  }
  // Step 2. Compute eigen values of A & solve "reduced" system for each node 
  double atol = 1.E-10;
  double rtol = 1.E-8;
  for(int i=0; i<nNodes; i++) {
    // Comute eigen values vectors of 3x3 matrix A 
    double eigVals[3], eigVects[9], work[9];
    int info;
    _FORTRAN(dspev)('V', 'U', 3, matA[i], eigVals, eigVects, 3, work, info);
    if(info) printf(" *** ERROR from dspev at node %d, info = %d\n",i,info);
    //printf(" surface %d, node %6d, eigen values = %e %e %e\n",Id,i+1,eigVals[0],eigVals[1],eigVals[2]);
//#define DEBUG
#ifdef DEBUG
      FILE* fdebug = stderr;
      fprintf(fdebug," --- surface %d, node %6d\n",Id,i+1);
      fprintf(fdebug," * A = %12.6e  %12.6e  %12.6e\n",matA[i][0],matA[i][1],matA[i][2]);
      fprintf(fdebug,"       x             %12.6e  %12.6e\n",matA[i][3],matA[i][4]);
      fprintf(fdebug,"       x             x             %12.6e\n",matA[i][5]);
      fprintf(fdebug," * eig values = %e %e %e\n",eigVals[0],eigVals[1],eigVals[2]);
      fprintf(fdebug," * eig vectors= %e %e %e\n",eigVects[0],eigVects[1],eigVects[2]);
      fprintf(fdebug," *              %e %e %e\n",eigVects[3],eigVects[4],eigVects[5]);
      fprintf(fdebug," *              %e %e %e\n",eigVects[6],eigVects[7],eigVects[8]);
#endif
    // Solve "reduced" system -> direction of nodal normal
    double (*eigVs)[3] = reinterpret_cast<double (*)[3]>(eigVects);
    for(int k=2; k>=0; k--) { // dspev outputs eigen values in ascending order
      if(eigVals[k]<0) eigVals[k] = 0.0; // make sure all the eigen values are >=0
      if(eigVals[k]>rtol*eigVals[2]+atol) {
        double sdot = (eigVs[k][0]*rhsb[i][0]
                      +eigVs[k][1]*rhsb[i][1]
                      +eigVs[k][2]*rhsb[i][2])/eigVals[k];
        NdNormals[i][0] -= sdot*eigVs[k][0];
        NdNormals[i][1] -= sdot*eigVs[k][1];
        NdNormals[i][2] -= sdot*eigVs[k][2];
      }
    }
    vect3DNormalize(NdNormals[i]);
  }
  // Delete working arrays
  if(matA) delete [] matA;
  if(rhsb) delete [] rhsb;
  if(NdCard) delete [] NdCard;
  return(reinterpret_cast<double*>(NdNormals));
}

void
SurfaceEntity::PrintNodalNormals()
{
  filePrint(stderr," ------------------------------- \n");
  filePrint(stderr," * Nodal normals of surface %d\n", Id);
  filePrint(stderr," ------------------------------- \n");
  if(gNodeIds)
      for(int i=0; i<nNodes; i++)
        filePrint(stderr,"  global Id %6d --- local Id %6d, normal = %e %e %e\n",
                  gNodeIds[i]+1,i+1,NdNormals[i][0],NdNormals[i][1],NdNormals[i][2]);
  else
      for(int i=0; i<nNodes; i++)
        filePrint(stderr,"  local Id %6d, normal = %e %e %e\n",
                  i+1,NdNormals[i][0],NdNormals[i][1],NdNormals[i][2]);
}
#endif
