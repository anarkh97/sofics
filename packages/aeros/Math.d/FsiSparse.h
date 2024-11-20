#ifndef _FSISPARSE_H_
#define _FSISPARSE_H_

#include <Driver.d/Mpc.h>
#include <Utils.d/GlobalToLocalMap.h>
#ifdef HB_COUPLED_PRECOND
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h>
#endif

#ifdef HB_DOFBASE_8
#define DOFBASE 8
#define HELMDOF 7
#else
#define DOFBASE 7
#define HELMDOF 6
#endif

template<class Scalar>
class GenFsiSparse {
 private:
   ResizeArray<SubLMPCons<Scalar> *> *fsi;
   int numFSI;
   Scalar* wweight;
#ifdef HB_COUPLED_PRECOND
   int numLocalFSI;
   ResizeArray<SubLMPCons<Scalar> *> *local_fsi;
   GlobalToLocalMap& localDofMap;
#endif

 public:
   GenFsiSparse(ResizeArray<LMPCons *> &global_fsi, int glNumFSI, GlobalToLocalMap& local_map,int initSize=256);
   ~GenFsiSparse();
   void scale(double cscale_factor);
   void split(GlobalToLocalMap& local_map, Scalar *weight);
   Scalar csum(GlobalToLocalMap& local_map, int widof);
   void multAdd(const Scalar *rhs, Scalar *result,
                const GlobalToLocalMap& local_map,
                const GlobalToLocalMap& neighb_map, bool fluidOnly = false) const;
   void multAdd(const Scalar *rhs, Scalar *result,
                const GlobalToLocalMap& local_map, bool fluidOnly = false) const;
   void initialize();
#ifdef HB_COUPLED_PRECOND
   int getNumFsi() { return(numFSI); } 
   ResizeArray<SubLMPCons<Scalar> *>* viewFsi() { return(fsi); }
   int extractLocalFsi(ResizeArray<SubLMPCons<Scalar> *> &local_fsi,GlobalToLocalMap& local_map,int* glToLlNdMap=0, int* maxGlNdId=0);
   Connectivity* makeLocalFsiToNodes(GlobalToLocalMap& local_map, int* glToLlNdMap=0, int* maxGlNdId=0, int initSize=256);
   int addLocalFsiToMatrix(GenSparseMatrix<Scalar>* K, DofSetArray* dsa, int* glToLlNdMap=0, Scalar* kSumWI=0);
   void scaleLocalFsi(double cscale_factor);
   void splitLocalFsi(GlobalToLocalMap& local_map, Scalar *weight);
#endif
};

template<class Scalar>
void
GenFsiSparse<Scalar>::initialize()
{
  numFSI      = 0;
  fsi         = 0;
  wweight     = 0;
#ifdef HB_COUPLED_PRECOND
  numLocalFSI = 0;
  local_fsi   = 0;
#endif  
}

template<class Scalar>
GenFsiSparse<Scalar>::GenFsiSparse(ResizeArray<LMPCons *> &global_fsi, int glNumFSI, GlobalToLocalMap& local_map, int initSize)
#ifdef HB_COUPLED_PRECOND
:localDofMap(local_map)
#endif  
{
  initialize();

  fsi = new ResizeArray<SubLMPCons<Scalar> *>(0,initSize);
  numFSI = 0;
  for(int i=0; i<glNumFSI; ++i) {
    bool isLocal = false;  
    if(local_map[DOFBASE*global_fsi[i]->lmpcnum+HELMDOF] > -1) {  // C mult
      (*fsi)[numFSI] = new SubLMPCons<Scalar>(global_fsi[i]->lmpcnum, global_fsi[i]->template getRhs<Scalar>(), 
                                           global_fsi[i]->template getTerm<Scalar>(0), global_fsi[i]->nterms, 0);
      isLocal = true;
      for(int j=1; j<global_fsi[i]->nterms; ++j) (*fsi)[numFSI]->addterm(global_fsi[i]->template getTerm<Scalar>(j), j);
    }
    else {
      for(int j=0; j<global_fsi[i]->nterms; ++j) {  // C^transpose mult
        if(local_map[DOFBASE*global_fsi[i]->terms[j].nnum + global_fsi[i]->terms[j].dofnum] > -1) {
          if(!isLocal) {
            (*fsi)[numFSI] = new SubLMPCons<Scalar>(global_fsi[i]->lmpcnum, global_fsi[i]->template getRhs<Scalar>(), 
                                                    global_fsi[i]->template getTerm<Scalar>(j), global_fsi[i]->nterms, j);
            isLocal = true;
          }
          else (*fsi)[numFSI]->addterm(global_fsi[i]->template getTerm<Scalar>(j), j);
        }
      }
    }
    if(isLocal) numFSI++; 
  }
/*
  cerr<<" In GenFsiSparse(...): fsi array = "<<endl;
  for(i=0; i<numFSI; ++i) 
    (*fsi)[i]->print();
*/
}

template<class Scalar>
void
GenFsiSparse<Scalar>::scale(double cscale_factor)
{
  for(int i=0; i<numFSI; ++i) 
    for(int j=0; j<(*fsi)[i]->nterms; ++j) {
      (*fsi)[i]->terms[j].coef *= cscale_factor;
    }
/*
  cerr<<" In GenFsiSparse::scale(...): fsi array = "<<endl;
  for(int i=0; i<numFSI; ++i) 
    (*fsi)[i]->print();
*/
#ifdef HB_COUPLED_PRECOND
  if(local_fsi) scaleLocalFsi(cscale_factor);
#endif
}

template<class Scalar>
void
GenFsiSparse<Scalar>::split(GlobalToLocalMap& local_map, Scalar* weight)
{
  wweight = weight; // the splitting is now done "on-the-fly" in GenFsiSparse::multAdd 
  /*for(int i=0; i<numFSI; ++i) {
    int fluid_index = local_map[DOFBASE*(*fsi)[i]->lmpcnum+HELMDOF];
    Scalar fw = (fluid_index > -1) ? weight[fluid_index] : 1.0;
    for(int j=0; j<(*fsi)[i]->nterms; ++j) {
      int struct_index = local_map[DOFBASE*(*fsi)[i]->terms[j].nnum + (*fsi)[i]->terms[j].dofnum];
      Scalar sw = (struct_index > -1) ? weight[struct_index] : 1.0;
      Scalar sw = (structure_index > -1) ? weight[structure_index] : 1.0;
      (*fsi)[i]->terms[j].coef /= (fw*sw);
    }
  }*/
/*
  cerr<<" In GenFsiSparse::split(...): fsi array = "<<endl;
  for(int i=0; i<numFSI; ++i) 
    (*fsi)[i]->print();
*/
//#ifdef HB_COUPLED_PRECOND
//  if(local_fsi) splitLocalFsi(local_map, weight);
//#endif
}

template<class Scalar>
Scalar
GenFsiSparse<Scalar>::csum(GlobalToLocalMap& local_map, int widof)
{
  Scalar ret = 0.0;
  for(int i=0; i<numFSI; ++i) {
    int fluid_index = local_map[7*(*fsi)[i]->lmpcnum+6];
    for(int j=0; j<(*fsi)[i]->nterms; ++j) {
      int structure_index = local_map[7*(*fsi)[i]->terms[j].nnum + (*fsi)[i]->terms[j].dofnum];
      if(fluid_index == widof || structure_index ==widof) ret += (*fsi)[i]->terms[j].coef;
    }
  }
  return ret;
}

template<class Scalar>
GenFsiSparse<Scalar>::~GenFsiSparse()
{
  if(fsi) {
    for(int i=0; i<numFSI; ++i) if((*fsi)[i]) { delete (*fsi)[i]; (*fsi)[i] = 0; }
    delete fsi; fsi = 0;
  }
  numFSI = 0;
#ifdef HB_COUPLED_PRECOND
  if(local_fsi) {
    for(int i=0; i<numLocalFSI; ++i) if((*local_fsi)[i]) { delete (*local_fsi)[i]; (*local_fsi)[i] = 0; }
    delete local_fsi; local_fsi = 0;
  } 
  numLocalFSI = 0;
#endif
}

template<class Scalar>
void
GenFsiSparse<Scalar>::multAdd(const Scalar *rhs, Scalar *result,
                              const GlobalToLocalMap& local_map,
                              const GlobalToLocalMap& neighb_map, bool fluidOnly) const
{
  // fsi is a local array of fluid-structure interaction objects containing all interactions that involve either a 
  //   fluid or structure node in this subdomain
  // fsi[i]->lmpcnum is the global id of the fluid node in fluid-structure interaction i
  // fsi[i]->terms[j]->nnum is the global id of the structure node in term j of fluid-structure interaction i
  // local_map is a map from global dof id to local wet interface dof id of this subdomain 
  // result_map is a map from global dof id to local wet interface dof id of neighbor to whom result will be sent
  int i, j, fluid_index, structure_index;
  for(i=0; i<numFSI; ++i) {
    if((fluid_index = neighb_map[DOFBASE*(*fsi)[i]->lmpcnum+HELMDOF]) > -1) { // neighb has fluid node
      for(j=0; j<(*fsi)[i]->nterms; ++j) {
        if((structure_index = local_map[DOFBASE*(*fsi)[i]->terms[j].nnum + (*fsi)[i]->terms[j].dofnum]) > -1) {
          //result[fluid_index] += (*fsi)[i]->terms[j].coef * rhs[structure_index];
          result[fluid_index] += (*fsi)[i]->terms[j].coef * rhs[structure_index]/wweight[structure_index];
        }
      }
    }
    // else if(((fluid_index = local_map[DOFBASE*(*fsi)[i]->lmpcnum+HELMDOF]) > -1) && !fluidOnly) {
    if(((fluid_index = local_map[DOFBASE*(*fsi)[i]->lmpcnum+HELMDOF]) > -1) && !fluidOnly) {
      for(j=0; j<(*fsi)[i]->nterms; ++j) {
        if((structure_index = neighb_map[DOFBASE*(*fsi)[i]->terms[j].nnum + (*fsi)[i]->terms[j].dofnum]) > -1) {
          //result[structure_index] += (*fsi)[i]->terms[j].coef * rhs[fluid_index];
          result[structure_index] += (*fsi)[i]->terms[j].coef * rhs[fluid_index]/wweight[fluid_index];
        }
      }
    }
  }
}

template<class Scalar>
void
GenFsiSparse<Scalar>::multAdd(const Scalar *rhs, Scalar *result, const GlobalToLocalMap& local_map,bool fluidOnly) const
{
  // fsi is a local array of fluid-structure interaction objects containing all interactions that involve either a
  //   fluid or structure node in this subdomain
  // fsi[i]->lmpcnum is the global id of the fluid node in fluid-structure interaction i
  // fsi[i]->terms[j]->nnum is the global id of the structure node in term j of fluid-structure interaction i
  // local_map is a map from global dof id to local wet interface dof id of this subdomain
  // result_map is a map from global dof id to local wet interface dof id of neighbor to whom result will be sent
  int i, j, fluid_index, structure_index;
  for(i=0; i<numFSI; ++i) {
    if((fluid_index = local_map[DOFBASE*(*fsi)[i]->lmpcnum+HELMDOF]) > -1) { // has fluid node
      for(j=0; j<(*fsi)[i]->nterms; ++j) {
        if((structure_index = local_map[DOFBASE*(*fsi)[i]->terms[j].nnum + (*fsi)[i]->terms[j].dofnum]) > -1) {
          result[fluid_index] += (*fsi)[i]->terms[j].coef * rhs[structure_index]/(wweight[structure_index]*wweight[fluid_index]);
        }
      }
    }
    if(((fluid_index = local_map[DOFBASE*(*fsi)[i]->lmpcnum+HELMDOF]) > -1) && !fluidOnly) {
      for(j=0; j<(*fsi)[i]->nterms; ++j) {
        if((structure_index = local_map[DOFBASE*(*fsi)[i]->terms[j].nnum + (*fsi)[i]->terms[j].dofnum]) > -1) {
          result[structure_index] += (*fsi)[i]->terms[j].coef * rhs[fluid_index]/(wweight[structure_index]*wweight[fluid_index]);
        }
      }
    }
  }
}
                                                                                                                    
                                                                                                                    

#ifdef HB_COUPLED_PRECOND
//HB: extract the "local" fsi from all the sub fsi. A local fsi is defined as a fsi linking at least a pressure dof
//    and a displacement dofs of this sub ~> "local" fsi can only exist for "mixed subdomain"
//    The globalToLocal map "local_map" is dof based
template<class Scalar>
int
GenFsiSparse<Scalar>::extractLocalFsi(ResizeArray<SubLMPCons<Scalar> *> &localFsi, GlobalToLocalMap& local_map, int* glToLlNdMap, int* maxGlNdId)
{
  int fluid_index, structure_index;
  int numLocalFsi = 0;
  for(int i=0; i<numFSI; i++) {
    bool isLocal = false;
    if((fluid_index = local_map[DOFBASE*(*fsi)[i]->lmpcnum+HELMDOF]) > -1) { // pressure dof is "local"
#ifdef HB_COUPLED_PRECOND_CHECK
      if(maxGlNdId!=0 & (*fsi)[i]->lmpcnum>*maxGlNdId)
        fprintf(stderr," *** ERROR: GenFsiSparse::extractLocalFsi: pNode = %d > %d\n",(*local_fsi)[i]->lmpcnum,*maxGlNdId);
      if(glToLlNdMap!=0 & glToLlNdMap[(*fsi)[i]->lmpcnum]<0)
        fprintf(stderr," *** PROBLEM in GenFsiSparse::extractLocalFsi, pNode, glToLlNdMap[%d] = %d\n",
                (*fsi)[i]->lmpcnum,glToLlNdMap[(*fsi)[i]->lmpcnum]);
#endif
      int numLocTerms = 0;
      for(int j=0; j<(*fsi)[i]->nterms; ++j) { // loop over displacement dof
        if((structure_index = local_map[DOFBASE*(*fsi)[i]->terms[j].nnum + (*fsi)[i]->terms[j].dofnum]) > -1) { // displacement dof is "local"
#ifdef HB_COUPLED_PRECOND_CHECK
          if(maxGlNdId!=0 & (*fsi)[i]->terms[j].nnum>*maxGlNdId)
            fprintf(stderr," *** ERROR: GenFsiSparse::extractLocalFsi: uNode = %d > %d\n",(*local_fsi)[i]->terms[j].nnum,*maxGlNdId);
          if(glToLlNdMap!=0 & glToLlNdMap[(*fsi)[i]->terms[j].nnum]<0)
            fprintf(stderr," *** PROBLEM in GenFsiSparse::extractLocalFsi, uNode, glToLlNdMap[%d] = %d\n",
                    (*fsi)[i]->terms[j].nnum,glToLlNdMap[(*fsi)[i]->terms[j].nnum]);
#endif
          if(!isLocal) {
            localFsi[numLocalFsi] = new SubLMPCons<Scalar>((*fsi)[i]->lmpcnum, (*fsi)[i]->rhs,
                                                            (*fsi)[i]->terms[j], (*fsi)[i]->nterms, j);
            isLocal = true;
          }
          else localFsi[numLocalFsi]->addterm((*fsi)[i]->terms[j], j);
        }
      }
    }
    if(isLocal) numLocalFsi++;
  }
  fprintf(stderr," ... GenFsiSparse::extractLocalFsi: numLocalFsi = %d (%f %%)\n",numLocalFsi,100*numLocalFsi/double(numFSI));
  //cerr<<" ... GenFsiSparse::extractLocalFsi: localFsi array = "<<endl;
  //for(int i=0; i<numLocalFsi; ++i)
  //  localFsi[i]->print();

  return(numLocalFsi);
}

//HB: The globalToLocal map "local_map" is dof based
//    The globalToLocal map "glToLlNdMap" is node based
template<class Scalar>
Connectivity*
GenFsiSparse<Scalar>::makeLocalFsiToNodes(GlobalToLocalMap& local_map, int* glToLlNdMap, int* maxGlNdId, int initSize)
{
  // Step 1. Extract "local" fsi
  if(!local_fsi) {
    local_fsi   = new ResizeArray<SubLMPCons<Scalar> *> (0,initSize);
    numLocalFSI = extractLocalFsi(*local_fsi,local_map,glToLlNdMap,maxGlNdId);
  }
  // Step 1.1 Create a flag array to flag for each loca_fsi the "unique" fsi nodes
  // ~> to avoid re-do again the "search loop" when filling the target array 
  //    (at the expense of allocating a temporary array & overhead to fill it...)
  int* ptr0  = new int[numLocalFSI+1]; // ptr0[i] points to the starting of the flags for local_fsi i
  ptr0[0] = 0;
  int count = 0;
  for(int i=0; i<numLocalFSI; i++) {
    ptr0[i+1] = ptr0[i]+(*local_fsi)[i]->nterms;
    count += (*local_fsi)[i]->nterms;
  }
  bool* isUnique = new bool[count];
  for(int i=0; i<count; i++) isUnique[i] = false;

  // Step 2. Count number of "targets"
  int* ptr = new int[numLocalFSI+1];  
  ptr[0] = 0;
  int ntargets = 0; 
  for(int i=0; i<numLocalFSI; i++) {
    ntargets++; // fluid node
    ptr[i+1] = ptr[i]+1;
    for(int j=0; j<(*local_fsi)[i]->nterms; j++) { // loop over structure's nodes
      if((*local_fsi)[i]->terms[j].nnum==(*local_fsi)[i]->lmpcnum) continue; // case of same fluid & structure node
      bool found = false;
      for(int jj=0; jj<j; jj++) // look if node already found among structure nodes
        if((*local_fsi)[i]->terms[j].nnum==(*local_fsi)[i]->terms[jj].nnum) { found = true; break; }
      if(!found) { ntargets++; ptr[i+1]++; isUnique[ptr0[i]+j] = true; }
    }
  }
  // Step 3. Fill target array
  int* target = new int[ntargets];
  count = 0;
  for(int i=0; i<numLocalFSI; i++) {
#ifdef HB_COUPLED_PRECOND_CHECK
    if(maxGlNdId!=0 & (*local_fsi)[i]->lmpcnum>*maxGlNdId)
      fprintf(stderr," *** ERROR: GenFsiSparse::makeLocalFsiToNodes: pNode = %d > %d\n",(*local_fsi)[i]->lmpcnum,*maxGlNdId);
#endif
    target[count] = (glToLlNdMap) ? glToLlNdMap[(*local_fsi)[i]->lmpcnum] : (*local_fsi)[i]->lmpcnum; // pressure node 
#ifdef HB_COUPLED_PRECOND_CHECK
    if(target[count]<0) fprintf(stderr," *** ERROR: GenFsiSparse::makeLocalFsiToNodes: target[%d] = %d\n",count,target[count]);
#endif 
    count++;     
    for(int j=0; j<(*local_fsi)[i]->nterms; j++) { // loop over structure's nodes
      if((*local_fsi)[i]->terms[j].nnum==(*local_fsi)[i]->lmpcnum) continue; // case of same fluid & structure node
      if(!isUnique[ptr0[i]+j]) continue; // skip "dupplicate" fsi nodes
#ifdef HB_COUPLED_PRECOND_CHECK
      if(maxGlNdId!=0 & (*local_fsi)[i]->terms[j].nnum>*maxGlNdId)
        fprintf(stderr," *** ERROR: GenFsiSparse::makeLocalFsiToNodes: uNode = %d > %d\n",(*local_fsi)[i]->terms[j].nnum,*maxGlNdId);
#endif
      target[count] = (glToLlNdMap) ? glToLlNdMap[(*local_fsi)[i]->terms[j].nnum] : (*local_fsi)[i]->terms[j].nnum; 
#ifdef HB_COUPLED_PRECOND_CHECK
      if(target[count]<0) fprintf(stderr," *** ERROR: GenFsiSparse::makeLocalFsiToNodes: target[%d] = %d\n",count,target[count]);
#endif  
      count++; 
    }
  }
  if(count!=ntargets) fprintf(stderr," *** ERROR: GenFsiSparse::makeLocalFsiToNodes: count (%d) != ntargets (%d)\n",count,ntargets);
  Connectivity* localFsiToNodes = new Connectivity(numLocalFSI,ptr,target);
  if(isUnique) { delete [] isUnique; isUnique = 0; }
  if(ptr0)     { delete [] ptr0;     ptr0     = 0; }

#define HB_PAIR_FSI_TO_NODE
#ifdef HB_PAIR_FSI_TO_NODE 
  // Try to generate a localFsiToNodes connectivity which represents the connecttion between ONE pressure and ONE structure node 
  // at a time: localFsiToNodes will then have the size be the number of pairs (pressure node,structure node) representing the
  // local_fsi connections. For instance, if there are n loca_fsi with m terms each (i.e. pressure nodes), there will be n.m of 
  // such pairs.
  // This may create a "better" nodeToNode connectivity to build Kii as only small dense blocks will be created for each pairs
  // and not create some additional connections between the structure's nodes
  // It should be possible to make this connectivity directly instead of first making the "standard" localFsiToNodes connectivity
  // and then post-process it as it is actually implemented.
  int npairs= 0;
  ntargets = 0;
  
  // Count number of pairs & targets
  for(int i=0; i<numLocalFSI; i++) {
    if(localFsiToNodes->num(i)==1) { // case of a unique fluid/structure node
      ntargets++;
      npairs++;
    } else {  
      npairs  += localFsiToNodes->num(i)-1;
      ntargets+= 2*(localFsiToNodes->num(i)-1);
    }
  }

  // Fill ptr & target array
  int* ptr1 = new int[npairs+1];
  int* target1 = new int[ntargets];
  ptr1[0] = 0; 
  int count1=0;
  count    = 0;
  for(int i=0; i<numLocalFSI; i++) {
    if(localFsiToNodes->num(i)==1) { // case of a unique fluid/structure node
      ptr1[count1+1] = ptr1[count1]+1;
      target1[count++] = (*localFsiToNodes)[i][0];
      count1++;
    } else {
      for(int j=1; j<localFsiToNodes->num(i); j++) { // Use the fact that the fluid node is (*localFsiToNodes)[i][0]
        ptr1[count1+1] = ptr1[count1]+2;
        target1[count++] = (*localFsiToNodes)[i][0];
        target1[count++] = (*localFsiToNodes)[i][j];
        count1++;
      }
    }
  }  
  if(count!=ntargets) fprintf(stderr," *** ERROR: GenFsiSparse::makeLocalFsiToNodes: count (%d) != ntargets (%d)\n",count,ntargets);
  if(count1!=npairs) fprintf(stderr," *** ERROR: GenFsiSparse::makeLocalFsiToNodes: count1 (%d) != npairs(%d)\n",count1,npairs);
  if(localFsiToNodes) delete localFsiToNodes; 

  Connectivity* localFsiPairToNodes = new Connectivity(npairs,ptr1,target1);
  return(localFsiPairToNodes);
#endif
  //localFsiToNodes->print();
  return(localFsiToNodes); 
}

//HB: the addition of kSumWI to the diagonal terms should better be done somewhere else using directly the map wiInternalMap
//(see BaseSub::makeIMaps) The only difficulty" is then to write a add method into the matrix classes that directly deals with
//given unconstrained numbering, i.e. wiInternalMap[i] is "RCN" (RowColunmNumber) of the wetdof i=0,...,numWIdof-1 
//in the matrix Kii
template<class Scalar>
int
GenFsiSparse<Scalar>::addLocalFsiToMatrix(GenSparseMatrix<Scalar>* K, DofSetArray* dsa, int* glToLlNdMap, Scalar* kSumWI)
{
  int nadded = 0;
  double nrm = 0.0;
  double nrm0= 0.0;
  int pIndex, uIndex;
  bool* done = (kSumWI) ? (bool*)dbg_alloca(localDofMap.numKeys()*sizeof(bool)) : 0;//to avoid adding kSumWI multiple times to the same dof
  if(done) for(int i=0; i<localDofMap.numKeys(); i++) done[i] = false;

  for(int i=0; i<numLocalFSI; i++) {
    int pNode = (glToLlNdMap) ? glToLlNdMap[(*local_fsi)[i]->lmpcnum] : (*local_fsi)[i]->lmpcnum; // pressure node
    int pDof  = dsa->locate(pNode,DofSet::Helm);
    if(pDof>=0) {
      if(kSumWI)
         if((pIndex = localDofMap[DOFBASE*(*local_fsi)[i]->lmpcnum+HELMDOF])>-1 & !done[pIndex]) {
           K->add(pDof,pDof,kSumWI[pIndex]); 
           done[pIndex] = true;
           nrm0 += ScalarTypes::sqNorm(kSumWI[pIndex]); 
         }
      for(int j=0; j<(*local_fsi)[i]->nterms; j++) { // loop over structure's nodes
        int uNode = (glToLlNdMap) ? glToLlNdMap[(*local_fsi)[i]->terms[j].nnum] : (*local_fsi)[i]->terms[j].nnum;
        int uDof  = dsa->locate(uNode,1<<(*local_fsi)[i]->terms[j].dofnum);
        if(uDof>=0) {
          if(kSumWI)
            if((uIndex=localDofMap[DOFBASE*(*local_fsi)[i]->terms[j].nnum+(*local_fsi)[i]->terms[j].dofnum])>-1 & !done[uIndex]){ 
              K->add(uDof,uDof,kSumWI[uIndex]);
              done[uIndex] = true;
              nrm0 += ScalarTypes::sqNorm(kSumWI[uIndex]); 
            }
          if(uDof>pDof) // skyline only store upper triangular par      
            K->add(pDof,uDof,(*local_fsi)[i]->terms[j].coef);
          else
            K->add(uDof,pDof,(*local_fsi)[i]->terms[j].coef); //for sparse matrix, this calls ::addone that adds both to Kij and Kij 
          nrm += ScalarTypes::sqNorm((*local_fsi)[i]->terms[j].coef);
          nadded += 2;
        }
#ifdef HB_COUPLED_PRECOND_CHECK
        else
          fprintf(stderr," *** WARNING: GenFsiSparse::addLocalFsiToMatrix: uNode = %6d, uDof = %6d, gNode = %6d, dofnum = %d\n",uNode,uDof,(*local_fsi)[i]->terms[j].nnum,(*local_fsi)[i]->terms[j].dofnum);
#endif
      }
    }
#ifdef HB_COUPLED_PRECOND_CHECK
    else 
       fprintf(stderr," *** WARNING: GenFsiSparse::addLocalFsiToMatrix: pNode = %6d, pDof = %d, gNode = %6d\n",pNode,pDof,(*local_fsi)[i]->lmpcnum);
#endif
  }
#ifdef HB_COUPLED_PRECOND_CHECK
  fprintf(stderr," -> in GenFsiSparse::addLocalFsiToMatrix: added %d terms into the given matrix, norm = %f\n",nadded,sqrt(nrm));
  if(kSumWI) fprintf(stderr," -> in GenFsiSparse::addLocalFsiToMatrix: added diagonal terms into given matrix, norm = %f\n",sqrt(nrm0)); 
#endif
  return(nadded);
}

template<class Scalar>
void
GenFsiSparse<Scalar>::scaleLocalFsi(double cscale_factor)
{
 if(local_fsi)
   for(int i=0; i<numLocalFSI; i++)
     for(int j=0; j<(*local_fsi)[i]->nterms; ++j)
       (*local_fsi)[i]->terms[j].coef *= cscale_factor;
}
                                                                                                                                
template<class Scalar>
void
GenFsiSparse<Scalar>::splitLocalFsi(GlobalToLocalMap& local_map, Scalar* weight)
{
 if(local_fsi)
   for(int i=0; i<numLocalFSI; i++) {
     int fluid_index = local_map[DOFBASE*(*local_fsi)[i]->lmpcnum+HELMDOF];
     Scalar fw = (fluid_index > -1) ? weight[fluid_index] : 1.0;
     for(int j=0; j<(*local_fsi)[i]->nterms; ++j) {
       int struct_index = local_map[DOFBASE*(*local_fsi)[i]->terms[j].nnum + (*local_fsi)[i]->terms[j].dofnum];
       Scalar sw = (struct_index > -1) ? weight[struct_index] : 1.0;
       (*local_fsi)[i]->terms[j].coef /= (fw*sw);
       //(*local_fsi)[i]->terms[j].coef /= fw; // just for test ...
//#ifdef HB_COUPLED_PRECOND_CHECK
//       if(ScalarTypes::norm(1./(fw*sw))>1.)  
//         cerr<<" *** PROBLEM in GenFsiSparse:splitLocalFsi: |1/fw.sw| = "<<ScalarTypes::norm(1./(fw*sw))<<", pNd = "<<(*local_fsi)[i]->lmpcnum<<", |fw| = "<<ScalarTypes::norm(fw)<<", |sw| = "<<ScalarTypes::norm(sw)<<endl;
//#endif
     }
   }
}

#endif
#endif

