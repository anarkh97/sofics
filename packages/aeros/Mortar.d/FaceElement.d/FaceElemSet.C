// ----------------------------------------------------------------
// HB - 06/25/03
// ----------------------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <unistd.h>

// STL
#include <map>

// FEM headers
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>

#ifdef SOWER_SURFS
#include <Utils.d/BinFileHandler.h>
#endif

// Last returns the true last defined element
int 
FaceElemSet::last() const
{
/*
 int last = size();
 while(--last >= 0)
    if(elem[last] != 0) break;

 return last+1;
*/
 return _last;
}

typedef FaceElement* ElemP;

FaceElemSet::FaceElemSet(int initsize) : ba(4096, initsize)
{
  emax = initsize;
  if(initsize > 0) {
    elem = new ElemP[initsize];
    int i;
    for(i = 0; i < emax; ++i)
      elem[i] = 0;
  }
  else elem = 0;
  _last = 0;
}

void
FaceElemSet::elemadd(int num, FaceElement *el)
{
  if(num >= emax) // resize elem[]
   {
    int newsize = ((num+1)*3)/2;
    ElemP *np = new ElemP[newsize];
    int i;
    for(i= 0; i < emax; ++i)
     np[i] = elem[i];
    for(emax = newsize; i < emax; ++i)
     np[i] = 0;
    delete[] elem;
    elem = np;
   }
  elem[num] = el;
  _last++;
}

int 
FaceElemSet::nElems() { return last(); }

void
FaceElemSet::Renumber(std::map<int,int>& OldToNewNodeIds)
{
  int nEl = last(); 
  for(int iel=0; iel<nEl; iel++)
    elem[iel]->Renumber(OldToNewNodeIds); 
}

void
FaceElemSet::print()
{
   int nEl = last(); 
   //fprintf(stderr," -> face el. set contains %d face els.\n",nEl);
   filePrint(stderr," -> face el. set contains %d face els.\n",nEl);
   //pause();
   //fprintf(stderr," (... wait 2s ...)\n"); sleep(2);
   for(int iel=0; iel<nEl; iel++){
   	elem[iel]->print();
   }
}

#ifdef SOWER_SURFS
void 
FaceElemSet::WriteSower(BinFileHandler& file, int* ndMap)
{
  int nEls = last();
  int nodes[64];
  for(int iel=0; iel<nEls; iel++){
    file.write(&iel,1);
    int elType = elem[iel]->GetFaceElemType();
    int nNds   = elem[iel]->nNodes();
    elem[iel]->GetNodes(nodes,ndMap); 
    file.write(&elType, 1);
    file.write(&nNds, 1);
    file.write(nodes, nNds);
    //elem[iel]->WriteSower(BinFileHandler& file, int* ndMap);
  }
}
#endif

std::map<int,locoord> FaceElemSet::computeNodeLocalCoords(int* fnId, int size) 
{
  std::map<int,locoord> exy;  //Node Id -> (iElem, (x,y))
  std::map<int,locoord>::iterator it;

  for(int iel=0; iel<last(); iel++) {
    double* coords = elem[iel]->ViewRefCoords();
    for(int k=0; k<elem[iel]->nNodes(); k++) {
      int gId = elem[iel]->GetNode(k);
      it = exy.find(gId);
      if(it==exy.end()) //not found before 
        exy[gId] = locoord(iel,std::pair<double,double>(coords[2*k],coords[2*k+1]));
    }
  }

  return exy;
}

void FaceElemSet::remove(int num)
{
  elem[num] = 0;
}

void FaceElemSet::repack()
{
  int nEls = last();
  int count = 0;
  for(int i = 0; i < nEls; ++i) {
    if(elem[i] == 0) count++;
    else if(count > 0) { elem[i-count] = elem[i]; elem[i] = 0; }
  }
  _last = nEls - count;
}
