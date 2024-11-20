#include <iostream>
#include <cstdio>
#include <Utils.d/BinFileHandler.h>
#include "ConnectivityT.h"

template<typename IndexType, typename DataType, typename Map>
ConnectivityT<IndexType,DataType,Map>::ConnectivityT(IndexType _size, IndexType *_pointer, DataType *_target,
                                                 bool _removeable, float *_weight) 
	:
	  pointer{_pointer, _pointer+_size+1},
	  target{_target,_target+_pointer[_size]}, 
	  weight{_weight, _weight+_pointer[_size]}
{
	if(_removeable) {
		delete [] _pointer;
		delete [] _target;
		delete [] _weight;
	}
		
}

template<typename IndexType, typename DataType, typename Map>
ConnectivityT<IndexType,DataType,Map>::ConnectivityT(BinFileHandler& f, bool oldSower)
{
	size_t size;
	f.read(&size,1);
	pointer.resize(size);
	f.read(pointer.data(),size);
	size_t numtarget = pointer.back();
	f.read(&numtarget,1);
	target.resize(numtarget);
	f.read(target.data(),numtarget);
}

template<typename IndexType, typename DataType, typename Map>
size_t ConnectivityT<IndexType,DataType,Map>::write(BinFileHandler& f) const
{
	size_t _size = pointer.size();
	f.write(&_size,1);
	f.write(pointer.data(),_size);
	size_t numtarget = getNumTarget();
	f.write(&numtarget,1);
	f.write(target.data(),numtarget);
	return 0;
}

//template<typename IndexType, typename DataType, typename Map>
//ConnectivityT<IndexType,DataType,Map>::ConnectivityT(IndexType _size, IndexType *_count)
//{
//	size    = _size;
//	pointer.resize(size+1);
//	pointer[0] = 0;
//	IndexType i;
//	for(i=0; i < _size; ++i)
//		pointer[i+1] = pointer[i] + _count[i];
//	target.resize(pointer[size]);
//}
//
//template<typename IndexType, typename DataType, typename Map>
//ConnectivityT<IndexType,DataType,Map>::ConnectivityT(IndexType _size, IndexType count)
//{
//	size    = _size;
//	pointer = new IndexType[size+1];
//	pointer[0] = 0;
//	IndexType i;
//	for(i=0; i < _size; ++i)
//		pointer[i+1] = pointer[i] + count;
//	target.resize(pointer[size]);
//	for(IndexType i=0; i < pointer[size]; ++i)
//		target[i] = i;
//}

/*ConnectivityT::ConnectivityT(BinFileHandler &file)
{
  removeable = true;
  file.read(&size, 1);
  file.read(&numtarget, 1);

  pointer = new IndexType[size+1];
  target = new IndexType[numtarget];
  weight = 0;

  file.read(pointer, size+1);
  file.read(target, numtarget);
}*/

template<typename IndexType, typename DataType, typename Map>
IndexType
ConnectivityT<IndexType,DataType,Map>::offset(IndexType i, DataType j)
{
	for(auto ii = pointer[i]; ii < pointer[i+1]; ++ii)
		if(target[ii] == j) return ii;

	return -1; // We didn't find a connection between i and j
}

template<typename IndexType, typename DataType, typename Map>
IndexType
ConnectivityT<IndexType,DataType,Map>::cOffset(IndexType i, DataType j)
{
	for(auto ii = pointer[i]; ii < pointer[i+1]; ++ii)
		if(target[ii] == j) return (ii - pointer[i]);

	return -1; // We didn't find a connection between i and j
}

#if 0
#include <Utils.d/resize_array.h>
#include <Utils.d/dbg_alloca.h>
#include <algorithm>
using std::stable_sort;

ConnectivityT*
ConnectivityT::transconOne( ConnectivityT* tc)
// D. Rixen :for every pointer
// of tc only one entry in target of tc is considered
// (used to associate a mpc term for a interface node
//  to only one sub, with a preference for already included sub)
{
 IndexType i,j,k;

 // First find the biggest target so we can size arrays correctly
 IndexType tgmax=-1;

 for(i =0; i < tc->numtarget; ++i)
   if(tc->target[i] > tgmax) tgmax=tc->target[i];
 tgmax++; // Important adjustment

 // Now we can size the array that flags if a target has been visited
 IndexType *flags = (IndexType *) dbg_alloca(sizeof(IndexType)*tgmax);

 // For every pointer, build the number of occurrence of targets
 // and choose one target per tc->pointer (max occurrence)
 // At the same time build new pointers np

 IndexType cp = 0;
 IndexType *np = new IndexType[size+1];
 IndexType ii;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   for(ii = 0; ii < tgmax; ++ii)
    flags[ii] = 0;
   //-- store number of occurence in flag
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      IndexType intermed = target[j];
      for (k = 0; k < tc->num(intermed); ++k)
          flags[(*tc)[intermed][k]]++;
   }
   //-- set pointer and size target
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      IndexType intermed = target[j];
      //-- target with max occ.
      IndexType targMaxOcc;
      IndexType maxOcc = 0;
      for (k = 0; k < tc->num(intermed); ++k){
          if(flags[(*tc)[intermed][k]]==-1){
            maxOcc=0;
            break;
          }
          else if(flags[(*tc)[intermed][k]]> maxOcc){
            maxOcc = flags[(*tc)[intermed][k]];
            targMaxOcc = (*tc)[intermed][k];
          }
      }
      if(maxOcc!=0){cp++;flags[targMaxOcc]=-1;}
    }
 }
 np[size] = cp;
 // Now allocate and fill new target
 IndexType *ntg = new IndexType[cp];
 cp = 0;
 for(i = 0; i < size; ++i) {
   for(ii = 0; ii < tgmax; ++ii)
    flags[ii] = 0;
   //-- store number of occurence in flag
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      IndexType intermed = target[j];
      for (k = 0; k < tc->num(intermed); ++k)
          flags[(*tc)[intermed][k]]++;
   }
   //-- fill target
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      IndexType intermed = target[j];
      //-- target with max occ.
      IndexType targMaxOcc;
      IndexType maxOcc = 0;
      for (k = 0; k < tc->num(intermed); ++k){
          if(flags[(*tc)[intermed][k]]==-1){
            maxOcc=0;
            break;
          }
          else if(flags[(*tc)[intermed][k]]> maxOcc){
            maxOcc = flags[(*tc)[intermed][k]];
            targMaxOcc = (*tc)[intermed][k];
          }
      }
      if(maxOcc!=0){
         ntg[cp]=targMaxOcc;cp++;flags[targMaxOcc]=-1;}
    }
 }

 ConnectivityT *res = new ConnectivityT();
 res->size      = size;
 res->pointer   = np;
 res->numtarget = np[size];
 res->target    = ntg;
 res->weight    = 0;
 return res;
}

IndexType
ConnectivityT::num(IndexType nd, IndexType *mask)
{
 IndexType res=0;
 IndexType jstrt = pointer[nd];
 IndexType jstop = pointer[nd+1];
 IndexType j;
 for(j = jstrt; j < jstop; ++j)
    if(mask[target[j]]) res++;
 return res;
}

void
ConnectivityT::findPseudoDiam(IndexType *s, IndexType *e, IndexType *mask)
{
 IndexType i,k,nw;
 // Select the node with the lowest connectivity
 IndexType cmin = numtarget+1;
 IndexType cmax = 0;
 for(i = 0; i < size; ++i) {
  if(mask[i]) {
     if((nw = num(i,mask)) < cmin) {
       cmin=nw;
       *s = i;
     }
    if(nw > cmax) cmax = nw;
  }
 }

 IndexType *ls  =  new IndexType[numtarget];
 IndexType *xls =  new IndexType[numtarget+1];

 // Created rooted level structure
 IndexType w;
 IndexType h = rootLS(*s, xls, ls, w, mask);
 IndexType subconsize = xls[h];
 IndexType *sorted = (IndexType *) dbg_alloca((cmax+1)*sizeof(IndexType));
 *e = ls[xls[h-1]]; // Give a default end point in case h == subconsize.
 while(h < subconsize) {
   for(k=0; k <= cmax; ++k) sorted[k] = -1;
   // Find in the last level the node with the minimum connectivity
   //IndexType maxweight = subconsize;
   IndexType kstrt = xls[h-1];
   IndexType kstop = xls[h];
   for(k=kstrt; k < kstop; ++k)
      {
        sorted[num(ls[k],mask)] = ls[k];
      }
   IndexType w_e = subconsize;
   for(k = 0; k <= cmax; ++k)
    if(sorted[k] >= 0) {
      IndexType nh = rootLS(sorted[k], xls, ls, w, mask);
      if(w < w_e) {
         if(nh > h) {
           *s = sorted[k];
           h = nh;
           break;
         }
         *e = sorted[k];
         w_e = w;
      }
    }
   if(k > cmax) break;
 }
 delete [] xls;
 delete [] ls;
 return;
}

IndexType
ConnectivityT::rootLS(IndexType root, IndexType *xls, IndexType *ls, IndexType &w, IndexType *mask)
{
 IndexType i, j;
 w = 0;

 IndexType *locMask = new IndexType[size];

 if(mask)
   for(i = 0; i < size; ++i) locMask[i] = mask[i];
 else
   for(i = 0; i < size; ++i) locMask[i] = 1;

 locMask[root] = 0;
 ls[0] = root;
 xls[0] = 0;
 IndexType nlvl = 1;
 IndexType nf = 1;
 while(nf > xls[nlvl-1]) {
   xls[nlvl] = nf;
   IndexType lbegin = xls[nlvl-1];
   IndexType lend   = xls[nlvl];
   for (i = lbegin; i <lend; ++i) {
     IndexType n1 = ls[i];
     IndexType jstart = pointer[n1];
     IndexType jstop  = pointer[n1+1];
     for(j=jstart; j < jstop; ++j) {
       IndexType n2 = target[j];
       if(locMask[n2]) {
          locMask[n2] = 0;
          ls[nf++] = n2;
       }
     }
   }
   if(nf-xls[nlvl] > w) w = nf-xls[nlvl];
   nlvl++;
 }

 delete [] locMask;

 return nlvl-1;
}

compStructT
ConnectivityT::renumByComponent(IndexType renumAlg)
{
  // size = total number of nodes
  IndexType *globalRenum = new IndexType[size];
  IndexType *mark = new IndexType[size];
  IndexType *ls   = new IndexType[size];
  IndexType *xls  = new IndexType[size+1];
  ResizeArray<IndexType> xcomp(0,2);
  // Initialize mark to zero, accounting for missing node #s
  // Initialize globalMask

  IndexType inode;
  for(inode = 0; inode < size; ++inode) {
    mark[inode] = (num(inode) != 0) ? 1 : 0;
    globalRenum[inode] = -1;
  }

  // Loop over nodes checking which ones are marked
  // and belong to the same component.
  IndexType j, k, nextNum = 0, count = 0;
  IndexType *locMask = new IndexType[size];
  for(inode = 0; inode < size; ++inode)
     locMask[inode] = 0;

  IndexType currentNumber = 0;
  xcomp[0] = currentNumber;

  IndexType *lrenum = 0;
  if(renumAlg>0) { lrenum = new IndexType[size]; }
 
  for(inode = 0; inode < size; ++inode)
  {
    if(mark[inode] == 1) {
      // Find all neighbors of inode
 
      IndexType w;
      IndexType h = rootLS(inode, xls, ls, w);

      // Declare and set local mask

      for(j=0; j<xls[h]; ++j) {
          locMask[ls[j]] = 1;
          mark[ls[j]] = 0;
      }

      // call renumbering for local mask
      switch(renumAlg) {
           case 1:
              renumSloan(locMask, nextNum, lrenum);
              break;
           case 2:
              renumRCM(locMask, nextNum, lrenum);
              break;
           default:
              break;
      }

      // Assemble local mask into global mask
      if(renumAlg>0) { 
        for(j=0; j<xls[h]; ++j) {
            k = ls[j];
            globalRenum[k] = lrenum[k];
        }
      }
      else {
#ifdef HB_RENUMBYCOMP_OPT
        if(renumAlg==-1){      
          //HB: an optimization would be to loop over ls (0...xls[h])
          // but to recover the same answer as below, one has to sort
          // the ls array. so this may be not worthy ...
          // anyway, do we really need to sort the ls array ??? Not sorting
          // it will give us an output with a different ordering but we can 
          // probably take care of it by providing an associated globalRenum array ???  
          std::sort(ls,ls+xls[h]);
          for(j=0; j<xls[h]; ++j) {
            k = ls[j];
            globalRenum[k] = currentNumber + j;
          }
        }else{
#endif
          IndexType cnt = 0;
          for(j=0; j<size; ++j)
            if(locMask[j])
              globalRenum[j] = currentNumber+cnt++;
#ifdef HB_RENUMBYCOMP_OPT
        } 
#endif
      }
      currentNumber += xls[h];
      count += 1;
      xcomp[count] = currentNumber;
      // reset locMask to zero
      for(j=0; j<xls[h]; ++j)
        locMask[ls[j]] = 0;
    }
  }

  compStructT ret;
  ret.numComp = count;
  ret.xcomp   = xcomp.yield(); 
  ret.renum   = globalRenum;
  delete [] mark; delete [] ls; delete [] xls;
  delete [] locMask;
  if(lrenum) { delete [] lrenum; }
  return ret;
}

IndexType *
ConnectivityT::renumRCM(IndexType *mask, IndexType &nextNum, IndexType *renum)
{
 IndexType i,j,k;
 if(mask == 0)
  {
    mask = (IndexType *) dbg_alloca(sizeof(IndexType)*size);
    for(i=0; i < size; ++i)
       mask[i] = (num(i)) ? 1 : 0;
  }

 IndexType s_node, e_node;

 findPseudoDiam( &s_node, &e_node, mask);

 if(renum == 0) renum = new IndexType[size];
 IndexType *order  = (IndexType *) dbg_alloca(sizeof(IndexType)*size);
 IndexType *degree = (IndexType *) dbg_alloca(sizeof(IndexType)*size);

 // get the degree of all the nodes
 for(i =0; i < size; ++i)
   degree[i] = num(i,mask);

 order[0] =e_node;
 mask[e_node] = -mask[e_node]; // mark we have seen this node
 IndexType lastNode=1; // number of nodes which have been assigned a number
 for(i = 0; i < lastNode; ++i) {
   IndexType curNode = order[i];
   IndexType firstNeighb = lastNode;
   // Look at the neighbors of this node and add that to the list
   for(j = pointer[curNode]; j < pointer[curNode+1]; ++j) {
     IndexType neighbNode = target[j];
     if(mask[neighbNode] > 0) {
       order[lastNode] = neighbNode;
       mask[neighbNode] = -mask[neighbNode];
       lastNode += 1;
     }
   }
   // now sort the added nodes by degree
   if(firstNeighb < lastNode-1) {
     for(j = firstNeighb; j < lastNode-1; ++j)
      for(k = j+1; k < lastNode; ++k) 
       if(degree[order[k]] < degree[order[j]]) {
         IndexType tmp = order[k];
         order[k] = order[j];
         order[k] = tmp;
       }
   }
 }
 // now reverse the order
 for(i = lastNode; i--; ) {
   renum[order[i]] = (nextNum++);
 }
 return renum;
}


IndexType *
ConnectivityT::renumSloan(IndexType *mask, IndexType &nextNum, IndexType *renum)
{
 IndexType i,j,k;
 IndexType s_node, e_node;
 float w1=1.0, w2=2.0;

 if(mask == 0)
  {
    mask = (IndexType *) dbg_alloca(sizeof(IndexType)*size);
    for(i=0; i < size; ++i)
       mask[i] = (num(i)) ? 1 : 0;
  }

 findPseudoDiam( &s_node, &e_node, mask);

 IndexType *ls  = new IndexType[numtarget];
 IndexType *xls = new IndexType[numtarget+1];

 IndexType w;
 IndexType h = rootLS(e_node, xls,ls,w,mask);
 // now give a distance to each point

 if(renum == 0) renum = new IndexType[size];
 IndexType *distance = renum;
 IndexType *status = distance;

 for(i = 0; i < size; ++i)
   distance[i] = -1;

 for(i = 0; i < h; ++i) {
    for(j = xls[i]; j < xls[i+1]; ++j)
      distance[ls[j]] = i;
  }
 delete [] xls;
 delete [] ls;

 float *priority = (float *) dbg_alloca(sizeof(float)*size);
 // initialize the priority values
 for(i = 0; i < size; ++i) {
   if(mask[i]) {
     priority[i] = w1*distance[i] - w2*num(i);
     distance[i] = -2; // status and distance are stored in the same array
   }
 }

 // initalize the queue with the starting point
 IndexType *q = (IndexType *) dbg_alloca(sizeof(IndexType)*size);
       // maximum size is all the points in the q
 q[0] = s_node;
 IndexType nn=1;
 status[s_node] = -1;

 // While the queue is not empty
 while(nn > 0) {
   // Find in the queue the point with maximum priority
   IndexType next = 0;
   float maxpri = priority[q[next]];
   for(i = 1; i < nn; ++i) {
      if(priority[q[i]] > maxpri) {
         maxpri = priority[q[i]];
         next = i;
      }
   }
   // Remove the next node numbered from the queue
   IndexType nextnode = q[next];
   q[next] = q[nn-1];
   nn = nn-1;
   IndexType istart = pointer[nextnode], istop = pointer[nextnode+1];
   if(status[nextnode] == -1) {
     status[nextnode] = 0; // So we don't step on our feet.
     // Preactive node. Examine its neighbors
     for(i = istart; i < istop; ++i) {
        IndexType neighbor = target[i];
        priority[neighbor] += w2; // update the priority
        if(status[neighbor] == -2) { // if neighbor was inactive, it becomes pre-active
           // NB the next loop will set the same nodes as active
           q[nn] = neighbor;
           nn++;
           status[neighbor] = -1;
        }
     }
   }
   status[nextnode] = nextNum;
   nextNum++;
   // Now scan the preactive neighbors, make them active and their neighbors become
   // preactive if they were inactive
   for(i = istart; i < istop; ++i) {
      IndexType neighbor = target[i];
      if(status[neighbor] == -1) {
        status[neighbor] = 0;
        priority[neighbor] += w2;
        IndexType kstart = pointer[neighbor], kstop = pointer[neighbor+1];
        for(k = kstart; k < kstop ; ++k) {
            IndexType kn = target[k];
            priority[kn] += w2;
            if(status[kn] == -2) { // This node is inactive. must become preactive
               status[kn] = -1;
               q[nn] = kn;
               nn++;
            }
        }
      }
    }
 }
 return status;
}

void
ConnectivityT::print(FILE *f, IndexType node)
{
 if(node == -1) {
   IndexType i;
   IndexType maxdist = 0;
   fprintf(f, "Size %d\n",size);
   fflush(f);
   for(i = 0; i < size; ++i) {
     fprintf(f, "%d ->", i+1);
     IndexType j;
     for(j = pointer[i]; j < pointer[i+1]; ++j) {
       fprintf(f, " %d", target[j]+1);
       if(i-target[j] > maxdist ) maxdist = i-target[j];
     }
     fprintf(f,"\n");
   }
   fprintf(f, "Max dist %d\n", maxdist);
   fflush(f);
 } else {
   fprintf(f,"%d ->",node+1);
   IndexType j;
   for(j=pointer[node]; j<pointer[node+1]; ++j) {
     fprintf(f," %d", target[j]+1);
   }
   fprintf(f,"\n");
 }
}

IndexType
ConnectivityT::findMaxDist(IndexType *renum)
{
 IndexType count = numNonZeroP();
 IndexType i;
 IndexType maxdist = 0;
 IndexType localDist;
 double avgDist=0;
 for(i = 0; i < size; ++i) {
   localDist = 0;
   IndexType j;
   for(j = pointer[i]; j < pointer[i+1]; ++j) {
     if(renum[i]-renum[target[j]] > maxdist )
	 maxdist = renum[i]-renum[target[j]];
     if(renum[i]-renum[target[j]] > localDist)
         localDist = renum[i]-renum[target[j]];
   }
   avgDist += localDist;
 }
 fprintf(stderr, "Total   distance: %e\n", avgDist);
 fprintf(stderr, "Average distance: %e\n", avgDist/count);
 fprintf(stderr, "Maximum distance: %d\n", maxdist);
 return maxdist;
}

ConnectivityT*
ConnectivityT::merge(ConnectivityT *con2)
{
 IndexType size1 = csize();
 IndexType size2 = con2->csize();

 IndexType i;
 IndexType *cp = new IndexType[size1 + size2+1];
 IndexType fp = 0;
 for(i = 0; i < size1; ++i) {
   cp[i] = fp;
   fp += num(i);
 }
 for(i = 0; i < size2; ++i) {
   cp[i+size1] = fp;
   fp += con2->num(i);
 }
 cp[size1 + size2] = fp;

 IndexType *ct = new IndexType[fp];
 fp = 0;
 for(i=0; i<size1; ++i) {
   IndexType j;
   for(j =0; j < num(i); ++j)
     ct[fp++] = (*this)[i][j];
 }
 for(i = 0; i < size2; ++i) {
   IndexType j;
   for(j =0; j < con2->num(i); ++j)
     ct[fp++] = (*con2)[i][j];
 }

 return new ConnectivityT(size1+size2, cp, ct);
}

void
ConnectivityT::combine(ConnectivityT *con2, IndexType *&cmap, IndexType *cmap2)
{
 // add con2 to this and make combined cmap
 IndexType i, j;
 IndexType size1 = csize();
 IndexType size2 = con2->csize();
 IndexType *tmp1 = (IndexType *) dbg_alloca(sizeof(IndexType)*size1);
 for(i=0; i<size1; ++i) tmp1[i] = -1;
 IndexType *tmp2 = (IndexType *) dbg_alloca(sizeof(IndexType)*size2);
 for(i=0; i<size2; ++i) tmp2[i] = -1;
  
 IndexType *cp = new IndexType[size1 + size2+1];
 IndexType fp = 0;
 for(i = 0; i < size1; ++i) {
   cp[i] = fp;
   fp += num(i);
   for(j = 0; j < size2; j++) {
     if(cmap2[j] == cmap[i]) {
       tmp2[j] = i;
       tmp1[i] = j;
       fp += con2->num(j);
       break;
     }
   }
 }
 IndexType count = 0;
 for(i = 0; i < size2; ++i) {
   if(tmp2[i] == -1) {
     cp[count+size1] = fp;
     fp += con2->num(i);
     count++;
   }
 }
 cp[size1 + count] = fp;
 IndexType *new_cmap = new IndexType[size1 + count];

 IndexType *ct = new IndexType[fp];
 fp = 0;
 for(i=0; i<size1; ++i) {
   new_cmap[i] = cmap[i]; 
   for(j =0; j < num(i); ++j) 
     ct[fp++] = (*this)[i][j];
   if(tmp1[i] != -1) {
     for(j =0; j < con2->num(tmp1[i]); ++j) 
       ct[fp++] = -1 - (*con2)[tmp1[i]][j];  // merged are -ve
   }
 }
 count = 0;
 for(i = 0; i < size2; ++i) {
   if(tmp2[i] == -1) {
     new_cmap[size1+count] = cmap2[i];
     count++;
     for(j =0; j < con2->num(i); ++j) 
       ct[fp++] = -1 - (*con2)[i][j]; // merged are -ve
   }
 }
 delete [] cmap;
 cmap = new_cmap;
 
 if(removeable) {
   delete [] pointer;
   delete [] target;
 }
 size      = size1+count;
 pointer   = cp;
 target    = ct;
 numtarget = cp[size1+count];
 removeable = true;
}

IndexType
ConnectivityT::numNonZeroP()
{
 IndexType count = 0;
 IndexType i;
 for(i=0; i < size; ++i)
   if(pointer[i+1] != pointer[i]) ++count;
 return count;
}

// Collapse creates a connecivity that represents the grouping
// of a reflexive connectivity into subconnected groups 
ConnectivityT *
ConnectivityT::collapse() {
  IndexType i,j;
  IndexType count = 0;
  bool *flag = new bool[size];
  IndexType grIndex;
  IndexType *group = new IndexType[size];
  for(i=0; i < size; ++i)
    flag[i] = false;
  grIndex = 0;
  for(i = 0; i < size; ++i)
    if(flag[i] == false) {
      count++;
      IndexType follow = grIndex;
      flag[i] = true;
      group[grIndex++] = i;
      while(follow < grIndex) {
	IndexType s = group[follow];
	for(j = 0; j < num(s); ++j)
	  if((*this)[s][j] >= 0 && flag[ (*this)[s][j] ] == false) {
	    flag[ (*this)[s][j] ] = true;
	    group[grIndex++] = (*this)[s][j];
	  }
        follow++;
      }
    }
  IndexType *ptr = new IndexType[count+1];
  count = 0;
  grIndex = 0;
  for(i=0; i < size; ++i)
    flag[i] = false;
  for(i = 0; i < size; ++i)
    if(flag[i] == false) {
      ptr[count++] = grIndex;
      IndexType follow = grIndex;
      flag[i] = true;
      group[grIndex++] = i;
      while(follow < grIndex) {
	IndexType s = group[follow];
	for(j = 0; j < num(s); ++j)
	  if((*this)[s][j] >= 0 && (*this)[s][j] <size && flag[ (*this)[s][j] ] == false) {
	    flag[ (*this)[s][j] ] = true;
	    group[grIndex++] = (*this)[s][j];
	  }
          else if((*this)[s][j] >= size)
           fprintf(stderr, "error connectivity is not reflexive\n");
        follow++;
      }
    }
 ptr[count] = grIndex;
 delete [] flag;
 return new ConnectivityT(count, ptr, group);      
}

ConnectivityT *
ConnectivityT::copy()
{
  if(size==0) return new ConnectivityT();
  IndexType *nptr = new IndexType[size+1];
  IndexType *ntrg = new IndexType[numConnect()];
  IndexType i;
  for(i = 0; i <= size; ++i)
    nptr[i] = pointer[i];
  for(i = 0; i < numConnect(); ++i)
    ntrg[i] = target[i];
  return new ConnectivityT(size,nptr,ntrg);
}

// Trim a connectivity to remove targets that have only one corresponding
// target in the marker connectivity
ConnectivityT *
ConnectivityT::trim(ConnectivityT *marker)
{
  IndexType i,j;
  IndexType count;
  count = 0;
  for(i = 0; i < numConnect(); ++i)
    if(marker->num( target[i] ) > 1)
      count++;
  IndexType *ntrg = new IndexType[count];
  IndexType *nptr = new IndexType[size+1];
  count = 0;
  for(i = 0; i < size; ++i) {
    nptr[i] = count;
    for(j = 0; j < num(i); ++j)
      if( marker->num( (*this)[i][j] ) > 1)
	ntrg[count++] = (*this)[i][j];
  }
  nptr[size] = count;
  return new ConnectivityT(size, nptr, ntrg);  
}

// PJSA: modify a connectivity to insert connection between each target and itself 
// needed to make bodyToBody, since some bodys may have no MPCs
ConnectivityT *
ConnectivityT::modify()
{
  IndexType i,j;
  IndexType count;
  count = 0;
  for(i = 0; i < size; ++i)
    if(num(i) == 0) count++;
  IndexType *ntrg = new IndexType[numtarget+count];
  IndexType *nptr = new IndexType[size+1];
  count = 0;
  for(i = 0; i < size; ++i) {
    nptr[i] = count;
    if(num(i) == 0) ntrg[count++] = i;
    else {  
      for(j = 0; j < num(i); ++j)
        ntrg[count++] = (*this)[i][j];
    }
  }
  nptr[size] = count;
  return new ConnectivityT(size, nptr, ntrg);
}

// this one is to remove connection with self
ConnectivityT *
ConnectivityT::modifyAlt()
{
  IndexType i,j;
  IndexType count;
  count = 0;
  for(i = 0; i < size; ++i)
    for(j = 0; j < num(i); ++j)
      if((*this)[i][j] == i) count++;
  IndexType *ntrg = new IndexType[numtarget-count];
  IndexType *nptr = new IndexType[size+1];
  count = 0;
  for(i = 0; i < size; ++i) {
    nptr[i] = count;
    for(j = 0; j < num(i); ++j)
      if((*this)[i][j] != i) ntrg[count++] = (*this)[i][j];
  }
  nptr[size] = count;
  return new ConnectivityT(size, nptr, ntrg);
}

ConnectivityT *
ConnectivityT::combineAll(IndexType addSize, IndexType* cmap)
{
  IndexType i,j,k;
  IndexType count;
  bool foundi, foundj;
  //count = 0;
  //for(i = 0; i < size; ++i)
    //if(num(i) == 0) count++;
  IndexType *ntrg = new IndexType[numtarget+(addSize)*(addSize)];
  IndexType *nptr = new IndexType[size+1];

  count = 0;
  for(i = 0; i < size; ++i) {
    nptr[i] = count;
    //if(num(i) == 0) ntrg[count++] = i;
    //else {  
    for(j = 0; j<num(i); ++j)
      ntrg[count++] = (*this)[i][j];
    //}
    IndexType oldNumI = num(i);
    //fprintf(stderr, " ... i = %d, cmap[addCount] = %d ...\n", i,cmap[addCount]);
    foundi = false;
    for(k = 0; k < addSize; ++k)
      if(cmap[k] == i) { foundi = true; break; }
    if(foundi) {
      for(k = 0; k < addSize; ++k) {
        foundj = false;
        for(j = 0; j < oldNumI; ++j) {
          //fprintf(stderr, " ... i = %d, j = %d, k = %d, cmap[k] = %d, (*this)[i][j] = %d ...\n", i,j,k,cmap[k],(*this)[i][j]);
          if (cmap[k] == (*this)[i][j]) {
            foundj = true;
            //fprintf(stderr, " ... found!!! : i = %d, j = %d, k = %d ...\n", i,j,k);
            break;
          }
        }
        if (!foundj) {
            //fprintf(stderr, " ... NOT found!!! : i = %d, j = %d, k = %d ...\n", i,j,k);
          ntrg[count++] = cmap[k];
        }
      }
    }
  }
  nptr[size] = count;
  return new ConnectivityT(size, nptr, ntrg);
}

void
ConnectivityT::sortTargets()
{
  for(IndexType i = 0; i < size; ++i) 
    stable_sort(target+pointer[i], target+pointer[i+1]);
}

void ConnectivityT::renumberTargets(IndexType *map)  {

  for (IndexType i = 0; i < numtarget; i++)  {
    if (map[target[i]] < 0)
      fprintf(stderr, "target(i) is neg: %d\n", target[i]);
    target[i] = map[target[i]];
  }
}

void ConnectivityT::renumberTargets(map<IndexType, IndexType> &map)  {
 
  for (IndexType i = 0; i < numtarget; i++)  {
    if (map.find(target[i]) == map.end())
      fprintf(stderr, "target(i) does not exist: %d\n", target[i]);
    target[i] = map[target[i]];
  }
}

ConnectivityT *
ConnectivityT::subSection(bool *select)
{
  IndexType i, j;
  IndexType tgcnt = 0;
  for(i = 0; i < size; ++i)
    if(select[i])
      tgcnt += num(i);
  IndexType *newPtr  = new IndexType[size+1];
  IndexType *newTrg = new IndexType[tgcnt];
  // fprintf(stderr, "Size is %d\n", tgcnt);
  tgcnt = 0;
  for(i = 0; i < size; ++i) {
    newPtr[i] = tgcnt;
    if(select[i])
      for(j = 0; j < num(i); ++j)
        newTrg[tgcnt++] = (*this)[i][j];
  }
  newPtr[size] = tgcnt;
  return new ConnectivityT(size, newPtr, newTrg);
}

CountedConnectivityT::CountedConnectivityT(IndexType ns) : ConnectivityT(ns) {
 cnt = new IndexType[csize()];
}

CountedConnectivityT::~CountedConnectivityT()
{
 delete [] cnt;
}

void
CountedConnectivityT::end_count()
{
 IndexType i;
 for(i=0; i < csize(); ++i)
   cnt[i] = pointer[i];
 ConnectivityT::end_count();
}
void
ConnectivityT::end_count()
{
 for(IndexType i=0; i < size; ++i)
  pointer[i+1] += pointer[i];
 numtarget = pointer[size];
 target = new IndexType[numtarget];
}

void
ConnectivityT::countlink(IndexType from, IndexType)
{
 pointer[from]++;
}

void
ConnectivityT::addlink(IndexType from, IndexType to)
{
 IndexType t = --pointer[from];
 target[t] = to;
}

void
CountedConnectivityT::remove(IndexType from, IndexType to)
{
 cnt[from]--;
 IndexType i, ib, ie;
 ib = pointer[from];
 ie = ib+cnt[from];
 for(i =ib; i < ie; ++i)
  if(target[i] == to) break;
 for(; i < ie; ++i)
  target[i] = target[i+1];
}

ConnectivityT::ConnectivityT(IndexType ns)
{
 removeable = true;
 size = ns;
 pointer = new IndexType[size+1];
 IndexType i;
 for(i=0; i < size+1; ++i)
   pointer[i] = 0;
 numtarget = 0;
 target = 0;
 weight = 0;
}

ConnectivityT::ConnectivityT(FILE *f, IndexType n)
{
  removeable = true;
  weight = (float *) 0;
  numtarget = n;
  target = new IndexType[numtarget];

  IndexType error = fscanf(f,"%d",&size);

  pointer = new IndexType[size+1];

  IndexType k = 0, m;
  for (IndexType i = 0; i < size; ++i) {
    fscanf(f,"%d",&m);
    pointer[i] = k;
    if(k + m > numtarget) {
      fprintf(stderr," *** ERROR: ConnectivityT has too many connections\n");
      exit(1);
    }
    for (IndexType iele = 0; iele < m; ++iele) {
      fscanf(f,"%d",target+k);
      target[k++]--;
    }
  }
  pointer[size] = k;
}
#endif
