#include <cstdio>
#include <iostream>
#include <Utils.d/dofset.h>
#include <Element.d/Element.h>
#include "dofset.h"


// PJSA: note changed lagrange dofs to start at 1 (instead of 6) for rigid element corotators

int
      DofSet::Xdisp = 1 << 0,
      DofSet::Ydisp = 1 << 1,
      DofSet::Zdisp = 1 << 2,
      DofSet::XYZdisp = (1 << 0) + (1 << 1) + (1 << 2),
      DofSet::Xrot = 1 << 3,
      DofSet::Yrot = 1 << 4,
      DofSet::Zrot = 1 << 5,
      DofSet::XYZrot = (1 << 3) + (1 << 4) + (1 << 5),
      DofSet::Temp = 1 << 6,
      DofSet::Helm = 1 << 7,
      DofSet::IntPress = 1 << 8,
      DofSet::Contact = 1 << 9,
      DofSet::Potential = 1 << 10,
      DofSet::LagrangeE = 1 << 11, // Lagrange multiplier for equality constraint
      DofSet::LagrangeI = 1 << 12; // Lagrange multiplier for inequality constraint

int DofSet::nonL_dof = (1 << 0) + (1 << 1) + (1 << 2) + (1 << 3) + (1 << 4) + (1 << 5) + (1 << 6) + (1 << 11) + (1 << 12);
    
DofSet DofSet::nullDofset(-1);

void
DofSet::print(char* msg) const
{
  if(msg) std::cerr<<" --- DofSet "<<msg<<" contains "<<count()<<" dofs: ";
  else    std::cerr<<" --- DofSet contains "<<count()<<" dofs: ";
  int dmap[10];
  DofSet sevenDofs;
  sevenDofs.mark(DofSet::XYZdisp | DofSet::XYZrot | DofSet::Temp | DofSet::Helm | DofSet::Contact | DofSet::IntPress);
  sevenDofs.number(*this, dmap);
  for(int i=0; i<count(); i++)
    std::cerr<<dmap[i]<<" ";
  std::cerr<<std::endl;
}

int
DofSet::count() const
{
  int cdof = flags;
  int num = 0;
  int dofn;
  for(dofn = max_known_dof; dofn--; cdof >>= 1)
    if(cdof & 1) num++;
  return num;
}

int 
DofSet::locate(int dof) const
{
  int cdof = flags;
  int seekeddof = (int) dof;
  int num = 0;
  int dofn;
  for(dofn = max_known_dof; dofn--; ) {
    if(seekeddof & 1)  { 
       return (cdof & 1) ? num : -1;
    }
    if(cdof & 1) num++;
    cdof >>= 1;
    seekeddof >>=1;
  }
  return -1; // should normally never happen if well used
}

int
DofSet::number(DofSet r, int *list) const
{
 int count = 0;
 int cdof = flags;
 int seekeddof = r.flags;
 int num = 0;
 int dofn;
 for(dofn = max_known_dof; dofn--; ) {
    if(seekeddof & 1)
       list[count++] = (cdof & 1) ? num : -1;
    if(cdof & 1) num++;
    cdof >>= 1;
    seekeddof >>=1;
  }
 return count;
}

//for DEC
DofSetArray::DofSetArray(int nnode, int *dofsPerNode, int *renumtable)
{
 numnodes = nnode;
 dofs.resize(numnodes);
 node_num_dofs.resize(numnodes);
 for(int inode=0; inode < nnode; ++inode)
   {
	dofs[inode].mark(DofSet::XYZdisp);
	if(dofsPerNode[inode] == 6)
		dofs[inode].mark(DofSet::XYZrot);
		   
   }
   if (renumtable != nullptr)
       renummap.assign(renumtable, renumtable + nnode);

 makeOffset();
}

DofSetArray::DofSetArray(int nnode, Elemset &eles, int *renumtable, int _myMap)
{
 numnodes = nnode;
 dofs.resize(numnodes);
 node_num_dofs.resize(numnodes);
 int nele = eles.size();
 int iele;
 for(iele=0; iele < nele; ++iele) {
   Element *c_ele = eles[iele];
   if(c_ele)
     c_ele->markDofs(*this);
 }
 if(renumtable != nullptr)
     renummap.assign(renumtable, renumtable+nnode);
 if(_myMap)
 	delete [] renumtable;
 makeOffset();
}

DofSetArray::DofSetArray(Element *ele)
{
 // make a DofSetArray for one element, using the element's local node numbering
 numnodes = ele->numNodes();
 dofs.resize(numnodes);
 node_num_dofs.resize(numnodes);

 int *orig_nn = ele->nodes();
 int *new_nn = new int[numnodes];
 for(int i=0; i<numnodes; ++i) new_nn[i] = i;
 ele->renum(new_nn);
 ele->markDofs(*this);
 ele->renum(orig_nn);
 delete [] orig_nn;
 delete [] new_nn;

 makeOffset();
}

void
DofSetArray::clean_up()
{
	dofs.clear();
	node_num_dofs.clear();
	rowcolnum.clear();
	invrowcol.clear();
}

DofSetArray::DofSetArray(int nnode, int *renumtable, int _myMap)
{
 numnodes = nnode;
 dofs.resize(numnodes);
 node_num_dofs.resize(numnodes);
 int i;
 for(i=0; i<numnodes; ++i)
   node_num_dofs[i] = 0;
 if(renumtable != nullptr)
     renummap.assign(renumtable, renumtable+nnode);
}

int *
DofSetArray::makeDofTypeArray()
{
 if(dofType.size()) return dofType.data();

 dofType.assign(size(), 0);

 // mark the rotational degrees of freedom as 1
 for(int i = 0; i < numnodes; ++i) {
   int xr = locate(i, DofSet::Xrot);
   if(xr >= 0) dofType[xr] = 1;
   int yr = locate(i, DofSet::Yrot);
   if(yr >= 0) dofType[yr] = 1;
   int zr = locate(i, DofSet::Zrot);
   if(zr >= 0) dofType[zr] = 1;
 }

 return dofType.data();
}

void
DofSetArray::makeOffset()
{
 std::vector<int> remapedoffset(numnodes+1);

 int inode;
 for(inode = 0; inode < numnodes; ++inode)
    remapedoffset[inode] = 0;

 for(inode = 0; inode < numnodes; ++inode)
  {
   int node = (renummap.size()) ? renummap[inode] : inode;
   if(node >= 0) {
      node_num_dofs[inode] = dofs[inode].count();
      remapedoffset[node]  = node_num_dofs[inode];
   } else
      node_num_dofs[inode] = 0;
  }

 int offset = 0;
 for(inode = 0; inode < numnodes; ++inode) {
    int tmp = remapedoffset[inode];
    remapedoffset[inode] = offset;
    offset += tmp;
 }

 node_offset.clear();
 if(renummap.size()) {
   node_offset.resize(numnodes+1);
   for(inode = 0; inode < numnodes; ++inode)
    {
     int jnode = renummap[inode];
     node_offset[inode] = (jnode >=0 && node_num_dofs[inode] > 0) ? 
                                           remapedoffset[jnode] : -1;
    }
 } else 
   node_offset = std::move(remapedoffset);
 node_offset[numnodes] = offset;
}

ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, const ConstrainedDSA &cdsa, int ns, const int *sing)
{
  int i, inode;
  numnodes = dsa.numnodes;

  dofs         .resize(numnodes);
  node_num_dofs.resize(numnodes);
  node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;
  for(i = 0; i < numnodes; ++i) 
     dofs[i] = dsa.dofs[i];

  // set weight for augmented virtual nodes, assuming no real corner nodes have zero active dofs
  for(inode = 0; inode < numnodes; ++inode) 
    if(dofs[inode].count()==0) node_num_dofs[inode] = dsa.getWeight(inode);
    else node_num_dofs[inode] = 0;  // these are set in makeModifiedOffset()

  int ndof  = dsa.size();
  rowcolnum.resize(ndof);
  invrowcol.resize(ndof);

  // 1. initialize rowcolnum to the same constraints as in cdsa
  for(i=0; i<ndof; ++i)
    rowcolnum[i] = (cdsa.rowcolnum[i] >= 0) ? 0 : 1;
  // 2. loop over the number of singularities   
  for(i=0; i<ns; ++i)
    if(sing[i] >= 0)
      rowcolnum[sing[i]] = 1;

  // 3. final loop where you start numbering.
  int cnt  = 0;
  int cnt2 = 0;
  for(i=0; i<ndof; ++i) {
    if(rowcolnum[i]) {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;
    } else {
      rowcolnum[i] = cnt++;   // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
  }
  invrowcolmax = cnt2;

  // loop over the nodes and unmark singular dofs
  for(i = 0; i < numnodes; ++i) {
    int dofNums[DofSet::max_known_dof];
    DofSet dofSet = dsa.dofs[i];
    int ndofs = dsa.number(i, dofSet, dofNums);
    int j;
    int flag = 1;
    for(j = 0; j < ndofs; ++j) {
      while((flag & dofSet.list()) == 0)
        flag <<= 1;
      if(rowcolnum[dofNums[j]] == -1)
        dofs[i].unmark(flag);
      flag <<= 1;
    }
  }

  makeModifiedOffset();
}


ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, int ns, const int *sing)
{
  int i, inode;
  numnodes = dsa.numnodes;

  dofs          .resize(numnodes);
  node_num_dofs .resize(numnodes);
  node_offset   .resize(numnodes+1);

  renummap = dsa.renummap;
  for(i = 0; i < numnodes; ++i) 
     dofs[i] = dsa.dofs[i];

  // set weight for augmented virtual nodes, assuming no real corner nodes have zero active dofs
  for(inode = 0; inode < numnodes; ++inode) 
    if(dofs[inode].count()==0) node_num_dofs[inode] = dsa.getWeight(inode);
    else node_num_dofs[inode] = 0;  // these are set in makeModifiedOffset()

  int ndof  = dsa.size();
  rowcolnum.resize(ndof);
  invrowcol.resize(ndof);

  // 1. initialize rowcolnum to zero
  for(i=0; i<ndof; ++i)
    rowcolnum[i] = 0;
  // 2. loop over the number of singularities   
  for(i=0; i<ns; ++i)
    if(sing[i] >= 0)
      rowcolnum[sing[i]] = 1;
  // 3. final loop where you start numbering.
  int cnt  = 0;
  int cnt2 = 0;
  for(i=0; i<ndof; ++i) {
    if(rowcolnum[i]) {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;
    } else {
      rowcolnum[i] = cnt++;   // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
  }
  invrowcolmax = cnt2;

  // loop over the nodes and unmark singular dofs
  for(i = 0; i < numnodes; ++i) {
    int dofNums[DofSet::max_known_dof];
    DofSet dofSet = dsa.dofs[i];
    int ndofs = dsa.number(i, dofSet, dofNums);
    int j;
    int flag = 1;
    for(j = 0; j < ndofs; ++j) {
      while((flag & dofSet.list()) == 0)
        flag <<= 1;
      if(rowcolnum[dofNums[j]] == -1)
        dofs[i].unmark(flag);
      flag <<= 1;
    }
  }

  makeModifiedOffset();
}


ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, int nbc, const BCond *bcd, int ncbc, const ComplexBCond *cbcd)
{
  int i, inode;
  numnodes = dsa.numnodes;

  dofs         .resize(numnodes);
  node_num_dofs.resize(numnodes);
  node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;

  for(i = 0; i < numnodes; ++i)
     dofs[i] = dsa.dofs[i];

  for(i = 0; i < nbc; ++i) 
  {
    int flag = (1 << bcd[i].dofnum);
    if(bcd[i].nnum < numnodes) dofs[bcd[i].nnum].unmark(flag);
  }

  for(i = 0; i < ncbc; ++i)
  {
    int flag = (1 << cbcd[i].dofnum);
    if(cbcd[i].nnum < numnodes) dofs[cbcd[i].nnum].unmark(flag);
  }

  for(inode = 0; inode < numnodes; ++inode)
    node_num_dofs[inode] = dofs[inode].count();

  makeOffset();

  int ndof  = dsa.size();
  rowcolnum.resize(ndof);
  invrowcol.resize(ndof);

  // 1. initialize rowcolnum to zero
  for(i=0; i<ndof; ++i)
    rowcolnum[i] = 0;

  // 2. loop over the nodes
  for(i = 0; i < numnodes; ++i) {
    int dofNums[DofSet::max_known_dof];
    DofSet fixedDofs = dofs[i] ^ dsa.dofs[i];
    int ndofs = dsa.number(i, fixedDofs, dofNums);
    int j;
    for(j = 0; j < ndofs; ++j)
       rowcolnum[dofNums[j]] = 1;
  }

  // 3. final loop where you start numbering.

  int cnt  = 0;
  int cnt2 = 0;
  for(i=0; i<ndof; ++i) {
    if(rowcolnum[i]) {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;
    } else {
      rowcolnum[i] = cnt++;   // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
  }
  invrowcolmax = cnt2;

}

//new constructor for the updating taking in account both bound cond :Dirichlet and measured dof

ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, const DofSetArray &c_dsa, int nbc, const BCond *bcd, int info)
{
  int i, inode;
  numnodes = dsa.numnodes;

  // The dofs are created equal to zero (no dof marked)
	dofs         .resize(numnodes);
	node_num_dofs.resize(numnodes);
	node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;

  if(info==1) {    // for Z11
    for(i = 0; i < nbc; ++i) {
      // Should check that  bcd[i].nnum < numnodes
      int flag = (1 << bcd[i].dofnum);
      int thisDofInCDSA = c_dsa.locate(bcd[i].nnum, flag);
      if(thisDofInCDSA < 0) {
        fprintf(stderr," *** WARNING: A measured DOF is constrained\n");
        continue;
      }
      dofs[bcd[i].nnum].mark(flag);
    }
  }
  else if (info==0) {            // for Z12    and Z22
    for(i = 0; i < numnodes; ++i)
      dofs[i] = c_dsa.dofs[i];
 
    for(i = 0; i < nbc; ++i) {
      int flag = (1 << bcd[i].dofnum);
      if(bcd[i].nnum < numnodes) dofs[bcd[i].nnum].unmark(flag);
    }
  }

  for(inode = 0; inode < numnodes; ++inode) {
    node_num_dofs[inode] = dofs[inode].count();
  }
  makeOffset();
}

ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, int nbc, const BCond *bcd, const int *bc)
{
  int i, inode;
  numnodes = dsa.numnodes;

	dofs         .resize(numnodes);
	node_num_dofs.resize(numnodes);
	node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;

  for(i = 0; i < numnodes; ++i)
     dofs[i] = dsa.dofs[i];

  for(i = 0; i < nbc; ++i) {
    int flag = (1 << bcd[i].dofnum);
    if(bcd[i].nnum < numnodes) dofs[bcd[i].nnum].unmark(flag);
  }

  for(inode = 0; inode < numnodes; ++inode)
    node_num_dofs[inode] = dofs[inode].count();

  makeOffset();

  rowcolnum.resize(dsa.size());
  invrowcol.resize(dsa.size());

  int cnt = 0, cnt2 = 0;

  for(i=0; i<dsa.size(); ++i)
  {
    if( bc[i] != BCFIXED ) 
    {  
      rowcolnum[i] = cnt++; // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
    else
    {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;  // To renumber constrained dofs.
    } 
  }
  invrowcolmax = cnt2;

}

// New constructor for Experimental version of FETI
ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, int nbc, const BCond *bcond,
                               int numCornerNodes, const std::vector<int> &cornerNodes,
                               const std::vector<DofSet> &cornerDofs,
                               int ncbc, const ComplexBCond *cbcond, int numWetInterfaceNodes,
                               const int *wetInterfaceNodes, const DofSet *wetInterfaceDofs)
{
  int i, inode;

  numnodes = dsa.numnodes;

	dofs         .resize(numnodes);
	node_num_dofs.resize(numnodes);
	node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;

  for(i = 0; i < numnodes; ++i)
     dofs[i] = dsa.dofs[i];

  // work on the real applied boundary conditions
  for(i = 0; i < nbc; ++i) {
    int flag = (1 << bcond[i].dofnum);
    dofs[bcond[i].nnum].unmark(flag);
  }

  // work on the complex applied boundary conditions
  for(i = 0; i < ncbc; ++i) {
    int flag = (1 << cbcond[i].dofnum);
    dofs[cbcond[i].nnum].unmark(flag);
  }

  // work on the corner dofs now
  for(i = 0; i < numCornerNodes; ++i) {
    dofs[cornerNodes[i]].unmark(cornerDofs[i].list());
  }

  // work on the wet interface nodes now
  for(i = 0; i < numWetInterfaceNodes; ++i) {
    dofs[wetInterfaceNodes[i]].unmark(wetInterfaceDofs[i].list());  
  }

  for(inode = 0; inode < numnodes; ++inode)
    node_num_dofs[inode] = dofs[inode].count();

  makeOffset();

  int ndof  = dsa.size();
	rowcolnum.resize(ndof);
	invrowcol.resize(ndof);

  // 1. initialize rowcolnum to zero
  for(i=0; i<ndof; ++i)
    rowcolnum[i] = 0;

  // 2. loop over the nodes
  for(i = 0; i < numnodes; ++i) {
    int dofNums[DofSet::max_known_dof];
    DofSet fixedDofs = dofs[i] ^ dsa.dofs[i];
    int ndofs = dsa.number(i, fixedDofs, dofNums);
    int j;
    for(j = 0; j < ndofs; ++j)
       rowcolnum[dofNums[j]] = 1;
  }

  // 3. final loop where you start numbering.

  int cnt  = 0;
  int cnt2 = 0;
  for(i=0; i<ndof; ++i) {
    if(rowcolnum[i]) {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;
    } else {
      rowcolnum[i] = cnt++;   // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
  }
}


void
DofSetArray::setWeight(int node, int w)
{
 node_num_dofs[node] = w;
}

int 
DofSetArray::getWeight(int node) const
{
 return node_num_dofs[node];
}


void
DofSetArray::makeModifiedOffset()
{
 std::vector<int> remapedoffset(numnodes+1);

 int inode;
 for(inode = 0; inode < numnodes; ++inode)
    remapedoffset[inode] = 0;

 for(inode = 0; inode < numnodes; ++inode)
  {
   int node = (renummap.size() != 0) ? renummap[inode] : inode;
   if(node >= 0) {
      if(node_num_dofs[inode] == 0)
        node_num_dofs[inode] = dofs[inode].count();
      remapedoffset[node]  = node_num_dofs[inode];
   } else
      node_num_dofs[inode] = 0;
  }

 int offset = 0;
 for(inode = 0; inode < numnodes; ++inode) {
    int tmp = remapedoffset[inode];
    remapedoffset[inode] = offset;
    offset += tmp;
 }

 if(renummap.size() != 0) {
   node_offset.resize(numnodes+1);
   for(inode = 0; inode < numnodes; ++inode)
    {
     int jnode = renummap[inode];
     node_offset[inode] = (jnode >=0 && node_num_dofs[inode] > 0) ?
                                           remapedoffset[jnode] : -1;
    }
 } else
   node_offset = std::move(remapedoffset);
 node_offset[numnodes] = offset;
}

ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, int nbc, const ComplexBCond *bcd,
			       const int *bc)
{
  int i, inode;

  numnodes = dsa.numnodes;

	dofs         .resize(numnodes);
	node_num_dofs.resize(numnodes);
	node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;

  for(i = 0; i < numnodes; ++i)
     dofs[i] = dsa.dofs[i];

  for(i = 0; i < nbc; ++i)
  {
    int flag = (1 << bcd[i].dofnum);
    if(bcd[i].nnum < numnodes) dofs[bcd[i].nnum].unmark(flag);
  }

  for(inode = 0; inode < numnodes; ++inode)
    node_num_dofs[inode] = dofs[inode].count();

  makeOffset();

	rowcolnum.resize(dsa.size());
	invrowcol.resize(dsa.size());

  int cnt  = 0;
  int cnt2 = 0;

  for(i=0; i<dsa.size(); ++i) {
    if( bc[i] != BCFIXED ) {  
      rowcolnum[i] = cnt++; // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
    else {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;  // To renumber constrained dofs.
    } 
  }
  invrowcolmax = cnt2;
}

ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, const BCond *bcond, int nbc,
			       const ComplexBCond *cbcond, int ncbc, const int *bc)
{
  // bcond   = boundary condition
  // nbc     = number boundary conditions

  // cbcond  = complex boundary condition
  // ncbc    = number complex boundary conditions 
  // dsa     = dof set array

  int i, inode;

  numnodes = dsa.numnodes;

	dofs         .resize(numnodes);
	node_num_dofs.resize(numnodes);
	node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;

  for(i = 0; i < numnodes; ++i)
     dofs[i] = dsa.dofs[i];

  for(i = 0; i < nbc; ++i)
  {
    int flag = (1 << bcond[i].dofnum);
    if(bcond[i].nnum < numnodes) dofs[bcond[i].nnum].unmark(flag);
  }

  for(i = 0; i < ncbc; ++i)
  {
    int flag = (1 << cbcond[i].dofnum);
    if(cbcond[i].nnum < numnodes) dofs[cbcond[i].nnum].unmark(flag);
  }

  for(inode = 0; inode < numnodes; ++inode)
    node_num_dofs[inode] = dofs[inode].count();

  makeOffset();

  int ndof = dsa.size();

	rowcolnum.resize(dsa.size());
	invrowcol.resize(dsa.size());

  int cnt  = 0;
  int cnt2 = 0;

  for(i=0; i<ndof; ++i) {
    if( bc[i] != BCFIXED ) { 
      rowcolnum[i] = cnt++;   // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
    else {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;  // To renumber constrained dofs.
    }
  }
  invrowcolmax = cnt2;
}

ConstrainedDSA::ConstrainedDSA(const DofSetArray &dsa, const ConstrainedDSA &cdsa)
{
  // make a ConstrainedDSA from another ConstrainedDSA and in addition
  // constrained all of the Lagrange multiplier dofs
  int i, inode;
  numnodes = dsa.numnodes;

	dofs         .resize(numnodes);
	node_num_dofs.resize(numnodes);
	node_offset  .resize(numnodes+1);

  renummap = dsa.renummap;

  for(i = 0; i < numnodes; ++i) {
     dofs[i] = cdsa.dofs[i];
     dofs[i].unmark(DofSet::LagrangeE);
     dofs[i].unmark(DofSet::LagrangeI);
  }

  for(inode = 0; inode < numnodes; ++inode)
    node_num_dofs[inode] = dofs[inode].count();

  makeOffset();

  int ndof  = dsa.size();
  rowcolnum.resize(ndof);
  invrowcol.resize(ndof);

  // 1. initialize rowcolnum to zero
  for(i=0; i<ndof; ++i)
    rowcolnum[i] = 0;

  // 2. loop over the nodes
  for(i = 0; i < numnodes; ++i) {
    int dofNums[DofSet::max_known_dof];
    DofSet fixedDofs = dofs[i] ^ dsa.dofs[i];
    int ndofs = dsa.number(i, fixedDofs, dofNums);
    int j;
    for(j = 0; j < ndofs; ++j)
       rowcolnum[dofNums[j]] = 1;
  }

  // 3. final loop where you start numbering.

  int cnt  = 0;
  int cnt2 = 0;
  for(i=0; i<ndof; ++i) {
    if(rowcolnum[i]) {
      rowcolnum[i] = -1;
      invrowcol[i] = cnt2++;
    } else {
      rowcolnum[i] = cnt++;   // To renumber unconstrained dofs.
      invrowcol[i] = -1;
    }
  }
  invrowcolmax = cnt2;

}

void
DofSetArray::mark(int node, int ds)
{
 dofs[node].mark(ds);
}

void
DofSetArray::mark(const int *node, int numNodes, int ds)
{
 int iNode;
 for(iNode=0; iNode<numNodes; ++iNode)
   dofs[node[iNode]].mark(ds);
}

int
DofSetArray::locate(int node, int ds) const
{
 if(node >= numnodes) return -1;
 int dnum = dofs[node].locate(ds);
 if(dnum < 0) return dnum;
 return  node_offset[node] + dnum;
}

int
DofSetArray::number(int node, DofSet ds, int *p) const
{
 if(node >= numnodes) return 0;
 int nn = dofs[node].number(ds,p);
 int i;
 for(i=0; i < nn; ++i) {
   if(p[i] >= 0) p[i] += node_offset[node];
 }
 return nn;
}

DofSetArray::~DofSetArray()
{
}

SimpleNumberer::SimpleNumberer(int nnodes, int *renum, int _myMap)
{
	numnodes = nnodes;
	node_offset.resize(numnodes+1);
	node_num_dofs.resize(numnodes);
	if(renum != nullptr)
        renummap.assign(renum, renum+nnodes);

	if(_myMap) {
		delete []renum;
	}
}

void
SimpleNumberer::setWeight(int node, int w)
{
 node_num_dofs[node] = w;
}

void
SimpleNumberer::makeOffset()
{
	std::vector<int> remapedoffset(numnodes+1);

	int inode;
	for(inode = 0; inode < numnodes; ++inode)
		remapedoffset[inode] = 0;

	for(inode = 0; inode < numnodes; ++inode) {
		int node = (renummap.size()) ? renummap[inode] : inode;
		if(node >= 0)
			remapedoffset[node] = node_num_dofs[inode];
	}

	int offset = 0;
	for(inode = 0; inode < numnodes; ++inode) {
		int tmp = remapedoffset[inode];
		remapedoffset[inode] = offset;
		offset += tmp;
	}

	node_offset.clear();
	if(renummap.size()) {
		node_offset.resize(numnodes+1);
		for(inode = 0; inode < numnodes; ++inode) {
			int jnode = renummap[inode];
			node_offset[inode] = (jnode >=0 && node_num_dofs[inode] > 0) ?
			                     remapedoffset[jnode] : -1;
		}
	}
	else
		node_offset = remapedoffset;
	node_offset[numnodes] = offset;
}

void
EqNumberer::print()
{
 fprintf(stderr,"--- Equation Weights ---\n");
 int i;
 for(i=0; i<numnodes; ++i)
   fprintf(stderr,"Equation Block %d has %d equations\n",i+1,node_num_dofs[i]);
 fprintf(stderr,"--- Offsets ---\n");
 for(i=0; i<=numnodes; ++i)
   fprintf(stderr,"node_offset[%d] = %d\n",i+1,node_offset[i]);
}


