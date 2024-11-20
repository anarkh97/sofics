#ifndef ROM_DISTRDOMAINUTILS_H
#define ROM_DISTRDOMAINUTILS_H

#include <Driver.d/SubDomain.h>
#include <Driver.d/DecDomain.h>
#include <Dist.d/DistDom.h>

#include <vector>
#include <algorithm>

namespace Rom {

// Decomposed domain creation
template <typename Scalar>
inline
GenDecDomain<Scalar> *
createDecDomain(Domain *domain) {
  GenDecDomain<Scalar> *result;
#ifdef DISTRIBUTED
  result = new GenDistrDomain<Scalar>(domain);
#else
  result = new GenDecDomain<Scalar>(domain);
#endif
  return result;
}

// Master nodes of the subdomain (in local indexing)
template <typename BoolOutIt>
BoolOutIt
master_node_flags(const SubDomain &subDom, BoolOutIt result) {
  SubDomain &sd = const_cast<SubDomain &>(subDom); // Fix constness
  const int mySubId = sd.subNum();

  const int nodeCount = sd.numNode();
  std::vector<bool> work(nodeCount, true); // vector<bool> specialization acceptable here

  auto &sharedNodes = sd.getSComm()->sharedNodes;
  
  const int neighborCount = sd.getSComm()->numNeighb;
  for (int iNeighbor = 0; iNeighbor < neighborCount; ++iNeighbor) {
    const int neighborId = sd.getSComm()->subNums[iNeighbor];
    if (mySubId > neighborId) {
      const int sharedNodeCount = sharedNodes->num(iNeighbor);
      for (int iNode = 0; iNode < sharedNodeCount; ++iNode) {
        const int nodeId = (*sharedNodes)[iNeighbor][iNode];
        work[nodeId] = false;
      }
    }
  }

  // remove nodes without any displacement/temperature dofs (e.g. internal Lagrange multiplier nodes)
  for(int iNode = 0; iNode < nodeCount; ++iNode) {
    if(! (*sd.getDSA())[iNode].contains(DofSet::XYZdisp | DofSet::XYZrot | DofSet::Temp) ) work[iNode] = false;
  }

#if defined(HACK_INTEL_COMPILER_ITS_CPP11) && (__GLIBC__ == 2) && (__GLIBC_MINOR__ == 12)
  // workaround issue in older version of glibc when compiling with -std=c++11
  std::vector<bool>::iterator first = work.begin(), last = work.end();
  while (first!=last) {
    bool firstval = *first;
    *result = firstval;
    ++result; ++first;
  }
  return result;
#else
  return std::copy(work.begin(), work.end(), result);
#endif
}

// Internal nodes of the subdomain (in local indexing)
template <typename BoolOutIt>
BoolOutIt
internal_node_flags(const SubDomain &subDom, BoolOutIt result) {
  SubDomain &sd = const_cast<SubDomain &>(subDom); // Fix constness
  const int mySubId = sd.subNum();

  const int nodeCount = sd.numNode();
  std::vector<bool> work(nodeCount, false); // vector<bool> specialization acceptable here

  // add nodes without any displacement/temperature dofs (e.g. internal Lagrange multiplier nodes)
  for(int iNode = 0; iNode < nodeCount; ++iNode) {
    if(! (*sd.getDSA())[iNode].contains(DofSet::XYZdisp | DofSet::XYZrot | DofSet::Temp) ) work[iNode] = true;
  }

#if defined(HACK_INTEL_COMPILER_ITS_CPP11) && (__GLIBC__ == 2) && (__GLIBC_MINOR__ <= 12)
  // workaround issue in older version of glibc when compiling with -std=c++11
  std::vector<bool>::iterator first = work.begin(), last = work.end();
  while (first!=last) {
    bool firstval = *first;
    *result = firstval;
    ++result; ++first;
  }
  return result;
#else
  return std::copy(work.begin(), work.end(), result);
#endif
}

} /* end namespace Rom */

#endif /* ROM_DISTRDOMAINUTILS_H */
