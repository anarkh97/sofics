#ifndef ROM_RENUMBERINGUTILS_H
#define ROM_RENUMBERINGUTILS_H

#include <Driver.d/Domain.h>

#include <utility>
#include <map>
#include <vector>

class Connectivity;

namespace Rom {

class MeshRenumbering {
public:
  typedef std::vector<int> Restriction;
  typedef std::map<int, int> Extension;

  const Restriction &reducedNodeIds() const { return reducedNodeIds_; }
  const Extension &nodeRenumbering() const { return nodeRenumbering_; }

  const Restriction &reducedElemIds() const { return reducedElemIds_; }
  const Extension &elemRenumbering() const { return elemRenumbering_; }

  template <typename IdxInIt>
  MeshRenumbering(IdxInIt firstSampleElem, IdxInIt lastSampleElem,
                  const Connectivity &elemToNode, bool verboseFlag = true);
  
  // additional constructor for initializing when surface topologies are specified
  template <typename IdxInIt>
  MeshRenumbering(IdxInIt firstSampleElem, IdxInIt lastSampleElem,
                  const Connectivity &elemToNode, Domain * domain, 
                  bool verboseFlag = true);

protected:
  MeshRenumbering() {}

  Restriction reducedNodeIds_;
  Extension nodeRenumbering_;
  
  Restriction reducedElemIds_;
  Extension elemRenumbering_;

private:
  void init(const Connectivity &, bool verboseFlag);
  void init(const Connectivity &, Domain * domain, bool verboseFlag);

  // Disallow copy & assignment
  MeshRenumbering(const MeshRenumbering &);
  MeshRenumbering &operator=(const MeshRenumbering &);
};

template <typename IdxInIt>
MeshRenumbering::MeshRenumbering(IdxInIt firstSampleElem, IdxInIt lastSampleElem,
                                 const Connectivity &elemToNode, bool verboseFlag) :
  reducedElemIds_(firstSampleElem, lastSampleElem)
{
  init(elemToNode, verboseFlag);
}

template <typename IdxInIt>
MeshRenumbering::MeshRenumbering(IdxInIt firstSampleElem, IdxInIt lastSampleElem,
                                 const Connectivity &elemToNode, Domain * domain, bool verboseFlag) :
  reducedElemIds_(firstSampleElem, lastSampleElem)
{
  init(elemToNode, domain, verboseFlag);
}

class SampledMeshRenumbering : public MeshRenumbering {
public:
  const Restriction &sampleNodeIds() const { return sampleNodeIds_; }
  const Restriction &reducedSampleNodeIds() const { return reducedSampleNodeIds_; }

  template <typename IdxInIt>
  SampledMeshRenumbering(IdxInIt firstSampleNode, IdxInIt lastSampleNode,
                         const Connectivity &nodeToNode,
                         const Connectivity &nodeToElem);

private:
  Restriction sampleNodeIds_;
  Restriction reducedSampleNodeIds_;

  void init(const Connectivity &, const Connectivity &);
};

template <typename IdxInIt>
SampledMeshRenumbering::SampledMeshRenumbering(IdxInIt firstSampleNode,
                                               IdxInIt lastSampleNode,
                                               const Connectivity &nodeToNode,
                                               const Connectivity &nodeToElem) :
  sampleNodeIds_(firstSampleNode, lastSampleNode)
{
  init(nodeToNode, nodeToElem);
}

// Precondition: The range [first, last) must not contain duplicate elements
template <typename InputIterator, typename OutputIterator>
OutputIterator
inverse_numbering(InputIterator first, InputIterator last,
                  OutputIterator result) {
  int index = 0;

  for (InputIterator it = first; it != last; ++it) {
    *result++ = std::pair<const int, int>(*it, index++);
  }

  return result;
}

template <typename InputIterator, typename OutputIterator>
OutputIterator
renumber(const std::map<int, int> &renum, InputIterator first, InputIterator last, OutputIterator result) {
  for (InputIterator it = first; it != last; ++it) {
    std::map<int, int>::const_iterator r = renum.find(*it);
    if (r != renum.end()) {
      *result++ = r->second;
    }
  }

  return result;
}

} /* end namespace Rom */

#endif /* ROM_RENUMBERINGUTILS_H */
