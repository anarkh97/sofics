#include "MeshDesc.h"

#include "MeshOutput.h"
#include "RenumberingUtils.h"
#include "SimpleBuffer.h"

#include <Driver.d/EFrameData.h>

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>

#include <Mortar.d/MortarDriver.d/MortarHandler.h>

#include <functional>
#include <algorithm>
#include <deque>

namespace Rom {

namespace Detail {

template <typename InputIterator>
const CoordSet &
reduce(InputIterator first, InputIterator last, const CoordSet &origin, CoordSet &target) {
  int index = 0;
  for (InputIterator it = first; it != last; ++it) {
    target.nodeadd(index++, const_cast<CoordSet &>(origin).getNode(*it));
  }
  return target;
}

const Elemset &
renumber_nodes(const std::map<int, int> &nodeIndices, Elemset &target);

template <typename InputIterator>
const Elemset &
reduce(InputIterator first, InputIterator last, const Elemset &origin, Elemset &target) {
  int index = 0;
  for (InputIterator it = first; it != last; ++it) {
    target.elemadd(index++, origin[*it]);
  }
  return target;
}

template <typename InputIterator>
const Elemset &
reduce(InputIterator first, InputIterator last, const std::map<int, int> &nodeIndices, const Elemset &origin, Elemset &target) {
  reduce(first, last, origin, target);
  return renumber_nodes(nodeIndices, target);
}

template <typename ValueType>
struct MeshSectionTraits {
  static int ValueType::* const renumbered_slot;
};

template <typename BCondInputIterator, typename BCondOutputIterator>
BCondOutputIterator
reduce(const std::map<int, int> &indices, BCondInputIterator first, BCondInputIterator last, BCondOutputIterator result) {
  typedef typename std::iterator_traits<BCondInputIterator>::value_type ValueType;
  int ValueType::* const slot = MeshSectionTraits<ValueType>::renumbered_slot;
  
  for (BCondInputIterator source = first; source != last; ++source) {
    const ValueType &value = *source;
    const std::map<int, int>::const_iterator it = indices.find(value.*slot);
    if (it != indices.end()) {
      ValueType newValue(value);
      newValue.*slot = it->second;
      *result++ = newValue;
    }
  }

  return result;
}

template <typename BCondInputIterator, typename BCondOutputIterator>
BCondOutputIterator
reduce(const std::map<int, int> &indices, BCondInputIterator first, BCondInputIterator last, BCondOutputIterator result, bool temp) {
  typedef typename std::iterator_traits<BCondInputIterator>::value_type ValueType;
  int ValueType::* const slot = MeshSectionTraits<ValueType>::renumbered_slot;

  for (BCondInputIterator source = first; source != last; ++source) {
    const ValueType &value = *source;
    if((!temp && value.dofnum == 6) || (temp && value.dofnum != 6)) continue;
    const std::map<int, int>::const_iterator it = indices.find(value.*slot);
    if (it != indices.end()) {
      ValueType newValue(value);
      newValue.*slot = it->second;
      *result++ = newValue;
    }
  }

  return result;
}

template <typename V, typename PairOutputIterator>
PairOutputIterator
reduce(const std::map<int, int> &indices, const std::map<int, V> &input, PairOutputIterator result) {
  typedef typename std::map<int, V>::const_iterator InputIt;
  for (InputIt source = input.begin(), last = input.end();
       source != last;
       ++source) {
    const typename std::map<int, V>::value_type &p = *source;
    const std::map<int, int>::const_iterator it = indices.find(p.first);
    if (it != indices.end()) {
      *result++ = std::make_pair(it->second, p.second);
    }
  }

  return result;
}

/*
// Pressure on individual elements
template <typename PairOutputIterator>
PairOutputIterator
reduce(const std::map<int, int> &indices, const std::vector<std::pair<int, std::pair<double,int> > > &input, PairOutputIterator result) {
  typedef std::vector<std::pair<int, std::pair<double,int> > >::const_iterator InputIt;
  for (InputIt source = input.begin(), last = input.end();
       source != last;
       ++source) {
    const std::vector<std::pair<int, std::pair<double,int> > >::value_type &p = *source;
    const std::map<int, int>::const_iterator it = indices.find(p.first);
    if (it != indices.end()) {
      *result++ = std::make_pair(it->second, p.second);
    }
  }

  return result;
}
*/

// EFrames from Elements
template <typename EFrameDataOutputIterator>
EFrameDataOutputIterator
copy_eframes(const Elemset & elemSet, EFrameDataOutputIterator result) {
  const int iElemEnd = elemSet.last();
  for (int iElem = 0; iElem != iElemEnd; ++iElem) {
    Element *elem = elemSet[iElem];
    if (!elem) {
      continue;
    } 
    const EFrame *eframe = elem->getFrame();
    if (eframe) {
      EFrameData data;
      data.elnum = iElem;
      std::memcpy(&data.frame, eframe, sizeof(data.frame));
      data.next = NULL;
      
      *result++ = data;
    }
  }

  return result;
}

void
copy_cframes(const std::vector<Attrib> &attributes, std::map<int,FrameData> &cframes) {
  for(std::vector<Attrib>::const_iterator it = attributes.begin(); it != attributes.end(); ++it) {
    if(it->cmp_frm > -1) {
      if(cframes.find(it->cmp_frm) == cframes.end()) {
        FrameData frame;
        frame.num = it->cmp_frm;
        for(int i=0; i<9; ++i) frame.d[i] = geoSource->getCframes()[it->cmp_frm][i];
        cframes[it->cmp_frm] = frame;
      }
    }
  }
}

void
copy_coef(const std::vector<Attrib> &attributes, std::map<int,CoefData> &coefdata) {
  for(std::vector<Attrib>::const_iterator it = attributes.begin(); it != attributes.end(); ++it) {
    if(it->cmp_attr > -1) {
      if(geoSource->getCoefData(it->cmp_attr) == NULL) {
        std::cerr << " *** WARNING: only COEF-type composites are copied to reduced mesh file\n";
      }
      else if(coefdata.find(it->cmp_attr) == coefdata.end()) {
        CoefData coef;
        for(int i=0; i<7; ++i)
          for(int j=0; j<6; ++j)
            coef.c[i][j] = geoSource->getCoefData(it->cmp_attr)->c[i][j];
        coef.coefFlag = geoSource->getCoefData(it->cmp_attr)->coefFlag;
        coefdata[it->cmp_attr] = coef;
      }
    }
  }
}

// Map to sequence conversion

template <typename MapIterator>
class value_iterator {
public:
  typedef typename std::bidirectional_iterator_tag iterator_category;
  typedef typename MapIterator::value_type::second_type value_type;
  typedef const value_type *pointer;
  typedef const value_type &reference;
  typedef typename MapIterator::difference_type difference_type;

  bool operator==(const value_iterator& other) const { return mapIt_ == other.mapIt_; }
  bool operator!=(const value_iterator& other) const { return !(*this == other); }

  const value_type &operator*() const { return mapIt_->second; }
  const value_type *operator->() const { return &mapIt_->second; }

  value_iterator &operator++() { ++mapIt_; return *this; }
  value_iterator operator++(int) { value_iterator temp(*this); ++(*this); return temp; }
  
  value_iterator &operator--() { --mapIt_; return *this; }
  value_iterator operator--(int) { value_iterator temp(*this); --(*this); return temp; }

  explicit value_iterator(MapIterator mapIt) : mapIt_(mapIt) {}

private:
  MapIterator mapIt_;
};

template <typename MapIterator>
value_iterator<MapIterator>
make_value_iterator(MapIterator mapIt) {
  return value_iterator<MapIterator>(mapIt);
}

const Elemset &
renumber_nodes(const std::map<int, int> &nodeIndices, Elemset &target) {
  const int indexEnd = nodeIndices.size() > 0 ? nodeIndices.rbegin()->first + 1 : 0;
  SimpleBuffer<int> buffer(indexEnd);
 
  for (std::map<int, int>::const_iterator it = nodeIndices.begin(); it != nodeIndices.end(); ++it) {
    buffer[it->first] = it->second;
  }

  const int iElemEnd = target.last();
  for (int iElem = 0; iElem < iElemEnd; ++iElem) {
    if (Element * e = target[iElem]) {
      e->renum(buffer.array());
    }
  }

  return target;
}

template <>
int BCond::* const
MeshSectionTraits<BCond>::renumbered_slot = &BCond::nnum;

template <>
int Attrib::* const
MeshSectionTraits<Attrib>::renumbered_slot = &Attrib::nele;

template <>
int EFrameData::* const
MeshSectionTraits<EFrameData>::renumbered_slot = &EFrameData::elnum;

template <>
int PressureBCond::* const
MeshSectionTraits<PressureBCond>::renumbered_slot = &PressureBCond::elnum;

} // end namespace Detail


using namespace Detail;

MeshDesc::MeshDesc(Domain *domain, GeoSource *geoSource, const MeshRenumbering &ren, const std::map<int, double> &weights) :
  properties_(&geoSource->getStructProps()),
  materialLaws_(&geoSource->getMaterialLaws()),
  elemWeights_(1),
  ActiveSurfaces(0,0),
  domain_(domain)
{
  init(domain, geoSource, ren);
  reduce(ren.elemRenumbering(), weights, std::inserter(elemWeights_[0], elemWeights_[0].end())); 
}

MeshDesc::MeshDesc(Domain *domain, GeoSource *geoSource, const MeshRenumbering &ren, const std::vector<std::map<int, double> > &weights) :
  properties_(&geoSource->getStructProps()),
  materialLaws_(&geoSource->getMaterialLaws()),
  elemWeights_(weights.size()),
  ActiveSurfaces(0,0),
  domain_(domain)
{
  init(domain, geoSource, ren);
  for(int j=0; j<weights.size(); ++j) {
    reduce(ren.elemRenumbering(), weights[j], std::inserter(elemWeights_[j], elemWeights_[j].end()));
  }
}
                                                                                      
MeshDesc::MeshDesc(Domain *domain, GeoSource *geoSource, const SampledMeshRenumbering &ren) :
  properties_(&geoSource->getStructProps()),
  materialLaws_(&geoSource->getMaterialLaws()),
  sampleNodeIds_(ren.reducedSampleNodeIds().begin(), ren.reducedSampleNodeIds().end()),
  ActiveSurfaces(0,0),
  domain_(domain)
{
  init(domain, geoSource, ren);
}

void
MeshDesc::init(Domain *domain, GeoSource *geoSource, const MeshRenumbering &ren) {

  // Nodes
  reduce(ren.reducedNodeIds().begin(), ren.reducedNodeIds().end(), domain->getNodes(), nodes_);

  // Elements & element frames
  reduce(ren.reducedElemIds().begin(), ren.reducedElemIds().end(), ren.nodeRenumbering(), *geoSource->getElemSet(), elements_);
  copy_eframes(elements_, std::back_inserter(elemFrames_));

  // Attributes
  const AttribContainer &attrib = geoSource->getAttributes(); 
  reduce(ren.elemRenumbering(), make_value_iterator(attrib.begin()), make_value_iterator(attrib.end()), std::back_inserter(attributes_)); 

    // Composites
  copy_cframes(attributes_, compositeFrames_);
  copy_coef(attributes_, coefData_);

  // Material laws
  const std::map<int, int> &matLawMap = geoSource->getMaterialLawMapping();
  reduce(ren.elemRenumbering(), geoSource->getMaterialLawMapping(), std::inserter(materialLawMapping_, materialLawMapping_.end()));

  // Boundary conditions
  reduce(ren.nodeRenumbering(), domain->getDBC(), domain->getDBC() + domain->nDirichlet(), std::back_inserter(dirichletBConds_), false);
/* Now both constant and time-dependent forces are outputted in reduced coordinates, see ElementSamplingDriver::postProcess
   TODO configuration-dependent nodal forces and moments still should be outputted here, though
  reduce(ren.nodeRenumbering(), domain->getNBC(), domain->getNBC() + domain->nNeumann(), std::back_inserter(neumannBConds_));
*/
  // Element pressures
  reduce(ren.elemRenumbering(), geoSource->getElementPressure().begin(), geoSource->getElementPressure().end(),
         std::inserter(elemPressures_, elemPressures_.end()));

/* Now the initial conditions are outputted in the reduced coordinates, see ElementSamplingDriver::postProcess
  // Initial conditions
  reduce(ren.nodeRenumbering(), domain->getInitDisp(), domain->getInitDisp() + domain->numInitDisp(), std::back_inserter(initDisp_));
  reduce(ren.nodeRenumbering(), domain->getInitDisp6(), domain->getInitDisp6() + domain->numInitDisp6(), std::back_inserter(initDisp_));
  reduce(ren.nodeRenumbering(), domain->getInitVelocity(), domain->getInitVelocity() + domain->numInitVelocity(), std::back_inserter(initVel_));
*/
  // Nodal temperatures
  reduce(ren.nodeRenumbering(), domain->getDBC(), domain->getDBC() + domain->nDirichlet(), std::back_inserter(temperatures_), true);
 
  // Contact Surfaces
  int numMC = domain->GetnMortarConds();
  if(numMC > 0) { // check if any contact conditions are specified
    for(int mc = 0; mc<numMC; ++mc) { // if so, first loop over mortar conditions
      // if the mortar condition is CTC type, store for output to reduced mesh
      if(domain->GetMortarCond(mc)->GetInteractionType() == MortarHandler::CTC){
        ContactContainer CC(GetMortarCond(mc)->GetMasterEntityId(),
                            GetMortarCond(mc)->GetSlaveEntityId(),
                            0,
                            GetMortarCond(mc)->GetNormalTol(),
                            GetMortarCond(mc)->GetTangentialTol());
        ContactPairs_.push_back(CC);
        activeSurfs_.insert(domain->GetMortarCond(mc)->GetMasterEntityId());
        activeSurfs_.insert(domain->GetMortarCond(mc)->GetSlaveEntityId());
      }
    }

    int nCS = 0;
    std::map<int,int> renumMap = ren.nodeRenumbering(); 
    // loop over all active surfaces
    for(std::set<int>::iterator it = activeSurfs_.begin(); it != activeSurfs_.end(); ++it){ 
      // then loop over all surface topologies to find matching ID 
      for(int surf = 0; surf < domain->getNumSurfs(); ++surf) { 
        if(*it == domain->GetSurfaceEntity(surf)->GetId()) {
          ActiveSurfaces[nCS++] = domain->GetSurfaceEntity(surf);// save pointer to SurfaceEntity class
          // then apply new global numbering to node map stored in surface entity
          int *glRenumber = ActiveSurfaces[nCS-1]->GetPtrGlNodeIds();  
          int numNodes = ActiveSurfaces[nCS-1]->GetnNodes();
          for(int lnn = 0; lnn < numNodes; ++lnn) {
            int intBuffer = renumMap[glRenumber[lnn]];
            glRenumber[lnn] = intBuffer;
          }
        }
      }
    }
  }
 
}

std::ostream &
operator<<(std::ostream &out, const MeshDesc &mesh) {
  out << mesh.nodes();
  out << mesh.elements();
  if (!mesh.elemFrames().empty()) {
    out << make_section(mesh.elemFrames().begin(), mesh.elemFrames().end());
  }

  if(mesh.GetnContactSurfacePairs() != 0)
    out << mesh.ContactSurfaces();

  out << make_section(mesh.attributes().begin(), mesh.attributes().end());
  out << mesh.properties();

  if (!mesh.compositeFrames().empty())
    out << make_section(mesh.compositeFrames().begin(), mesh.compositeFrames().end(), CFrameTag());
  if(!mesh.coefData().empty())
    out << make_section(mesh.coefData().begin(), mesh.coefData().end());
 
  if(!mesh.materialLawMapping().empty())
    out << make_section(mesh.materialLawMapping().begin(), mesh.materialLawMapping().end(), MatUsageTag());
  if(!mesh.materialLaws().empty())
    out << make_section(mesh.materialLaws().begin(), mesh.materialLaws().end(), MatLawTag());

  if(!mesh.dirichletBConds().empty())
    out << make_section(mesh.dirichletBConds().begin(), mesh.dirichletBConds().end(), BCond::Displacements);
  if(!mesh.neumannBConds().empty()) 
    out << make_section(mesh.neumannBConds().begin(), mesh.neumannBConds().end(), BCond::Forces);

  if(!mesh.initDisp().empty())
    out << make_section(mesh.initDisp().begin(), mesh.initDisp().end(), BCond::Idisplacements);
  if(!mesh.initVel().empty())
    out << make_section(mesh.initVel().begin(), mesh.initVel().end(), BCond::Ivelocities);

  if(!mesh.elemPressures().empty())
    out << make_section(mesh.elemPressures().begin(), mesh.elemPressures().end(), ElementPressureTag());

  if(!mesh.activeSurfaces().empty())
    out << mesh.GetSurfaceEntities();

  if(!mesh.temperatures().empty())
    out << make_section(mesh.temperatures().begin(), mesh.temperatures().end(), BCond::Temperatures);

  if (!mesh.sampleNodeIds().empty())
    out << make_section(mesh.sampleNodeIds().begin(), mesh.sampleNodeIds().end(), SampleNodeTag());

  for(int j=0; j<mesh.numLocal(); ++j) {
    if (!mesh.elemWeights(j).empty()) {
      out.precision(std::numeric_limits<double>::digits10+1);
      out << make_section(mesh.elemWeights(j).begin(), mesh.elemWeights(j).end(), ElementWeightTag(), j+1);
    }
  }

  return out;
}

} /* end namespace Rom */
