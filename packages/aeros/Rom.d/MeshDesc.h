#ifndef ROM_MESHDESC_H
#define ROM_MESHDESC_H

#include "RenumberingUtils.h"

#include <Element.d/Element.h>
#include <Driver.d/Attrib.h>
#include <Parser.d/AuxDefs.h>
#include <Utils.d/CompositeInfo.h>

#include <Driver.d/Domain.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>

#include <vector>
#include <map>
#include <iterator>
#include <cstring>
#include <utility>

//class Domain;
class GeoSource;

class EFrameData;
class Attrib;
class NLMaterial;
class BCond;

namespace Rom {

struct ContactContainer {
  
  int MasterId;
  int SlaveId;
  int MortarType;

  double NormalTol;
  double TangentTol;

  ContactContainer(int m, int s, int t, double n, double tol): 
    MasterId(m),
    SlaveId(s),
    MortarType(t),
    NormalTol(n),
    TangentTol(tol)
  {}

};

class MeshDesc {
public:
  const CoordSet &nodes() const { return nodes_; }
  const Elemset &elements() const { return elements_; }

  const std::vector<EFrameData> &elemFrames() const { return elemFrames_; }
  const std::map<int,FrameData> &compositeFrames() const { return compositeFrames_; }
  const SPropContainer &properties() const { return *properties_; } 
  const std::map<int,CoefData> &coefData() const { return coefData_; }
  const std::vector<Attrib> &attributes() const { return attributes_; }
  const std::map<int, NLMaterial *> &materialLaws() const { return *materialLaws_; }
  const std::map<int, int> &materialLawMapping() const { return materialLawMapping_; }

  const std::vector<BCond> &dirichletBConds() const { return dirichletBConds_; }
  const std::vector<BCond> &neumannBConds() const { return neumannBConds_; }

  const std::vector<BCond> &initDisp() const { return initDisp_; }
  const std::vector<BCond> &initVel() const { return initVel_; }

  const std::vector<BCond> &temperatures() const { return temperatures_; }

  const std::vector<PressureBCond> &elemPressures() const { return elemPressures_; }
  const std::set<int>              &activeSurfaces() const { return activeSurfs_; }
  
  // Reduced mesh only
  const std::vector<int> &sampleNodeIds() const { return sampleNodeIds_; }
  const std::map<int, double> &elemWeights(int j) const { return elemWeights_[j]; }
  int numLocal() const { return elemWeights_.size(); } // local bases

  // Create element-based reduced mesh
  MeshDesc(Domain *, GeoSource *, const MeshRenumbering &, const std::map<int, double> &weights);

  // Create collection of element-based reduced meshes for local bases method
  MeshDesc(Domain *, GeoSource *, const MeshRenumbering &, const std::vector<std::map<int, double> > &weights);
  
  // Create node-based reduced mesh 
  MeshDesc(Domain *, GeoSource *, const SampledMeshRenumbering &);

  // Contact Surface Utilities
  int GetnContactSurfacePairs() const { return domain_->GetnContactSurfacePairs(); } 
  int GetnMortarConds() const { return domain_->GetnMortarConds(); }
  MortarHandler* GetMortarCond(int i) const { return domain_->GetMortarCond(i); }
  const ResizeArray<SurfaceEntity*>* GetSurfaceEntities() const { return (&ActiveSurfaces); }
  const std::vector<ContactContainer> &ContactSurfaces() const { return ContactPairs_; }
  
private:
  void init(Domain *domain, GeoSource *geoSource, const MeshRenumbering &ren);

  typedef std::map<int, Attrib> AttribContainer;
  
  CoordSet nodes_; 
  Elemset elements_;

  Domain *domain_;

  std::vector<EFrameData> elemFrames_;
  std::map<int,FrameData> compositeFrames_;
  std::vector<Attrib> attributes_;
  const SPropContainer *properties_;
  std::map<int,CoefData> coefData_;
  std::map<int, int> materialLawMapping_;
  const std::map<int, NLMaterial *> *materialLaws_;

  std::vector<BCond> dirichletBConds_;
  std::vector<BCond> neumannBConds_;
  std::vector<BCond> initDisp_; 
  std::vector<BCond> initVel_;
  std::vector<BCond> temperatures_;
  std::vector<PressureBCond> elemPressures_;
  std::vector<ContactContainer> ContactPairs_; 
  std::set<int> activeSurfs_;

  ResizeArray<SurfaceEntity*> ActiveSurfaces;  

  const std::vector<int> sampleNodeIds_;
  std::vector<std::map<int, double> > elemWeights_;

  bool mfttFlag_;

  // Disallow copy & assignment
  MeshDesc(const MeshDesc &);
  MeshDesc &operator=(const MeshDesc &);
};

std::ostream &
operator<<(std::ostream &, const MeshDesc &);

} /* end namespace Rom */

#endif /* ROM_MESHDESC_H */
