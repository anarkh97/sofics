#include "MeshOutput.h"

#include "SimpleBuffer.h"

#include <Driver.d/Domain.h>

#include <Element.d/NonLinearity.d/NLMaterial.h>

#include <algorithm>

namespace Rom {

// Atoms

std::ostream &
operator<<(std::ostream &out, const Attrib &source) {
  out << source.nele + 1 << " "
      << source.attr + 1;
  if(source.cmp_attr > -1 && source.cmp_frm > -1) {
    out << " " << source.cmp_attr + 1 << " " << source.cmp_frm + 1;
  }
  return out;
}

std::ostream &
operator<<(std::ostream &out, const BCond &source) {
  if(source.dofnum == 6) { // temperature
    out << source.nnum   + 1 << " "
        << source.val;
  }
  else {
    out << source.nnum   + 1 << " "
        << source.dofnum + 1 << " "
        << source.val;
  }
  return out;
}

std::ostream &
operator<<(std::ostream &out, const PressureBCond &source) {
  out << source.elnum + 1 << " ";
  if(source.face > -1) out << "face " << source.face + 1 << " ";
  out << source.val;
  if(!source.conwepswitch) out << " off";
  return out;
}

std::ostream &
operator<<(std::ostream &out, const EFrameData &source) {
  out << source.elnum + 1 << " ";
  
  const EFrame & eframe = source.frame;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      out << eframe[i][j];
      if (i != 3 || j != 3) {
        out << " ";
      }
    }
  }

  return out;
}

std::ostream &
operator<<(std::ostream &out, const FrameData &source) {
  
  for (int i = 0; i < 9; ++i) {
    out << source.d[i];
    if (i != 8) {
        out << " ";
    }
  }

  return out;
}

std::ostream &
operator<<(std::ostream &out, const std::pair<int,CoefData> &source) {

  out << "COEF " << source.first + 1;
  for(int j=0; j<6; ++j)
    out << " " << source.second.c[6][j];
  if (source.second.coefFlag)
    out << " On";    
  out << std::endl;
  for(int i=0; i<6; ++i)
    for(int j=0; j<6; ++j) {
      out << i+1 << " " << j+1 << " " << source.second.c[i][j];
      if(!(i==5 && j==5)) out << std::endl;
    }
  return out;
}

std::ostream &
operator<<(std::ostream &out, const NLMaterial &source) {
  out.setf(std::ios_base::scientific, std::ios_base::floatfield);
  source.print(out);
  source.print2(out);
  return out;
}

// Sections

std::ostream &
operator<<(std::ostream &out, const CoordSet &source) {
  out << "NODES\n";

  const int size = const_cast<CoordSet &>(source).size();
  for (int i = 0; i < size; ++i) {
    const Node &node = const_cast<CoordSet &>(source).getNode(i);
    out << i + 1
        << " " << node.x
        << " " << node.y
        << " " << node.z;
    out << "\n";
  }

  return out;
}

std::ostream &
operator<<(std::ostream &out, const Elemset &source) {
  out << "*\nTOPOLOGY\n";

  const int size = source.last();
  for (int i = 0; i < size; ++i) {
    Element &ele = *source[i];

    out << i + 1 << " " << ele.getElementType();
    
    SimpleBuffer<int> nodes(ele.numNodes());
    ele.nodes(nodes.array());
    for (int j = 0; j < nodes.size(); ++j) {
      out << " " << nodes[j] + 1;
    }
    out << "\n";
  }

  return out;
}

std::ostream &
operator<<(std::ostream &out, const std::vector<ContactContainer> &source) {

  out << "*\nCONTACTSURFACES\n";
    for(int mc = 0; mc<source.size(); ++mc){
      out << mc + 1 << " " << source[mc].MasterId   << " "
                           << source[mc].SlaveId    << " "
                           << source[mc].MortarType << " "
                           << source[mc].NormalTol  << " "
                           << source[mc].TangentTol;
    }
    out << "\n";

  return out;
}

// this function is a f***ing work of art
std::ostream &
operator<<(std::ostream &out, const ResizeArray<SurfaceEntity*>* source) {
  int numSurf = const_cast<ResizeArray<SurfaceEntity*>* >(source)->max_size();

  int glEleNum = 1; 
  for (int isurf = 0; isurf < numSurf; ++isurf) { // loop over each surface
    // output surface number and thickness
    int surfId = (*const_cast<ResizeArray<SurfaceEntity*>* >(source))[isurf]->GetId();
    out << "*\nSURFACETOPO " << surfId;
    int divideBy = 1;
    if ((*const_cast<ResizeArray<SurfaceEntity*>* >(source))[isurf]->GetIsShellFace()){
      out << " surface_thickness " << (*const_cast<ResizeArray<SurfaceEntity*>* >(source))[isurf]->GetShellThickness();
      divideBy = 2;// specifying a surface thickness causes elements to be listed twice 
    }
    out << "\n";
    // get number of elements in this face
    int numEle = (*const_cast<ResizeArray<SurfaceEntity*>* >(source))[isurf]->GetnFaceElems()/divideBy; 
    FaceElemSet * fEleSet = (*const_cast<ResizeArray<SurfaceEntity*>* >(source))[isurf]->GetPtrFaceElemSet();
    int *glNodeNum = (*const_cast<ResizeArray<SurfaceEntity*>* >(source))[isurf]->GetPtrGlNodeIds(); 
    for(int iEle = 0; iEle < numEle; ++iEle) { // loop over each element in that surface
      int ftype = (*fEleSet)[iEle]->GetFaceElemType();
      int nn = (*fEleSet)[iEle]->nNodes();
      out << glEleNum << " " << ftype << " "; 
      for(int iNode = 0; iNode < nn; ++iNode) { // loop over each node in that element
         int localNodeId = (*fEleSet)[iEle]->GetNode(iNode); 
         out << glNodeNum[localNodeId] + 1 << " "; // map from numbering local to surface to global reduced mesh numbering
      } // win
      glEleNum++;
      out << "\n"; 
    }
  } 

  return out;
}

std::ostream &
operator<<(std::ostream &out, const SPropContainer &source) {
  out << "*\nMATERIALS\n";

  const SPropContainer::const_iterator itEnd = source.end();
  for (SPropContainer::const_iterator it = source.begin(); it != itEnd; ++it) {
    out << it->first + 1 << " ";
    const StructProp &sp = it->second;
    out << std::scientific
        << sp.A    << " "
        << sp.E    << " "
        << sp.nu   << " "
        << sp.rho  << " "
        << sp.c    << " "
        << sp.k    << " "
        << sp.eh   << " "
        << sp.P    << " "
        << sp.Ta   << " "
        << sp.Q    << " "
        << sp.W    << " "
        << sp.Ixx  << " "
        << sp.Iyy  << " "
        << sp.Izz  << " "
        << sp.ymin << " "
        << sp.ymax << " "
        << sp.zmin << " "
        << sp.zmax << " "
        << sp.betaDamp << " "
        << sp.alphaDamp << " "
        << sp.lagrangeMult << " "
        << sp.penalty << " "
        << sp.initialPenalty << " "
        << sp.funtype << " "
        << sp.type << " "
        << sp.k1 << " "
        << sp.k2 << " "
        << sp.k3 << "\n";
  }

  return out;
}

// Headers

template <>
const std::string &
InputFileSectionHelper<int, SampleNodeTag>::header(SampleNodeTag) {
  static const std::string result("SAMPLENODES");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<EFrameData, EmptyTag>::header(EmptyTag) {
  static const std::string result("EFRAMES");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<std::pair<const int,FrameData>, CFrameTag>::header(CFrameTag) {
  static const std::string result("CFRAMES");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<std::pair<const int,CoefData>, EmptyTag>::header(EmptyTag) {
  static const std::string result("COMPOSITE");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<Attrib, EmptyTag>::header(EmptyTag) {
  static const std::string result("ATTRIBUTES");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<BCond, BCond::BCType>::header(BCond::BCType tag) {
  static const std::string result[] = { "FORCES", "DISPLACEMENTS", "IDISPLACEMENTS", "IVELOCITIES", "TEMPERATURES" };
  switch (tag) {
    case BCond::Forces:
      return result[0];
    case BCond::Displacements:
      return result[1];
    case BCond::Idisplacements:
      return result[2];
    case BCond::Ivelocities:
      return result[3];
    case BCond::Temperatures:
      return result[4];
    default:
      throw std::logic_error("Unknown section tag");
  }
}

template <>
const std::string &
InputFileSectionHelper<PressureBCond, ElementPressureTag>::header(ElementPressureTag) {
  static const std::string result("PRESSURE");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<std::pair<const int, ElementWeightTag::SecondType>, ElementWeightTag>::header(ElementWeightTag) {
  static const std::string result("ELLUMP");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<std::pair<const int, MatUsageTag::SecondType>, MatUsageTag>::header(MatUsageTag) {
  static const std::string result("MATUSAGE");
  return result;
}

template <>
const std::string &
InputFileSectionHelper<std::pair<const int, MatLawTag::SecondType>, MatLawTag>::header(MatLawTag) {
  static const std::string result("MATLAW");
  return result;
}

} /* end namespace Rom */
