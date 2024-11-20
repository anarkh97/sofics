#ifndef ROM_MESHOUTPUT_H
#define ROM_MESHOUTPUT_H

#include <Driver.d/Attrib.h>
#include <Driver.d/EFrameData.h>
#include <Element.d/Element.h>
#include <Parser.d/AuxDefs.h>
#include <Utils.d/CompositeInfo.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include "MeshDesc.h"


class NLMaterial;

#include <ostream>
#include <iterator>
#include <string>
#include <sstream>
#include <utility>
#include <stdexcept>
#include <limits>

namespace Rom {

// Specific sections

std::ostream &
operator<<(std::ostream &, const CoordSet &);

std::ostream &
operator<<(std::ostream &, const Elemset &);

std::ostream &
operator<<(std::ostream &, const std::vector<ContactContainer> &);

std::ostream &
operator<<(std::ostream &, const ResizeArray<SurfaceEntity*>*);

std::ostream &
operator<<(std::ostream &, const SPropContainer &);

// Atoms

std::ostream &
operator<<(std::ostream &, const Attrib &);

std::ostream &
operator<<(std::ostream &, const BCond &);

std::ostream &
operator<<(std::ostream &, const EFrameData &);

std::ostream &
operator<<(std::ostream &, const FrameData &);

std::ostream &
operator<<(std::ostream &, const std::pair<int,CoefData> &);

std::ostream &
operator<<(std::ostream &, const NLMaterial &);

std::ostream &
operator<<(std::ostream &, const PressureBCond &);

// Sections from atoms

template <typename ValueType, typename TagType>
struct InputFileSectionHelper {
  static const std::string &header(TagType);
  static ValueType transformation(const ValueType &v) { return v; }
};

template <typename InputIterator, typename TagType>
class InputFileSection {
public:
  typedef typename std::iterator_traits<InputIterator>::value_type ValueType;

  const std::string &header() const {
    return InputFileSectionHelper<ValueType, TagType>::header(tag_);
  }

  InputIterator begin() const { return first_; }
  InputIterator end()   const { return last_;  }

  int j() const { return j_; }

  InputFileSection(InputIterator first, InputIterator last, TagType tag, int j) :
    first_(first), last_(last), tag_(tag), j_(j)
  {}

private:
  InputIterator first_, last_;
  TagType tag_;
  int j_;
};

template <typename InputIterator, typename TagType>
std::ostream &
operator<<(std::ostream &out, const InputFileSection<InputIterator, TagType> &source) {
  if(source.j() == 0) out << "*\n" << source.header() << "\n";
  else out << "*\n" << source.header() << " " << source.j() << "\n";
  InputIterator itEnd = source.end();
  for (InputIterator it = source.begin(); it != itEnd; ++it) {
    typedef typename InputFileSection<InputIterator, TagType>::ValueType ValueType;
    try {
      out << InputFileSectionHelper<ValueType, TagType>::transformation(*it) << "\n";
    }
    catch(std::exception& e) {
      std::cerr << "caught exception: " << e.what() << std::endl;
    }
  }

  return out;
}


struct EmptyTag {};

struct SampleNodeTag {};

template <>
inline
int
InputFileSectionHelper<int, SampleNodeTag>::transformation(const int &v) {
  return v + 1;
}

template <typename TagType>
struct InputFileSectionHelper<std::pair<const int, typename TagType::SecondType>, TagType> {
  typedef std::pair<const int, typename TagType::SecondType> ValueType;
  static const std::string &header(TagType);
  static std::string transformation(const ValueType &);
};

template <typename TagType>
std::string
InputFileSectionHelper<std::pair<const int, typename TagType::SecondType>, TagType>::transformation(const ValueType &p) {
  std::ostringstream result;
  result.precision(std::numeric_limits<double>::digits10+1);
  result << p.first + 1 << " " << TagType::valueTransformation(p.second);
  return result.str();
}

struct ElementPressureTag {
  typedef std::pair<double,int> SecondType;
  static double valueTransformation(std::pair<double,int> x) { return x.first; }
};

struct ElementWeightTag {
  typedef double SecondType;
  static double valueTransformation(double x) { return x; }
};

struct MatUsageTag {
  typedef int SecondType;
  static int valueTransformation(int i) { return i + 1; }
};

struct MatLawTag {
  typedef NLMaterial* SecondType;
  static const NLMaterial &valueTransformation(const NLMaterial *m) { return *m; }
};

struct CFrameTag {
  typedef FrameData SecondType;
  static const FrameData& valueTransformation(const FrameData &x) { return x; }
};

struct CTSurfaceTag {
  typedef FrameData SecondType; 
  static const FrameData& valueTransformation(const FrameData &x) { return x; }
};

// Convenience functions

template <typename InputIterator>
InputFileSection<InputIterator, EmptyTag>
make_section(InputIterator first, InputIterator last) {
  return InputFileSection<InputIterator, EmptyTag>(first, last, EmptyTag(), 0);
}

template <typename InputIterator, typename TagType>
InputFileSection<InputIterator, TagType>
make_section(InputIterator first, InputIterator last, TagType tag) {
  return InputFileSection<InputIterator, TagType>(first, last, tag, 0);
}

template <typename InputIterator, typename TagType>
InputFileSection<InputIterator, TagType>
make_section(InputIterator first, InputIterator last, TagType tag, int j) {
  return InputFileSection<InputIterator, TagType>(first, last, tag, j);
}

} /* end namespace Rom */

#endif /* ROM_MESHOUTPUT_H */
