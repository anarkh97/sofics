#include "BasisFileStream.h"

#include <stdexcept>

namespace Rom {

template <int DOFS_PER_NODE>
BasisInputStream<DOFS_PER_NODE>::BasisInputStream(const std::string &fileName, const VecNodeDofConversion<DOFS_PER_NODE> &converter) :
  file_(fileName),
  isValid_(true),
  converter_(converter),
  buffer_(file_.nodeIdBegin(), file_.nodeIdEnd())
{}

template <int DOFS_PER_NODE>
BasisOutputStream<DOFS_PER_NODE>::BasisOutputStream(const std::string &fileName, const VecNodeDofConversion<DOFS_PER_NODE> &converter, bool restart) :
  file_(fileName, converter.dofSetNodeCount(), restart, DOFS_PER_NODE), // Conservative, potentially overallocating
  converter_(converter),
  buffer_(file_.nodeCount())
{}

template
BasisInputStream<6>::BasisInputStream(const std::string &fileName, const VecNodeDofConversion<6> &converter);

template
BasisInputStream<1>::BasisInputStream(const std::string &fileName, const VecNodeDofConversion<1> &converter);

template
BasisOutputStream<6>::BasisOutputStream(const std::string &fileName, const VecNodeDofConversion<6> &converter, bool restart);

template
BasisOutputStream<1>::BasisOutputStream(const std::string &fileName, const VecNodeDofConversion<1> &converter, bool restart);

} /* end namespace Rom */
