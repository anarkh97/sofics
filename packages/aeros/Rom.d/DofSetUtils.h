#ifndef ROM_NODEDOF_H
#define ROM_NODEDOF_H

#include <cstddef>

namespace Rom {

struct NodeDof {
  typedef int DofType;

  int nodeRank;
  DofType dofId;

  NodeDof(int n, DofType d) :
    nodeRank(n),
    dofId(d)
  {}

  NodeDof() {}
};
 
const size_t DOF_ID_COUNT = 6;
extern NodeDof::DofType DOF_ID(int iDof);

} // end namespace Rom

#endif /* ROM_NODEDOF_H */
