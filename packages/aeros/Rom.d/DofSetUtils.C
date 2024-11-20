#include "DofSetUtils.h"

#include <Utils.d/dofset.h>
#include <Utils.d/SolverInfo.h>

extern SolverInfo &solInfo;

namespace Rom {

NodeDof::DofType DOF_ID(int iDof) {
  switch(iDof) {
    default :
    case 0 :
      return (solInfo.soltyp == 2) ? DofSet::Temp : DofSet::Xdisp;
    case 1 :
      return DofSet::Ydisp;
    case 2 :
      return DofSet::Zdisp;
    case 3 :
      return DofSet::Xrot;
    case 4 :
      return DofSet::Yrot;
    case 5 :
      return DofSet::Zrot;
  }
}

} // end namespace Rom
