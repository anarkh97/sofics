//
// Created by Michel Lesoinne on 3/30/18.
//

#ifndef FEM_DOFMANIP_H
#define FEM_DOFMANIP_H

#include <Utils.d/dofset.h>
#include <FetiLib/DOFInfo.h>

namespace FetiLib {
namespace Details {

/** \brief Create a DofSetArray representing the information in the DOFInfo. */
DofSetArray buildDSA(const DOFInfo &dofInfo);

/** \brief Create a DofSetArray representing the information in the DOFInfo. */
DofSetArray buildCDSA(const DOFInfo &dofInfo);

std::pair<DofSetArray, ConstrainedDSA> buildNodalDOFInfo(const DOFInfo &dofInfo);

}
}
#endif //FEM_DOFMANIP_H
