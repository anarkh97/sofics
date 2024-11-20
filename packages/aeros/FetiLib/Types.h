//
// Created by Michel Lesoinne on 3/28/18.
//

#ifndef FETILIB_TYPES_H
#define FETILIB_TYPES_H

namespace FetiLib {

/// \brief An integer type for node indices in numbering within a single subdomain.
using local_node_index = int;
/// \brief An integer type for DOF indices in numbering within a single subdomain.
using local_dof_index = int;
/// \brief An integer type for node indices in a global numbering for the complete problem.
using global_node_index = int;
/// \brief An integer type for subdomain indices in a global numbering for the complete problem.
using global_subdomain_index = int;
/// \brief An integer type for subdomain indices in a global numbering for the complete problem.
using local_subdomain_index = int;

} // namespace FetiLib

#endif //FETILIB_TYPES_H
