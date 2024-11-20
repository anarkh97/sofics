//
// Created by Michel Lesoinne on 2/1/18.
//

#ifndef FEM_DOFINFO_H
#define FEM_DOFINFO_H

#include <vector>
#include <FetiLib/Types.h>

namespace FetiLib {

//!< Enum representing types of Degrees of Freedom.
enum class DOFType  : int {
	XDisp = 1, //!< \brief X Displacement DOF type.
	YDisp = 2, //!< \brief Y Displacement DOF type.
	ZDisp = 4, //!< \brief Z Displacement DOF type.
	XRot =  8,  //!< \brief X rotation DOF type.
	YRot = 16,  //!< \brief Y rotation DOF type.
	ZRot = 32,  //!< \brief Z rotation DOF type.
	P    = 64,  //!< \brief Pressure DOF type.
};



/** \brief DOFInfo gives for each DOF the local node number to which the DOF belongs and the type of DOF it is. */
using DOFInfo = std::vector<std::pair<local_node_index, DOFType>>;

}

#endif //FEM_DOFINFO_H
