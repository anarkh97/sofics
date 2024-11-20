//
// Created by Michel Lesoinne on 4/3/18.
//

#include "SubImpl.h"

namespace FetiLib {
SubImpl::SubImpl(VectorReference<const global_node_index> globalNodeIndices, const SharedNodes &sharedNodes)
		: glNodes({globalNodeIndices.begin(), globalNodeIndices.end()})
{

}

}