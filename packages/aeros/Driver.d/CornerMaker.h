#ifndef _CORNER_MAKER_H_
#define _CORNER_MAKER_H_

#include <vector>
#include <Feti.d/CornerSelector.h>
#include <Types.h>

class Elemset;

class SubCornerHandler : public FetiSubCornerHandler {
public:
	SubCornerHandler(gl_sub_idx sub, int nn, CoordSet &n, Elemset &ele, Connectivity &nTn, DofSetArray &d,
	                 Connectivity &sh, const std::vector<lc_node_idx> &nsb, ConstrainedDSA *c_dsa, FetiBaseSub *_subPre);
};

#endif
