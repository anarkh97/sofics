#include <Driver.d/CornerMaker.h>
#include <Element.d/Element.h>

SubCornerHandler::SubCornerHandler(gl_sub_idx sub, int nn, CoordSet &n, Elemset &eles, Connectivity &nTn, DofSetArray &d,
                                   Connectivity &sh, const std::vector<lc_node_idx> &nsb, ConstrainedDSA *c_dsa,
                                   FetiBaseSub *_subPre)
		: FetiSubCornerHandler(sub, nn, n, nTn, d, sh, nsb, c_dsa,
		                       _subPre) {
	for(int i=0; i<eles.last(); ++i) {
		int *nodes = eles[i]->nodes();
		for(int j=0; j<eles[i]->numNodes(); ++j)
			if(eles[i]->isRotMidSideNode(j)) isRotMidSideNode[nodes[j]] = true;
		delete [] nodes;  // Element::nodes() should always return a copy to be deleted by caller
		if(eles[i]->dim() > dim) dim = eles[i]->dim();
	}

	glSafe.assign(nnodes, false);
	// mark the safe nodes by setting glSafe to true if node is touched by at least one safe element
	for(int iEle = 0; iEle < eles.last(); ++iEle) {
		if((eles[iEle]->isSafe() && !eles[iEle]->isPhantomElement()) || allSafe) {
			int nnd = eles[iEle]->numNodes();
			int *nodes = new int[nnd];
			eles[iEle]->nodes(nodes);
			for(int iNode = 0; iNode < nnd; ++iNode) {
				int n = nodes[iNode];
				glSafe[n] = true;
			}
			delete [] nodes;
		}
	}
}
