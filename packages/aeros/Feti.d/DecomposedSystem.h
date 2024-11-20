//
// Created by Michel Lesoinne on 6/25/18.
//

#ifndef FEM_DECOMPOSEDSYSTEM_H
#define FEM_DECOMPOSEDSYSTEM_H


#include "FetiSub.h"

class DecomposedSystem {
public:
	void preprocess();

private:
	std::vector<FetiBaseSub *> subdomains;
};


#endif //FEM_DECOMPOSEDSYSTEM_H
