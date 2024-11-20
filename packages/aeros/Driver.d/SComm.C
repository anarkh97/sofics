#include <iostream>
#include <Driver.d/SComm.h>

SComm::SComm(int nN, std::vector<gl_sub_idx> subIds, std::vector<lc_sub_idx> ids, std::unique_ptr<Connectivity> con) :
	subNums(std::move(subIds)), remoteId(std::move(ids)), sharedNodes{std::move(con)}
{
	numNeighb = nN;

	// type specific lists
	numDofType = 8;
}

void
SComm::deleteTypeSpecificList(DofType type)
{
	if(SubNums.size() >= type)
		SubNums[type].clear();
	TypeMap.clear();
	sharedDOFsPerType.clear();
	NumNeighb[type] = 0;
}

void
SComm::setEdgeNeighb(int _numEdgeNeighb, std::vector<bool> _isEdgeNeighb)
{
	numEdgeNeighb = _numEdgeNeighb;
	isEdgeNeighb = std::move(_isEdgeNeighb);
}

void
SComm::setTypeSpecificList(DofType type,
                           std::vector<gl_sub_idx> _subNums,
                           std::unique_ptr<Connectivity> _sharedDOFs)
{
	if (sharedDOFsPerType.size() == 0) {
		NumNeighb.assign(numDofType, 0);
		SubNums.assign(numDofType, std::vector<gl_sub_idx>{});
		sharedDOFsPerType.assign(numDofType, Connectivity{});
	}
	NumNeighb[type] = _sharedDOFs->csize();
	SubNums[type] = _subNums;
	sharedDOFsPerType[type] = std::move(*_sharedDOFs);
}

std::vector<int>
SComm::mergeTypeSpecificLists()
{
	int i, j, k, l;
	// step 1. initialize combined list to type 0
	std::unique_ptr<Connectivity> allSharedDOFs = std::make_unique<Connectivity>(sharedDOFsPerType[0]);
	std::vector<int> allSubNums(SubNums[0].begin(), SubNums[0].begin() + NumNeighb[0]);
	TypeMap.resize(numDofType);
	TypeMap[0].resize(sharedDOFsPerType[0].numConnect());

	// step 2. loop over all other types and add to combined list
	for (i = 1; i < int(all); ++i) {
		if (NumNeighb[i] > 0) {
			allSharedDOFs->combine(&sharedDOFsPerType[i], allSubNums, SubNums[i]);
			TypeMap[i].resize(sharedDOFsPerType[i].numConnect());
		} else {
			TypeMap[i].clear();
		}
	}
	//int allNumNeighb = allSharedDOFs->csize();

	// step 3. make boundDofFlag and TypeMaps
	std::vector<int> boundDofFlag(allSharedDOFs->numConnect());
	int count = 0;
	for (i = 0; i < allSharedDOFs->csize(); ++i) {
		int neighb = allSubNums[i];
		for (j = 0; j < int(all); ++j) {
			for (k = 0; k < sharedDOFsPerType[j].csize(); ++k) {
				if (SubNums[j][k] == neighb) {
					for (l = 0; l < sharedDOFsPerType[j].num(k); ++l) {
						TypeMap[j][sharedDOFsPerType[j].offset(k) + l] = count;
						boundDofFlag[count++] = j;
					}
					break;
				}
			}
		}
	}
	setTypeSpecificList(all, allSubNums, std::move(allSharedDOFs) );

	for (i = int(all); i < numDofType; ++i)
		TypeMap[i].clear();

	return boundDofFlag;
}

void SComm::print(DofType t)
{
	std::cerr << "NumNeighb = " << NumNeighb[t] << std::endl;
	std::cerr << "SubNums = ";
	for (int i = 0; i < NumNeighb[t]; ++i) std::cerr << SubNums[t][i] << " ";
	std::cerr << std::endl;
	std::cerr << "sharedDOFsPerType = \n";
	sharedDOFsPerType[t].print();
}