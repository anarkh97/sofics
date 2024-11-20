#include <Driver.d/Sower.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/DistHelper.h>
#include <algorithm>

#ifdef SOWER_SURFS

#include <Utils.d/resize_array.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>

#endif

extern GeoSource *geoSource;

//#define SOWER_DEBUG
//#define SOWER_DEBUG_SURFS

Sower::Sower(Connectivity *subToElem, Elemset &eset, int nClus, ResizeArray<SurfaceEntity *> *Surfs,
             Connectivity *cpuToSub)
{
//Warning, subToElem is packed and eset is not ...
	subToClus = 0;
	tocRead = false;
#ifdef SOWER_DEBUG
	std::cerr << " ** Sower Constructor : " << std::endl;
	std::cerr << " ** subdomain to elements connectivity : " << std::endl;
	subToElem->print();
#endif
	Connectivity *eTos = subToElem->alloc_reverse();
#ifdef SOWER_DEBUG
	std::cerr << " ** elements to subdomains connectivity : " << std::endl;
	eTos->print();
#endif
	Connectivity *sTos = subToElem->transcon(eTos);
#ifdef SOWER_DEBUG
	std::cerr << " ** subdomains to subdomains connectivity : " << std::endl;
	sTos->print();
#endif
	Connectivity *eToN = new Connectivity(eset.asSet());
#ifdef SOWER_DEBUG
	std::cerr << " ** elements to nodes connectivity : " << std::endl;
	eToN->print();
#endif
	Connectivity *sToN = subToElem->transcon(eToN);
#ifdef SOWER_DEBUG
	std::cerr << " ** subdomains to nodes connectivity : " << std::endl;
	sToN->print();
#endif
	Connectivity *nToS = sToN->alloc_reverse();
#ifdef SOWER_DEBUG
	std::cerr << " ** nodes to subdomain connectivity : " << std::endl;
	nToS->print();
#endif
	int nSub = subToElem->csize();
	typedef long long superlong;

	// read clusToSub from file if possible   JAT 021215
	const char *clusName = "CLUSMAP";
	FILE *f = fopen(clusName, "r");
	if (f) {
		int numCLUS;
		int r1 = fscanf(f, "%d", &numCLUS);
		fprintf(stderr, " ... Reading Cluster Map from file %s, numCLUS = %d ... \n", clusName, numCLUS);
		if (numCLUS != nClus) {
			fprintf(stderr, " *** ERROR: CLUSMAP file is for %d clusters\n", numCLUS);
			exit(-1);
		}
		std::vector<size_t> clusp(nClus + 2);
		std::vector<int> subs(nSub + 1);
		int nSubFile = 0;
		int i, j, k = 0, m, n;
		clusp[0] = k;
		for (i = 0; i < nClus; i++) {
			int r2 = fscanf(f, "%d", &n);
			for (j = 0; j < n; j++) {
				int r3 = fscanf(f, "%d", &m);
				if (m > nSubFile) nSubFile = m;
				if (k == nSub) {
					fprintf(stderr, " *** ERROR: CLUSMAP file has too many subdomains\n");
					exit(-1);
				}
				subs[k++] = m - 1;
			}
			clusp[i + 1] = k;
		}
		fclose(f);
		if (nSubFile != nSub) {
			fprintf(stderr, " *** ERROR: CLUSMAP file is for %d subdomains\n", nSubFile);
			exit(-1);
		}
		clusToSub = std::make_unique<Connectivity>(nClus, std::move(clusp), std::move(subs));
	} else {
		// for new clustering method, all the subdomains on a cpu must be in the same cluster
		if (nClus == 1) {
			clusToSub = std::make_unique<Connectivity>(nClus, nSub);
		} else {
			int nCpu = cpuToSub->csize();
			if (nCpu == nClus) {
				clusToSub = std::make_unique<Connectivity>(*cpuToSub);
			} else {
				if (nClus > nCpu) {
					filePrint(stderr, "ERROR: number of clusters greater than number of CPUs is not supported\n");
					exit(-1);
				}
				Connectivity *subToCpu = cpuToSub->alloc_reverse();
				Connectivity *cpuToCpu = cpuToSub->transcon(subToCpu);
				long long *sizes = new superlong[nCpu];
				for (int i = 0; i < nCpu; i++) {
					sizes[i] = 0;
					for (int j = 0; j < cpuToSub->num(i); ++j) {
						sizes[i] += subToElem->num((*cpuToSub)[i][j]);
					}
				}
				Connectivity *clusToCpu = greedy(cpuToCpu, subToCpu, cpuToSub, sizes, nCpu, nClus);
				delete[] sizes;
				delete cpuToCpu;
				delete subToCpu;
				clusToSub = std::make_unique<Connectivity>(clusToCpu->transcon(*cpuToSub));
				delete clusToCpu;
			}
		}
	}
	//addParentToChildData<Elemset*,Connectivity*>(ELEMENTS,CLUSTER,0,&eset,clusToSub);

#ifdef SOWER_DEBUG
	std::cerr << " ** cluster to subdomain connectivity : " << std::endl;
	clusToSub->print();
#endif
	// PJSA: write clusToSub file
	BinFileHandler file(subdomains_.c_str(), "w");
	subToClus = std::make_unique<Connectivity>(clusToSub->reverse());
	subToClus->write(file);
	// PJSA: write subToElem & subToNode files for each cluster
	// new method: write enough to reconstruct implict (sparse) connectivities
	double t1 = getTime();
	int csize = std::max(eTos->csize(), nToS->csize());
	int *allptr = new int[csize + 1];
	for (int i = 0; i < csize + 1; ++i) allptr[i] = i;
	//Connectivity *clusToElem = clusToSub->transcon(subToElem);
	Connectivity *clusToNode = clusToSub->transcon(sToN);
	Connectivity *clusToSub2 = clusToNode->transcon(nToS);
	for (int i = 0; i < nClus; ++i) {
		std::ostringstream oss;
		oss << decomposition_ << i + 1;
		BinFileHandler file2(oss.str().c_str(), "w");
		Connectivity *cnodeToNode = new Connectivity(clusToNode->num(i), allptr, (*clusToNode)[i].data(), 0);
		Connectivity *cnodeToSub = cnodeToNode->transcon(nToS);
		Connectivity *subToCnode = cnodeToSub->alloc_reverse();
		Connectivity *subToNode_i = subToCnode->transcon(cnodeToNode);
		Connectivity *csubToSub2 = new Connectivity(clusToSub2->num(i), allptr, (*clusToSub2)[i].data(), 0);
		Connectivity *csubToNode = csubToSub2->transcon(subToNode_i);
		csubToNode->sortTargets();
		csubToSub2->write(file2);
		csubToNode->writeg(file2);
		cnodeToNode->writeg(file2);
		cnodeToSub->write(file2);
		delete cnodeToNode;
		delete cnodeToSub;
		delete subToCnode;
		delete subToNode_i;
		delete csubToSub2;
		delete csubToNode;
	}
	delete[] allptr;
	//delete clusToElem;
	delete clusToNode;
	delete clusToSub2;
	//std::cerr << "elapsed time = " << (getTime()-t1)/1000 << std::endl;

	// PJSA: write connectivity file
	BinFileHandler file3(connectivity_.c_str(), "w");
	clusToSub->write(file3);
#ifdef SUBTOSUBINFILE
	Connectivity *subToSub = sToN->transcon(nToS);  // JAT 110116
	subToSub->write(file3);
	delete subToSub;
#else
	file3.write(&nSub, 1); // JAT 110116
#endif
	int nGlobNodes = nToS->csize();
	file3.write(&nGlobNodes, 1);
	int nGlobElems = eTos->csize();
	file3.write(&nGlobElems, 1);

	nCluster = clusToSub->csize(); /* true number of cluster */
	entries[ELEMENTS_TYPE] = new GenDataStruct<Elemset *, ElemsetIO>(&eset,/*clusToSub->transcon(*/*subToElem/*)*/);

	numSubdomains = subToElem->csize();

	nSurfaces = 0;
	Surfaces = Surfs;
	if (Surfaces) { // count real number of surfaces
		for (int i = 0; i < Surfaces->max_size(); i++)
			if ((*Surfaces)[i]) nSurfaces++;
	}
	if (nSurfaces) std::cerr << " ... " << nSurfaces << " surfaces have been defined in input file" << std::endl;
	// Write (global) number of surfaces into connectivity files
	file3.write(&nSurfaces, 1);

	delete eTos;
	delete sTos;
	delete eToN;
	delete sToN;
	delete nToS;
}

Sower::~Sower()
{
	for (std::map<TypeTag, DataStruct *>::iterator it = entries.begin(); it != entries.end(); it++)
		if ((*it).second != 0)
			delete ((*it).second);
}

void Sower::printDebug() const
{
	std::cerr << " def : END=0,NODES=1,ELEMENTS=2,ATTRIBUTES=3,FORCES=4,MATERIALS=5,DISPLACEMENTS=6,\
	     EFRAMES=7,COMPOSITE=8,CFRAMES=9,CLUSTER=10 etc " << std::endl;
	for ( const auto &it : entries ) {
		std::cerr << std::endl << " Cluster to " << it.first << " connectivity : " << std::endl;
		it.second->getClusterToData(clusToSub.get()).print();
		std::cerr << std::endl << " SubDomain to " << it.first << " connectivity : " << std::endl;
		it.second->getSubToData().print();
	}
}

void Sower::nomask(Connectivity *mn, int *maskno, int nelem)
{
	int i;
	for (i = 0; i < mn->num(nelem); i++) maskno[(*mn)[nelem][i]]--;
}

Connectivity *
Sower::greedy(Connectivity *etoe, Connectivity *ntoe,
              Connectivity *eton, long long *sizes,
              int decompNumSubdomains, int expectedNumSubdomains)
{
	int i, minw, minw2;
	int node, node2, ntot, nelem, ptr1, ptr2;
	int *maskno, *elsub, *me, *color;
	int *newelsub;
	int nsub = expectedNumSubdomains;

	// Compute total profile size
	int totalSize = 0;
	for (i = 0; i < decompNumSubdomains; ++i)
		totalSize += sizes[i];

	double percent = 0.95;
	int exact_num = etoe->numNonZeroP();
	int numele = etoe->csize();
	int numnod = ntoe->csize();
	elsub = new int[nsub + 1];
	newelsub = new int[nsub + 2];
	me = new int[exact_num + 1];
	color = new int[numele];
	maskno = new int[numnod];

	//  initialize newelsub as a pointer into me
	newelsub[0] = 0;

	int div = totalSize / nsub;
	int rem = totalSize % nsub;
	// fprintf(stderr,"Per Cluster %d Remainder %d\n",div,rem);

	int div1 = exact_num / nsub;
	int rem1 = exact_num % nsub;

	// Compute an a priori balanced distribution of element numbers
	for (i = 0; i < nsub; i++) {
		elsub[i] = (i < rem1) ? div1 + 1 : div1;
	}

	for (i = 0; i < nsub; i++) {
		newelsub[i + 1] = newelsub[i] + elsub[i];
	}

	for (i = 0; i < nsub; i++) {
		elsub[i] = (i < rem) ? div + 1 : div;
	}

	for (i = 0; i < numele; i++) color[i] = -1;

	for (i = 0; i < numnod; i++) maskno[i] = ntoe->num(i);

	ptr1 = 0;
	ptr2 = 0;
	int numsub = 0;
	ntot = 0;

	while (numsub < nsub && ntot <= percent * elsub[numsub]) {

		// Locate a starting (interface if possible) node
		// The part of code is in O(Numnod). It could be done in
		// O(I) (I = interface size) by handling properly the list of
		// interface nodes.

		minw = minw2 = 32000000;
		node = -1;
		for (i = 0; i < numnod; i++) {
			if (maskno[i] == 0) continue;
			if (maskno[i] == ntoe->num(i)) {
				if (maskno[i] < minw2) {
					minw2 = maskno[i];
					node2 = i;
				}
			} else {
				if (maskno[i] < minw) {
					minw = maskno[i];
					node = i;
				}
			}
		}

		if (node < 0) node = node2;

		// Initialize the list of elements atached to the starting node
		ptr1 = ptr2;

		for (i = 0; i < ntoe->num(node) && ntot < percent * elsub[numsub]; i++) {
			// find an element attached to "node"
			nelem = (*ntoe)[node][i];
			if (color[nelem] >= 0) continue;
			color[nelem] = numsub;
			me[ptr2] = nelem;
			ptr2++;
			ntot += sizes[nelem]; // MODIFICATION
			// reduce mask value for all nodes connected to this element
			nomask(eton, maskno, nelem);
		}

		/* RECURSIVELY ADD TO LIST NEW ADJACENT ELEMENTS */

		while (ptr2 > ptr1 && ntot < percent * elsub[numsub]) {
			for (i = 0; i < etoe->num(me[ptr1]) && ntot < percent * elsub[numsub]; i++) {
				nelem = (*etoe)[me[ptr1]][i];
				if (color[nelem] >= 0) continue;
				color[nelem] = numsub;
				me[ptr2] = nelem;
				ntot += sizes[nelem]; // MODIFICATION
				ptr2++;
				nomask(eton, maskno, nelem);
			}
			ptr1++;
		}

		// MODIFICATION
		if (ntot >= percent * elsub[numsub] && ntot <= elsub[numsub] / percent) {
			numsub++;
			ntot = 0;
		}
	}

	for (i = 0; i < numele; i++) {
		if (etoe->num(i) && color[i] < 0) {
			me[ptr2++] = i;
			color[i] = nsub;
		}
	}

	// delete [] maskno;
	// delete [] color;
	// delete [] elsub;

	Connectivity *cToS = new Connectivity(nsub, newelsub, me);

	return cToS;
}

bool
ObjectOrdering::operator()(int n1, int n2)
{
	int nt1 = 0, nt2 = 0, t1, t2;
	do {
		while (nt1 < objToSub->num(n1) && subIsInClus[(*objToSub)[n1][nt1]] == false)
			nt1++;
		while (nt2 < objToSub->num(n2) && subIsInClus[(*objToSub)[n2][nt2]] == false)
			nt2++;
		t1 = (nt1 < objToSub->num(n1)) ? (*objToSub)[n1][nt1] : -1;
		t2 = (nt2 < objToSub->num(n2)) ? (*objToSub)[n2][nt2] : -1;
		if (t1 < t2)
			return (true);
		if (t2 < t1)
			return (false);
		nt1++;
		nt2++;
	} while (t1 >= 0);
	return (n1 < n2);
}

struct DataEntry {
	TypeTag dataType;
	INT64BIT offset;
};

void Sower::write()
{
	for (int currentClusNum = 1; currentClusNum <= nCluster; currentClusNum++) {
		// opening cluster file
		char filename[FILENAME_LENGTH];
		char clusterNumStr[11] = {'\0'};
		sprintf(clusterNumStr, "%u", currentClusNum);
		strcpy(filename, clusterData_.c_str());
		strcat(filename, clusterNumStr);
#ifdef SOWER_DEBUG
		std::cout << " ** Writing to file " << filename << std::endl;
#endif
		BinFileHandler file(filename, "wb");

		// counting number of data types in cluster
		int dataTypesNumber = 0;
		std::map<TypeTag, DataStruct *>::iterator it = entries.begin();
		for (; it != entries.end(); it++) {
			if ((*it).second->getClusterToData(clusToSub.get()).num(currentClusNum - 1) > 0) {
				dataTypesNumber++;
			}
		}
		if (nSurfaces) dataTypesNumber++; //HB: make room for surfaces header/toc
#ifdef SOWER_DEBUG
		std::cout << " ** " << dataTypesNumber << " different data types will be written to this file." << std::endl;
#endif

		// saving space for the table of content
		// Type of data ID / pointer to beginning of data / pointer to associated rangeset
		int tocSize = (sizeof(TypeTag) + sizeof(INT64BIT) + sizeof(INT64BIT)) * (dataTypesNumber + 1);
		INT64BIT nextEntryOffset = file.tell() + tocSize;       /** offset in the data **/
		INT64BIT tocCurrentOffset = file.tell();  /** offset in the table of content **/
		INT64BIT curRangeSetLocation = 0;
		// now writing each datatype
		it = entries.begin();
		for (; it != entries.end(); it++) {
			// oh no ! this datatype is not in this cluster
			if ((*it).second->getClusterToData(clusToSub.get()).num(currentClusNum - 1) <= 0)
				continue;
			// updating toc with entry for this datatype
			file.seek(tocCurrentOffset);

			file.write(&(*it).first, 1);
			file.write(&nextEntryOffset, 1);
			tocCurrentOffset = file.tell();

			// preparing to write all the data of this datatype
			file.seek(nextEntryOffset);
			(*it).second->write(currentClusNum, clusToSub.get(), numSubdomains, file, curRangeSetLocation);
			nextEntryOffset = file.tell();

			// adding location of offset of the Rangesets for this datatype to T O C
			file.seek(tocCurrentOffset);
			file.write(&curRangeSetLocation, 1);
			tocCurrentOffset = file.tell();

		}

		//HB: create & write toc, range set & data for surfaces
		if (nSurfaces) { writeSurfaces(file, tocCurrentOffset, nextEntryOffset, curRangeSetLocation); }

		// adding 0 entry to TOC to mark its end
		INT64BIT zero = 0;

		int zeroi = END_TYPE;
		file.seek(tocCurrentOffset);
		file.write(&zeroi, 1);
		file.write(&zero, 1);
		file.write(&zero, 1);
		// cluster file closed by BinFileHandler destructor

	} // end cluster loop
}

size_t
DataStruct::write(int clusNumber, Connectivity *clusToSub, int numSubdomains,
                  BinFileHandler &file, INT64BIT &curRangeSetLocation)
{
	std::map<int, RangeSet *> rangeSet;  // where is each subdomain's data written ?

	if (clusterToData.csize() == 0) clusterToData = clusToSub->transcon(subToData);
	Connectivity *objToSub = subToData.alloc_reverse();
	bool *subIsInClus = new bool[numSubdomains];
	for (int i = 0; i < numSubdomains; i++)
		subIsInClus[i] = false;
	for (int j = 0; j < clusToSub->num(clusNumber - 1); j++) {
		subIsInClus[((*clusToSub)[clusNumber - 1][j])] = true;
		rangeSet[((*clusToSub)[clusNumber - 1][j])] = new RangeSet();
	}

	// sorting
	// the oject ids will be in an optimum order to be written
	objToSub->sortTargets();
	ObjectOrdering order(objToSub, subIsInClus);
	auto ctd = clusterToData[clusNumber - 1];
	std::sort(ctd.begin(), ctd.end(), order);
	// now let's write each object using its objectID
	std::list<int> currentSubs; // list of subdomains currently been written -- for use with rangeset
	int k;
	for (k = 0; k < clusterToData.num(clusNumber - 1); k++) {
		int curObjID = clusterToData[clusNumber - 1][k];
		// range set
		std::list<int> nextCurrentSubs;
		auto listOfSubsObjIsIn = (*objToSub)[curObjID];
		int numOfSubsObjIsIn = objToSub->num(curObjID);

		for (int i = 0; i < numOfSubsObjIsIn; ++i) {
			if (subIsInClus[listOfSubsObjIsIn[i]]) { // sub is in cluster
				// std::list<int>::iterator it = find( currentSubs.begin(), currentSubs.end(), listOfSubsObjIsIn[i]);
				std::list<int>::iterator it = std::find(currentSubs.begin(), currentSubs.end(),
				                                        listOfSubsObjIsIn[i]); // PJSA: for sgi intel
				if (currentSubs.empty() || it == currentSubs.end()) { // this subdomain just appeared
					rangeSet[listOfSubsObjIsIn[i]]->start(file.tell(), curObjID);
					nextCurrentSubs.push_back(listOfSubsObjIsIn[i]); // this range might have to be ended at the next k
				} else { // this subdomain stays
					currentSubs.remove(listOfSubsObjIsIn[i]);
					nextCurrentSubs.push_back(listOfSubsObjIsIn[i]);
				}
			}
		}
		// all subdomains that are still in currentSubs have disapeared ! time to end them
		for (std::list<int>::iterator it = currentSubs.begin(); it != currentSubs.end(); ++it) {
			rangeSet[(*it)]->end(clusterToData[clusNumber - 1][k - 1]);
		} // can never be called for k = 0
		currentSubs.clear();
		currentSubs.merge(nextCurrentSubs);
		// end of range set

		// actual writing of data
		//file.write(&curObjID, 1);
		writeObjectData(clusterToData[clusNumber - 1][k], file, curObjID);
	}
	// all subdomains that are still in currentSubs have disapeared ! time to end them
	for (std::list<int>::iterator it = currentSubs.begin(); it != currentSubs.end(); ++it)
		rangeSet[(*it)]->end(clusterToData[clusNumber - 1][k - 1]);
	curRangeSetLocation = file.tell();
	file.seek(file.tell() + sizeof(int));
	// writing range set for this datatype
	int realNumOfSubs = 0;
	for (std::map<int, RangeSet *>::iterator it = rangeSet.begin(); it != rangeSet.end(); ++it) {
		if (!(*it).second->empty()) {
			++realNumOfSubs;
			int subNum = (*it).first;
			file.write(&subNum, 1);
#ifdef SOWER_DEBUG
			(*it).second->print();
#endif
			(*it).second->write(file);
		}
	}
	INT64BIT tmp = file.tell();
	file.seek(curRangeSetLocation);
	file.write(&realNumOfSubs, 1);
	file.seek(tmp);
	return 0;
}

#include <Driver.d/GeoSource.h>

size_t ElemsetIO::write(Elemset *eset, int index, BinFileHandler &file, int curObjID)
{
	if ((*eset)[index]) {//test JF for gap in elemset
		int *nodes = (*eset)[index]->allNodes();
		int numNodes = (*eset)[index]->numNodes();
		int elType = (*eset)[index]->getElementType(); // PJSA
		int globalID = ((*geoSource).localToGlobalElementsNumber)[curObjID];          // translation to global coords
#ifdef SOWER_DEBUG
		std::cerr << "writing element index = " << index << " curObjID = " << curObjID << " globalID = " << globalID << std::endl;
#endif
		file.write(&globalID, 1);
		file.write(&elType, 1);
		file.write(&numNodes, 1);
		file.write(nodes, numNodes);
		double pressure = ((*eset)[index]->getPressure() != NULL) ? (*eset)[index]->getPressure()->val : 0.0;
		file.write(&pressure, 1);
		std::vector<double> preload = (*eset)[index]->getPreLoad();
		int preload_size = preload.size();
		file.write(&preload_size, 1);
		file.write(&(preload[0]), preload_size);
	} else {
		std::cerr << "Sower.h, void element in eset. index = " << index << std::endl;
		int elType = -1;
		file.write(&elType, 1);
	}
	return 0;
}

//HB
void Sower::writeSurfaces(BinFileHandler &file, INT64BIT &tocCurrentOffset, INT64BIT &nextEntryOffset,
                          INT64BIT &curRangeSetLocation)
{
#ifdef SOWER_DEBUG_SURFS
	std::cerr<<" ** Write surfaces"<<std::endl;
#endif
	// updating toc with entry for this datatype
	file.seek(tocCurrentOffset);
	int surfTagType = SURFACES;

	file.write(&surfTagType, 1); // i.e. DataType
	file.write(&nextEntryOffset, 1);
	tocCurrentOffset = file.tell();

	// preparing to write all the data of this datatype
	file.seek(nextEntryOffset);
	std::map<int, RangeSet *> surfRangeSet;
	for (int isurf = 0; isurf < nSurfaces; isurf++) { //assume the surface are "packed" in the Surfaces array
#ifdef SOWER_DEBUG_SURFS
		std::cerr<<" ** Write surface "<< isurf<<std::endl;
#endif
		surfRangeSet[isurf] = new RangeSet();
		surfRangeSet[isurf]->start(file.tell(), isurf);
		file.write(&isurf, 1);
		(*Surfaces)[isurf]->WriteSower(file);
		surfRangeSet[isurf]->end(isurf);
	}
	curRangeSetLocation = file.tell();
	file.seek(file.tell() + sizeof(int));
	// writing range set for this datatype
	int realNumOfSurfs = 0;
	for (std::map<int, RangeSet *>::iterator it = surfRangeSet.begin(); it != surfRangeSet.end(); ++it) {
		if (!(*it).second->empty()) {
			++realNumOfSurfs;
			int surfNum = (*it).first;
			file.write(&surfNum, 1);
#ifdef SOWER_DEBUG_SURFS
			(*it).second->print();
#endif
			(*it).second->write(file);
		}
	}
	INT64BIT tmp = file.tell();
	file.seek(curRangeSetLocation);
	file.write(&realNumOfSurfs, 1);
	file.seek(tmp);

	nextEntryOffset = file.tell();
	// adding location of offset of the Rangesets for this datatype to T O C
	file.seek(tocCurrentOffset);
	file.write(&curRangeSetLocation, 1);
	tocCurrentOffset = file.tell();
}
