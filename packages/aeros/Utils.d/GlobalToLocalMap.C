#include <Utils.d/GlobalToLocalMap.h>
#include <Driver.d/Communicator.h>
#include <iostream>

#ifdef MAP_MIN_MEMORY

GlobalToLocalMap::GlobalToLocalMap() : globalToLocalArray() {
}

GlobalToLocalMap::GlobalToLocalMap(int localSize, int *localToGlobalArray) {
	initialize(localSize, localToGlobalArray);
}

GlobalToLocalMap::GlobalToLocalMap(gsl::span<const int> localToGlobalArray) {
	initialize(localToGlobalArray.size(), localToGlobalArray.data());
}

GlobalToLocalMap::~GlobalToLocalMap() {
}

void
GlobalToLocalMap::initialize(int localSize, const int *localToGlobalArray) {
	if (~globalToLocalArray.empty()) globalToLocalArray.clear();
	for (int i = 0; i < localSize; ++i) {
		if (localToGlobalArray[i] < 0) continue; //HB: skip <0 values (for instance, -1 used as flag)
		globalToLocalArray[localToGlobalArray[i]] = i;
	}
}

int
GlobalToLocalMap::operator[](int glNum) const {
	auto I = globalToLocalArray.find(glNum);
	if (I != globalToLocalArray.end())
		return I->second;
	else return -1;
}

void
GlobalToLocalMap::print(bool skip) {
	std::cerr << "GlobalToLocalMap of size " << size() << ": ";
	std::map<int, int>::iterator I = globalToLocalArray.begin();
	while (I != globalToLocalArray.end()) {
		std::cerr << "(" << (*I).first << "," << (*I).second << ") ; ";
		I++;
	}
	std::cerr << std::endl;
}

int
GlobalToLocalMap::size() {
	return globalToLocalArray.size();
}

int
GlobalToLocalMap::numKeys() {
	return size();
}

template<class Scalar>
void
GlobalToLocalMap::setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID) {
	pat->setLen(localID, remoteID, 2 * size() + 3);
}

template<class Scalar>
void
GlobalToLocalMap::pack(FSCommPattern<Scalar> *pat, int localID, int remoteID) {
	FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(localID, remoteID);
	int l = size();
	sInfo.data[0] = 0;
	sInfo.data[1] = 0;
	sInfo.data[2] = l;
	int i = 0;
	std::map<int, int>::iterator I = globalToLocalArray.begin();
	while (I != globalToLocalArray.end()) {
		sInfo.data[i + 3] = (*I).first;
		sInfo.data[i + l + 3] = (*I).second;
		I++;
		i++;
	}
}

template<class Scalar>
void
GlobalToLocalMap::unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID) {
	if (~globalToLocalArray.empty()) globalToLocalArray.clear();
	FSSubRecInfo<Scalar> rInfo = pat->recData(remoteID, localID);
	int l = rInfo.data[2];
	for (int i = 0; i < l; i++)
		globalToLocalArray[rInfo.data[i + 3]] = rInfo.data[i + l + 3];
}

#else
GlobalToLocalMap::GlobalToLocalMap()
{
  min = 0;
  max = -1;
  globalToLocalArray = 0;
  nKeys = 0;
}

GlobalToLocalMap::GlobalToLocalMap(int localSize, int *localToGlobalArray)
{
 initialize(localSize, localToGlobalArray);
}

GlobalToLocalMap::~GlobalToLocalMap()
{
  if(globalToLocalArray) delete [] globalToLocalArray;
}

void
GlobalToLocalMap::initialize(int localSize, int *localToGlobalArray)
{
  if(globalToLocalArray) delete [] globalToLocalArray;
  if(localSize<=0) { min = 0; max = -1; globalToLocalArray = 0; return; }
  min = max = localToGlobalArray[0];
  if(min<0) { min = max = 0; }
  for(int i=0; i<localSize; ++i) {
	if(localToGlobalArray[i]<0) continue; //HB: skip <0 values (for instance, -1 used as flag)
	if(localToGlobalArray[i] < min) min = localToGlobalArray[i];
	if(localToGlobalArray[i] > max) max = localToGlobalArray[i];
  }
  if(min<0)
	std::cerr<<" ### PB in GlobalToLocalMap::initialize(...): min ("<<min<<") < 0"<<std::endl;
  globalToLocalArray = new int[max-min+1];
  for(int i=0; i < (max-min+1); ++i) globalToLocalArray[i] = -1;
  nKeys = 0;
  for(int i=0; i<localSize; ++i)
	if(localToGlobalArray[i]>=0) {
	  globalToLocalArray[localToGlobalArray[i]-min] = i;
	  nKeys++;
	}
}

int
GlobalToLocalMap::operator[](int glNum)
{
  if((glNum < min) || (glNum > max)) return -1;
  else return globalToLocalArray[glNum-min];
}

void
GlobalToLocalMap::print(bool skip)
{
  std::cerr << "GlobalToLocalMap of size "<<max-min+1<<": ";
  for(int i=min; i<=max; ++i){
	if(skip && (*this)[i]<0) continue;
	std::cerr <<"("<<i<<","<<(*this)[i]<<") ; ";
  }
  std::cerr << std::endl;
}

int
GlobalToLocalMap::size()
{
  return (max-min+1>0) ? max-min+1 : 0;
} 

int
GlobalToLocalMap::numKeys()
{
  return nKeys;
}

long long
GlobalToLocalMap::mem()
{
   return size()+3;
}

template<class Scalar>
void
GlobalToLocalMap::setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID)
{
  pat->setLen(localID, remoteID, size()+3);
}

template<class Scalar>
void
GlobalToLocalMap::pack(FSCommPattern<Scalar> *pat, int localID, int remoteID)
{
  FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(localID, remoteID);
  sInfo.data[0] = min;
  sInfo.data[1] = max;
  sInfo.data[2] = nKeys;
#ifdef HB_USE_MEMCPY
  memcpy(&sInfo.data[3], globalToLocalArray, size());
#else
  for(int i=0; i<size(); ++i) sInfo.data[i+3] = globalToLocalArray[i];
#endif
}

template<class Scalar>
void
GlobalToLocalMap::unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID)
{
  if(globalToLocalArray) delete [] globalToLocalArray;
  FSSubRecInfo<Scalar> rInfo = pat->recData(remoteID, localID);
  min  = rInfo.data[0];
  max  = rInfo.data[1];
  nKeys= rInfo.data[2];
  if(max-min+1<=0) { min = 0; max = -1; globalToLocalArray = 0; nKeys = 0; return; }
  globalToLocalArray = new int[max-min+1];
#ifdef HB_USE_MEMCPY
  memcpy(globalToLocalArray,&rInfo.data[3], size());
#else
  for(int i=0; i<size(); ++i)  globalToLocalArray[i] = rInfo.data[i+3];
#endif
}
#endif

template
void
GlobalToLocalMap::setCommSize(FSCommPattern<int> *pat, int localID, int remoteID);

template
void
GlobalToLocalMap::pack(FSCommPattern<int> *pat, int localID, int remoteID);

template
void
GlobalToLocalMap::unpack(FSCommPattern<int> *pat, int remoteID, int localID);

