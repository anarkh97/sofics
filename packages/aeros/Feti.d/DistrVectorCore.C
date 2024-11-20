#include <Feti.d/DistrVector.h>

void
DistrInfo::initialize()
{
  domLen = nullptr;
  subLen = nullptr;
  subOffset = nullptr;
  threadOffset = nullptr;
  threadLen = nullptr;
  masterFlag = nullptr;
}

DistrInfo::DistrInfo(int ns)
{
 initialize();
 numDom = ns;
 domLen = new int[numDom];
}

DistrInfo::~DistrInfo()
{
	delete [] domLen;
	delete [] masterFlag;
	delete [] subOffset;
	delete [] threadOffset;
	delete [] threadLen;
}

void DistrInfo::setMasterFlag()
{
  // this version for internal DistrInfo's or if interface values are weighted
  masterFlag = new bool[len];
  int i;
  for(i = 0; i < len; ++i)
    masterFlag[i] = true;
}

void
DistrInfo::computeOffsets()
{
  if(subOffset)
    return;
  int iSub, iThread;

  subOffset = new int[numLocSub];
  int cOffset = 0;
  for(iSub = 0; iSub < numLocSub; ++iSub) {
    subOffset[iSub] = cOffset;
    cOffset += subLen[iSub];
  }
  len = cOffset;
 
  // XML NEED TO UPDATE FOR OPENMP
  numLocThreads = threadManager->numThr();
  threadOffset = new int[numLocThreads+1];
  threadLen = new int[numLocThreads+1]; // PJSA: for DistVec class in Math.d
  int nRemain = len%numLocThreads;
  int npt = len/numLocThreads;
  for(iThread = 0; iThread < numLocThreads+1; ++iThread) {
    threadOffset[iThread] = (iThread < nRemain ? iThread : nRemain)
       +iThread*npt;
    threadLen[iThread] = (iThread < nRemain ? 1 : 0) + npt; // PJSA: for DistVec class in Math.d
  }
  threadOffset[0] = 0;
}

void
DistrInfo::recomputeOffsets()
{
  if(subOffset) { delete [] subOffset; subOffset = nullptr; }
  if(threadOffset) { delete [] threadOffset; threadOffset = nullptr; }
  if(threadLen) { delete [] threadLen; threadLen = nullptr; }

  computeOffsets();
}

int
DistrInfo::masterLen() const
{
  GenDistrVector<int> toto(*this);
  toto = 1;
  return toto.sqNorm();
}

bool
DistrInfo::operator==(const DistrInfo& other) const
{
  if((len != other.len) || (numLocSub != other.numLocSub)) return false;
  for(int i=0; i<numLocSub; ++i)
    if(subLen[i] != other.subLen[i]) return false;
  return true;
}

bool
DistrInfo::operator!=(const DistrInfo& other) const
{
  return !(*this == other);
}
