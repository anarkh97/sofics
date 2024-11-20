#include <Feti.d/DistrVector.h>
#include <Driver.d/SubDomain.h>
#include <iostream>

template<class Scalar>
GenBasicAssembler<Scalar>::GenBasicAssembler(int _numSub, GenSubDomain<Scalar> **_sd, FSCommPattern<Scalar> *_pat)
{
  numSub = _numSub;
  sd = _sd;
  pat = _pat;
}

template<class Scalar>
void
GenBasicAssembler<Scalar>::assemble(GenDistrVector<Scalar> &v, int flag)
{
  execParal(numSub, this, &GenBasicAssembler<Scalar>::sendSubVec, v);
  pat->exchange();
  execParal(numSub, this, &GenBasicAssembler<Scalar>::assembleSubVec, v, flag);
}

template<class Scalar>
void
GenBasicAssembler<Scalar>::sendSubVec(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subV = v.subData(sd[iSub]->localSubNum());
  sd[iSub]->extractAndSendInterf(subV, pat);
}

template<class Scalar>
void
GenBasicAssembler<Scalar>::assembleSubVec(int iSub, GenDistrVector<Scalar> &v, int flag)
{
  Scalar *subV = v.subData(sd[iSub]->localSubNum());
  if(flag == 0) sd[iSub]->assembleInterf(subV, pat);
  else if(flag == 1) sd[iSub]->assembleInterfInvert(subV, pat);
  else std::cerr << "BasicAssembler::assembleSubVec bad flag\n";
}

template<class Scalar>
void
GenBasicAssembler<Scalar>::split(GenDistrVector<Scalar> &v)
{
  execParal(numSub, this, &GenBasicAssembler<Scalar>::splitSubVec, v);
}

template<class Scalar>
void
GenBasicAssembler<Scalar>::splitSubVec(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subV = v.subData(sd[iSub]->localSubNum());
  sd[iSub]->splitInterf(subV);
}

