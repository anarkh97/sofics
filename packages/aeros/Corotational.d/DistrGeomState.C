#include <Corotational.d/DistrGeomState.h>
#include <Feti.d/DistrVector.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/DecDomain.h>
#include <Threads.d/PHelper.h>
#include <Corotational.d/TemperatureState.h>
#ifdef USE_MPI
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

DistrGeomState::DistrGeomState(DecDomain* domain)
{
 // Number of subdomains
 numSub = domain->getNumSub();

 // array of pointers to the subdomain's geometry states
 gs = new GeomState*[numSub];

 // array of data structures to store lagrange multipliers 
 mu = new std::map<std::pair<int,int>,double>[numSub];
 lambda = new std::vector<double>[numSub];

 // parallel execution of subdomain geometry state construction
 execParal(numSub,this,&DistrGeomState::makeSubGeomStates,domain); 
}

// Subdomain geom state construction
void
DistrGeomState::makeSubGeomStates(int isub, DecDomain *domain)
{
 SubDomain *sd = domain->getSubDomain(isub);
 if(sd->solInfo().soltyp == 2)
   gs[isub] = new TemperatureState(*sd->getDSA(), *sd->getCDSA(), sd->getNodes());
 else
   gs[isub] = new GeomState(*sd->getDSA(), *sd->getCDSA(), sd->getNodes(), &sd->getElementSet(),
                            sd->getNodalTemperatures());
}

DistrGeomState::DistrGeomState(const DistrGeomState &g2)
{
  numSub = g2.getNumSub();
  gs = new GeomState*[numSub];
  mu = new std::map<std::pair<int,int>,double>[numSub];
  lambda = new std::vector<double>[numSub];
  execParal(numSub,this,&DistrGeomState::subCopyConstructor,g2);
}

void
DistrGeomState::subCopyConstructor(int isub, const DistrGeomState &g2)
{
  TemperatureState* ts;
  if((ts = dynamic_cast<TemperatureState*>(g2[isub]))) gs[isub] = new TemperatureState(*ts);
  else gs[isub] = new GeomState(*(g2[isub]));

  mu[isub] = g2.mu[isub];
  lambda[isub] = g2.lambda[isub];
}

DistrGeomState::~DistrGeomState()
{
  for(int i=0; i<numSub; ++ i) delete gs[i];
  delete [] gs;
  delete [] mu;
  delete [] lambda;
}

// Subdomain update
void
DistrGeomState::subUpdate(int isub, DistrVector &v, int SO3param)
{
  StackVector vec(v.subData(isub), v.subLen(isub));
  gs[isub]->update(vec, SO3param);
}

void
DistrGeomState::subUpdateRef(int isub, DistrGeomState &ref, DistrVector &v, int SO3param)
{
  StackVector vec(v.subData(isub), v.subLen(isub));
  gs[isub]->update(*ref[isub], vec, SO3param);
}

void
DistrGeomState::subExplicitUpdate(int isub, DistrVector &v, GenDecDomain<double> *decDomain)
{
  StackVector vec(v.subData(isub), v.subLen(isub));
  gs[isub]->explicitUpdate(decDomain->getSubDomain(isub)->getNodes(), vec); 
}

void
DistrGeomState::subSetVelocity(int isub, DistrVector &v, int SO3param)
{
  StackVector vsub(v.subData(isub), v.subLen(isub));
  gs[isub]->setVelocity(vsub, SO3param);
}

void
DistrGeomState::subSetAcceleration(int isub, DistrVector &a, int SO3param)
{
  StackVector asub(a.subData(isub), a.subLen(isub));
  gs[isub]->setAcceleration(asub, SO3param);
}

void
DistrGeomState::subSetVelocityAndAcceleration(int isub, DistrVector &v, DistrVector &a)
{
  StackVector vsub(v.subData(isub), v.subLen(isub));
  StackVector asub(a.subData(isub), a.subLen(isub));
  gs[isub]->setVelocityAndAcceleration(vsub, asub);
}

void
DistrGeomState::subSetNodalTemperatures(int isub, DistrVector &temps)
{
  gs[isub]->setNodalTemperatures(temps.subData(isub));
}

void
DistrGeomState::subStep_update(int isub, DistrVector &v_n, DistrVector &a_n, 
                               double &delta, DistrGeomState &ss,
                               double beta, double gamma, double alphaf, double alpham, bool zeroRot)
{
 StackVector vel(v_n.subData(isub),v_n.subLen(isub));
 StackVector acc(a_n.subData(isub),a_n.subLen(isub));
 GeomState &step = *ss[isub];
 gs[isub]->midpoint_step_update(vel,acc,delta,step,beta,gamma,alphaf,alpham,zeroRot);
}

void
DistrGeomState::midpoint_step_update(DistrVector &veloc_n, DistrVector &accel_n,
                                     double &delta, DistrGeomState &ss,
                                     double beta, double gamma, double alphaf, double alpham, bool zeroRot)
{
 execParal(numSub,this,&DistrGeomState::subStep_update,veloc_n,accel_n,delta,ss,beta,gamma,alphaf,alpham,zeroRot);
}

void
DistrGeomState::subInc_get(int isub,DistrVector &inc_vec, DistrGeomState &ss, bool zeroRot)
{
 StackVector v(inc_vec.subData(isub),inc_vec.subLen(isub));
 gs[isub]->get_inc_displacement(v,*ss[isub],zeroRot);
}

void
DistrGeomState::get_inc_displacement(DistrVector &inc_vec, DistrGeomState &ss, bool zeroRot)
{
 execParal(numSub,this,&DistrGeomState::subInc_get, inc_vec, ss, zeroRot);
}

void
DistrGeomState::subPushForward(int isub, DistrVector &f)
{
 StackVector subf(f.subData(isub), f.subLen(isub));
 gs[isub]->push_forward(subf);
}

void
DistrGeomState::push_forward(DistrVector &f)
{
 execParal(numSub, this, &DistrGeomState::subPushForward, f);
}

void
DistrGeomState::subPullBack(int isub, DistrVector &f)
{
 StackVector subf(f.subData(isub), f.subLen(isub));
 gs[isub]->pull_back(subf);
}

void
DistrGeomState::pull_back(DistrVector &f)
{
 execParal(numSub, this, &DistrGeomState::subPullBack, f);
}

void
DistrGeomState::subTransform(int isub, DistrVector &f, int type, bool unscaled)
{
 StackVector subf(f.subData(isub), f.subLen(isub));
 gs[isub]->transform(subf, type, unscaled);
}

void
DistrGeomState::transform(DistrVector &f, int type, bool unscaled)
{
 execParal(numSub, this, &DistrGeomState::subTransform, f, type, unscaled);
}


void
DistrGeomState::subTot_get(int isub, DistrVector &tot_vec, bool rescaled)
{
 StackVector v(tot_vec.subData(isub), tot_vec.subLen(isub));
 gs[isub]->get_tot_displacement(v, rescaled);
}

void
DistrGeomState::get_tot_displacement(DistrVector &tot_vec, bool rescaled)
{
	execParal(numSub, this, &DistrGeomState::subTot_get, tot_vec, rescaled);
}

void
DistrGeomState::subInterp(int isub, double &alpha, DistrGeomState &u,
                          DistrGeomState &un)
{
 GeomState &uR  = *u[isub];
 GeomState &unR = *un[isub];
 gs[isub]->interp(alpha, uR, unR );
}

void
DistrGeomState::interp(double alpha, DistrGeomState &u, DistrGeomState &un)
{
  execParal(numSub, this, &DistrGeomState::subInterp, alpha, u, un);
}

void
DistrGeomState::subDiff(int isub, DistrGeomState &unp, DistrVector &un)
{
 StackVector u(un.subData(isub), un.subLen(isub));
 GeomState &unpR = *unp[isub];
 gs[isub]->diff(unpR, u);
}

void
DistrGeomState::diff(DistrGeomState &unp, DistrVector &un)
{
  execParal(numSub, this, &DistrGeomState::subDiff, unp, un);
}

void
DistrGeomState::update(DistrVector &v, int SO3param)
{
  execParal(numSub, this, &DistrGeomState::subUpdate, v, SO3param);
}

void
DistrGeomState::update(DistrGeomState &ref, DistrVector &v, int SO3param)
{
  execParal(numSub, this, &DistrGeomState::subUpdateRef, ref, v, SO3param);
}

void
DistrGeomState::explicitUpdate(GenDecDomain<double> *decDomain, DistrVector &v)
{
  execParal(numSub, this, &DistrGeomState::subExplicitUpdate, v, decDomain);
}

void
DistrGeomState::setVelocity(DistrVector &v, int SO3param)
{
  execParal(numSub, this, &DistrGeomState::subSetVelocity, v, SO3param);
}

void
DistrGeomState::setAcceleration(DistrVector &a, int SO3param)
{
  execParal(numSub, this, &DistrGeomState::subSetAcceleration, a, SO3param);
}

void
DistrGeomState::setVelocityAndAcceleration(DistrVector &v, DistrVector &a)
{
  execParal(numSub, this, &DistrGeomState::subSetVelocityAndAcceleration, v, a);
}

void
DistrGeomState::setNodalTemperatures(DistrVector &temps)
{
  execParal(numSub, this, &DistrGeomState::subSetNodalTemperatures, temps);
}

DistrGeomState &
DistrGeomState::operator=(DistrGeomState &unp)
{
  execParal(numSub, this, &DistrGeomState::subCopy, unp);
  return *this;
}

void
DistrGeomState::subCopy(int isub, DistrGeomState &unp)
{
  GeomState &unpR = *unp[isub];
  *(gs[isub]) = unpR; 

  mu[isub] = unp.mu[isub];
  lambda[isub] = unp.lambda[isub];
}

int 
DistrGeomState::getTotalNumElemStates()
{
  int ret = 0;
  for(int i=0; i<numSub; ++i) ret += gs[i]->getTotalNumElemStates();
#ifdef USE_MPI
  ret = structCom->globalSum(ret);
#endif
  return ret;
}

bool
DistrGeomState::getHaveRot()
{
  int ret = 0;
  for(int i=0; i<numSub; ++i) ret += int(gs[i]->getHaveRot());
#ifdef USE_MPI
  ret = structCom->globalSum(ret);
#endif
  return bool(ret);
}

void
DistrGeomState::resize(DecDomain* domain, std::map<std::pair<int,int>,double> *mu)
{
  execParal(numSub, this, &DistrGeomState::subResize, domain);
  if(mu) setMultipliers(*mu);
  domain->exchangeInterfaceGeomState(this);
}

void
DistrGeomState::subResize(int isub, DecDomain* domain)
{
 SubDomain *sd = domain->getSubDomain(isub);
 gs[isub]->resize( *sd->getDSA(), *sd->getCDSA(), sd->getNodes(), &sd->getElementSet() );
}

void
DistrGeomState::print()
{
  for(int i=0; i<numSub; ++i) {
    gs[i]->print();
  }
}

void
DistrGeomState::getMultipliers(std::map<std::pair<int,int>,double> &mu)
{
  for(std::map<std::pair<int,int>,double>::iterator it = mu.begin(); it != mu.end(); it++) {
    it->second = 0;
    for(int isub=0; isub<numSub; ++isub) {
      double val = gs[isub]->getMultiplier(it->first);
      if(val != 0) { it->second = val; break; }
    }
  }
#ifdef USE_MPI
  std::vector<double> values;
  for(std::map<std::pair<int,int>,double>::iterator it = mu.begin(); it != mu.end(); it++) values.push_back(it->second);
  structCom->globalMax(values.size(), values.data()); // XXX inequlity constraint multipliers always non-negative?
  std::vector<double>::iterator it2 = values.begin();
  for(std::map<std::pair<int,int>,double>::iterator it = mu.begin(); it != mu.end(); it++) { it->second = *it2; it2++; }
#endif
}

void
DistrGeomState::setMultipliers(std::map<std::pair<int,int>,double> &mu)
{
  execParal(numSub, this, &DistrGeomState::subSetMultipliers, mu);
}

void
DistrGeomState::subSetMultipliers(int isub, std::map<std::pair<int,int>,double> &mu)
{
  gs[isub]->setMultipliers(mu);
}
