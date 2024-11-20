#ifndef __SUBDOMAIN_FACTORY_H__
#define __SUBDOMAIN_FACTORY_H__

#include <memory>
#include <Driver.d/SubDomain.h>

template<class Scalar>
class GenSubDomainFactory
{
 public:
  virtual ~GenSubDomainFactory() {}
  virtual
    GenSubDomain<Scalar>* createSubDomain(Domain& d, int sn, int nNodes, int *nds,
					  int nElems, int *elems, int gn) const
      { return new GenSubDomain<Scalar>(d, sn, nNodes, nds, nElems, elems, gn); }

  virtual
    GenSubDomain<Scalar>* createSubDomain(Domain & d, int sn, Connectivity &con, 
					  Connectivity &nds, int gn) const
      { return new GenSubDomain<Scalar>(d, sn, con, nds, gn); }

  virtual
    GenSubDomain<Scalar>* createSubDomain(Domain& d, int sn, CoordSet* nodes, 
					  Elemset* elems, int *glNodeNums, 
					  int *glElemNums, int gn) const
      { return new GenSubDomain<Scalar>(d, sn, nodes, elems, glNodeNums, glElemNums, gn); }


  static GenSubDomainFactory* getFactory();
};

extern std::unique_ptr<GenSubDomainFactory<double> >   subDomainFactory;
extern std::unique_ptr<GenSubDomainFactory<DComplex> > subDomainFactoryC;

#endif
