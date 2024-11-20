#include <Driver.d/SubDomain.h>
#include <Driver.d/SubDomainFactory.h>

template<>
GenSubDomainFactory<double>* GenSubDomainFactory<double>::getFactory() { return subDomainFactory.get(); }

template<>
GenSubDomainFactory<DComplex>* GenSubDomainFactory<DComplex>::getFactory() { return subDomainFactoryC.get(); }
