#ifndef FWK_MANAGER_H
#define FWK_MANAGER_H

#include "NamedInterface.h"
#include "Exception.h"
#include "Macros.h"
#include <map>

namespace Fwk {

/* Generic Manager Interface
** PtrT => Managed pointer type
** K => Key */
template <typename PtrT, typename K>
class GenManagerInterface : public Fwk::PtrInterface<GenManagerInterface<PtrT, K> > {
public:
  EXPORT_PTRINTERFACE_TYPES(GenManagerInterface);

  // Read accessors
  virtual PtrT instance(const K & key) const = 0;
  virtual size_t instanceCount() const = 0;

  // Mutators
  virtual PtrT instanceNew(const K & key) = 0;
  virtual void instanceDel(const K & key) = 0;
};

template <typename PtrT, typename K, typename Base>
class GenManagerSubInterface : public Base {
public:
  EXPORT_PTRINTERFACE_TYPES(GenManagerSubInterface);
  
  // Read accessors
  virtual PtrT instance(const K & key) const = 0;
  virtual size_t instanceCount() const = 0;

  // Mutators
  virtual PtrT instanceNew(const K & key) = 0;
  virtual void instanceDel(const K & key) = 0;
};

/* Generic Manager Implementation
** T => Managed type
** K => Key                       */
template <typename T, typename K, typename CreatePtrT = T*>
class GenManagerImpl {
public:
  typedef std::map<K, Ptr<T> > InstanceMap;
  typedef typename InstanceMap::size_type InstanceCount;
  typedef typename InstanceMap::const_iterator IteratorConst;
  typedef typename InstanceMap::iterator Iterator;
  typedef typename InstanceMap::value_type Pair;
  
  T * instance(const K & key) const;
  InstanceCount instanceCount() const { return instance_.size(); }
  
  T * instanceNew(const K & key);
  void instanceDel(const K & key) { instance_.erase(key); } 
  
  IteratorConst instanceBegin() const { return instance_.begin(); }
  IteratorConst instanceEnd() const { return instance_.end(); }

  Iterator instanceBegin() { return instance_.begin(); }
  Iterator instanceEnd() { return instance_.end(); }

protected:
  virtual ~GenManagerImpl() {}

  virtual CreatePtrT createNewInstance(const K & key) = 0;

private:
  InstanceMap instance_;
};

template <typename T, typename K, typename CreatePtrT>
T *
GenManagerImpl<T, K, CreatePtrT>::instance(const K & key) const {
  typename InstanceMap::const_iterator it = instance_.find(key);
  return (it != instance_.end()) ? it->second.ptr() : NULL;
}

template <typename T, typename K, typename CreatePtrT>
T *
GenManagerImpl<T, K, CreatePtrT>::instanceNew(const K & key) {
  // Find insertion point
  typename InstanceMap::iterator it = instance_.lower_bound(key);
  if (it != instance_.end() && it->first == key)
    throw NameInUseException();
 
  // Build and add new instance
  Ptr<T> newInstance = createNewInstance(key);
  instance_.insert(it, std::make_pair(key, newInstance));

  return newInstance.ptr();
}

/* Factory-based generic manager implementation */
template <typename T, typename K>
struct InstanceFactory {
  typedef T InstanceType;
  typedef K KeyType;

  /* Derived classes should implement [T / Ptr<T>] operator()(const K & key) [const] */
};

template <typename FactoryType> 
class FactoryManagerImpl {
public:
  typedef typename FactoryType::InstanceType InstanceType;
  typedef typename FactoryType::KeyType KeyType;
  typedef std::map<KeyType, Ptr<InstanceType> > InstanceMap;
  typedef typename InstanceMap::size_type InstanceCount;
  typedef typename InstanceMap::const_iterator IteratorConst;
  typedef typename InstanceMap::iterator Iterator;
  typedef typename InstanceMap::value_type Pair;

  InstanceType * instance(const KeyType & key) const;
  InstanceCount instanceCount() const { return instance_.size(); }
  
  InstanceType * instanceNew(const KeyType & key);
  void instanceDel(const KeyType & key) { instance_.erase(key); } 
  
  IteratorConst instanceBegin() const { return instance_.begin(); }
  IteratorConst instanceEnd() const { return instance_.end(); }

  Iterator instanceBegin() { return instance_.begin(); }
  Iterator instanceEnd() { return instance_.end(); }

  explicit FactoryManagerImpl(FactoryType factory) :
    factory_(factory),
    instance_()
  {}

private:
  FactoryType factory_;
  InstanceMap instance_;
};

template <typename FactoryType>
typename FactoryManagerImpl<FactoryType>::InstanceType *
FactoryManagerImpl<FactoryType>::instance(const typename FactoryManagerImpl<FactoryType>::KeyType & key) const {
  typename InstanceMap::const_iterator it = instance_.find(key);
  return (it != instance_.end()) ? it->second.ptr() : NULL;
}

template <typename FactoryType>
typename FactoryManagerImpl<FactoryType>::InstanceType *
FactoryManagerImpl<FactoryType>::instanceNew(const typename FactoryManagerImpl<FactoryType>::KeyType & key) {
  // Find insertion point
  typename InstanceMap::iterator it = instance_.lower_bound(key);
  if (it != instance_.end() && it->first == key)
    throw NameInUseException();
 
  // Build and add new instance
  Ptr<InstanceType> newInstance = factory_(key);
  instance_.insert(it, std::make_pair(key, newInstance));

  return newInstance.ptr();
}

/* Template Manager */
template <typename T, typename K, typename CreatePtrT = Ptr<T> >
class GenManager : public PtrInterface<GenManager<T, K, CreatePtrT> >, private GenManagerImpl<T, K, CreatePtrT> {
public:
  EXPORT_PTRINTERFACE_TYPES(GenManager);
 
  typedef typename GenManagerImpl<T, K, CreatePtrT>::InstanceCount InstanceCount;

  T * instance(const K & key) const { return GenManagerImpl<T, K, CreatePtrT>::instance(key); } 
  InstanceCount instanceCount() const { return GenManagerImpl<T, K, CreatePtrT>::instanceCount(); }
  
  T * instanceNew(const K & key) { return GenManagerImpl<T, K, CreatePtrT>::instanceNew(key); }
  void instanceDel(const K & key) { GenManagerImpl<T, K, CreatePtrT>::instanceDel(key); }

protected:
  GenManager() {}
  ~GenManager() {}
 
  typedef K KeyType;
  typedef GenManagerImpl<T, K, CreatePtrT> Impl;
  typedef typename Impl::IteratorConst IteratorConst;
  typedef typename Impl::Iterator Iterator;
  
  IteratorConst instanceBegin() const { return Impl::instanceBegin(); }
  IteratorConst instanceEnd() const { return Impl::instanceEnd(); }

  Iterator instanceBegin() { return Impl::instanceBegin(); }
  Iterator instanceEnd() { return Impl::instanceEnd(); }

  virtual CreatePtrT createNewInstance(const K & key) = 0;

  DISALLOW_COPY_AND_ASSIGN(GenManager);
};


/* Specialized Manager for NamedInterface */
typedef GenManagerImpl<NamedInterface, String> NamedInterfaceManagerImpl;

template <typename T>
class GenNamedInterfaceManager : public PtrInterface<GenNamedInterfaceManager<T> >, private NamedInterfaceManagerImpl {
public:
  EXPORT_PTRINTERFACE_TYPES(GenNamedInterfaceManager<T>);
 
  typedef NamedInterfaceManagerImpl::InstanceCount InstanceCount;

  T * instance(const String & name) const { return static_cast<T *>(NamedInterfaceManagerImpl::instance(name)); } 
  InstanceCount instanceCount() const { return NamedInterfaceManagerImpl::instanceCount(); }
  
  T * instanceNew(const String & name) { return static_cast<T *>(NamedInterfaceManagerImpl::instanceNew(name)); }
  void instanceDel(const String & name) { NamedInterfaceManagerImpl::instanceDel(name); }

protected:
  GenNamedInterfaceManager() {}
  GenNamedInterfaceManager(const GenNamedInterfaceManager<T> &); // No implementation
  GenNamedInterfaceManager<T> & operator=(const GenNamedInterfaceManager<T> &); // No implementation

  virtual T * createNewInstance(const String & key) = 0;
};

typedef GenNamedInterfaceManager<NamedInterface> NamedInterfaceManager;

} // end namespace Fwk

#endif /* FWK_MANAGER_H */
