#ifndef FWK_MACROS_H
#define FWK_MACROS_H

// Standard C assert
#include <cassert>

// A macro to disallow the copy constructor and operator= functions
#define DISALLOW_COPY_AND_ASSIGN(TypeName)   \
  TypeName(const TypeName &);                 \
  const TypeName & operator=(const TypeName &)

// A macro to export the smart pointer typedefs Ptr/PtrConst
#define EXPORT_PTRINTERFACE_TYPES(Typename) \
  typedef Fwk::Ptr< Typename > Ptr;           \
  typedef Fwk::Ptr< const Typename > PtrConst  

#endif /* FWK_MACROS_H */
