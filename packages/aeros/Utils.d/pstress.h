#ifndef _PSTRESS_H_
#define _PSTRESS_H_

template<class Scalar>
  void pstress(Scalar stress[6], Scalar pstress[3], Scalar *pdirections = 0);

#ifdef _TEMPLATE_FIX_
  #include <Utils.d/pstress.C>
#endif

#endif
