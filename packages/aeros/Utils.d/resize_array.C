template <class Type>
ResizeArray<Type>::ResizeArray(Type iv, int is) 
{
  init_v = iv;
  csize  = is;
  if(csize > 0) {
    d = new Type[csize];
    int i;
    for(i=0; i < csize; ++i)
      d[i] = init_v;
  }
  else d = 0;

  myData = true;
}

template <class Type>
void
ResizeArray<Type>::resize(int ns)
{
  int nsize = (3*csize)/2;
  if(ns >= nsize) nsize=ns+1;
  Type *nd = new Type[nsize];
  int j;
  for(j=0; j < csize; ++j)
     nd[j] = d[j];
  for(; j <nsize; ++j)
    nd[j] = init_v;
#if !defined(DEBUG_OPENMP)
  if(d) delete [] d;
#endif
  csize=nsize;
  d=nd;
}

template <class Type>
ResizeArray<Type>::~ResizeArray()
{
  if(d && myData) { delete [] d; d = 0; }
}

template  <class Type>
Type *
ResizeArray<Type>::yield()
{
  Type *old = d;
  d = (Type *)0;
  return old;
}
