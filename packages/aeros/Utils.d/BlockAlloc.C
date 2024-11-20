#include <Utils.d/BlockAlloc.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>

void *
BlockAlloc::getMem(size_t nbyte)
{
  if(int(nbyte) > blLen) {
      auto u_p = std::make_unique<char[]>(nbyte);
      char *d = u_p.get();
      if(allBlocks.size() > 0) {
          std::swap(allBlocks.back(), u_p);
      } else
          index = blLen; // Force not to use the tail block for any other allocation.
      allBlocks.push_back(std::move(u_p));
      return d;
  }
  if(allBlocks.size() == 0 || (index+int(nbyte)) > blLen) {
     allBlocks.push_back(std::make_unique<char[]>(blLen) );
     index = 0;
  }
  if(nbyte & 0x7) {
    nbyte = (nbyte+8)-(nbyte & 0x7);
  }
  void *p = allBlocks.back().get()+index;
  index += nbyte;
  return p;
}

void * operator new(size_t nbyte, BlockAlloc &block)
{
 return block.getMem(nbyte);
}

void * operator new[](size_t nbyte, BlockAlloc &block)
{
 return block.getMem(nbyte);
}

BlockAlloc global_ba;

void operator delete(void *p, BlockAlloc &block)
{
  // this is only used for exception unwinding
  free(p);
}

void operator delete[](void *p, BlockAlloc &block)
{
  // this is only used for exception unwinding
  free(p);
}
