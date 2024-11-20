#ifndef _BLOCK_ALLOC_H_
#define _BLOCK_ALLOC_H_

#include <cstdlib>
#include <cstddef>
#include <new>
#include <vector>
#include <memory>
#include <sys/types.h>
#include <Utils.d/resize_array.h>

class BlockAlloc {
	int index = 0;
	int blLen;
	std::vector<std::unique_ptr<char[]>> allBlocks;
public:
	BlockAlloc(int l = 4096, int initnb = 16) : blLen(l)
	{ }
	BlockAlloc(BlockAlloc &&) = default;
	~BlockAlloc() = default;
	void *getMem(size_t nb);
};

void * operator new(size_t nbyte, BlockAlloc &block);
void * operator new[](size_t nbyte, BlockAlloc &block);
void operator delete(void *p, BlockAlloc &block);
void operator delete[](void *p, BlockAlloc &block);


extern BlockAlloc global_ba;

template<typename T>
class block_allocator
{
 public:
  typedef size_t     size_type;
  typedef ptrdiff_t  difference_type;
  typedef T*         pointer;
  typedef const T*   const_pointer;
  typedef T&         reference;
  typedef const T&   const_reference;
  typedef T          value_type;

  template<typename T1> struct rebind 
  { typedef block_allocator<T1> other; };

  block_allocator() throw() { }
  
  template<typename T1> block_allocator(const block_allocator<T1>&) throw() { }
  block_allocator(const block_allocator&) throw() { }
  ~block_allocator() throw() { }
  
  pointer address(reference x) const { return &x; }
  const_pointer address(const_reference x) const { return &x; }
  // NB: __n is permitted to be 0.  The C++ standard says nothing
  // about what the return value is when __n == 0.
  pointer allocate(size_type n, const void* = 0)
	{
	  if(n > this->max_size())
	{ throw std::bad_alloc(); }

	  pointer ret = static_cast<T*>(global_ba.getMem(n*sizeof(T)));
	  if (!ret)
	throw std::bad_alloc();
	  return ret;
	}
  
  // do nothing
  void deallocate(pointer p, size_type) {}

  size_type max_size() const throw() { return size_t(-1) / sizeof(T); }

  void construct(pointer p, const  T& val) { ::new(p) value_type(val); }
  void destroy(pointer p) { p->~T(); } 
};

template<typename T>
inline bool operator==(const block_allocator<T>&, const block_allocator<T>&)
{ return true; }
  
template<typename T>
inline bool operator!=(const block_allocator<T>&, const block_allocator<T>&)
{ return false; }


#endif
