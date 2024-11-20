#ifndef __DBG_ALLOCA_H__
#define __DBG_ALLOCA_H__

#ifdef __GNUC__
#include <new>
void print_trace_handler (void);
#endif

#include <alloca.h>
unsigned dbg_alloca_assert(unsigned size);

//#define USE_DBG_ALLOCA_ASSERT
#define _MAX_ALLOCA_SIZE (1<<16)

#ifdef USE_DBG_ALLOCA_ASSERT
 #define dbg_alloca(size) \
  alloca(dbg_alloca_assert(size))
#else
 #define dbg_alloca(size) \
  alloca(size)
#endif

#endif

