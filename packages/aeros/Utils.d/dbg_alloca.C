#include <Utils.d/dbg_alloca.h>

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
//--- UH ---
  #if defined(MACOSX) || defined(ANDROID)
    void print_trace (void)
    {
    }

    void print_trace_handler (void)
    {
    }

  #else
    // IMPORTANT: add -Xlinker --export-dynamic to use this functionality!!!
    #include <execinfo.h>
    #include <cstdlib>
    #include <cstdio>
    #define NFRAMES 64
    void print_trace (void)
    {
      void *array[NFRAMES];
      size_t size;
      char **strings;
      size_t i;
  
      size = backtrace (array, NFRAMES);
      strings = backtrace_symbols (array, size);
  
      fprintf (stderr, "Obtained %zd stack frames.\n", size);
  
      for (i = 0; i < size; i++)
        printf ("%s\n", strings[i]);
  
      free (strings);
    }

    void print_trace_handler (void)
    {
      print_trace();
      abort();
    }

  #endif
//--- UH ---

#endif

#include <cassert>
#include <iostream>
#include <cstdlib>
using namespace std;
// moved to dbg_alloca.h #define _MAX_ALLOCA_SIZE (1<<16)
unsigned dbg_alloca_assert(unsigned size)
{
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
//--- UH ---
  #if defined(MACOSX) || defined(ANDROID)
    assert(size <= _MAX_ALLOCA_SIZE);
    return size;
  #else
    if(size > _MAX_ALLOCA_SIZE)
      { /*char * mp_stack_size = getenv("MP_STACK_SIZE");
        char * stacksize = getenv("STACKSIZE");
        char * mpstkz = getenv("MPSTKZ");*/
        print_trace(); }
  #endif
//--- UH ---
#endif
  assert(size <= _MAX_ALLOCA_SIZE);
  return size;
}
