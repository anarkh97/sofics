/*
  holds dec initializers if set from parser
*/
#ifndef DECINIT_H
#define DECINIT_H
class DecInit
{
  public : // FEM class design style
  int nsubs;
  int nthreads;
  int nproc;
  char* file;
  bool weight;
  bool memory;
  bool exitAfterDec, skip;
  bool nosa;
  bool trivial;
  bool fsgl;
  DecInit() :
    nsubs(1),
    nthreads(1),
    nproc(1),
    file(0),
    weight(false),
    memory(false),
    exitAfterDec(false),
    skip(false),
    nosa(false),
    trivial(false),
    fsgl(false)
    {}
};
#endif
