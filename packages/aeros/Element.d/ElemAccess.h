#ifndef _ELEMACCESS_H_
#define _ELEMACCESS_H_

#include <Element.d/Element.h>

#include <iostream>
#include <utility>

// This file contains Accessors for all of the types defined in Element.d

class BCDataAccessor {
 public:
  static int getNum(std::pair<int,BCond *> &, int )
    { return 1; }
  static int getSize(const std::pair<int,BCond *> &o)
    { return o.first; }
  static int *getData(std::pair<int,BCond *> &o, int i, int *nd)
    {
      if(i>=o.first)
        std::cout << "BCDataAccessor" << std::endl;
      if(nd) {
        nd[0] = o.second[i].nnum;
        return nd;
      }
      else return &o.second[i].nnum;
    }
};

class ComplexBCDataAccessor {
 public:
  static int getNum(std::pair<int,ComplexBCond *> &, int )
    { return 1; }
  static int getSize(const std::pair<int,ComplexBCond *> &o)
    { return o.first; }
  static int *getData(std::pair<int,ComplexBCond *> &o, int i, int *nd)
    {
      if(i>=o.first)
        std::cout << "ComplexBCDataAccessor" << std::endl;
      if(nd) {
        nd[0] = o.second[i].nnum;
        return nd;
      }
      else return &o.second[i].nnum;
    }
};

class CurveMatAccessor {
  // protect StructProp->__SOWER_TMP for multi-threading TG
 public:
  static int getNum(std::pair<int, SPropContainer*> &o, int i)
    {
      if(((int) - (*(o.second))[i].W) > 0)
        return 1;
      else
        return 0;
    }
  static int getSize(const std::pair<int, SPropContainer*> &o)
    { return o.first; }
  static int *getData(std::pair<int, SPropContainer*> &o, int i, int *nd)
    {
      if(((int) - (*(o.second))[i].W) > 0) {
        if(nd) {
          nd[0] = (int) - (*(o.second))[i].W;
          return nd;
        }
        else {
          (*(o.second))[i].__SOWER_TMP = (int) - (*(o.second))[i].W;
          return &((*(o.second))[i].__SOWER_TMP);
        }
      }
      std::cerr << "warning getting data for size = 0 !" << std::endl;
      return 0;
    }
};

class CurveYoungMatAccessor {
 // protect StructProp->__SOWER_TMP for multi-threading TG
 public:
  static int getNum(std::pair<int, SPropContainer *> &o, int i)
    {
      if(((int) - (*(o.second))[i].E) > 0)
        return 1;
      else
        return 0;
    }
  static int getSize(const std::pair<int,SPropContainer*> &o)
    { return o.first; }
  static int *getData(std::pair<int,SPropContainer*> &o, int i, int *nd)
    {
      if(((int) - (*(o.second))[i].E) > 0) {
        if(nd) {
          nd[0] = (int) - (*(o.second))[i].E;
          return nd;
        }
        else {
          (*(o.second))[i].__SOWER_TMP = (int) - (*(o.second))[i].E;
          return &((*(o.second))[i].__SOWER_TMP);
        }
      }
      std::cerr << "warning getting data for size = 0 !" << std::endl;
      return 0;
    }
};

class EsetGeomAccessor {
 public:
  static int getNum(Elemset &eSet, int elNum)
    { return eSet[elNum] ? eSet[elNum]->numNodes() : 0; }
  static int getSize(const Elemset &eSet)
    { return eSet.last(); }
  static int *getData(Elemset &eSet, int elNum, int *nd)
    { return eSet[elNum] ? eSet[elNum]->allNodes(nd) : 0; }
};

class EsetDataAccessor {
 public:
  static int getNum(Elemset &eSet, int elNum)
    { return 1; }
  static int getSize(const Elemset &eSet)
    { return eSet.last(); }
  static int *getData(Elemset &eSet, int elNum, int *nd)
    {
      if(eSet[elNum]) {
        if(nd) { nd[0] = elNum; return nd; }
        else { std::cerr << "(EE) error : unefficient use of EsetDataAccessor" << std::endl; return 0; }
      }
      else {
        return 0;
      }
    }
};

#endif
