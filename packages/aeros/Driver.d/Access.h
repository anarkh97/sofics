#ifndef _ACCESS_H_
#define _ACCESS_H_

#include <Driver.d/OffsetData.h>
#include <Driver.d/DMassData.h>
#include <Driver.d/EFrameData.h>
#include <Driver.d/Attrib.h>
#include <Driver.d/Mpc.h>

#include <iostream>
#include <utility>
#include <vector>

// This file contains Accessors for all of the types defined in Driver.d

class BoffsetAccessor {
 public:
  static int getNum(std::vector<OffsetData> & l, int i) 
    { return ((l[i]).last - (l[i]).first + 1); }
  static int getSize(const std::vector<OffsetData> & l)
    { return l.size(); }
  static int *getData(std::vector<OffsetData> & l, int i, int *nd)
    { 
      if(i < int(l.size())) {
        int j = 0;
        int *ndr;
        if(nd != 0)
           ndr = nd;
        else {
          ndr = new int(getNum(l,i));
          std::cerr << "(W) warning : inefficient use of BoffsetAccessor" << std::endl;
        }
        for(int k = (l[i]).first; k <= (l[i]).last; ++k,++j)
          ndr[j] = k;
        return ndr;
      }
      else {
        std::cerr << "BoffsetAccessor : offset out of range" << std::endl;
      }
      return(0);
    }
};

class DimassAccessor
{
 public:
  static int getNum(std::vector<DMassData*> &, int)
    { return 1; }
  static int getSize(const std::vector<DMassData*> &o)
    { return o.size(); }
  static int *getData(std::vector<DMassData*> &o, int i, int *nd)
    {
      if(nd != 0) {
        nd[0] = (*(o[i])).node;
        return nd;
      }
      else
        return &(*(o[i])).node;
    }
};

class EFrameDataAccessor {
 public:
  static int getNum(std::pair<int,ResizeArray<EFrameData>* > &, int)
   { return 1; }
  static int getSize(const std::pair<int,ResizeArray<EFrameData>* > &o)
   { return o.first; }
  static int *getData(std::pair<int,ResizeArray<EFrameData>* > &o, int i, int *nd)
   {
     if(i >= o.first)
       std::cout << "EFrameDataAccessor corruption" << std::endl;
     if(nd) {
       nd[0] = ((*(o.second))[i]).elnum;
       return nd;
     }
     else return &(*o.second)[i].elnum;
   }
};

class ElemAttrAccessor {
 public:
  static int getNum(std::pair<int,std::map<int,Attrib> *> &, int)
    { return 1; }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
    {
      if(nd) {
        nd[0] = (*(o.second))[i].nele;
        return nd;
      }
      else return &(*(o.second))[i].nele;
    }
};

class MatAttrAccessor {
 public:
  static int getNum(std::pair<int,std::map<int,Attrib> *> &, int)
    { return 1; }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
    {
      if(nd) {
        nd[0] = (*(o.second))[i].attr;
        return nd;
      }
      else return &(*(o.second))[i].attr;
    }
};

class CmpAttrAccessor {
 public:
  static int getNum(std::pair<int,std::map<int,Attrib> *> &o, int i)
    {
      if((*(o.second))[i].cmp_attr != -1)
        return 1;
      else
        return 0;
    }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
    {
      if((*(o.second))[i].cmp_attr != -1) {
        if(nd) {
          nd[0] = (*(o.second))[i].cmp_attr;
          return nd;
        }
        else return &(*(o.second))[i].cmp_attr;
      }
      std::cerr << "warning getting data for size = 0 !" << std::endl;
      return 0;
    }
};

class CmpFrAttrAccessor{
 public:
  static int getNum(std::pair<int,std::map<int,Attrib> *> &o, int i)
    {
      if((*(o.second))[i].cmp_frm != -1)
        return 1;
      else
        return 0;
    }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
    {
      if((*(o.second))[i].cmp_frm != -1) {
        if(nd) {
          nd[0] = (*(o.second))[i].cmp_frm;
          return nd;
        }
        else return &(*(o.second))[i].cmp_frm;
      }
      std::cerr << "warning getting data for size = 0 !" << std::endl;
      return 0;
    }
};

class LMPCAccessor {
 public:
  static int getNum(std::pair<int,ResizeArray<LMPCons*> *> &o, int i)
    {
      if((*(o.second))[i]->nterms > 0)
        return (*(o.second))[i]->nterms;
      else
        return 0;
    }
  static int getSize(const std::pair<int,ResizeArray<LMPCons*> *> &o)
    { return o.first; }
  static int *getData(std::pair<int,ResizeArray<LMPCons*> *> &o, int i, int *nd)
    {
      if((*(o.second))[i]->nterms > 0) {
        if(nd) {
          for(int j = 0; j < (*(o.second))[i]->nterms; ++j)
            nd[j] = (*(o.second))[i]->terms[j].nnum;
          return nd;
        }
        else {
          nd = new int[(*(o.second))[i]->nterms];
          for(int j = 0; j < (*(o.second))[i]->nterms; ++j)
            nd[j] = (*(o.second))[i]->terms[j].nnum;
          return nd;
        }
      }
      std::cerr << "warning getting data for size = 0 !" << std::endl;
      return 0;
    }
};

#endif
