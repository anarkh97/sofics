#ifndef _GLOBALTOLOCALMAP_H_
#define _GLOBALTOLOCALMAP_H_

//#define HB_USE_MEMCPY
#ifdef HB_USE_MEMCPY
#include <cstring>
#endif

//#define MAP_MIN_MEMORY
#ifdef MAP_MIN_MEMORY
#include <map>
#include <Extermal.d/include/gsl/span>

#endif

template<class Scalar> class FSCommPattern;

class CommunicatableObject
{
 protected:
  template<class Scalar>
    void setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID);
  template<class Scalar>
    void pack(FSCommPattern<Scalar> *pat, int localID, int remoteID);
  template<class Scalar>
    void unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID);
};

class GlobalToLocalMap : public CommunicatableObject
{
  // HB: encapsulate a STL map<int,int>
  private:
    std::map<int,int> globalToLocalArray;

  public:  
    GlobalToLocalMap();
    GlobalToLocalMap(int localSize, int *localToGlobalArray);
    GlobalToLocalMap(gsl::span<const int> localToGlobalArray);
    ~GlobalToLocalMap();

    void initialize(int localSize, const int *localToGlobalArray);
    /** \brief Obtain the localnumber or -1 if it was not found. */
    int operator[](int glNum) const;
    /** \brief Obtain the localnumber or throw if it was not found. */
    int at(int glNum) const { return globalToLocalArray.at(glNum); }

    void print(bool skip=false);
    int size();
    int numKeys(); 

    template<class Scalar>
      void setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID);
    template<class Scalar>
      void pack(FSCommPattern<Scalar> *pat, int localID, int remoteID);
    template<class Scalar>
      void unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID);
};

#endif
