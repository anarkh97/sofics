#ifndef _SPARSECONNECTIVITYT_H_
#define _SPARSECONNECTIVITYT_H_

#include <Utils.d/ConnectivityT.h>
#include <map>

template<typename IndexType, typename DataType>
class SparseConnectivityAccessorT {
 public:
  static IndexType getNum(std::pair<std::map<IndexType,int>, ConnectivityT<size_t,DataType>*> &o, IndexType i)
    { typename std::map<IndexType,int>::const_iterator it = o.first.find(i);
      return (it != o.first.end()) ? o.second->num(it->second) : 0; }
  static IndexType getSize(const std::pair<std::map<IndexType,int>,ConnectivityT<size_t,DataType>*> &o)
    { return o.first.rbegin()->first+1; }
  static gsl::span<const DataType> getData(std::pair<std::map<IndexType,int>,ConnectivityT<size_t,DataType>*> &o,
      IndexType i, gsl::span<DataType> nd)
    {
      typename std::map<IndexType,int>::const_iterator it = o.first.find(i);
      if(it == o.first.end()) {
        return (nd.data()) ? nd : gsl::span<const DataType>{};
      }
      else {
        if(nd.data()) {
          for(int j=0; j<o.second->num(it->second); ++j)
            nd[j] = (*o.second)[it->second][j];
          return nd;
        }
        else return (*o.second)[it->second];
      }
    }
};

#include <Types.h>

typedef std::pair<std::map<int,int>,ConnectivityT<size_t,gl_num_t>*> SparsePairType1;
typedef ImplicitConnectivityT<SparsePairType1*, SparseConnectivityAccessorT<int,gl_num_t> > SparseConnectivityType1;

typedef std::pair<std::map<gl_num_t,int>,ConnectivityT<size_t,int>*> SparsePairType2;
typedef ImplicitConnectivityT<SparsePairType2*, SparseConnectivityAccessorT<gl_num_t,int> > SparseConnectivityType2;

template <typename _IndexType, typename _DataType>
struct base_traits<SparseConnectivityAccessorT<_IndexType,_DataType> > {
  typedef _IndexType IndexType;
  typedef _DataType DataType;
};

#endif
