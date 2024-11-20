// ---------------------------------------------------------------- 
// HB - 06/25/03
// ---------------------------------------------------------------- 
// WARNING: IT IS CURRENTLY IMPLICITLY ASSUMED THAT NO GAP EXIST IN
//          THE NUMBERING OF THE FACE ELEMENTS IN THE FaceElemSet,
//          AND THE NUMBERING START FROM 0 TO N-1 (WHERE N IS THE 
//          NUMBER OF FACE ELEMENTS REALLY INSTANTIATED), SO THAT
//          FaceElemSet::nElems() RETURN N.
// ---------------------------------------------------------------- 
#ifndef _FACEELEMSET_H_
#define _FACEELEMSET_H_

// STL
#include <utility>
#include <map>
typedef std::pair<int, std::pair<double, double> > locoord;
//                elem id        xi1     xi2

// FEM headers
#include <Utils.d/BlockAlloc.h>

#ifdef SOWER_SURFS
#include <Utils.d/BinFileHandler.h>
#endif
                                                                                                              
class FaceElement;

class FaceElemSet {
  protected:
    FaceElement **elem;
    int _last;
    int emax;
    BlockAlloc ba;
  public:
    FaceElemSet(int = 256);
    //~FaceElemSet() { ba.~BlockAlloc(); deleteElems(); }
    ~FaceElemSet() { deleteElems(); }
    int size() const { return emax; }
    int last() const;
    FaceElement *operator[] (int i) const { return elem[i]; }
    void elemadd(int num, FaceElement *);
    void elemadd(int num, int type, int nnodes, int *nodes);

    int nElems();
    void print();

    void Renumber(std::map<int,int>& OldToNewNodeIds);

    void deleteElems() { if(elem) { delete [] elem; elem = 0; } emax = 0; _last = 0; }
    void remove(int num);
    void repack();

    std::map<int,locoord> computeNodeLocalCoords(int*, int);

#ifdef SOWER_SURFS
    void WriteSower(BinFileHandler& file, int* ndMap=0);
#endif
};

#endif
