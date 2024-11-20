#ifndef IPOINTS_H
#define IPOINTS_H

#include <Utils.d/resize_array.h>
#include <sys/types.h>
#include <Utils.d/BlockAlloc.h>
//#include <Utils.d/Connectivity.h>

class IPoint {  // Interface point
  protected:
    int from,to ;
  public:
    IPoint(int f, int t) { from = f; to = t ; }
    friend class IPointset ;
} ;
 
class IPointLink : public IPoint {
  protected:
    int num ;
    int attrib;
    IPointLink * next ;
  public:
    IPointLink (int f, int t, int n, IPointLink *nx, int attr = -1)
       : IPoint(f,t)
       {  num = n; next = nx; attrib = attr; } ;
    friend class IPointset ;
} ;

class CountedConnectivity;


class IPointset {
     int np ;    // Current number of diffent points
     int nodemax;
     ResizeArray<IPointLink*> ipl;
     BlockAlloc blkal;
  public:
     IPointset(int = 16) ;
     ~IPointset() ;
     int size() { return np ; }
     int num(int,int, int = -1) ;
     int num(IPoint &) ;    // Find the number corresponding to an interface
                            // point
     CountedConnectivity* two_way_connect();
     float (* coord(float *,float *,float *))[3] ;
           // Create the coordinates fo a given scalar value
} ;
#endif 
