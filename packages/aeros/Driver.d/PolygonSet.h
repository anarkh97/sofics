#ifndef _POLYGON_SET_H_
#define _POLYGON_SET_H_

#include <cstdio>
#include <Element.d/Element.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Utils.d/resize_array.h>

class PolyLine2;
class PolyTri3;
class PolyQuad4;
class PolygonSet;
class PolygonSetMessage;

#define POLYGON_LINE2 1
#define POLYGON_TRI3 2
#define POLYGON_QUAD4 3

#define POLYGON_TYPE_MASK 31

#define POLYGON_SOLID 32


class Polygon {
  public:
    Polygon *next;
    int flag;
    Element *el;

    Polygon() { next = 0; }
    virtual ~Polygon() {/*TODO*/}
    void link(Polygon *p) { next = p; }
    virtual int min_node() {return 0;}
    virtual int isLikeLine(PolyLine2 *);
    virtual int isLikeTri(PolyTri3 *);
    virtual int isLikeQuad(PolyQuad4 *);
    virtual int findSelf(PolygonSet *); // Is this polygon present in a set?
    virtual int findSelf(PolygonSetMessage *);
    virtual SommerElement *getElem(int *map)=0;
    virtual SommerElement *getAxiElem(int *map);
    virtual int dataSize() { return 0; }
    virtual int copyData(int *p) { return 0; }
};

/*
class PolygonSet {
    int nmax;
    float *v;
    float *nrm;
    int numNrm;
    int *table;

  protected:
    ResizeArray<PolyMesher *> pma;
    int npm;
  public:
    ResizeArray<Polygon *> polys;
    int cur_attrib;
    int curElemNum;
    int isID;

    BlockAlloc loc;

    PolygonSet();
    ~PolygonSet();
    int addLine2(int,int); // OK
    int addTri3(int,int,int);
    int addQuad4(int,int,int,int);
    int addTri6(int,int,int,int,int,int);
    int addQuad9(int,int,int,int,int,int,int,int,int);
    int addQuad8(int,int,int,int,int,int,int,int);

    void normbuild();
    int *indexing() { return table; }
    float *normals() { return nrm; }
    int numNormals() { return numNrm; }
    void setElemNum(int i) { curElemNum = i; }
} ;*/

class PolygonSet {
      int *ndToInterfMap;
      ResizeArray<Polygon *> polys;
      int size;
   public:
      // constructor with the number of interface nodes and an array
      // that maps subdomain number to interface node number or -1 if it
      // is not an interface node. This array must be unmodified until all
      // polygons have been added
      int isID; // dec
      PolygonSet(); //dec
      PolygonSet(int *);
      ~PolygonSet();
      void addLine(Element*,int, int);
      void addTri(Element*,int, int, int);
      void addQuad(Element*,int, int, int, int);
      void addLine(int, int);
      void addTri(int, int, int);
      int findLine(PolyLine2 *);
      int findTri(PolyTri3 *);
      int findQuad(PolyQuad4 *);

      // compute the intersection with another set
      // return the size of the intersection
      int selfIntersect(PolygonSet *neighb);
      int selfIntersect(PolygonSet *neighb, int* solid_flag);
      int selfIntersect(PolygonSetMessage *neighb_m, int* solid_flag);

      // return the corresponding sommer elements
      void getSommerElems(int *map, SommerElement **);
      int getSize() {return size;};
      void getAxiSommerElem(int *map, SommerElement **);

      friend class PolygonSetMessage;
};

class PolygonSetMessage {
      int *message;
public:
      PolygonSetMessage(PolygonSet&);
      PolygonSetMessage(int* _message) { message = _message; }
      int getSizeInInt();
      int * getMessagePtr() { return message; }
      int findLine(PolyLine2 *);
      int findTri(PolyTri3 *);
      int findQuad(PolyQuad4 *);
};


class PolyLine2 : public Polygon {
    int n[2];
  public:
    PolyLine2(Element*,int,int);
    PolyLine2(int,int);
    int isLikeLine(PolyLine2 *);
    int min_node();
    SommerElement *getElem(int *map);
    int findSelf(PolygonSet *);
    int findSelf(PolygonSetMessage *);
    int dataSize() { return 3; }
    int copyData(int *p);
    SommerElement *getAxiElem(int *map);
};


class PolyTri3 : public Polygon {
  protected:
    int n[3];
  public:
    PolyTri3(Element*,int,int,int);
    PolyTri3(int,int,int);
    int isLikeTri(PolyTri3 *);
    int min_node();
    SommerElement *getElem(int *map);
    int findSelf(PolygonSet *);
    int findSelf(PolygonSetMessage *);
    int dataSize() { return 4; }
    int copyData(int *p);
};


class PolyQuad4 : public Polygon {
  protected:
    int n[4];
  public:
    PolyQuad4(Element*,int,int,int,int);
    int isLikeQuad(PolyQuad4 *);
    int min_node();
    SommerElement *getElem(int *map);
    int findSelf(PolygonSet *);
    int findSelf(PolygonSetMessage *);
    int dataSize() { return 5; }
    int copyData(int *p);
};


#endif
