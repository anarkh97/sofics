#include <Driver.d/PolygonSet.h>
#include <Element.d/Sommerfeld.d/QuadSommerBC.h>
#include <Element.d/Sommerfeld.d/TriangleSommerBC.h>
#include <Element.d/Sommerfeld.d/LineSommerBC.h>
#include <Utils.d/dbg_alloca.h>
#include <cstdlib>


PolyLine2::PolyLine2(Element *_el,int n1, int n2)
{
 el = _el;
 flag = 0;
 n[0] = n1;
 n[1] = n2;
}

PolyTri3::PolyTri3(Element *_el,int n1,int n2, int n3)
{
 el = _el;
 flag = 0;
 n[0] = n1; n[1] = n2; n[2] = n3;
}

PolyQuad4::PolyQuad4(Element *_el,int n1,int n2, int n3, int n4)
{
 el = _el;
 flag = 0;
 n[0] = n1; n[1] = n2; n[2] = n3; n[3] = n4;
}

int
Polygon::findSelf(PolygonSet*)
{
 fprintf(stderr,"Error... wrong findSelf\n");
 return 0;
}

PolygonSet::PolygonSet() : polys(0)
{
 isID = 0;
}

/*
SommerElement *
Polygon::getElem(int *)
*/

int
Polygon::findSelf(PolygonSetMessage*)
{
 fprintf(stderr,"Error... wrong findSelf\n");
 return 0;
}


SommerElement *
Polygon::getAxiElem(int *)
{
 return 0;
}


int minof(int nn, int*n)
{
 int ret = n[0];
 int i;
 for(i=1; i < nn; ++i)
  if(n[i] < ret) ret = n[i] ;
 return ret;
}

int Polygon::isLikeLine(PolyLine2 *) { return 0; }
int Polygon::isLikeTri(PolyTri3 *) { return 0; }
int Polygon::isLikeQuad(PolyQuad4 *) { return 0; }


int
PolyLine2::min_node()
{
 return minof(2,n);
}

int PolyLine2::isLikeLine(PolyLine2 *p)
{
 if((n[0] == p->n[0] && n[1] == p->n[1]) ||
    (n[0] == p->n[1] && n[1] == p->n[0])) return 1;
 return 0;
}

SommerElement *
PolyLine2::getElem(int *map)
{
 HelmElement *he = dynamic_cast<HelmElement *>(el);
 if (he==0) {
   fprintf(stderr,"A non-Helmholtz element found.\n");
   exit(-1);
 }
 return new LineSommerBC(map[n[0]], map[n[1]],el);
}

int PolyLine2::findSelf(PolygonSet *ps)
{
 return ps->findLine(this);
}

int PolyLine2::findSelf(PolygonSetMessage *psm)
{
 return psm->findLine(this);
}

int
PolyLine2::copyData(int *p) {
 p[0] = POLYGON_LINE2;
 HelmElement *he = dynamic_cast<HelmElement *>(el);
 if (he==0) {
   fprintf(stderr,"A non-Helmholtz element found.\n");
   exit(-1);
 }
 if (!he->isFluid()) p[0] |= POLYGON_SOLID;
 p[1] = n[0];
 p[2] = n[1];
 return 3;
}

int PolyTri3::findSelf(PolygonSet *ps)
{
 return ps->findTri(this);
}

int PolyTri3::findSelf(PolygonSetMessage *psm)
{
 return psm->findTri(this);
}

SommerElement *
PolyTri3::getElem(int *map)
{
 HelmElement *he = dynamic_cast<HelmElement *>(el);
 if (he==0) {
   fprintf(stderr,"A non-Helmholtz element found.\n");
   exit(-1);
 }
 else {
   return new TriangleSommerBC(map[n[0]], map[n[1]], map[n[2]],el); 
 }
 return 0;
}

int
PolyTri3::isLikeTri(PolyTri3 *p)
{
 if((n[0] == p->n[0] && n[1] == p->n[1] && n[2] == p->n[2]) ||
    (n[1] == p->n[0] && n[2] == p->n[1] && n[0] == p->n[2]) ||
    (n[2] == p->n[0] && n[0] == p->n[1] && n[1] == p->n[2])) return 1;
 if((n[0] == p->n[2] && n[1] == p->n[1] && n[2] == p->n[0]) ||
    (n[1] == p->n[2] && n[2] == p->n[1] && n[0] == p->n[0]) ||
    (n[2] == p->n[2] && n[0] == p->n[1] && n[1] == p->n[0])) return -1;
 return 0;
}

int
PolyTri3::min_node()
{
 return minof(3,n);
}

int
PolyTri3::copyData(int *p) {
 p[0] = POLYGON_TRI3;
 HelmElement *he = dynamic_cast<HelmElement *>(el);
 if (he==0) {
   fprintf(stderr,"PolyTri3::copyData: A non-Helmholtz element found.\n");
   exit(-1);
 }
 if (!he->isFluid()) p[0] |= POLYGON_SOLID;
 p[1] = n[0];
 p[2] = n[1];
 p[3] = n[2];
 return 4;
}

int PolyQuad4::findSelf(PolygonSet *ps)
{
 return ps->findQuad(this);
}

int PolyQuad4::findSelf(PolygonSetMessage *psm)
{
 return psm->findQuad(this);
}

SommerElement *
PolyQuad4::getElem(int *map)
{
 HelmElement *he = dynamic_cast<HelmElement *>(el);
 if (he==0) {
   fprintf(stderr,"PolyQuad4::getElem: A non-Helmholtz element found.\n");
   exit(-1);
 }
 if (!he->isFluid() && (flag&POLYGON_SOLID)) { 
   fprintf(stderr,"PolyQuad4::getElem: Solid Interface element not implemented.\n");
   return 0;
 }
 else if (he->isFluid() && (!(flag&POLYGON_SOLID))) {
   return new QuadSommerBC(map[n[0]], map[n[1]], map[n[2]], map[n[3]],el);
 }
 else {
   fprintf(stderr,"PolyQuad4::getElem: Fluid-Solid Interface element not implemented.\n");
   return 0;
 }
}

int
PolyQuad4::isLikeQuad(PolyQuad4 *p)
{
 if((n[0] == p->n[0] && n[1] == p->n[1] && n[2] == p->n[2] && n[3] == p->n[3]) ||
    (n[1] == p->n[0] && n[2] == p->n[1] && n[3] == p->n[2] && n[0] == p->n[3]) ||
    (n[2] == p->n[0] && n[3] == p->n[1] && n[0] == p->n[2] && n[1] == p->n[3]) ||
    (n[3] == p->n[0] && n[0] == p->n[1] && n[1] == p->n[2] && n[2] == p->n[3]))
 {
   //fprintf (stderr, "\n\n ISLIKEQUAD 1 \n\n");
   return 1;
 }

 if((n[0] == p->n[3] && n[1] == p->n[2] && n[2] == p->n[1] && n[3] == p->n[0]) ||
    (n[1] == p->n[3] && n[2] == p->n[2] && n[3] == p->n[1] && n[0] == p->n[0]) ||
    (n[2] == p->n[3] && n[3] == p->n[2] && n[0] == p->n[1] && n[1] == p->n[0]) ||
    (n[3] == p->n[3] && n[0] == p->n[2] && n[1] == p->n[1] && n[2] == p->n[0]))
 {
   //fprintf (stderr, "\n\n ISLIKEQUAD -1 \n\n");
   return -1;
 }

 return 0;
}

int
PolyQuad4::min_node()
{
 return minof(4,n);
}

int
PolyQuad4::copyData(int *p) {
 p[0] = POLYGON_QUAD4;
 HelmElement *he = dynamic_cast<HelmElement *>(el);
 if (he==0) {
   fprintf(stderr,"A non-Helmholtz element found.\n");
   exit(-1);
 }
 if (!he->isFluid()) p[0] |= POLYGON_SOLID;
 p[1] = n[0];
 p[2] = n[1];
 p[3] = n[2];
 p[4] = n[3];
 return 5;
}

PolygonSet::PolygonSet(int *map) : polys(0)
{
 ndToInterfMap = map;
}

void
PolygonSet::addLine(Element* el,int n1, int n2)
{
 n1 = ndToInterfMap[n1];
 if(n1 < 0) return;
 n2 = ndToInterfMap[n2];
 if(n2 < 0) return;

 PolyLine2 l(el,n1,n2);

 int id = l.min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeLine(&l) != 0)
     return;
   p = p->next;
 }
 Polygon *newp = new PolyLine2(l);
 newp->link(polys[id]);
 polys[id] = newp;
}

void
PolygonSet::addTri(Element* el,int n1, int n2, int n3)
{
 n1 = ndToInterfMap[n1];
 if(n1 < 0) return;
 n2 = ndToInterfMap[n2];
 if(n2 < 0) return;
 n3 = ndToInterfMap[n3];
 if(n3 < 0) return;

 PolyTri3 l(el,n1,n2,n3);
 int id = l.min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeTri(&l) != 0)
     return;
   p = p->next;
 }
 Polygon *newp = new PolyTri3(l);
 newp->link(polys[id]);
 polys[id] = newp;
}

void
PolygonSet::addQuad(Element* el,int n1, int n2, int n3, int n4)
{
 n1 = ndToInterfMap[n1];
 if(n1 < 0) return;
 n2 = ndToInterfMap[n2];
 if(n2 < 0) return;
 n3 = ndToInterfMap[n3];
 if(n3 < 0) return;
 n4 = ndToInterfMap[n4];
 if(n4 < 0) return;

 PolyQuad4 l(el,n1,n2,n3,n4);
 int id = l.min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeQuad(&l) != 0)
     return;
   p = p->next;
 }
 Polygon *newp = new PolyQuad4(l);
 newp->link(polys[id]);
 polys[id] = newp;
}

int
PolygonSet::findLine(PolyLine2 *l)
{
 int id = l->min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeLine(l) != 0) {
     int ret = POLYGON_LINE2;
     HelmElement *he = dynamic_cast<HelmElement *>(p->el);
     if (he==0) {
       fprintf(stderr,"A non-Helmholtz element found.\n");
       exit(-1);
     }
     if (!he->isFluid()) ret |= POLYGON_SOLID;
     return ret;
   }
   p = p->next;
 }
 return 0;
}

int
PolygonSetMessage::findLine(PolyLine2 *l) {
 int len = message[0];
 int *offset = message + 1;
 int id = l->min_node();
 if (id >= len) return 0;
 int i;
 for(i=offset[id]+len+2;i<offset[id+1]+len+2;) {
   if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_LINE2) {
     PolyLine2 p(0,message[i+1],message[i+2]);
     if(p.isLikeLine(l) != 0) { 
       int ret = POLYGON_LINE2;
       if (message[i] & POLYGON_SOLID) ret |= POLYGON_SOLID;
       return ret;
     }
     i += p.dataSize();
   } else if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_TRI3) {
     PolyTri3 p(0,message[i+1],message[i+2],message[i+3]);
     i += p.dataSize();
   } else if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_QUAD4) {
     PolyQuad4 p(0,message[i+1],message[i+2],message[i+3],message[i+4]);
     i += p.dataSize();
   } else fprintf(stderr,"Error in PolygonSetMessage::findLine.\n");
 }
 return 0;
}

int
PolygonSet::findTri(PolyTri3 *l)
{
 int id = l->min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeTri(l)) {
     int ret = POLYGON_TRI3;
     HelmElement *he = dynamic_cast<HelmElement *>(p->el);
     if (he==0) {
       fprintf(stderr,"A non-Helmholtz element found.\n");
       exit(-1);
     }
     if (!he->isFluid()) ret |= POLYGON_SOLID;
     return ret;
   }
   p = p->next;
 }
 return 0;
}

int
PolygonSetMessage::findTri(PolyTri3 *l) {
 int sz=getSizeInInt();
 int len = message[0];
 int *offset = message + 1;
 int id = l->min_node();
 int i;
 for(i=offset[id]+len+2;i<offset[id+1]+len+2;) {
   if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_LINE2) {
     PolyLine2 p(0,message[i+1],message[i+2]);
     i += p.dataSize();
   } else if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_TRI3) {
     PolyTri3 p(0,message[i+1],message[i+2],message[i+3]);
     if(p.isLikeTri(l) != 0) { 
       int ret = POLYGON_TRI3;
       if (message[i] & POLYGON_SOLID) ret |= POLYGON_SOLID;
       return ret;
     }
     i += p.dataSize();
   } else if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_QUAD4) {
     PolyQuad4 p(0,message[i+1],message[i+2],message[i+3],message[i+4]);
     i += p.dataSize();
   } else fprintf(stderr,"Error in PolygonSetMessage::findTri %d.\n",sz);
 }
 return 0;
}

int
PolygonSet::findQuad(PolyQuad4 *l)
{
 int id = l->min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeQuad(l)) {
     int ret = POLYGON_QUAD4;
     HelmElement *he = dynamic_cast<HelmElement *>(p->el);
     if (he==0) {
       fprintf(stderr,"A non-Helmholtz element found.\n");
       exit(-1);
     }
     if (!he->isFluid()) ret |= POLYGON_SOLID;
     return ret;
   }
   p = p->next;
 }
 return 0;
}

int
PolygonSetMessage::findQuad(PolyQuad4 *l) {
 int len = message[0];
 int *offset = message + 1;
 int id = l->min_node();
 int i;
 for(i=offset[id]+len+2;i<offset[id+1]+len+2;) {
   if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_LINE2) {
     PolyLine2 p(0,message[i+1],message[i+2]);
     i += p.dataSize();
   } else if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_TRI3) {
     PolyTri3 p(0,message[i+1],message[i+2],message[i+3]);
     i += p.dataSize();
   } else if ((message[i] & POLYGON_TYPE_MASK) == POLYGON_QUAD4) {
     PolyQuad4 p(0,message[i+1],message[i+2],message[i+3],message[i+4]);
     if(p.isLikeQuad(l) != 0) { 
       int ret = POLYGON_QUAD4;
       if (message[i] & POLYGON_SOLID) ret |= POLYGON_SOLID;
       return ret;
     }
     i += p.dataSize();
   } else fprintf(stderr,"Error in PolygonSetMessage::findQuad.\n");
 }
 return 0;
}


int
PolygonSet::selfIntersect(PolygonSet *neighb, int *solid_flag)
{
 *solid_flag = 0;
 int nnp = 0;
 size = 0;
 int maxNode = polys.max_size();
 int i;
 for(i = 0; i < maxNode; ++i) {
   Polygon *p = polys[i];
   while(p) {
     nnp += 1;
     p->flag = p->findSelf(neighb);
     if (p->flag & POLYGON_TYPE_MASK) size++;
     if (p->flag & POLYGON_SOLID) *solid_flag = 1;
     p = p->next;
   }
 }
 return size;
}


int
PolygonSet::selfIntersect(PolygonSetMessage *neighb_m, int *solid_flag)
{
 *solid_flag = 0;
 int nnp = 0;
 size = 0;
 int maxNode = polys.max_size();
 for(int i = 0; i < maxNode; ++i) {
   Polygon *p = polys[i];
   while(p) {
     nnp += 1;
     p->flag = p->findSelf(neighb_m);
     if (p->flag & POLYGON_TYPE_MASK) size++;
     if (p->flag & POLYGON_SOLID) *solid_flag = 1;
     p = p->next;
   }
 }
 //fprintf(stderr, "My size is %d. Intersect size is %d\n", nnp,size);
 return size;
}


void
PolygonSet::getSommerElems(int *map, SommerElement **se) {
 int maxNode = polys.max_size(); 
 int index =0;
 int i;
 for(i = 0; i < maxNode; ++i) {
   Polygon *p = polys[i];
   while(p) {
     if(p->flag & POLYGON_TYPE_MASK) {
       se[index++] = p->getElem(map);
     }
     p = p->next;
   }
 }
}

PolygonSet::~PolygonSet() {
 int i;
 for(i=0;i<polys.max_size();i++) {
   Polygon *p = polys[i];
   while(p) {
     Polygon *pp = p;
     p = p->next;
     delete pp;
   }
 }
}


PolygonSetMessage::PolygonSetMessage(PolygonSet &ps) {
 // First find the size of the set
 int size = 0;
 int *offset = (int*)dbg_alloca(sizeof(int)*(ps.polys.max_size()+1));
 int i;
 for(i=0;i<ps.polys.max_size();i++) {
   Polygon *p = ps.polys[i];
   offset[i] = size;
   while(p) {
     size += p->dataSize();
     p = p->next;
   }
 }
 offset[ps.polys.max_size()] = size;

 int * msg = new int[1+ps.polys.max_size()+1+size];
 msg[0] = ps.polys.max_size();
 for(i=0;i<ps.polys.max_size()+1;i++) msg[1+i] = offset[i];

 int cur_offset = 1+ps.polys.max_size()+1;
 for(i=0;i<ps.polys.max_size();i++) {
   Polygon *p = ps.polys[i];
   while(p) {
     int s = p->copyData(msg+cur_offset);
     cur_offset += s;
     p = p->next;
   }
 }

 message = msg;

}

int
PolygonSetMessage::getSizeInInt() {
 int len = message[0];
 int last_offset = message[len+1];
 return len+2+last_offset;
}



void
PolygonSet::getAxiSommerElem(int *map, SommerElement **se) {
 int maxNode = polys.max_size();
 int index =0;
 int i;
 for(i = 0; i < maxNode; ++i) {
   Polygon *p = polys[i];
   while(p) {
     if(p->flag) 
       se[index++] = p->getAxiElem(map);
     p = p->next;
   }
 }
}


SommerElement *
PolyLine2::getAxiElem(int *map)
{
    throw std::logic_error("Axi Helm is not supported anymore");
}

PolyLine2::PolyLine2(int n1, int n2)
{
 n[0] = n1;
 n[1] = n2;
}

PolyTri3::PolyTri3(int n1,int n2, int n3)
{
 n[0] = n1; n[1] = n2; n[2] = n3;
}

void
PolygonSet::addLine(int n1, int n2)
{
 n1 = ndToInterfMap[n1];
 if(n1 < 0) return;
 n2 = ndToInterfMap[n2];
 if(n2 < 0) return;

 PolyLine2 l(n1,n2);

 int id = l.min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeLine(&l) != 0)
     return;
   p = p->next;
 }
 Polygon *newp = new PolyLine2(l);
 newp->link(polys[id]);
 polys[id] = newp;
}

void
PolygonSet::addTri(int n1, int n2, int n3)
{
 n1 = ndToInterfMap[n1];
 if(n1 < 0) return;
 n2 = ndToInterfMap[n2];
 if(n2 < 0) return;
 n3 = ndToInterfMap[n3];
 if(n3 < 0) return;

 PolyTri3 l(n1,n2,n3);
 int id = l.min_node();
 Polygon *p = polys[id];
 while(p) {
   if(p->isLikeTri(&l) != 0)
     return;
   p = p->next;
 }
 Polygon *newp = new PolyTri3(l);
 newp->link(polys[id]);
 polys[id] = newp;
}

int
PolygonSet::selfIntersect(PolygonSet *neighb)
{
 int nnp = 0;
 size = 0;
 int maxNode = polys.max_size();
 int i;
 for(i = 0; i < maxNode; ++i) {
   Polygon *p = polys[i];
   while(p) {
     nnp += 1;
     p->flag = p->findSelf(neighb);
     if(p->flag) size++;
     p = p->next;
   }
 }
 return size;
}

