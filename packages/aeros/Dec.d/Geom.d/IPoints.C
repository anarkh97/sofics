#include <cstdio>
#include <Dec.d/Geom.d/IPoints.h>
#include <Utils.d/Connectivity.h>

 
int IPointset::num(int n1, int n2, int attr)
{
 int smallest,biggest ; // smallest node number of the 2
 if(n1 < n2)
   {
    smallest = n1 ;
    biggest = n2 ;
   }
 else
   {
    smallest = n2 ;
    biggest = n1 ;
   }
 if(biggest >= nodemax) nodemax = biggest+1;
 
 // try to find it among known points
 for(IPointLink * link = ipl[smallest] ;
  link != 0;
 link = link->next)
  if(link->to == biggest && link->attrib == attr)
    return link->num ;
 
 // When not found, add a new one
 IPointLink *newlink = new (blkal)
      IPointLink(smallest, biggest,np++,ipl[smallest], attr);
 ipl[smallest] = newlink ;
 return np-1 ;
}

int IPointset::num(IPoint &ip)
{
 return num(ip.from,ip.to) ;
}
 
IPointset::IPointset(int n) : ipl(0,n+1)
 {
 np = 0 ;
 for(int i =0; i < n+1;++i)
   ipl[i] = 0 ;
 nodemax = 0;
 }
 
IPointset::~IPointset()
{
}
 
float (*
IPointset::coord(float *rv, float *nval, float *disp) )[3]
{
 float (*crd)[3] = new float [np][3] ;
 for(int fn =0; fn < ipl.max_size(); ++ fn)
   {
    for(IPointLink *cipl = ipl[fn]; cipl ; cipl = cipl->next)
      {
       float *cxyz = crd[cipl->num] ;
       float alpha = (rv[cipl->attrib] - nval[fn])/(nval[cipl->to]-nval[fn]) ;
       cxyz[0] = (1-alpha)* disp[fn*3] + alpha*disp[cipl->to*3] ;
       cxyz[1] = (1-alpha)* disp[fn*3+1] + alpha*disp[cipl->to*3+1] ;
       cxyz[2] = (1-alpha)* disp[fn*3+2] + alpha*disp[cipl->to*3+2] ;
      }
   }
 return crd ;
}

CountedConnectivity *
IPointset::two_way_connect()
{
 CountedConnectivity *cc = new CountedConnectivity(nodemax);

 int i;
 for(i=0; i < nodemax; ++i) {
  IPointLink *ip;
  ip = ipl[i];
  while(ip) {
   cc->countlink(ip->from,ip->to);
   cc->countlink(ip->to,ip->from);
   ip = ip->next;
  }
 }

 cc->end_count();

 for(i=0; i < nodemax; ++i) {
  IPointLink *ip;
  ip = ipl[i];
  while(ip) {
   cc->addlink(ip->from,ip->to);
   cc->addlink(ip->to,ip->from);
   ip = ip->next;
  }
 }

 return cc;
}
