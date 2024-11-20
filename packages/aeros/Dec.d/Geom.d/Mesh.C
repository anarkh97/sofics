#include <cstdio>
#include <Dec.d/Geom.d/Mesh.h>

// static int cheat = 0;

void
Mesh::info()
{
 printf("Numpieces: %d Vertices: %d\n", numpieces, cv);
}

Mesh::Mesh(unsigned char sc[3]) : msz(0), vertices(0), norms(0), swaps(0)
{
 numpieces =0;
 cv = 0;
 step = 0;
 for(int i=0; i < 3 ; ++i)
   scolor[i] = float(sc[i])/float(0xff);
 meshChoice = 3;
}

Mesh::Mesh(float sc[3]) : msz(0), vertices(0), norms(0), swaps(0)
{
 numpieces =0;
 cv = 0;
 step = 0;
 for(int i=0; i < 3 ; ++i)
   scolor[i] = sc[i];
 meshChoice = 3;
}

Mesh::Mesh() : msz(0), vertices(0), norms(0), swaps(0)
{
 numpieces =0;
 cv = 0;
 step = 0;
 for(int i=0; i < 3 ; ++i)
   scolor[i] = 1.0;

 rgb[0] = rgb[1] = rgb[2] = 1.0;
 meshChoice = 0;
}

void
Mesh::set_rgb(float r, float g, float b)
{
  rgb[0] = r; rgb[1] = g; rgb[2] = b;
}

void
Mesh::setUpDateChoice(int choice) 
{ 
  meshChoice = choice; 
}

void
Mesh::add(int vn, int sw)
{
 norms[cv] = vn;
 swaps[cv] = sw;
 cv++;
}

void
Mesh::end_submesh()
{
 msz[numpieces] = cv;
 numpieces++;
}

void
Mesh::apply_index(int *table)
{
 for(int i = 0; i < cv; ++i)
   vertices[i] = table[norms[i]];
}

#ifdef USE_GL
#include "gl.h"
void 
Mesh::gl_draw()
{
 int icv =0;
 int *sw = swaps+0;
  float *v = vtx+0;
  float *n = nrm+0;
  int *vtcs = vertices+0;
  int *nrms = norms+0;
  int *msize = msz+0;

 icv = (cheat == 0) ? 0 : msize[cheat-1];
 int *color = (step > 0) ? sr.getstep(step-1) : 0;
 if(color) {
  lmcolor(LMC_AD);
 }
 else lmcolor(LMC_COLOR);

 for(int mn = cheat; mn < numpieces; ++mn)
  {
   bgntmesh();
   while(icv < msize[mn])
    {
     if(sw[icv]) swaptmesh();
     if(color) cpack(*((unsigned long *)(pal + 4*color[vtcs[icv]])));
     n3f(n+3*nrms[icv]);
     v3f(v+3*vtcs[icv]);
     icv++;
    }
     
   endtmesh();
  }
}
#endif
#ifdef USE_OPENGL
#include <GL/gl.h>
int start = 287;
void
Mesh::GL_draw()
{
 int icv =0;
 int l1,l2;
 int *sw = swaps+0;
  float *v = vtx+0;
  float *n = nrm+0;
  int *vtcs = vertices+0;
  int *nrms = norms+0;
  int *msize = msz+0;

 icv = (cheat == 0) ? 0 : msize[cheat-1];
 int *color = 0;

 if(step > 0) 
   color = (step <= sr.get_num_steps()) ? sr.getstep(step-1) : 
                                               sr.getstep(sr.get_num_steps()-1);

 glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
 glShadeModel(GL_SMOOTH);
 glEnable(GL_LIGHTING);
 glEnable(GL_COLOR_MATERIAL);

 if(meshChoice < 2) glColor4f(rgb[0], rgb[1], rgb[2], 1.0);
 else if(!color) glColor3fv(scolor);

 for(int mn = cheat; mn < numpieces; ++mn)
 {
   glBegin(GL_TRIANGLE_STRIP);
   while(icv < msize[mn])
   {
     if(sw[icv]) {
       l1=l2;
       if(color) glColor3ubv(pal + 4*color[vtcs[l1]]);
       glNormal3fv(n+3*nrms[l1]);
       glVertex3fv(v+3*vtcs[l1]);
     }

     if(color) glColor3ubv(pal + 4*color[vtcs[icv]]);
     glNormal3fv(n+3*nrms[icv]);
     glVertex3fv(v+3*vtcs[icv]);
     l2 = l1;
     l1 = icv;
     icv++;
   }
   glEnd();
 }
 glDisable(GL_COLOR_MATERIAL);
}

#endif

/*
void 
Mesh::set_v(float *v)
{
 ref=vtx=v;
}
*/
void 
Mesh::set_nrm(float *n, int nmax)
{
 nrm = n;
 // 
#ifdef USE_OPENGL
 snrm = new GLshort[3*nmax];
 for(int i=0; i < 3*nmax; ++i)
   snrm[i] = 0x7fff*n[i];
#endif
}

/*
void
Mesh::set_step(int s)
{
 step = s;
 if(s == 0) { vtx = ref; return; }
 if(vr.get_num_steps() != 0) {
   if(s > vr.get_num_steps())
       vtx = vr.getstep(vr.get_num_steps()-1);
   else vtx = vr.getstep(s-1);
 } 
 else
   vtx = ref;
}
*/

/*
void
Mesh::set_sres(ScalRes _sr)
{ 
 sr = _sr;
 pal = sr.get_pal();
}

void
Mesh::set_vres(VecRes _vr)
{
 vr = _vr;
 set_step(step); // make sure vtx is updated
}
*/

/*
int Mesh::get_num_steps()
{
 if(vr.get_num_steps() == 0) return sr.get_num_steps();
  else return vr.get_num_steps();
}

float
Mesh::getTime(int v)
{
 if(v <= 0)  return 0.0;
 if(vr.get_num_steps() == 0){
     if(v > sr.get_num_steps())
       return 0.0;
     else
       return sr.getTimes()[v-1];
 } else
     if(v > vr.get_num_steps())
       return 0.0;
     else
       return vr.getTimes()[v-1];
}

*/
