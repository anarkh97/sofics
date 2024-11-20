#ifndef _MESH_H_
#define _MESH_H_

#include <Utils.d/resize_array.h>

class Mesh {
  int numpieces, cv;
  ResizeArray<int> msz;
  ResizeArray<int> vertices;
  ResizeArray<int> norms;
  ResizeArray<int> swaps;
  float *vtx, *nrm;
  float *ref;
  int step;
  unsigned char *pal;
  float scolor[3];

  float rgb[3];
  int meshChoice;
 
 public:
  Mesh();
  Mesh(unsigned char[3]);
  Mesh(float[3]);
  int nlast() { return norms[cv-1]; }
  void add(int,int);
  void apply_index(int *);
  void end_submesh();
  void info();
  void set_rgb(float, float, float);
  void setUpDateChoice(int);
  void setRes_Choice(int);
  void set_nrm(float *,int);
  void set_step(int s);
  int get_num_steps();
  float getTime(int s);
};
#endif
