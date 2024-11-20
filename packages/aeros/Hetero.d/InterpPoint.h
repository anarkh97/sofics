#ifndef _INTERP_POINT_H_
#define _INTERP_POINT_H_
// This is an inerpolation point. The element elemNum uses x and y to
// compute the interpolated displacements/velocities
 
struct InterpPoint {
    int subNumber;
    int elemNum;
    double xy[2];
    double gap[3];
    int *dofs;
};

#endif
