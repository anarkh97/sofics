#ifndef _DMASSDATA_H_
#define _DMASSDATA_H_

struct DMassData {
 DMassData *next;
 int    node;	// node number of discrete mass
 int    dof;	// dof  number of discrete mass
 double diMass; // value of discrete mass
 int    jdof;   // use for off-diagonal mass, eg. products of inertia (I21, I22, I23)
 DMassData() { jdof = -1; }
};

#endif
