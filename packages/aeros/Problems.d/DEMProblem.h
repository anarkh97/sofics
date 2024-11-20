#include <Element.d/DEM.d/DEMElement.h>

class LMInput {
public:
 int en;
 int fi;
 int type;
};

class DEM {
public:
 DEM() { nlmi = 0; }
 int nlmi;
 LMInput *lmi;
 DEMLM* createLM(int type);
 void run(Domain *d, GeoSource *g);
};

class FaceInfo {
public:
 FaceInfo() { off = bcf = en2 = fi2 = -1; }
 int off;
 int bcf;
 int en2;
 int fi2;
};

class DEMFace {
public:
 DEMFace() { cornernodes = 0; }
 int en;
 int fi;
 int nc;
 int *cornernodes; 
 int isFluid;
};
