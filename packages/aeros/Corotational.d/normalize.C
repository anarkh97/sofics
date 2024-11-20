#include <cmath>
#include <Corotational.d/utilities.h>

// Coded by: Kendall H. Pierson

#ifdef NO_INLINE_COROT_UTILS
// This routine normalizes a 3 dimensional vector
void normalize(double v[3])
{
	double length = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0] /= length;
	v[1] /= length;
	v[2] /= length;
}

double magnitude(double v[3])
{
	return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}
#endif

void dnormalize(double v[3],double dv[3])
{ 
  double length = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  double dlength;

  if (length > 0) {

    dlength =(dv[0]*v[0]+dv[1]*v[1]+dv[2]*v[2])/length;

    dv[0]= (dv[0]*length-v[0]*dlength)/(length*length);
    dv[1]= (dv[1]*length-v[1]*dlength)/(length*length);
    dv[2]= (dv[2]*length-v[2]*dlength)/(length*length);
  
    v[0] /= length;
    v[1] /= length;
    v[2] /= length;
  }
  else {

    dv[0] = 0;
    dv[1] = 0;
    dv[2] = 0;

    //dlength =(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
  
    //if (dlength > 0) {
    //  dv[0] /= dlength; 
    //  dv[1] /= dlength;
    //  dv[2] /= dlength;
    //}
  }
}
