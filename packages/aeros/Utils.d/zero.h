#ifndef FETI_ZERO_H_
#define FETI_ZERO_H_

/* this contains routines to set vectors to zero. In many cases they can
   be used to set matrix elements to zero. All have the same format - they
   are called with overloaded calls.
*/

#include <string>

inline double* zero( double *data, int num ){
  memset( data, 0, num*sizeof(double) );
  return data; }

inline float* zero( float *data, int num ){
  memset( data, 0, num*sizeof(float) );
  return data; }

inline int* zero( int *data, int num ) {
  memset( data, 0, num*sizeof(int) );
  return data; }

#endif

