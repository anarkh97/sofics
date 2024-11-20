#include <Utils.d/MFTT.h>
#include <iostream>

template<typename DataType>
GenMFTTData<DataType>::GenMFTTData()
{
 maxval = 16;
 np     = 0;
 curp   = 0;
 time   = new double[maxval];
 value  = new DataType[maxval];
}

template<typename DataType>
GenMFTTData<DataType>::GenMFTTData(int _id)
{
 maxval = 16;
 np     = 0;
 curp   = 0;
 time   = new double[maxval];
 value  = new DataType[maxval];
 id = _id;
}

template<typename DataType>
void
GenMFTTData<DataType>::add(double t, DataType v)
{
 if(np == maxval) {
   int n_maxval = 2*maxval;
   double *n_time  = new double[n_maxval];
   DataType *n_value = new DataType[n_maxval];
   int i;
   for(i=0; i < maxval; ++i) {
     n_time[i] = time[i];
     n_value[i] = value[i];
   }
  delete[] time;
  delete[] value;
  time   = n_time;
  value  = n_value;
  maxval = n_maxval;
 }
 time[np] = t;
 value[np] = v;
 np++;
}

template<typename DataType>
DataType
GenMFTTData<DataType>::getVal(double t)
{
 // This function returns zero if t < t_min or t > t_max
 
 // np = total number of points
 if (np) {

 // interpolate and make sure we deal with the case 
 // curp = current point
 // curp==np-1 and curp=-1 correctly

 if(time[curp] > t) // Reverse the search
   while(curp >= 0 && time[curp] > t) curp--;
 else
   while(curp < np-1 && time[curp+1] < t) curp++;

   if ((curp < 0) || (curp == np - 1)) {
      return zero;
   }
   else {
      DataType v1 = value[curp], v2 = value[curp+1];
      double t1 =  time[curp], t2 =  time[curp+1];
      return v1 + (v2-v1) / (t2 - t1) * ( t - t1);
   }
 }
 else {
   return one;
 }
}

template<typename DataType>
DataType
GenMFTTData<DataType>::getValAlt(double t)
{
 // This function returns the value corresponding to t_min if t < t_min,
 // and the value corresponding to t_max if t > t_max

 // np = total number of points

 if (np) {
   curp = np/2; // starting point of search

   // interpolate and make sure we deal with the case
   // curp = current point
   // curp==np-1 and curp=-1 correctly

   if(time[curp] > t) // Reverse the search
     while(curp >= 0 && time[curp] > t) curp--;
   else
     while(curp < np-1 && time[curp+1] < t) curp++;

   if(curp < 0) 
     return value[0];
   else if (curp == np - 1) 
     return value[curp];
   else {
     DataType v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     if(t2 == t1) return v1;
     else return  v1 + (v2-v1) / (t2 - t1) * ( t - t1);
   }
 }
 else {
   return zero;
 }
}

template<typename DataType>
void 
GenMFTTData<DataType>::getValAndSlopeAlt(double t, DataType *v, DataType *s)
{
 // This function returns the value corresponding to t_min if t < t_min,
 // and the value corresponding to t_max if t > t_max
 // The slope is set to zero if t < t_min or t > t_max

 // np = total number of points

 if (np) {
   curp = np/2; // starting point of search

   // interpolate and make sure we deal with the case
   // curp = current point
   // curp==np-1 and curp=-1 correctly

   if(time[curp] > t) // Reverse the search
     while(curp >= 0 && time[curp] > t) curp--;
   else
     while(curp < np-1 && time[curp+1] < t) curp++;

   if(curp < 0 ) {
     *v = value[0];
     *s = zero;
   }
   else if (curp == np - 1) {
     *v = value[curp];
     *s = zero;
   }
   else {
     DataType v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     *v = v1 + (v2-v1) / (t2 - t1) * (t - t1);
     *s = (v2-v1) / (t2 - t1);
   }
 }
 else { *v = zero; *s = zero; }
}

template<typename DataType>
void
GenMFTTData<DataType>::getValAndSlopeAlt2(double t, DataType *v, DataType *s)
{
 // This function returns the value corresponding to t_min if t < t_min,
 // and extrapolates from the value corresponding to t_max if t > t_max
 // The slope is set to zero if t < t_min, and set to the slope at
 // t_max if t > t_max

 // np = total number of points

 if (np) {
   curp = np/2; // starting point of search

   // interpolate and make sure we deal with the case
   // curp = current point
   // curp==np-1 and curp=-1 correctly

   if(time[curp] > t) // Reverse the search
     while(curp >= 0 && time[curp] > t) curp--;
   else
     while(curp < np-1 && time[curp+1] < t) curp++;

   if(curp < 0 ) {
     *v = value[0];
     *s = zero;
   }
   else if (curp == np - 1) {
     DataType v1 = value[curp-1], v2 = value[curp];
     double t1 =  time[curp-1], t2 =  time[curp];
     *v = v2 + (v2-v1) / (t2 - t1) * (t - t2);
     *s = (v2-v1) / (t2 - t1);
   }
   else {
     DataType v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     *v = v1 + (v2-v1) / (t2 - t1) * (t - t1);
     *s = (v2-v1) / (t2 - t1);
   }
 }
 else { *v = zero; *s = zero; }
}

