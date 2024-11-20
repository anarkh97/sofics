#include <Utils.d/MFTT.h>
#include <Driver.d/GeoSource.h>

#include <algorithm>
#include <iostream>

#ifdef USE_EIGEN3
#include <Eigen/Core>
template<> const Eigen::Vector4d GenMFTTData<Eigen::Vector4d>::zero = Eigen::Vector4d::Zero();
template<> const Eigen::Vector4d GenMFTTData<Eigen::Vector4d>::one  = Eigen::Vector4d::Ones();
#endif

MFTTData::MFTTData()
{
 maxval = 16;
 np     = 0;
 curp   = 0;
 time   = new double[maxval];
 value  = new double[maxval];
}

MFTTData::MFTTData(int _id)
{
 maxval = 16;
 np     = 0;
 curp   = 0;
 time   = new double[maxval];
 value  = new double[maxval];
 id = _id;
}

void
MFTTData::add(double t, double v)
{
 if(np == maxval) {
   int n_maxval = 2*maxval;
   double *n_time  = new double[n_maxval];
   double *n_value = new double[n_maxval];
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

double
MFTTData::getVal(double t)
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

   if (curp < 0 )
      return 0;
   else if (curp == np - 1)
      return 0;
   else {
      double v1 = value[curp], v2 = value[curp+1];
      double t1 =  time[curp], t2 =  time[curp+1];
      return  v1 + (v2-v1) / (t2 - t1) * ( t - t1);
   }
 }
 else
   return 1;
}

double
MFTTData::getValAlt(double t)
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

   if(curp < 0 ) 
     return value[0];
   else if (curp == np - 1) 
     return value[curp];
   else {
     double v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     if(t2 == t1) return v1;
     else return  v1 + (v2-v1) / (t2 - t1) * ( t - t1);
   }
 }
 else 
   return 0.0;
}

void 
MFTTData::getValAndSlopeAlt(double t, double *v, double *s)
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
     *s = 0.0;
   }
   else if (curp == np - 1) {
     *v = value[curp];
     *s = 0.0;
   }
   else {
     double v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     *v = v1 + (v2-v1) / (t2 - t1) * ( t - t1);
     *s = (v2-v1) / (t2 - t1);
   }
 }
 else { *v = 0.0; *s = 0.0; }
}

void
MFTTData::getValAndSlopeAlt2(double t, double *v, double *s)
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
     *s = 0.0;
   }
   else if (curp == np - 1) {
     double v1 = value[curp-1], v2 = value[curp];
     double t1 =  time[curp-1], t2 =  time[curp];
     *v = v2 + (v2-v1) / (t2 - t1) * ( t - t2);
     *s = (v2-v1) / (t2 - t1);
   }
   else {
     double v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     *v = v1 + (v2-v1) / (t2 - t1) * ( t - t1);
     *s = (v2-v1) / (t2 - t1);
   }
 }
 else { *v = 0.0; *s = 0.0; }
}

SS2DTData::SS2DTData(int _id, DoubleList& _y, bool _engineeringFlag)
{
  id = _id;
  engineeringFlag = _engineeringFlag;
  for(int j = 0; j < _y.nval; ++j) y.push_back(_y.v[j]);
  // TODO: check if y is sorted
}

void
SS2DTData::add(double _x, DoubleList& _value)
{
  if(_value.nval != y.size()) { std::cerr << " *** Error reading SS2DT data.\n"; exit(-1); }
  x.push_back(_x); // TODO check if x is sorted
  value.resize(value.size()+1);
  for(int j = 0; j < _value.nval; ++j) value.back().push_back(_value.v[j]);
}

double
SS2DTData::getValAlt(double xP, double yP)
{
  std::vector<double>::iterator it1 = std::upper_bound(x.begin(), x.end(), xP),
                                it2 = std::upper_bound(y.begin(), y.end(), yP);
  if(it1 == x.begin()) {
    if(it2 == y.begin()) return value.front().front();
    else if(it2 == y.end()) return value.front().back();
    else {
      int j = std::distance(y.begin(), it2);
      const double& yA = *(it2-1);
      const double& yB = *it2;
      const double eta = (yP - yA) / (yB - yA);
      return (1-eta)*value.front()[j-1] + eta*value.front()[j];
    }
  }
  else if(it1 == x.end()) {
    if(it2 == y.begin()) return value.back().front();
    else if(it2 == y.end()) return value.back().back();
    else {
      int j = std::distance(y.begin(), it2);
      const double& yA = *(it2-1);
      const double& yB = *it2;
      const double eta = (yP - yA) / (yB - yA);
      return (1-eta)*value.back()[j-1] + eta*value.back()[j];
    }
  }
  else {
    int i = std::distance(x.begin(), it1);
    const double& xA = *(it1-1);
    const double& xB = *it1;
    const double xi = (xP - xA) / (xB - xA);
    if(it2 == y.begin()) {
      return (1-xi)*value[i-1].front() + xi*value[i].front();
    }
    else if(it2 == y.end()) {
      return (1-xi)*value[i-1].back() + xi*value[i].back();
    }
    else {
      int j = std::distance(y.begin(), it2);
      const double& yA = *(it2-1);
      const double& yB = *it2;
      const double eta = (yP - yA) / (yB - yA);
      return (1-xi)*(1-eta)*value[i-1][j-1] + xi*(1-eta)*value[i][j-1] + xi*eta*value[i][j] + (1-xi)*eta*value[i-1][j];
    }
  }
}

void
SS2DTData::getValAndSlopeAlt(double xP, double yP, double *v, double *sx, double *sy)
{
  std::vector<double>::iterator it1 = std::upper_bound(x.begin(), x.end(), xP),
                                it2 = std::upper_bound(y.begin(), y.end(), yP);
  if(it1 == x.begin()) {
    if(it2 == y.begin()) { *v = value.front().front(); *sx = *sy = 0; }
    else if(it2 == y.end()) { *v = value.front().back(); *sx = *sy = 0; }
    else {
      int j = std::distance(y.begin(), it2);
      const double& yA = *(it2-1);
      const double& yB = *it2;
      const double eta = (yP - yA) / (yB - yA);
      *v = (1-eta)*value.front()[j-1] + eta*value.front()[j];
      *sx = 0;
      *sy = (-value.front()[j-1] + value.front()[j]) / (yB - yA);
    }
  }
  else if(it1 == x.end()) {
    if(it2 == y.begin()) { *v = value.back().front(); *sx = *sy = 0; }
    else if(it2 == y.end()) { *v = value.back().back(); *sx = *sy = 0; }
    else {
      int j = std::distance(y.begin(), it2);
      const double& yA = *(it2-1);
      const double& yB = *it2;
      const double eta = (yP - yA) / (yB - yA);
      *v = (1-eta)*value.back()[j-1] + eta*value.back()[j];
      *sx = 0;
      *sy = (-value.back()[j-1] + value.back()[j]) / (yB - yA);
    }
  }
  else {
    int i = std::distance(x.begin(), it1);
    const double& xA = *(it1-1);
    const double& xB = *it1;
    const double xi = (xP - xA) / (xB - xA);
    if(it2 == y.begin()) {
      *v = (1-xi)*value[i-1].front() + xi*value[i].front();
      *sx = (-value[i-1].front() + value[i].front()) / (xB - xA);
      *sy = 0;
    }
    else if(it2 == y.end()) {
      *v = (1-xi)*value[i-1].back() + xi*value[i].back();
      *sx = (-value[i-1].back() + value[i].back()) / (xB - xA);
      *sy = 0;
    }
    else {
      int j = std::distance(y.begin(), it2);
      const double& yA = *(it2-1);
      const double& yB = *it2;
      const double eta = (yP - yA) / (yB - yA);
      *v = (1-xi)*(1-eta)*value[i-1][j-1] + xi*(1-eta)*value[i][j-1] + xi*eta*value[i][j] + (1-xi)*eta*value[i-1][j];
      *sx = (-(1-eta)*value[i-1][j-1] + (1-eta)*value[i][j-1] + eta*value[i][j] - eta*value[i-1][j]) / (xB - xA);
      *sy = (-(1-xi)*value[i-1][j-1] - xi*value[i][j-1] + xi*value[i][j] + (1-xi)*value[i-1][j]) / (yB - yA);
    }
  }
}
