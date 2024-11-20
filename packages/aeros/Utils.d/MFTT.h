#ifndef _MFTT_H_
#define _MFTT_H_

#include <vector>

struct DoubleList;

// Mechanical Force Time Table data
// also used for young's modulus vs temperature table, etc.
class MFTTData {
 public:
     int     id;
     int     np;  // number of points
     double* time;
     double* value;
     int maxval;
     int curp; // Because calls to getVal will usually be
               // for neighboring values of t
   public:
     MFTTData();
     MFTTData(int);
     void add(double, double);
     double getVal(double);
     double getValAlt(double);
     void getValAndSlopeAlt(double t, double *v, double *s);
     void getValAndSlopeAlt2(double t, double *v, double *s);
     int getID() { return id; }
     int getNumPoints() { return np; }
     double getT(int i) { return time[i]; }
     double getV(int i) { return value[i]; }
};

template<typename DataType>
class GenMFTTData {
 public:
     int     id;
     int     np;  // number of points
     double* time;
     DataType* value;
     int maxval;
     int curp; // Because calls to getVal will usually be
               // for neighboring values of t
     static const DataType zero;
     static const DataType one;

   public:
     GenMFTTData();
     GenMFTTData(int);
     void add(double, DataType);
     DataType getVal(double);
     DataType getValAlt(double);
     void getValAndSlopeAlt(double t, DataType *v, DataType *s);
     void getValAndSlopeAlt2(double t, DataType *v, DataType *s);
     int getID() { return id; }
     int getNumPoints() { return np; }
     double getT(int i) { return time[i]; }
     DataType getV(int i) { return value[i]; }
};

class SS2DTData {
 public:
     int id;
     bool engineeringFlag;
     std::vector<double> x, y;
     std::vector<std::vector<double> > value;
   public:
     SS2DTData(int, DoubleList&, bool);
     void add(double, DoubleList&);
     double getValAlt(double x, double y);
     void getValAndSlopeAlt(double x, double y, double *v, double *sx, double *sy);
     int getID() { return id; }
     int getNumPoints() { return x.size(); }
};

#ifdef _TEMPLATE_FIX_
#include <Utils.d/GenMFTT.C>
#endif

#endif
