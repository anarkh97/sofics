#ifndef _ELEMENTARY_FUNCTION_
#define _ELEMENTARY_FUNCTION_

#include <cmath>

class ElementaryFunction
{
  int type;
  double scale, shift;
  double a, b, c, d;

  public:
    ElementaryFunction(int _type, double _scale, double _shift, double _a, double _b, double _c = 0, double _d = 0)
      : type(_type), scale(_scale), shift(_shift), a(_a), b(_b), c(_c), d(_d) {}

    double operator()(double t) {
      using std::sin; 
      using std::pow;

      double f;
      switch(type) {
        case 0: { // sine
          f = sin(a*t + b);
        } break;
        case 1: { // bounded ramp
          if     (t <= a) f = 0;
          else if(t <= b) f = (t-a)/(b-a);
          else            f = 1;
        } break;
        case 2: { // triangle  
          if     (t <= a) f = 0;
          else if(t <= b) f = (t-a)/(b-a);
          else if(t <= c) f = (c-t)/(c-b);
          else            f = 0;
        } break;
        case 3: { // trapezoidal  
          if     (t <= a) f = 0;
          else if(t <= b) f = (t-a)/(b-a);
          else if(t <= c) f = 1;
          else if(t <= d) f = (d-t)/(d-c);
          else            f = 0;
        } break;
        case 4: { // S-shaped function
          if     (t <= a) f = 0;
          else if(t <= b) f = 0.5*pow((t-a)/(b-a),2);
          else if(t <= c) f = 1-0.5*pow((c-t)/(c-b),2);
          else            f = 1;
        } break;
      }

      return scale*f + shift;
    }

    double firstDerivative(double t) {
      using std::cos; 
      using std::pow;

      double dfdt;
      switch(type) {
        case 0: { // sine
          dfdt = a*cos(a*t + b);
        } break;
        case 1: { // bounded ramp
          if     (t <= a) dfdt = 0;
          else if(t <= b) dfdt = 1/(b-a);
          else            dfdt = 0;
        } break;
        case 2: { // triangle  
          if     (t <= a) dfdt = 0;
          else if(t <= b) dfdt = 1/(b-a);
          else if(t <= c) dfdt = -1/(c-b);
          else            dfdt = 0;
        } break;
        case 3: { // trapezoidal  
          if     (t <= a) dfdt = 0;
          else if(t <= b) dfdt = 1/(b-a);
          else if(t <= c) dfdt = 1;
          else if(t <= d) dfdt = -1/(d-c);
          else            dfdt = 0;
        } break;
        case 4: { // S-shaped function
          if     (t <= a) dfdt = 0;
          else if(t <= b) dfdt = (t-a)/pow(a-b,2);
          else if(t <= c) dfdt = (c-t)/pow(b-c,2);
          else            dfdt = 0;
        } break;
      }

      return scale*dfdt;
    }

    double secondDerivative(double t) {
      using std::sin;
      using std::pow;

      double d2fdt2;
      switch(type) {
        case 0: { // sine
          d2fdt2 = -pow(a,2)*sin(a*t + b);
        } break;
        case 1: { // bounded ramp
          d2fdt2 = 0;
        } break;
        case 2: { // triangle  
          d2fdt2 = 0;
        } break;
        case 3: { // trapezoidal  
          d2fdt2 = 0;
        } break;
        case 4: { // S-shaped function
          if     (t <= a) d2fdt2 = 0;
          else if(t <= b) d2fdt2 = 1/pow(a-b,2);
          else if(t <= c) d2fdt2 = -1/pow(b-c,2);
          else            d2fdt2 = 0;
        } break;
      }

      return scale*d2fdt2;
    }

};
#endif
