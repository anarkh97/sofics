// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University 
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.

    
#ifndef ContactBoundingBox_h_
#define ContactBoundingBox_h_
 
#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif
#include <algorithm>
#include <iostream>
#include <cmath>
#include "Contact_Defines.h"
typedef float RealBB;
#define BIGNUM_BB 1.E38

class ContactParOStream;

/* 
   Class for a simple axis aligned bounding box.  Operation to find if
   objects are withing the box, update one box to include another, etc. 
   These objects will be used at a very low level, so are made very
   lightweight
*/

class ContactBoundingBox {

 public:
  //
  //  Default box constructor.  This constructor will creat a box that encompasses the entire computational domain
  //
  ContactBoundingBox();
  //
  //  Create a box and explicitly initialize it's corners to the input triplets.
  //
  ContactBoundingBox(const RealBB &x1, const RealBB &y1, const RealBB &z1, 
                     const RealBB &x2, const RealBB &y2, const RealBB &z2);
  //
  //  Create a box intialized to an input point;
  //
  ContactBoundingBox(const Real *const point);
  //
  //  Create a box that is the interesection of 2 other bounding boxes.
  //
  ContactBoundingBox(const ContactBoundingBox &box1, const ContactBoundingBox &box2);
  ContactBoundingBox(const ContactBoundingBox *box1, const ContactBoundingBox *box2);
  //
  //  Reset the min/max.
  //
  void Reset();
  //
  //  Expand the bounding box to encompase the specified point.  Input is the points three coordinates.
  //
  void add_point(const Real *point);
  //
  //  Set the bounding box to equal the input point
  //
  void set_point(const Real *point);
  //
  //  Expand the bounding box to encompase a sphere, centered at the input point, with a given input radius
  //
  void add_sphere(const Real *center, const Real &radius);
  //
  //  Expand the bounding box to encompase another bounding box
  //
  void add_box(const ContactBoundingBox &box);
  void add_box(const ContactBoundingBox* box);
  //
  //  Expand the bounding box by ammount toler in all directions
  //
  void add_tolerance(const Real &toler);
  //
  //  Expand the bounding box by a different tolerance in each of the three directions
  //
  void add_tolerance(const Real *toler);
  //
  //  Expand the bounding box by ammount toler in all directions
  //
  void scale(const Real &toler);
  //
  //  Move the face box around keeping it's size constant, input is a vector to translate by.
  //
  ContactBoundingBox translate_negative(const Real* vector) const;
  //
  //  Determine the centroid of the box.  Used for box ordering routines
  //
  void calculate_centroid(Real *centroid) const;
  void calculate_centroid(float *centroid) const;
  void calculate_centroid2(Real *centroid) const;
  void calculate_centroid2(float *centroid) const;
  //
  //  Determine the largest axis aligned dimension of the bounding box
  //
  int longest_dimension();
  inline Real max_dimension() const;
  //
  //  Get the dimensions of the bounding box
  //
  inline void get_dimensions(Real* dimensions);
  //
  //  Determine if this box overlaps another input box or point
  //
  bool overlap(const ContactBoundingBox &box) const;
  bool overlap(const Real *point) const;
  bool overlap(const Real &x, const Real &y, const Real&z) const;
  //
  //  Explictly extract the corners of the bounding box.  Used for compatibility with some other
  //  ACME routines.
  //
  inline RealBB get_x_min() const {return x_min;};
  inline RealBB get_y_min() const {return y_min;};
  inline RealBB get_z_min() const {return z_min;};
  inline RealBB get_x_max() const {return x_max;};
  inline RealBB get_y_max() const {return y_max;};
  inline RealBB get_z_max() const {return z_max;};
  void get_min_max_points(Real *min, Real *max) const;
  //
  //  Explicitly set the corners of the bounding box.  Used for compatibility with some other 
  //  ACME routines
  //
  inline void set_x_min(const RealBB &value) {x_min = value;};
  inline void set_y_min(const RealBB &value) {y_min = value;};
  inline void set_z_min(const RealBB &value) {z_min = value;};
  inline void set_x_max(const RealBB &value) {x_max = value;};
  inline void set_y_max(const RealBB &value) {y_max = value;};
  inline void set_z_max(const RealBB &value) {z_max = value;};
  inline void set_box(Real *min, Real *max) {
    x_min = min[0];
    y_min = min[1];
    z_min = min[2];
    x_max = max[0];
    y_max = max[1];
    z_max = max[2];
  }
  inline void set_box(const ContactBoundingBox &box) {
    x_min = box.x_min;
    y_min = box.y_min;
    z_min = box.z_min;
    x_max = box.x_max;
    y_max = box.y_max;
    z_max = box.z_max;
  }
  inline void set_box(const ContactBoundingBox* box) {
    x_min = box->x_min;
    y_min = box->y_min;
    z_min = box->z_min;
    x_max = box->x_max;
    y_max = box->y_max;
    z_max = box->z_max;
  }
  
  inline bool is_valid() {
    if (x_min ==  BIGNUM_BB) return false;
    if (y_min ==  BIGNUM_BB) return false;
    if (z_min ==  BIGNUM_BB) return false;
    if (x_max == -BIGNUM_BB) return false;
    if (y_max == -BIGNUM_BB) return false;
    if (z_max == -BIGNUM_BB) return false;
    if (x_min > x_max) return false;
    if (y_min > y_max) return false;
    if (z_min > z_max) return false;
    return true;
  }
 
 friend ContactParOStream& operator<<( ContactParOStream& pos,
                                        const ContactBoundingBox &box );

  inline Real volume() {
    return std::max((x_max - x_min) * (y_max - y_min) * (z_max - z_min), 0.0f);
  }
  
 private:
  //
  //  The minimum and maximum corners of the bounding box.
  //
  RealBB x_min, y_min, z_min, x_max, y_max, z_max;
};

inline ContactBoundingBox::ContactBoundingBox() :
     x_min(BIGNUM_BB),
     y_min(BIGNUM_BB),
     z_min(BIGNUM_BB),
     x_max(-BIGNUM_BB),
     y_max(-BIGNUM_BB),
     z_max(-BIGNUM_BB)
{}

inline ContactBoundingBox::ContactBoundingBox(const Real *const point) :
     x_min(point[0]),
     y_min(point[1]),
     z_min(point[2]),
     x_max(point[0]),
     y_max(point[1]),
     z_max(point[2])
{}

inline ContactBoundingBox::ContactBoundingBox(const RealBB &x_min_, const RealBB &y_min_, const RealBB &z_min_, 
                                              const RealBB &x_max_, const RealBB &y_max_, const RealBB &z_max_) :
  x_min(x_min_),
  y_min(y_min_),
  z_min(z_min_),
  x_max(x_max_),
  y_max(y_max_),
  z_max(z_max_) 
{}

inline ContactBoundingBox::ContactBoundingBox(const ContactBoundingBox &box1, 
                                              const ContactBoundingBox &box2) :
     x_min(BIGNUM_BB),
     y_min(BIGNUM_BB),
     z_min(BIGNUM_BB),
     x_max(-BIGNUM_BB),
     y_max(-BIGNUM_BB),
     z_max(-BIGNUM_BB)
{
  if(box2.x_max < box1.x_min ||
     box2.y_max < box1.y_min ||
     box2.z_max < box1.z_min ||
     box2.x_min > box1.x_max ||
     box2.y_min > box1.y_max ||
     box2.z_min > box1.z_max) return;
  x_min = std::max(box2.x_min, box1.x_min);
  y_min = std::max(box2.y_min, box1.y_min);
  z_min = std::max(box2.z_min, box1.z_min);
  x_max = std::min(box2.x_max, box1.x_max);
  y_max = std::min(box2.y_max, box1.y_max);
  z_max = std::min(box2.z_max, box1.z_max);
}

inline ContactBoundingBox::ContactBoundingBox(const ContactBoundingBox *box1, 
                                              const ContactBoundingBox *box2) :
     x_min(BIGNUM_BB),
     y_min(BIGNUM_BB),
     z_min(BIGNUM_BB),
     x_max(-BIGNUM_BB),
     y_max(-BIGNUM_BB),
     z_max(-BIGNUM_BB)
{
  if(box2->x_max < box1->x_min ||
     box2->y_max < box1->y_min ||
     box2->z_max < box1->z_min ||
     box2->x_min > box1->x_max ||
     box2->y_min > box1->y_max ||
     box2->z_min > box1->z_max) return;
  x_min = std::max(box2->x_min, box1->x_min);
  y_min = std::max(box2->y_min, box1->y_min);
  z_min = std::max(box2->z_min, box1->z_min);
  x_max = std::min(box2->x_max, box1->x_max);
  y_max = std::min(box2->y_max, box1->y_max);
  z_max = std::min(box2->z_max, box1->z_max);
}

inline void ContactBoundingBox::Reset() {
  x_min =  BIGNUM_BB;
  y_min =  BIGNUM_BB;
  z_min =  BIGNUM_BB;
  x_max = -BIGNUM_BB;
  y_max = -BIGNUM_BB;
  z_max = -BIGNUM_BB;
}
                                              
inline void ContactBoundingBox::add_point(const Real *point) {
  x_min = std::min(x_min, (RealBB)point[0]);
  y_min = std::min(y_min, (RealBB)point[1]);
  z_min = std::min(z_min, (RealBB)point[2]);
  x_max = std::max(x_max, (RealBB)point[0]);
  y_max = std::max(y_max, (RealBB)point[1]);
  z_max = std::max(z_max, (RealBB)point[2]);
}

inline void ContactBoundingBox::set_point(const Real *point) {
  x_min = point[0];
  y_min = point[1];
  z_min = point[2];
  x_max = point[0];
  y_max = point[1];
  z_max = point[2];
}

inline void ContactBoundingBox::add_sphere(const Real *center, const Real &radius) {
  x_min = std::min(x_min, (RealBB)(center[0] - radius));
  y_min = std::min(y_min, (RealBB)(center[1] - radius));
  z_min = std::min(z_min, (RealBB)(center[2] - radius));
  x_max = std::max(x_max, (RealBB)(center[0] + radius));
  y_max = std::max(y_max, (RealBB)(center[1] + radius));
  z_max = std::max(z_max, (RealBB)(center[2] + radius));
}

inline void ContactBoundingBox::add_box(const ContactBoundingBox &box) {
  x_min = std::min(x_min, box.x_min);
  y_min = std::min(y_min, box.y_min);
  z_min = std::min(z_min, box.z_min);
  x_max = std::max(x_max, box.x_max);
  y_max = std::max(y_max, box.y_max);
  z_max = std::max(z_max, box.z_max);
}

inline void ContactBoundingBox::add_box(const ContactBoundingBox* box) {
  x_min = std::min(x_min, box->x_min);
  y_min = std::min(y_min, box->y_min);
  z_min = std::min(z_min, box->z_min);
  x_max = std::max(x_max, box->x_max);
  y_max = std::max(y_max, box->y_max);
  z_max = std::max(z_max, box->z_max);
}

inline void ContactBoundingBox::add_tolerance(const Real &toler) {
  Real abs_tol = std::fabs(toler);
  x_min -= abs_tol;
  y_min -= abs_tol;
  z_min -= abs_tol;
  x_max += abs_tol;
  y_max += abs_tol;
  z_max += abs_tol;
}

inline void ContactBoundingBox::add_tolerance(const Real *toler) {
  Real abs_tol_0 = std::fabs(toler[0]);
  Real abs_tol_1 = std::fabs(toler[1]);
  Real abs_tol_2 = std::fabs(toler[2]);

  x_min -= abs_tol_0;
  y_min -= abs_tol_1;
  z_min -= abs_tol_2;
  x_max += abs_tol_0;
  y_max += abs_tol_1;
  z_max += abs_tol_2;
}

inline void ContactBoundingBox::scale(const Real &toler) {
  Real xdelta = x_max - x_min;
  Real ydelta = y_max - y_min;
  Real zdelta = z_max - z_min;
  x_min -= toler*xdelta;
  y_min -= toler*ydelta;
  z_min -= toler*zdelta;
  x_max += toler*xdelta;
  y_max += toler*ydelta;
  z_max += toler*zdelta;
  if        (x_min==x_max) {
    x_min -= 0.0001*std::max(ydelta, zdelta);
    x_max += 0.0001*std::max(ydelta, zdelta);
  } else if (y_min==y_max) {
    y_min -= 0.0001*std::max(xdelta, zdelta);
    y_max += 0.0001*std::max(xdelta, zdelta);
  } else if (z_min==z_max) {
    z_min -= 0.0001*std::max(xdelta, ydelta);
    z_max += 0.0001*std::max(xdelta, ydelta);
  }
}

inline bool ContactBoundingBox::overlap(const ContactBoundingBox &box) const {
  if(x_max < box.x_min ||
     y_max < box.y_min ||
     z_max < box.z_min ||
     x_min > box.x_max ||
     y_min > box.y_max ||
     z_min > box.z_max) return false;
  return true;
}

inline bool ContactBoundingBox::overlap(const Real* point) const {
  if(point[0] < x_min ||
     point[1] < y_min ||
     point[2] < z_min ||
     point[0] > x_max ||
     point[1] > y_max ||
     point[2] > z_max) return false;
  return true;
}

inline bool ContactBoundingBox::overlap(const Real &x, 
                                        const Real &y,
                                        const Real &z) const {
  if(x < x_min ||
     y < y_min ||
     z < z_min ||
     x > x_max ||
     y > y_max ||
     z > z_max) return false;
  return true;
}

inline ContactBoundingBox ContactBoundingBox::translate_negative(const Real* vector) const{
  ContactBoundingBox return_box(*this);
  return_box.x_min -= vector[0];
  return_box.y_min -= vector[1];
  return_box.z_min -= vector[2];
  return_box.x_max -= vector[0];
  return_box.y_max -= vector[1];
  return_box.z_max -= vector[2];
  return return_box;
}

inline void ContactBoundingBox::get_min_max_points(Real *min, Real *max) const {
  min[0] = x_min;
  min[1] = y_min;
  min[2] = z_min;
  max[0] = x_max;
  max[1] = y_max;
  max[2] = z_max;
}

inline void ContactBoundingBox::calculate_centroid(Real *centroid) const {
  centroid[0] = (x_min + x_max) * 0.5;
  centroid[1] = (y_min + y_max) * 0.5;
  centroid[2] = (z_min + z_max) * 0.5;
}

inline void ContactBoundingBox::calculate_centroid2(Real *centroid) const {
  centroid[0] = (x_min + x_max);
  centroid[1] = (y_min + y_max);
  centroid[2] = (z_min + z_max);
}

inline void ContactBoundingBox::calculate_centroid(float *centroid) const {
  centroid[0] = (x_min + x_max) * 0.5;
  centroid[1] = (y_min + y_max) * 0.5;
  centroid[2] = (z_min + z_max) * 0.5;
}

inline void ContactBoundingBox::calculate_centroid2(float *centroid) const {
  centroid[0] = (x_min + x_max);
  centroid[1] = (y_min + y_max);
  centroid[2] = (z_min + z_max);
}

//
// Get maximum dimension of box
//
inline Real ContactBoundingBox::max_dimension() const{
  Real length_0 = x_max-x_min;
  Real length_1 = y_max-y_min;
  Real length_2 = z_max-z_min;
  if(length_0 > length_1 && length_0 > length_2) {
    return length_0;
  } else if(length_1 > length_2) {
    return length_1;
  } else {
    return length_2;
  }
}

//
// Get dimensions of box
//
inline void ContactBoundingBox::get_dimensions(Real* dimensions) {
  dimensions[0] = x_max-x_min;
  dimensions[1] = y_max-y_min;
  dimensions[2] = z_max-z_min;
}

//
//  Add a box to the current box
//
inline ContactBoundingBox operator+(const ContactBoundingBox &box, 
                                    const Real &toler) {
  ContactBoundingBox return_box(box);
  return_box.add_tolerance(toler);
  return return_box;
}



//
//  Add two boxes to produce a box that encompases both
//
inline ContactBoundingBox operator+(const ContactBoundingBox &box1, 
                                    const ContactBoundingBox &box2) {
  ContactBoundingBox return_box(box1);
  return_box.add_box(box2);
  return return_box;
}

inline std::ostream& operator<<(std::ostream &output, const ContactBoundingBox &box) {
  output<<"Min corner "<<box.get_x_min()<<" "<<box.get_y_min()<<" "<<box.get_z_min()<<std::endl;
  output<<"Max corner "<<box.get_x_max()<<" "<<box.get_y_max()<<" "<<box.get_z_max()<<std::endl;
  return output;
}


#define MIN_CONTACT_BOX_INT 0
#define MAX_CONTACT_BOX_INT 65535

class DomainBox : public ContactBoundingBox {
 public:
  DomainBox() {};
  void calculate_constants() {
    x_len  = this->get_x_max() - this->get_x_min();
    y_len  = this->get_y_max() - this->get_y_min();
    z_len  = this->get_z_max() - this->get_z_min();
    x_mult = MAX_CONTACT_BOX_INT / x_len;
    y_mult = MAX_CONTACT_BOX_INT / y_len;
    z_mult = MAX_CONTACT_BOX_INT / z_len;
  };
  inline RealBB get_x_len() const {return x_len;};
  inline RealBB get_y_len() const {return y_len;};
  inline RealBB get_z_len() const {return z_len;};
  inline RealBB get_x_mult() const {return x_mult;};
  inline RealBB get_y_mult() const {return y_mult;};
  inline RealBB get_z_mult() const {return z_mult;};

 private:
  RealBB x_len, y_len, z_len;
  RealBB x_mult, y_mult, z_mult;  
};

class ContactBoundingBox_Int {

 public:
  //
  //  Default box constructor.  This constructor will creat a box that encompasses the entire computational domain
  //
  ContactBoundingBox_Int();

  inline void set_box(const ContactBoundingBox &box,
                      const DomainBox &domain_box);
  inline bool overlap(const ContactBoundingBox_Int &box) const;
  inline unsigned short get_x_min() const {return x_min;};
  inline unsigned short get_y_min() const {return y_min;};
  inline unsigned short get_z_min() const {return z_min;};
  inline unsigned short get_x_max() const {return x_max;};
  inline unsigned short get_y_max() const {return y_max;};
  inline unsigned short get_z_max() const {return z_max;};

  inline void set_box(const ContactBoundingBox_Int &box) {
    x_min = box.x_min;
    y_min = box.y_min;
    z_min = box.z_min;
    x_max = box.x_max;
    y_max = box.y_max;
    z_max = box.z_max;
  }

  inline void add_box(const ContactBoundingBox_Int &box) {
    x_min = std::min(x_min, box.x_min);
    y_min = std::min(y_min, box.y_min);
    z_min = std::min(z_min, box.z_min);
    x_max = std::max(x_max, box.x_max);
    y_max = std::max(y_max, box.y_max);
    z_max = std::max(z_max, box.z_max);
  }

  inline void calculate_centroid(float *centroid) const {
    centroid[0] = 0.5 * (x_min + x_max);
    centroid[1] = 0.5 * (y_min + y_max);
    centroid[2] = 0.5 * (z_min + z_max);
  }

  inline void calculate_centroid2(float *centroid) const {
    centroid[0] = (x_min + x_max);
    centroid[1] = (y_min + y_max);
    centroid[2] = (z_min + z_max);
  }

 private:
  //
  //  The minimum and maximum corners of the bounding box.
  //
  unsigned short x_min, y_min, z_min, x_max, y_max, z_max;
};



inline ContactBoundingBox_Int::ContactBoundingBox_Int() :
     x_min(MAX_CONTACT_BOX_INT),
     y_min(MAX_CONTACT_BOX_INT),
     z_min(MAX_CONTACT_BOX_INT),
     x_max(MIN_CONTACT_BOX_INT),
     y_max(MIN_CONTACT_BOX_INT),
     z_max(MIN_CONTACT_BOX_INT)
{}

inline void ContactBoundingBox_Int::set_box(const ContactBoundingBox &box, 
                                            const DomainBox &domain_box) {
  //
  //  Convert the bounding box real limits into conservative integer bounds
  //
  RealBB real_pos_x_min = box.get_x_min() - domain_box.get_x_min();
  RealBB real_pos_y_min = box.get_y_min() - domain_box.get_y_min();
  RealBB real_pos_z_min = box.get_z_min() - domain_box.get_z_min();
  RealBB real_pos_x_max = box.get_x_max() - domain_box.get_x_min();
  RealBB real_pos_y_max = box.get_y_max() - domain_box.get_y_min();
  RealBB real_pos_z_max = box.get_z_max() - domain_box.get_z_min();

  RealBB x_len = domain_box.get_x_len();
  RealBB y_len = domain_box.get_y_len();
  RealBB z_len = domain_box.get_z_len();

  if(real_pos_x_min <= 0.0) {
    x_min = MIN_CONTACT_BOX_INT;
  } else if(real_pos_x_min >= x_len) {
    x_min = MAX_CONTACT_BOX_INT;
  } else {
    x_min = (unsigned short int)std::floor(real_pos_x_min * domain_box.get_x_mult());
  }
  if(real_pos_x_max <= 0.0) {
    x_max = MIN_CONTACT_BOX_INT;
  } else if(real_pos_x_max >= x_len) {
    x_max = MAX_CONTACT_BOX_INT;
  } else {
    x_max = (unsigned short int)std::floor(real_pos_x_max * domain_box.get_x_mult());
  }

  if(real_pos_y_min <= 0.0) {
    y_min = MIN_CONTACT_BOX_INT;
  } else if(real_pos_y_min >= y_len) {
    y_min = MAX_CONTACT_BOX_INT;
  } else {
    y_min = (unsigned short int)std::floor(real_pos_y_min * domain_box.get_y_mult());
  }
  if(real_pos_y_max <= 0.0) {
    y_max = MIN_CONTACT_BOX_INT;
  } else if(real_pos_y_max >= y_len) {
    y_max = MAX_CONTACT_BOX_INT;
  } else {
    y_max = (unsigned short int)std::floor(real_pos_y_max * domain_box.get_y_mult());
  }

  if(real_pos_z_min <= 0.0) {
    z_min = MIN_CONTACT_BOX_INT;
  } else if(real_pos_z_min >= z_len) {
    z_min = MAX_CONTACT_BOX_INT;
  } else {
    z_min = (unsigned short int)std::floor(real_pos_z_min * domain_box.get_z_mult());
  }
  if(real_pos_z_max <= 0.0) {
    z_max = MIN_CONTACT_BOX_INT;
  } else if(real_pos_z_max >= z_len) {
    z_max = MAX_CONTACT_BOX_INT;
  } else {
    z_max = (unsigned short int)std::floor(real_pos_z_max * domain_box.get_z_mult());
  }
}


inline bool ContactBoundingBox_Int::overlap(const ContactBoundingBox_Int &box) const {
  if(x_max < box.x_min ||
     y_max < box.y_min ||
     z_max < box.z_min ||
     x_min > box.x_max ||
     y_min > box.y_max ||
     z_min > box.z_max) return false;
  return true;
}

inline std::ostream& operator<<(std::ostream &output, const ContactBoundingBox_Int &box) {
  output<<"Min corner "<<box.get_x_min()<<" "<<box.get_y_min()<<" "<<box.get_z_min()<<std::endl;
  output<<"Max corner "<<box.get_x_max()<<" "<<box.get_y_max()<<" "<<box.get_z_max()<<std::endl;
  return output;
}

//
//  Bounding box of an object plus an index to that object.  This object is used to create bounding box arrays for
//  things like processors, faces, or nodes.
//
class ObjectBoundingBox : public ContactBoundingBox {
 public:
  ObjectBoundingBox(const ContactBoundingBox &box, const int obj_num_) :
    ContactBoundingBox(box),
    obj_num(obj_num_)
    {}


  //
  //  It is assumed that the object number refers to an index in a C array.  Thus a index of -1 is invalid and
  //  denotes an undefined object index
  //
  ObjectBoundingBox() {
    obj_num = -1;
  }

  //
  //  Explicity set or extract the object index for this bounding box
  //
  void set_object_number(const int &obj_num_);
  inline int get_object_number() const{return obj_num;};

  //
  //  Global box combine takes a set of boxes defined on mulitple processors.  It then uses global reductions to do
  //  box additions in such a way that after the operation completes each processor has an identical set of boxes.  
  //  Each box in the set encompases all boxes defined for each processor.  The primary use for this routine is to
  //  find bounding boxes for each processor of a parallel decomposition.
  //
  static void global_box_combine(ObjectBoundingBox *box_array, 
                                 const int &num_boxes, 
                                 MPI_Comm &SearchComm);
 private:
  int obj_num;
};

inline void ObjectBoundingBox::set_object_number(const int &obj_num_) {
  obj_num = obj_num_;
}



//
//  Bounding box of an object plus an index to that object.  This object is used to create bounding box arrays for
//  things like processors, faces, or nodes.
//
class ObjectBoundingBox_Int : public ContactBoundingBox_Int {
 public:
  //
  //  It is assumed that the object number refers to an index in a C array.  Thus a index of -1 is invalid and
  //  denotes an undefined object index
  //
  ObjectBoundingBox_Int() {
    obj_num = -1;
  }

  //
  //  Explicity set or extract the object index for this bounding box
  //
  void set_object_number(const int &obj_num_);
  inline int get_object_number() const{return obj_num;};

 private:
  int obj_num;
};

inline void ObjectBoundingBox_Int::set_object_number(const int &obj_num_) {
  obj_num = obj_num_;
}


inline std::ostream& operator<<(std::ostream &output, const ObjectBoundingBox &box) {
  output<<"Min corner "<<box.get_x_min()<<" "<<box.get_y_min()<<" "<<box.get_z_min()<<std::endl;
  output<<"Max corner "<<box.get_x_max()<<" "<<box.get_y_max()<<" "<<box.get_z_max()<<std::endl;
  output<<"object number "<<box.get_object_number()<<std::endl;
  return output;
}


#endif  // #ifdef ContactBoundingBox_h_
