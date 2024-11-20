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


#ifndef ContactVector_h_
#define ContactVector_h_

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <Contact_Defines.h>
#include <ContactBoundingBox.h>

template<typename DataType> class ContactNode;

namespace ACME {
  //
  //  For all practical purposes, ACME_Int_Vector is a equivalent to the simple STL vector std::vector<int>
  //  However ACME currently forbids templates so, here is a single limited instantiation of that class 
  //

  class Int_Vector {
  public:
    //
    //  Create a default empty vector
    //
    Int_Vector() :
      length(0),
      capacity(0),
      data(0){
    }
    //
    //  Create a vector of length length_
    //
    Int_Vector(const int length_) :
      length(length_),
      capacity(length_)
      {
	data = new int[length];
      }

    //
    //  Create a vector of length length_ with all elements having data init_data
    //
    Int_Vector(const int length_, const int init_data) :
      length(length_),
      capacity(length_),
      data(new int[length_]){
      for(int i = 0; i < length; ++i) {
        data[i] = init_data;
      }
    }
    ~Int_Vector() {
      if(data != 0) delete [] data;
    }
    //
    //  Subset of used STL vector operators...
    //
    inline int &operator[](const int index) const{
      return data[index];
    }
    //
    //  Add a single entry to the back of the list.  Use the standard 50% size increase, starting at size 16; 
    //   
    inline void push_back(const int new_entry) {
      if(length + 1> capacity) {
        if(capacity == 0) {
          capacity = 16;
	} else {
          capacity *= 3;
	  capacity /= 2;
	}
        int *new_data = new int[capacity];
	std::memcpy(new_data, data, length*sizeof(int));
        if(data != 0) delete [] data;
        data = new_data;
      }
      data[length++] = new_entry;
    }
    //
    //  Delete all data and reset length and capacity
    //
    inline void clear() {
      length = 0;
      capacity = 0;
      if(data!= 0) delete [] data;
      data = 0;
    }
    //
    //  Set the exact size for a vector
    //
    inline void resize(const int new_size) {
      clear();
      length = new_size;
      capacity = new_size;
      data = new int[capacity];
    }
    inline void resize(const int new_size, const int initial_val) {
      clear();
      length = new_size;
      capacity = new_size;
      data = new int[capacity];
      for(int i = 0; i < capacity; ++i) {
        data[i] = initial_val;
      }
    }
    //
    //  Return the current length of the list
    //
    inline int size() const{
      return length;
    }
    //
    //  Get a direct pointer to the data buffer
    //
    inline int* get_buffer() {
      return data;
    }
    //
    //  Make sure that the vector is at least size length_.
    //  If the vector is longer than length_, leave it unchanged
    //
    inline void reserve(const int required_capacity) {
      if(capacity < required_capacity) resize(required_capacity);
      empty();
    }
    //
    //  Reset length, but leave the vector capacity untouched
    //
    inline void empty() {
      length = 0;
    }
    //
    //  Does the vector contain a entry?
    //
    inline bool contain(const int test) {
      for(int i = 0; i < length; ++i) {
        if(test == data[i]) return true;
      }
      return false;
    }


  private:
    int length;
    int capacity;
    int *data;
  };


  //
  //  For all practical purposes, ACME_Real_Vector is a equivalent to the simple STL vector std::vector<Real>
  //  However ACME current forbids templates so, here is a single limited instantiation of that class 
  //

  class Real_Vector {
  public:
    Real_Vector() :
      length(0),
      capacity(0),
      data(0){
    }
    ~Real_Vector() {
      if(data != 0) delete [] data;
    }
    //
    //  Subset of used STL vector operators...
    //
    inline Real &operator[](const int index) const{
      return data[index];
    }
    //
    //  Add a single entry to the back of the list.  Use the standard 50% size increase, starting at size 16; 
    //   
    inline void push_back(const Real new_entry) {
      if(length + 1> capacity) {
        if(capacity == 0) {
          capacity = 16;
	} else {
          capacity *= 3;
	  capacity /= 2;
	}
        Real *new_data = new Real[capacity];
	std::memcpy(new_data, data, length*sizeof(Real));
        if(data != 0) delete [] data;
        data = new_data;
      }
      data[length++] = new_entry;
    }
    //
    //  Delete all data and reset length and capacity
    //
    inline void clear() {
      length = 0;
      capacity = 0;
      if(data!= 0) delete [] data;
      data = 0;
    }
    //
    //  Reset length, but leave the vector capacity untouched
    //
    inline void empty() {
      length = 0;
    }
    //
    //  Set the exact size for a vector
    //
    inline void resize(const int new_size) {
      clear();
      length = new_size;
      capacity = new_size;
      data = new Real[capacity];
    }
    //
    //  Return the current length of the list
    //
    inline int size() const{
      return length;
    }
    //
    //  Get a direct pointer to the data buffer
    //
    inline Real* get_buffer() {
      return data;
    }

  private:
    int length;
    int capacity;
    Real *data;
  };


  //
  //  For all practical purposes, ACME::ContactNode_Vector is a equivalent to the simple STL vector std::vector<ContactNode<Real>*>
  //  However ACME current forbids templates so, here is a single limited instantiation of that class 
  //

  class ContactNode_Vector {
  public:
    //
    //  Create a default empty vector
    //
    ContactNode_Vector() :
      length(0),
      capacity(0),
      data(0){
    }
    //
    //  Create a vector of length length_
    //
    ContactNode_Vector(const int length_) :
      length(length_),
      capacity(length_)
      {
	data = new ContactNode<Real>*[length];
      }

    //
    //  Create a vector of length length_ with all elements having data init_data
    //
    ContactNode_Vector(const int length_, ContactNode<Real> *const init_data) :
      length(length_),
      capacity(length_),
      data(new ContactNode<Real>*[length_]){
      for(int i = 0; i < length; ++i) {
        data[i] = init_data;
      }
    }
    ~ContactNode_Vector() {
      if(data != 0) delete [] data;
    }
    //
    //  Subset of used STL vector operators...
    //
    inline ContactNode<Real>* &operator[](const int index) const{
      return data[index];
    }
    //
    //  Add a single entry to the back of the list.  Use the standard 50% size increase, starting at size 16; 
    //   
    inline void push_back(ContactNode<Real> *const new_entry) {
      if(length + 1> capacity) {
        if(capacity == 0) {
          capacity = 16;
	} else {
          capacity *= 3;
	  capacity /= 2;
	}
        ContactNode<Real> **new_data = new ContactNode<Real>*[capacity];
	std::memcpy(new_data, data, length*sizeof(ContactNode<Real>*));
        if(data != 0) delete [] data;
        data = new_data;
      }
      data[length++] = new_entry;
    }
    //
    //  Delete all data and reset length and capacity
    //
    inline void clear() {
      length = 0;
      capacity = 0;
      if(data!= 0) delete [] data;
      data = 0;
    }
    //
    //  Reset length, but leave the vector capacity untouched
    //
    inline void empty() {
      length = 0;
    }
    //
    //  Set the exact size for a vector
    //
    inline void resize(const int new_size) {
      clear();
      length = new_size;
      capacity = new_size;
      data = new ContactNode<Real>*[capacity];
    }
    //
    //  Return the current length of the list
    //
    inline int size() const{
      return length;
    }
    //
    //  Get a direct pointer to the data buffer
    //
    inline ContactNode<Real>** get_buffer() {
      return data;
    }
    //
    //  Make sure that the vector is at least size length_.
    //  If the vector is longer than length_, leave it unchanged
    //
    inline void reserve(const int required_capacity) {
      if(capacity < required_capacity) resize(required_capacity);
      empty();
    }
    //
    //  Does the vector contain a entry?
    //
    inline bool contain(ContactNode<Real> *const test) {
      for(int i = 0; i < length; ++i) {
        if(test == data[i]) return true;
      }
      return false;
    }


  private:
    int length;
    int capacity;
    ContactNode<Real> **data;
  };


  //
  //  For all practical purposes, ACME_Int_Vector is a equivalent to the simple STL vector std::vector<ContactBoundingBox>
  //  However ACME current forbids templates so, here is a single limited instantiation of that class 
  //

  class ContactBoundingBox_Vector {
  public:
    //
    //  Create a default empty vector
    //
    ContactBoundingBox_Vector() :
      length(0),
      capacity(0),
      data(0){
    }
    //
    //  Create a vector of length length_
    //
    ContactBoundingBox_Vector(const int length_) :
      length(length_),
      capacity(length_)
      {
	data = new ContactBoundingBox[length];
      }

    //
    //  Create a vector of length length_ with all elements having data init_data
    //
    ContactBoundingBox_Vector(const int length_, const ContactBoundingBox init_data) :
      length(length_),
      capacity(length_),
      data(new ContactBoundingBox[length_]){
      for(int i = 0; i < length; ++i) {
        data[i] = init_data;
      }
    }
    ~ContactBoundingBox_Vector() {
      if(data != 0) delete [] data;
    }
    //
    //  Subset of used STL vector operators...
    //
    inline ContactBoundingBox &operator[](const int index) const{
      return data[index];
    }
    //
    //  Add a single entry to the back of the list.  Use the standard 50% size increase, starting at size 16; 
    //   
    inline void push_back(const ContactBoundingBox new_entry) {
      if(length + 1> capacity) {
        if(capacity == 0) {
          capacity = 16;
	} else {
          capacity *= 3;
	  capacity /= 2;
	}
        ContactBoundingBox *new_data = new ContactBoundingBox[capacity];
	std::memcpy(new_data, data, length*sizeof(ContactBoundingBox));
        if(data != 0) delete [] data;
        data = new_data;
      }
      data[length++] = new_entry;
    }
    //
    //  Delete all data and reset length and capacity
    //
    inline void clear() {
      length = 0;
      capacity = 0;
      if(data!= 0) delete [] data;
      data = 0;
    }
    //
    //  Reset length, but leave the vector capacity untouched
    //
    inline void empty() {
      length = 0;
    }
    //
    //  Set the exact size for a vector
    //
    inline void resize(const int new_size) {
      clear();
      length = new_size;
      capacity = new_size;
      data = new ContactBoundingBox[capacity];
    }
    //
    //  Return the current length of the list
    //
    inline int size() const{
      return length;
    }
    //
    //  Get a direct pointer to the data buffer
    //
    inline ContactBoundingBox* get_buffer() {
      return data;
    }
    //
    //  Make sure that the vector is at least size length_.
    //  If the vector is longer than length_, leave it unchanged
    //
    inline void reserve(const int required_capacity) {
      if(capacity < required_capacity) resize(required_capacity);
      empty();
    }
  private:
    int length;
    int capacity;
    ContactBoundingBox *data;
  };



  //
  //  For all practical purposes, ACME_Int_Vector is a equivalent to the simple STL vector std::vector<ObjectBoundingBox>
  //  However ACME current forbids templates so, here is a single limited instantiation of that class 
  //

  class ObjectBoundingBox_Vector {
  public:
    //
    //  Create a default empty vector
    //
    ObjectBoundingBox_Vector() :
      length(0),
      capacity(0),
      data(0){
    }
    //
    //  Create a vector of length length_
    //
    ObjectBoundingBox_Vector(const int length_) :
      length(length_),
      capacity(length_)
      {
	data = new ObjectBoundingBox[length];
      }

    //
    //  Create a vector of length length_ with all elements having data init_data
    //
    ObjectBoundingBox_Vector(const int length_, const ObjectBoundingBox init_data) :
      length(length_),
      capacity(length_),
      data(new ObjectBoundingBox[length_]){
      for(int i = 0; i < length; ++i) {
        data[i] = init_data;
      }
    }
    ~ObjectBoundingBox_Vector() {
      if(data != 0) delete [] data;
    }
    //
    //  Subset of used STL vector operators...
    //
    inline ObjectBoundingBox &operator[](const int index) const{
      return data[index];
    }
    //
    //  Add a single entry to the back of the list.  Use the standard 50% size increase, starting at size 16; 
    //   
    inline void push_back(const ObjectBoundingBox new_entry) {
      if(length + 1> capacity) {
        if(capacity == 0) {
          capacity = 16;
	} else {
          capacity *= 3;
	  capacity /= 2;
	}
        ObjectBoundingBox *new_data = new ObjectBoundingBox[capacity];
	std::memcpy(new_data, data, length*sizeof(ObjectBoundingBox));
        if(data != 0 ) delete [] data;
        data = new_data;
      }
      data[length++] = new_entry;
    }
    //
    //  Delete all data and reset length and capacity
    //
    inline void clear() {
      length = 0;
      capacity = 0;
      if(data!= 0) delete [] data;
      data = 0;
    }
    //
    //  Reset length, but leave the vector capacity untouched
    //
    inline void empty() {
      length = 0;
    }
    //
    //  Set the exact size for a vector
    //
    inline void resize(const int new_size) {
      clear();
      length = new_size;
      capacity = new_size;
      data = new ObjectBoundingBox[capacity];
    }
    //
    //  Return the current length of the list
    //
    inline int size() const{
      return length;
    }
    //
    //  Get a direct pointer to the data buffer
    //
    inline ObjectBoundingBox* get_buffer() {
      return data;
    }
    //
    //  Make sure that the vector is at least size length_.
    //  If the vector is longer than length_, leave it unchanged
    //
    inline void reserve(const int required_capacity) {
      if(capacity < required_capacity) resize(required_capacity);
      empty();
    }
  private:
    int length;
    int capacity;
    ObjectBoundingBox *data;
  };

} // end namespace ACME


#endif // #ifdef ContactVector_h_
