#ifndef _PRINT_DEBUG_H_
#define _PRINT_DEBUG_H_

#include <iostream>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>

template<class Scalar>
void print_debug(GenDistrVector<Scalar> *vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec->size(); ++i) std::cerr << i << "," << vec->data()[i] << " "; std::cerr << std::endl;
 }
 else {
   for(int i=0; i<vec->size(); ++i) std::cerr << vec->data()[i] << " "; std::cerr << std::endl;
 }
}

template<class Scalar>
void print_debug(GenDistrVector<Scalar> &vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec.size(); ++i) std::cerr << i << "," << vec.data()[i] << " "; std::cerr << std::endl;
 }
 else {
   for(int i=0; i<vec.size(); ++i) std::cerr << vec.data()[i] << " "; std::cerr << std::endl;
 }
}

template<class Scalar>
void print_debug(GenVector<Scalar> *vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec->size(); ++i) std::cerr << i << "," << vec->data()[i] << " "; std::cerr << std::endl; 
 }
 else {
   for(int i=0; i<vec->size(); ++i) std::cerr << vec->data()[i] << " "; std::cerr << std::endl;
 }
}

template<class Scalar>
void print_debug(GenVector<Scalar> &vec, bool print_index = false)
{
 if(print_index) {
   for(int i=0; i<vec.size(); ++i) std::cerr << i << "," << vec.data()[i] << " "; std::cerr << std::endl;
 }
 else {
   for(int i=0; i<vec.size(); ++i) std::cerr << vec.data()[i] << " "; std::cerr << std::endl;
 }
}

template<class Scalar>
void print_debug(Scalar *vec, int len)
{
 for(int i=0; i<len; ++i) std::cerr << vec[i] << " "; std::cerr << std::endl;
}

#endif
