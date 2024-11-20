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


#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "contact_assert.h"
#include "Contact_Defines.h"
#include "CString.h"

CString::CString(std::size_t n, char c)
{
  my_letter = (my_Letter*)(new char[sizeof(my_Letter)+n]);
  (std::size_t&)my_letter->capacity = my_letter->size = n;
  my_letter->rc = 1;
  std::memset(my_letter->data, c, n);
  my_letter->data[n] = 0;
}

CString::CString(const char* str)
{
  std::size_t sz = std::strlen(str);
  my_letter = (my_Letter*)(new char[sizeof(my_Letter)+sz]);
  (std::size_t&)my_letter->capacity = my_letter->size = sz;
  my_letter->rc = 1;
  std::strcpy(my_letter->data, str);
}

CString::CString(const char* str, std::size_t n)
{
  my_letter = (my_Letter*)(new char[sizeof(my_Letter)+n]);
  (std::size_t&)my_letter->capacity = n;
  my_letter->rc = 1;
  if (n==0) my_letter->size = 0;
  else my_letter->size = std::min(n, std::strlen(str));
  std::memcpy(my_letter->data, str, my_letter->size);
  my_letter->data[n] = 0;
}

CString::CString(const CString& str, std::size_t pos, std::size_t n)
{
  if (!pos && str.my_letter->size <= n){
    my_letter = str.my_letter;
    my_letter->rc++;
  }
  else
  {
    std::size_t size = std::min(str.my_letter->size - pos, n);
    my_letter = (my_Letter*)(new char[sizeof(my_Letter)+size]);
    (std::size_t&)my_letter->capacity = my_letter->size = size;
    my_letter->rc = 1;
    std::memcpy(my_letter->data, str.my_letter->data+pos, size);
    my_letter->data[size] = 0;
  }
}

CString::~CString(void)
{
  free();
}

CString& CString::append(const CString& s)
{
  std::size_t new_size = my_letter->size + s.my_letter->size;
  my_Letter* new_letter = (my_Letter*)( new char[ sizeof(my_Letter) + new_size ] );
  new_letter->capacity = my_letter->capacity;
  new_letter->size     = new_size;
  new_letter->rc       = 1;
  std::memcpy(new_letter->data, my_letter->data, size());
  std::memcpy(new_letter->data + size(), s.my_letter->data, s.my_letter->size);
  delete my_letter;
  my_letter = new_letter;
  return *this;
}

int operator==(const CString &b, const char *a)
{
  return std::strlen(a)==b.size() && !std::memcmp(a, b.data(), b.size());
}

void CString::free(void)
{
  if (!--my_letter->rc) delete[] (char*)my_letter;
}

std::ostream& operator<<(std::ostream &o, const CString &str)
{
  o << str.data();
  return o;
}

