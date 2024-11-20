/*
 * NeoHookean.cpp
 * DG++
 *
 * Created by Adrian Lew on 10/24/06.
 *  
 * Copyright (c) 2006 Adrian Lew
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */ 

#include <Material.d/Material.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>

bool SimpleMaterial::ConsistencyTest(const SimpleMaterial &SMat) 
{
  std::vector<double> strain(9);
  std::vector<double> stress(9);
  std::vector<double> tangents(81);

  std::vector<double> stressplus(9);
  std::vector<double> stressminus(9);
  std::vector<double> tangentsnum(81);
  
  double EPS = 1.e-6;
  double PERT = 1.e-1;

  srand(time(0));
  
  strain[0] = 1 + double(rand())/double(RAND_MAX)*PERT;
  strain[1] =     double(rand())/double(RAND_MAX)*PERT;
  strain[2] =     double(rand())/double(RAND_MAX)*PERT;
  strain[3] =     double(rand())/double(RAND_MAX)*PERT;
  strain[4] = 1 + double(rand())/double(RAND_MAX)*PERT;
  strain[5] =     double(rand())/double(RAND_MAX)*PERT;
  strain[6] =     double(rand())/double(RAND_MAX)*PERT;
  strain[7] =     double(rand())/double(RAND_MAX)*PERT;
  strain[8] = 1 + double(rand())/double(RAND_MAX)*PERT;

  for(unsigned int i=0; i<9; i++)
    {
      double Forig = strain[i];
      
      strain[i] = Forig + EPS;
      SMat.GetConstitutiveResponse(&strain, &stressplus);

      strain[i] = Forig - EPS;
      SMat.GetConstitutiveResponse(&strain, &stressminus);

      for(unsigned j=0; j<9; j++)
	tangentsnum[j*9 + i] = (stressplus[j] - stressminus[j])/(2*EPS);

      strain[i] = Forig;
    }

  SMat.GetConstitutiveResponse(&strain, &stress, &tangents);
  
  double error = 0;
  double norm = 0;
  for(int i=0; i<81; i++)
    {
      error += pow(tangents[i]-tangentsnum[i],2);
      norm += pow(tangents[i],2);
    }
  error = sqrt(error);
  norm = sqrt(norm);

  if(error/norm>EPS*100)
    {
      std::cerr << "SimpleMaterial::ConsistencyTest. Material not consistent\n";
      return false;
    }
  return true;
}
