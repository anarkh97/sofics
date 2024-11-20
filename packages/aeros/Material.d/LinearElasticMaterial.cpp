/*
 * LinearElasticMaterial.cpp
 * DG++
 *
 * Created by Adrian Lew on 11/19/06.
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
#include <cmath>
#include <iostream>





bool LinearElasticBase::GetConstitutiveResponse(const std::vector<double> * lstrain,
						std::vector<double> * lstress,
						std::vector<double> * ltangents) const
{
  const std::vector<double> &strain = *lstrain;
  std::vector<double> &stress = *lstress;
  
  const double Identity[] = {1,0,0,0,1,0,0,0,1};

  // Compute stress
  stress.resize(9);

  for(int i=0; i<3; i++)
    for(int J=0; J<3; J++)
      {
	stress[i*3+J] = 0;

	for(int k=0; k<3; k++)
	  for(int L=0; L<3; L++)
	    stress[i*3+J] += GetModuli(i,J,k,L)*(strain[k*3+L]-Identity[k*3+L]);
      }
  

  // Compute tangents, if needed
  if(ltangents)
    {
      std::vector<double> &tangents = *ltangents;
      
      tangents.resize(81);
      
      for(int i=0; i<3; i++)
	for(int J=0; J<3; J++)
	  for(int k=0; k<3; k++)
	    for(int L=0; L<3; L++)
	      tangents[i*27+J*9+k*3+L] = GetModuli(i,J,k,L);
    }

  return true;
}
