/*                                                                       
 * Copyright 10/12/2017 - Dr. Christopher H. S. Aylett                   
 *                                                                       
 * This program is free software; you can redistribute it and/or modify  
 * it under the terms of version 3 of the GNU General Public License as  
 * published by the Free Software Foundation.                            
 *                                                                       
 * This program is distributed in the hope that it will be useful,       
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
 * GNU General Public License for more details - YOU HAVE BEEN WARNED!   
 *                                                                       
 * Program: LAFTER V1.0                                                  
 *                                                                       
 * Authors: Chris Aylett                                                 
 *                                                                       
 */

// Library header inclusion for linking                                  
#include "lafter.h"

#define THREAD 0
#define STEP 1

// Updates out if in1/2 over noise - returns fractional recovery
double truncate_map(double *in1, double *in2, double *out, r_mrc *mask, list *node, double overfit_correction, int32_t size){
  int32_t i, max = size * size * size;
  double rcv = 0.0;
  double cur, cor, count = 0.0;
  double noise = 0.0;

  // Calculate max noise
  for (i = THREAD; i < max; i += STEP){
    if (mask->data[i] < 0.99){
      continue;
    }
    count += 1.0;
    cur = in1[i] - in2[i];
    cor = cur * cur;
    if (cor > noise){
      noise = cor;
    }
  }
  noise *= overfit_correction;
  // Pass through signal greater than noise
  for (i = THREAD; i < max; i += STEP){
    if (fabs(out[i]) > 0.0){
      rcv += 1;
      continue;
    }
    cur = in1[i] + in2[i];
    cor = cur * cur * mask->data[i];
    if (cor > noise){
      out[i] = cur;
      rcv += 1;
    }
  }
  return rcv / count;
}
